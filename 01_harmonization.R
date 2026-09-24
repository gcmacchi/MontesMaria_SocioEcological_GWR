# ============================================================================
# 01_harmonization.R
# Local ETM+ (Landsat 7) -> OLI (Landsat 8) EVI cross-sensor harmonization.
# ============================================================================
#
# PURPOSE:
#   Test empirically whether the CONUS-derived cross-sensor calibration of
#   Roy et al. (2016) is appropriate for this tropical dry forest study
#   area, by deriving a locally grounded EVI-to-EVI harmonization function
#   from five near-coincident (8-day) Landsat 7 / Landsat 8 scene pairs
#   acquired directly over the study area.
#
# REQUIRES: 00_config.R sourced first (directories, thresholds, bbox).
#
# OUTPUT: in addition to the usual raster/report outputs under
#   DIR_HARMONIZATION, this script writes a small results file to the
#   repository root, "01_harmonization_results.R", containing only the
#   final pooled calibration coefficients as R constants. Downstream
#   scripts (02_evi_epoch1.R) source that file rather than re-typing the
#   coefficients, so the two can never silently drift apart.
#
# DEPENDENCIES: rstac, terra, sf
# ============================================================================

source("00_config.R")

library(rstac)
library(terra)
library(sf)

cat("=== ETM+ / OLI LOCAL HARMONIZATION PIPELINE ===\n")
cat("Working directory:", getwd(), "\n")
cat("Study area shapefile:", STUDY_AREA_SHP, "\n")

study_area <- vect(STUDY_AREA_SHP)
aoi_wgs    <- st_as_sfc(st_bbox(BBOX_WGS84, crs = 4326))

# ----------------------------------------------------------------------------
# Pair registry: five near-coincident (8-day) Landsat 7 / Landsat 8 pairs,
# selected from an exhaustive STAC search over the study Path/Row
# (2013-2022) as the lowest-combined-cloud-cover candidates.
# ----------------------------------------------------------------------------

pair_registry <- data.frame(
  pair_id = c("pair_01", "pair_04", "pair_06", "pair_09", "pair_07"),
  l7_id = c(
    "LE07_L2SP_009053_20150103_02_T1",
    "LE07_L2SP_009053_20140217_02_T1",
    "LE07_L2SP_009053_20140217_02_T1",
    "LE07_L2SP_009053_20150308_02_T1",
    "LE07_L2SP_009053_20160122_02_T1"
  ),
  l8_id = c(
    "LC08_L2SP_009053_20150111_02_T1",
    "LC08_L2SP_009053_20140209_02_T1",
    "LC08_L2SP_009053_20140225_02_T1",
    "LC08_L2SP_009053_20150316_02_T1",
    "LC08_L2SP_009053_20160114_02_T1"
  ),
  stringsAsFactors = FALSE
)

cat("Pairs to process:", nrow(pair_registry), "\n\n")

# ============================================================================
# STAC / DOWNLOAD HELPERS
# ============================================================================

fetch_item_by_id <- function(scene_id) {
  res <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
    stac_search(collections = "landsat-c2-l2", ids = scene_id) |>
    get_request() |>
    items_sign(sign_fn = sign_planetary_computer())
  if (length(res$features) == 0) return(NULL)
  res$features[[1]]
}

get_asset_url <- function(item, name) {
  href <- item$assets[[name]]$href
  if (is.null(href) || !nzchar(href)) NULL else href
}

download_asset <- function(scene_id, asset_name, out_name, out_dir) {
  out_file <- file.path(out_dir, paste0(scene_id, "_", out_name, ".tif"))
  if (file.exists(out_file) && file.size(out_file) > 1000) {
    cat("    [skip]", out_name, "already exists\n")
    return(out_file)
  }
  item <- fetch_item_by_id(scene_id)
  if (is.null(item)) { cat("    ERROR: scene not found on STAC\n"); return(NULL) }
  url <- get_asset_url(item, asset_name)
  if (is.null(url)) { cat("    ERROR: asset '", asset_name, "' not available\n", sep = ""); return(NULL) }
  ok <- tryCatch({
    r <- rast(paste0("/vsicurl/", url))
    aoi_utm <- st_transform(aoi_wgs, crs(r))
    r_crop <- crop(r, vect(aoi_utm))
    writeRaster(r_crop, out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    cat("    OK:", out_name, "\n")
    TRUE
  }, error = function(e) { cat("    ERROR:", conditionMessage(e), "\n"); FALSE })
  if (ok) out_file else NULL
}

# ============================================================================
# STEP 1: DOWNLOAD (Blue, Red, NIR, QA_PIXEL, ATMOS_OPACITY/AEROSOL_QA)
# ============================================================================

cat("========================================\n")
cat("STEP 1: DOWNLOADING BANDS FOR ALL PAIRS\n")
cat("========================================\n\n")

for (i in seq_len(nrow(pair_registry))) {
  pair_id <- pair_registry$pair_id[i]
  l7_id   <- pair_registry$l7_id[i]
  l8_id   <- pair_registry$l8_id[i]

  pair_raw_dir <- file.path(DIR_HARM_RAW, pair_id)
  dir.create(pair_raw_dir, recursive = TRUE, showWarnings = FALSE)

  cat("--- ", pair_id, " ---\n", sep = "")
  cat("L7 (", l7_id, "):\n", sep = "")
  download_asset(l7_id, "blue", "BLUE", pair_raw_dir)
  download_asset(l7_id, "red", "RED", pair_raw_dir)
  download_asset(l7_id, "nir08", "NIR", pair_raw_dir)
  download_asset(l7_id, "qa_pixel", "QA", pair_raw_dir)
  download_asset(l7_id, "atmos_opacity", "ATMOS_OPACITY", pair_raw_dir)

  cat("L8 (", l8_id, "):\n", sep = "")
  download_asset(l8_id, "blue", "BLUE", pair_raw_dir)
  download_asset(l8_id, "red", "RED", pair_raw_dir)
  download_asset(l8_id, "nir08", "NIR", pair_raw_dir)
  download_asset(l8_id, "qa_pixel", "QA", pair_raw_dir)
  download_asset(l8_id, "qa_aerosol", "AEROSOL_QA", pair_raw_dir)
  cat("\n")
}

cat("Step 1 complete.\n\n")

# ============================================================================
# STEP 2: CROP TO STUDY AREA
# ============================================================================

cat("========================================\n")
cat("STEP 2: CROPPING TO STUDY AREA\n")
cat("========================================\n\n")

crop_and_save <- function(scene_id, band_name, pair_raw_dir, pair_crop_dir, study_area) {
  in_file  <- file.path(pair_raw_dir, paste0(scene_id, "_", band_name, ".tif"))
  out_file <- file.path(pair_crop_dir, paste0(scene_id, "_", band_name, "_cropped.tif"))
  if (!file.exists(in_file)) { cat("  MISSING:", basename(in_file), "\n"); return(NULL) }
  r <- rast(in_file)
  study_area_r <- if (crs(study_area) != crs(r)) project(study_area, crs(r)) else study_area
  r_crop <- crop(r, study_area_r)
  writeRaster(r_crop, out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("  OK:", basename(out_file), "\n")
  r_crop
}

for (i in seq_len(nrow(pair_registry))) {
  pair_id <- pair_registry$pair_id[i]
  l7_id   <- pair_registry$l7_id[i]
  l8_id   <- pair_registry$l8_id[i]

  pair_raw_dir  <- file.path(DIR_HARM_RAW, pair_id)
  pair_crop_dir <- file.path(DIR_HARM_CROP, pair_id)
  dir.create(pair_crop_dir, recursive = TRUE, showWarnings = FALSE)

  cat("--- ", pair_id, " ---\n", sep = "")
  bands_l7 <- c("BLUE", "RED", "NIR", "QA", "ATMOS_OPACITY")
  bands_l8 <- c("BLUE", "RED", "NIR", "QA", "AEROSOL_QA")

  for (b in bands_l7) crop_and_save(l7_id, b, pair_raw_dir, pair_crop_dir, study_area)
  for (b in bands_l8) crop_and_save(l8_id, b, pair_raw_dir, pair_crop_dir, study_area)

  l7_blue_c <- rast(file.path(pair_crop_dir, paste0(l7_id, "_BLUE_cropped.tif")))
  l8_blue_c <- rast(file.path(pair_crop_dir, paste0(l8_id, "_BLUE_cropped.tif")))
  geom_ok <- compareGeom(l7_blue_c, l8_blue_c, stopOnError = FALSE)
  cat("  Geometry match (L7 vs L8):", geom_ok, "\n\n")
}

cat("Step 2 complete.\n\n")

# ============================================================================
# STEP 3: BUILD HARMONIZATION MASK + COMPUTE EVI + EXTRACT PAIRED PIXELS
# ============================================================================

cat("========================================\n")
cat("STEP 3: MASKING, EVI COMPUTATION, PIXEL EXTRACTION\n")
cat("========================================\n\n")

scale_sr <- function(dn) clamp(dn * SR_SCALE + SR_OFFSET, 0, 1)

compute_evi <- function(nir, red, blue) {
  evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  clamp(evi, -1, 1)
}

process_pair <- function(pair_id, l7_id, l8_id, pair_crop_dir, pair_evi_dir) {
  cat("--- ", pair_id, " ---\n", sep = "")

  l7_blue <- rast(file.path(pair_crop_dir, paste0(l7_id, "_BLUE_cropped.tif")))
  l7_red  <- rast(file.path(pair_crop_dir, paste0(l7_id, "_RED_cropped.tif")))
  l7_nir  <- rast(file.path(pair_crop_dir, paste0(l7_id, "_NIR_cropped.tif")))
  l7_qa   <- rast(file.path(pair_crop_dir, paste0(l7_id, "_QA_cropped.tif")))
  l7_op   <- rast(file.path(pair_crop_dir, paste0(l7_id, "_ATMOS_OPACITY_cropped.tif")))

  l8_blue <- rast(file.path(pair_crop_dir, paste0(l8_id, "_BLUE_cropped.tif")))
  l8_red  <- rast(file.path(pair_crop_dir, paste0(l8_id, "_RED_cropped.tif")))
  l8_nir  <- rast(file.path(pair_crop_dir, paste0(l8_id, "_NIR_cropped.tif")))
  l8_qa   <- rast(file.path(pair_crop_dir, paste0(l8_id, "_QA_cropped.tif")))
  l8_aer  <- rast(file.path(pair_crop_dir, paste0(l8_id, "_AEROSOL_QA_cropped.tif")))

  cond_l7_qa <- l7_qa == QA_GOOD_VALUE_ETM
  cond_l8_qa <- l8_qa == QA_GOOD_VALUE_OLI
  cond_l7_op <- l7_op < L7_ATMOS_OPACITY_MAX
  cond_l8_ae <- (l8_aer == L8_AEROSOL_GOOD_VALUES[1]) |
                (l8_aer == L8_AEROSOL_GOOD_VALUES[2]) |
                (l8_aer == L8_AEROSOL_GOOD_VALUES[3])

  eligible <- cond_l7_qa & cond_l8_qa & cond_l7_op & cond_l8_ae
  n_eligible <- global(eligible, "sum", na.rm = TRUE)[1, 1]
  cat("  Eligible pixels:", n_eligible, "\n")

  if (is.na(n_eligible) || n_eligible == 0) {
    cat("  WARNING: zero eligible pixels -- pair excluded from pooled estimate.\n\n")
    return(NULL)
  }

  l7_blue_sr <- scale_sr(l7_blue); l7_red_sr <- scale_sr(l7_red); l7_nir_sr <- scale_sr(l7_nir)
  l8_blue_sr <- scale_sr(l8_blue); l8_red_sr <- scale_sr(l8_red); l8_nir_sr <- scale_sr(l8_nir)

  l7_blue_m <- ifel(eligible, l7_blue_sr, NA); l7_red_m <- ifel(eligible, l7_red_sr, NA); l7_nir_m <- ifel(eligible, l7_nir_sr, NA)
  l8_blue_m <- ifel(eligible, l8_blue_sr, NA); l8_red_m <- ifel(eligible, l8_red_sr, NA); l8_nir_m <- ifel(eligible, l8_nir_sr, NA)

  evi_l7 <- compute_evi(l7_nir_m, l7_red_m, l7_blue_m); names(evi_l7) <- paste0(l7_id, "_EVI")
  evi_l8 <- compute_evi(l8_nir_m, l8_red_m, l8_blue_m); names(evi_l8) <- paste0(l8_id, "_EVI")

  writeRaster(evi_l7, file.path(pair_evi_dir, paste0(l7_id, "_EVI_native.tif")), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(evi_l8, file.path(pair_evi_dir, paste0(l8_id, "_EVI_native.tif")), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(eligible, file.path(pair_evi_dir, paste0("HARMONIZATION_MASK_", pair_id, ".tif")), overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  v7 <- values(evi_l7, na.rm = FALSE)
  v8 <- values(evi_l8, na.rm = FALSE)
  ok <- !is.na(v7) & !is.na(v8)
  d <- data.frame(pair_id = pair_id, evi_l7 = as.numeric(v7[ok]), evi_l8 = as.numeric(v8[ok]))

  cat("  Paired pixels extracted:", nrow(d), "\n\n")
  d
}

all_pairs_data <- list()

for (i in seq_len(nrow(pair_registry))) {
  pair_id <- pair_registry$pair_id[i]
  pair_crop_dir <- file.path(DIR_HARM_CROP, pair_id)
  pair_evi_dir  <- file.path(DIR_HARM_EVI, pair_id)
  dir.create(pair_evi_dir, recursive = TRUE, showWarnings = FALSE)

  result <- process_pair(pair_id, pair_registry$l7_id[i], pair_registry$l8_id[i],
                          pair_crop_dir, pair_evi_dir)
  if (!is.null(result)) all_pairs_data[[pair_id]] <- result
}

cat("Step 3 complete.", length(all_pairs_data), "of", nrow(pair_registry), "pairs usable.\n\n")

# ============================================================================
# STEP 4: PER-PAIR STATISTICS, POOLED ESTIMATE, ROY ET AL. COMPARISON
# ============================================================================

cat("========================================\n")
cat("STEP 4: STATISTICAL SYNTHESIS\n")
cat("========================================\n\n")

pooled <- do.call(rbind, all_pairs_data)

per_pair_stats <- data.frame()
for (pid in names(all_pairs_data)) {
  d <- all_pairs_data[[pid]]
  m <- lm(evi_l8 ~ evi_l7, data = d)
  s <- summary(m)
  per_pair_stats <- rbind(per_pair_stats, data.frame(
    pair_id = pid,
    n_pixels = nrow(d),
    r2 = s$r.squared,
    slope = unname(coef(m)[2]),
    intercept = unname(coef(m)[1]),
    delta_mean = mean(d$evi_l8 - d$evi_l7)
  ))
}

cat("Per-pair results:\n")
print(per_pair_stats, row.names = FALSE)

w <- per_pair_stats$n_pixels
pooled_slope     <- sum(per_pair_stats$slope * w) / sum(w)
pooled_intercept <- sum(per_pair_stats$intercept * w) / sum(w)
pooled_delta     <- sum(per_pair_stats$delta_mean * w) / sum(w)
pooled_r2        <- sum(per_pair_stats$r2 * w) / sum(w)
pooled_n         <- sum(w)

pct_diff_slope     <- 100 * (pooled_slope - ROY_NIR_SLOPE) / ROY_NIR_SLOPE
pct_diff_intercept <- 100 * (pooled_intercept - ROY_NIR_INTERCEPT) / ROY_NIR_INTERCEPT

mean_slope <- mean(per_pair_stats$slope); sd_slope <- sd(per_pair_stats$slope)
mean_delta <- mean(per_pair_stats$delta_mean); sd_delta <- sd(per_pair_stats$delta_mean)
se_delta <- sd_delta / sqrt(nrow(per_pair_stats))
t_crit <- qt(0.975, df = nrow(per_pair_stats) - 1)
ci_low  <- mean_delta - t_crit * se_delta
ci_high <- mean_delta + t_crit * se_delta

# ============================================================================
# STEP 5: WRITE SYNTHESIS REPORT
# ============================================================================

cat("\n========================================\n")
cat("STEP 5: WRITING SYNTHESIS REPORT\n")
cat("========================================\n\n")

L <- character()
add <- function(...) L <<- c(L, paste0(...))

add("====================================================================")
add("   LOCAL ETM+ -> OLI EVI HARMONIZATION -- SYNTHESIS REPORT")
add("====================================================================")
add("")
add(sprintf("Pairs processed: %d of %d", length(all_pairs_data), nrow(pair_registry)))
add(sprintf("Total pooled pixels: %d", pooled_n))
add("")
add("--- PER-PAIR RESULTS ---")
for (i in seq_len(nrow(per_pair_stats))) {
  r <- per_pair_stats[i, ]
  add(sprintf("%-10s | n=%9d | R2=%.4f | slope=%.4f | intercept=%.4f | delta=%+.4f",
              r$pair_id, r$n_pixels, r$r2, r$slope, r$intercept, r$delta_mean))
}
add("")
add("--- PIXEL-WEIGHTED POOLED ESTIMATE ---")
add(sprintf("EVI_L8 = %.6f + %.6f * EVI_L7   (pooled R2 = %.4f)", pooled_intercept, pooled_slope, pooled_r2))
add("")
add("--- COMPARISON WITH ROY ET AL. (2016), NIR band, CONUS ---")
add(sprintf("Roy et al. slope = %.4f | Pooled slope = %.4f (%+.1f%%)", ROY_NIR_SLOPE, pooled_slope, pct_diff_slope))
add(sprintf("Roy et al. intercept = %.4f | Pooled intercept = %.4f (%+.1f%%)", ROY_NIR_INTERCEPT, pooled_intercept, pct_diff_intercept))
add("")
add("--- INTER-PAIR VARIABILITY ---")
add(sprintf("Slope: mean=%.4f, SD=%.4f, range=[%.4f, %.4f]", mean_slope, sd_slope, min(per_pair_stats$slope), max(per_pair_stats$slope)))
add(sprintf("Delta: mean=%.4f, SD=%.4f, 95%% CI=[%.4f, %.4f]", mean_delta, sd_delta, ci_low, ci_high))
add("")
add("====================================================================")
add("                         END REPORT")
add("====================================================================")

cat(paste(L, collapse = "\n"), "\n")

writeLines(L, file.path(DIR_HARM_SYNTH, "HARMONIZATION_SYNTHESIS_REPORT.txt"))
write.csv(per_pair_stats, file.path(DIR_HARM_SYNTH, "harmonization_per_pair_stats.csv"), row.names = FALSE)
write.csv(pooled, file.path(DIR_HARM_SYNTH, "harmonization_pooled_pixels.csv"), row.names = FALSE)

# ============================================================================
# STEP 6: WRITE RESULTS FILE FOR DOWNSTREAM SCRIPTS
# ============================================================================
# 02_evi_epoch1.R sources this file to obtain the calibration coefficients.
# Keeping this as a small, separate, auto-generated file (rather than
# re-typing the numbers in 02_evi_epoch1.R) guarantees the two scripts can
# never silently disagree on the calibration used.
# ============================================================================

results_file <- file.path(DIR_ROOT, "01_harmonization_results.R")

RL <- character()
addr <- function(...) RL <<- c(RL, paste0(...))

addr("# ============================================================================")
addr("# AUTO-GENERATED by 01_harmonization.R -- do not edit by hand.")
addr("# Local ETM+ -> OLI cross-sensor calibration, pixel-weighted pooled fit.")
addr(sprintf("# Generated: %s", Sys.time()))
addr("# ============================================================================")
addr("")
addr(sprintf("LOCAL_HARMONIZATION_INTERCEPT <- %.6f", pooled_intercept))
addr(sprintf("LOCAL_HARMONIZATION_SLOPE     <- %.6f", pooled_slope))
addr(sprintf("LOCAL_HARMONIZATION_R2        <- %.4f", pooled_r2))
addr(sprintf("LOCAL_HARMONIZATION_N         <- %d", pooled_n))

writeLines(RL, results_file)

cat("\nAll harmonization outputs saved to:", DIR_HARM_SYNTH, "\n")
cat("Results file for downstream scripts written to:", results_file, "\n")
cat("Pipeline complete.\n")