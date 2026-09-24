# ============================================================================
# 04_evi_checkpoint.R  (OPTIONAL)
# Intermediate EVI point (Landsat 8, Jan-Mar 2014-2015) for trajectory
# visualization only. NOT used as a predictor or response in the GWR model
# fitted in 09_analysis.R -- this script exists purely to support a
# descriptive figure showing the EVI trajectory between Epoch 1 and
# Epoch 2. It can be skipped entirely without affecting the regression
# results.
# ============================================================================
#
# REQUIRES: 00_config.R sourced first, and 02_evi_epoch1.R already run
#           (02_evi_epoch1_results.R present, for EPOCH1_N_TARGET -- the
#           same cumulative valid-pixel matching rule used in
#           03_evi_epoch2.R is applied here, for consistency).
#
# Scene list source: USGS EarthExplorer export, study Path/Row, 8 scenes.
#
#   LC08_L2SP_009053_20140108_02_T1   2014-01-08   cloud 18.05%
#   LC08_L2SP_009053_20140124_02_T1   2014-01-24   cloud  0.43%
#   LC08_L2SP_009053_20140209_02_T1   2014-02-09   cloud  5.23%
#   LC08_L2SP_009053_20140225_02_T1   2014-02-25   cloud  5.78%
#   LC08_L2SP_009053_20140313_02_T1   2014-03-13   cloud  6.56%
#   LC08_L2SP_009053_20140329_02_T1   2014-03-29   cloud  0.04%
#   LC08_L2SP_009053_20150111_02_T1   2015-01-11   cloud  0.02%
#   LC08_L2SP_009053_20150127_02_T1   2015-01-27   cloud 17.89%
#
# Note: this list deliberately omits LC08_L2SP_009053_20150212,
# _20150228, and _20150316 (present in the source USGS export) because
# they duplicate scenes already used as OLI members of the ETM+/OLI
# harmonization pairs (pair_04, pair_06, pair_09 in 01_harmonization.R).
# Reusing them here would mix calibration-window scenes into an
# independent descriptive composite.
#
# No harmonization is applied: OLI is the native reference sensor.
#
# DEPENDENCIES: rstac, terra, sf
# ============================================================================

source("00_config.R")
source("02_evi_epoch1_results.R")

library(rstac)
library(terra)
library(sf)

cat("Epoch 1 pixel-count target (from 02_evi_epoch1_results.R):", EPOCH1_N_TARGET, "\n\n")

# Chronologically ordered
checkpoint_scenes <- c(
  "LC08_L2SP_009053_20140108_02_T1",
  "LC08_L2SP_009053_20140124_02_T1",
  "LC08_L2SP_009053_20140209_02_T1",
  "LC08_L2SP_009053_20140225_02_T1",
  "LC08_L2SP_009053_20140313_02_T1",
  "LC08_L2SP_009053_20140329_02_T1",
  "LC08_L2SP_009053_20150111_02_T1",
  "LC08_L2SP_009053_20150127_02_T1"
)

study_area <- vect(STUDY_AREA_SHP)
aoi_wgs    <- st_as_sfc(st_bbox(BBOX_WGS84, crs = 4326))

# ----------------------------------------------------------------------------
# Download helpers
# ----------------------------------------------------------------------------

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
    cat("    [skip]", out_name, "\n"); return(out_file)
  }
  item <- fetch_item_by_id(scene_id)
  if (is.null(item)) { cat("    ERROR: scene not found\n"); return(NULL) }
  url <- get_asset_url(item, asset_name)
  if (is.null(url)) { cat("    ERROR: asset '", asset_name, "' not available\n", sep=""); return(NULL) }
  ok <- tryCatch({
    r <- rast(paste0("/vsicurl/", url))
    aoi_utm <- st_transform(aoi_wgs, crs(r))
    r_crop <- crop(r, vect(aoi_utm))
    writeRaster(r_crop, out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    cat("    OK:", out_name, "\n"); TRUE
  }, error = function(e) { cat("    ERROR:", conditionMessage(e), "\n"); FALSE })
  if (ok) out_file else NULL
}

cat("=== STEP 1: DOWNLOAD ===\n\n")
for (scene_id in checkpoint_scenes) {
  cat(scene_id, ":\n")
  download_asset(scene_id, "blue", "BLUE", DIR_EVI_RAW)
  download_asset(scene_id, "red", "RED", DIR_EVI_RAW)
  download_asset(scene_id, "nir08", "NIR", DIR_EVI_RAW)
  download_asset(scene_id, "qa_pixel", "QA", DIR_EVI_RAW)
  cat("\n")
}

# ----------------------------------------------------------------------------
# Crop, mask, EVI, cumulative pixel matching against EPOCH1_N_TARGET
# ----------------------------------------------------------------------------

scale_sr <- function(dn) clamp(dn * SR_SCALE + SR_OFFSET, 0, 1)

compute_evi <- function(nir, red, blue) {
  evi <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  clamp(evi, -1, 1)
}

cat("=== STEPS 2-5: CROP, MASK, EVI, CUMULATIVE PIXEL MATCHING ===\n\n")

evi_list <- list()
pixel_log <- data.frame()
cumulative_valid <- 0
target_reached <- FALSE

for (scene_id in checkpoint_scenes) {
  cat("Processing:", scene_id, "\n")

  blue_f <- file.path(DIR_EVI_RAW, paste0(scene_id, "_BLUE.tif"))
  red_f  <- file.path(DIR_EVI_RAW, paste0(scene_id, "_RED.tif"))
  nir_f  <- file.path(DIR_EVI_RAW, paste0(scene_id, "_NIR.tif"))
  qa_f   <- file.path(DIR_EVI_RAW, paste0(scene_id, "_QA.tif"))

  if (!all(file.exists(c(blue_f, red_f, nir_f, qa_f)))) {
    cat("  SKIP: missing downloaded files\n\n")
    next
  }

  blue <- rast(blue_f); red <- rast(red_f); nir <- rast(nir_f); qa <- rast(qa_f)

  study_area_r <- if (crs(study_area) != crs(qa)) project(study_area, crs(qa)) else study_area
  blue_c <- crop(blue, study_area_r, mask = TRUE)
  red_c  <- crop(red,  study_area_r, mask = TRUE)
  nir_c  <- crop(nir,  study_area_r, mask = TRUE)
  qa_c   <- crop(qa,   study_area_r, mask = TRUE)

  qa_mask <- qa_c == QA_GOOD_VALUE_OLI
  n_valid <- global(qa_mask, "sum", na.rm = TRUE)[1, 1]
  n_total <- global(!is.na(qa_c), "sum", na.rm = TRUE)[1, 1]

  cat("  Valid pixels (QA==", QA_GOOD_VALUE_OLI, "):", n_valid, "/", n_total, "\n", sep = "")

  blue_sr <- scale_sr(blue_c); red_sr <- scale_sr(red_c); nir_sr <- scale_sr(nir_c)
  blue_m <- ifel(qa_mask, blue_sr, NA)
  red_m  <- ifel(qa_mask, red_sr,  NA)
  nir_m  <- ifel(qa_mask, nir_sr,  NA)

  evi <- compute_evi(nir_m, red_m, blue_m)
  names(evi) <- scene_id
  writeRaster(evi, file.path(DIR_EVI_OUT, paste0(scene_id, "_EVI.tif")),
              overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  included_in_composite <- !target_reached
  if (included_in_composite) {
    evi_list[[scene_id]] <- evi
    cumulative_valid <- cumulative_valid + n_valid
  }

  pixel_log <- rbind(pixel_log, data.frame(
    scene_id = scene_id, n_valid = n_valid, n_total = n_total,
    cumulative_valid = ifelse(included_in_composite, cumulative_valid, NA),
    included_in_composite = included_in_composite
  ))

  cat("  Cumulative (composite scenes only):", cumulative_valid,
      "/ target", EPOCH1_N_TARGET, "\n")

  if (!target_reached && cumulative_valid >= EPOCH1_N_TARGET) {
    cat("  *** TARGET REACHED at this scene. Subsequent scenes excluded",
        "from composite. ***\n")
    target_reached <- TRUE
  }
  cat("\n")
}

cat("=== PIXEL LOG ===\n")
print(pixel_log, row.names = FALSE)
write.csv(pixel_log, file.path(DIR_EVI, "checkpoint_pixel_log.csv"), row.names = FALSE)

cat("\nScenes included in composite:", length(evi_list), "of", nrow(pixel_log), "\n")
cat("Final cumulative valid pixels:", cumulative_valid, "(target:", EPOCH1_N_TARGET, ")\n\n")

if (!target_reached) {
  cat("NOTE: target not reached with all 8 available scenes (this is",
      "expected -- fewer cloud-free Landsat 8 acquisitions exist in this",
      "short 2014-2015 window than in the full Epoch 1/Epoch 2 windows).",
      "Composite built from all", length(evi_list), "scenes.\n\n")
}

# ----------------------------------------------------------------------------
# Composite
# ----------------------------------------------------------------------------

ref <- evi_list[[1]]
evi_aligned <- lapply(evi_list, function(r) {
  if (!compareGeom(r, ref, stopOnError = FALSE)) resample(r, ref, method = "bilinear") else r
})
evi_stack <- rast(evi_aligned)
checkpoint_mean <- mean(evi_stack, na.rm = TRUE)
names(checkpoint_mean) <- "CHECKPOINT_EVI_MEAN"

writeRaster(checkpoint_mean, file.path(DIR_EVI_OUT, "CHECKPOINT_EVI_MEAN.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("Composite stats:\n")
print(global(checkpoint_mean, c("min","max","mean","sd"), na.rm = TRUE))

cat("\n04_evi_checkpoint.R complete (optional script -- not used in 09_analysis.R).\n")