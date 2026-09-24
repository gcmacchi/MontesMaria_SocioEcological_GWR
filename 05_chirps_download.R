# ============================================================================
# 05_chirps_download.R
# Downloads the full monthly CHIRPS v2.0 precipitation record and crops it
# to the study area using a fixed, geometrically explicit criterion.
# ============================================================================
#
# REQUIRES: 00_config.R sourced first.
#
# METHODOLOGICAL NOTE ON THE CROP CRITERION:
#   Earlier crop/mask approaches (bounding-box crop, terra::mask() with or
#   without touches=TRUE, terra::cells()) produced inconsistent pixel
#   counts (247-280 pixels for the same raster) specifically along the
#   western and southern edges of the study area, where the polygon
#   boundary falls very close to a native CHIRPS cell edge. To remove this
#   ambiguity, a pixel is retained here if its centroid lies within
#   CHIRPS_BUFFER_M meters of the study area boundary -- a purely
#   geometric rule, independent of cell-alignment edge cases, applied
#   identically to every monthly raster.
#
# DEPENDENCIES: terra, httr, R.utils
# ============================================================================

source("00_config.R")

library(terra)
library(httr)
library(R.utils)

cat("========================================\n")
cat("PHASE 1: DOWNLOAD + CROP\n")
cat("========================================\n\n")

study_area <- vect(STUDY_AREA_SHP)

study_area_utm          <- project(study_area, "EPSG:32618")
study_area_buffered_utm <- buffer(study_area_utm, width = CHIRPS_BUFFER_M)
study_area_buffered_wgs <- project(study_area_buffered_utm, "EPSG:4326")

cat("Buffered study area built (", CHIRPS_BUFFER_M / 1000, " km buffer, EPSG:4326).\n\n", sep = "")

download_and_crop <- function(year, month) {
  mm   <- sprintf("%02d", month)
  name <- sprintf("CHIRPS_%d_%02d", year, month)
  out  <- file.path(DIR_CHIRPS_MONTHLY, paste0(name, ".tif"))

  if (file.exists(out) && file.size(out) > 1000) {
    message("  [skip] ", name)
    return(out)
  }

  fname_gz <- paste0("chirps-v2.0.", year, ".", mm, ".tif.gz")
  dest_gz  <- file.path(DIR_CHIRPS_RAW, fname_gz)
  dest_tif <- file.path(DIR_CHIRPS_RAW, sub("\\.gz$", "", fname_gz))

  urls <- c(paste0(CHIRPS_URL_CAMER, fname_gz), paste0(CHIRPS_URL_GLOBAL, fname_gz))

  ok <- FALSE
  for (url in urls) {
    if (file.exists(dest_gz) && file.size(dest_gz) > 1000) { ok <- TRUE; break }
    message("  Downloading: ", basename(url))
    r <- try(GET(url, write_disk(dest_gz, overwrite = TRUE),
                 timeout(120), progress()), silent = TRUE)
    if (!inherits(r, "try-error") && status_code(r) == 200 &&
        file.exists(dest_gz) && file.size(dest_gz) > 1000) { ok <- TRUE; break }
  }

  if (!ok) { warning("Download failed for: ", name); return(NULL) }

  if (!file.exists(dest_tif) || file.size(dest_tif) < 1000) {
    gunzip(dest_gz, dest_tif, remove = FALSE, overwrite = TRUE)
  }

  r <- rast(dest_tif)
  NAflag(r) <- CHIRPS_NODATA

  r_cropped <- crop(r, study_area_buffered_wgs, mask = TRUE, touches = TRUE)
  names(r_cropped) <- name
  writeRaster(r_cropped, out, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message("  OK: ", name)
  return(out)
}

for (y in CHIRPS_YEARS) {
  for (m in CHIRPS_MONTHS) {
    message(sprintf("--- %d-%02d ---", y, m))
    download_and_crop(y, m)
  }
}

expected <- expand.grid(year = CHIRPS_YEARS, month = CHIRPS_MONTHS)
expected$file <- file.path(DIR_CHIRPS_MONTHLY, sprintf("CHIRPS_%d_%02d.tif", expected$year, expected$month))
expected$exists_ok <- file.exists(expected$file) &
  file.size(ifelse(file.exists(expected$file), expected$file, NA)) > 1000

cat("\n=== COMPLETENESS CHECK ===\n")
cat("Required:", nrow(expected), " | Present and valid:", sum(expected$exists_ok), "\n")
missing_p1 <- expected[!expected$exists_ok, ]
if (nrow(missing_p1) > 0) {
  cat("MISSING FILES:\n"); print(missing_p1[, c("year", "month")])
} else {
  cat("All monthly files present.\n")
}

# ----------------------------------------------------------------------------
# Monthly summary table (valid pixels only, physically implausible values
# excluded as data artefacts, not as valid extreme rainfall events)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("PHASE 2: MONTHLY SUMMARY TABLE\n")
cat("========================================\n\n")

files <- list.files(DIR_CHIRPS_MONTHLY, pattern = "^CHIRPS_\\d{4}_\\d{2}\\.tif$", full.names = TRUE)
cat("Cropped monthly files found:", length(files), "\n\n")

summary_results <- data.frame()

for (f in files) {
  fname <- basename(f)
  ym <- regmatches(fname, regexpr("\\d{4}_\\d{2}", fname))
  year  <- as.integer(substr(ym, 1, 4))
  month <- as.integer(substr(ym, 6, 7))

  r <- rast(f)
  v <- values(r, na.rm = TRUE)

  is_outlier <- v < CHIRPS_LOWER_BOUND | v > CHIRPS_UPPER_BOUND
  v_clean <- v[!is_outlier]

  summary_results <- rbind(summary_results, data.frame(
    year = year, month = month,
    n_valid_clean = length(v_clean),
    n_outliers_excluded = sum(is_outlier, na.rm = TRUE),
    mean_mm = if (length(v_clean) > 0) mean(v_clean) else NA,
    sd_mm   = if (length(v_clean) > 1) sd(v_clean) else NA
  ))
}

summary_results <- summary_results[order(summary_results$year, summary_results$month), ]
write.csv(summary_results, file.path(DIR_CHIRPS, "chirps_monthly_summary_table.csv"), row.names = FALSE)

cat("Months processed:", nrow(summary_results), "\n")
cat("Months with outliers excluded:", sum(summary_results$n_outliers_excluded > 0), "\n")
cat("Saved: chirps_monthly_summary_table.csv\n")

# ----------------------------------------------------------------------------
# Long-term trend test (context only; not used as a model predictor).
# Verifies that the baseline-vs-recent contrast used downstream does not
# conflate with an undetected secular trend.
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("PHASE 3: TREND ANALYSIS (", min(CHIRPS_YEARS), "-", max(CHIRPS_YEARS), ")\n", sep = "")
cat("========================================\n\n")

annual_total <- aggregate(mean_mm ~ year, data = summary_results, FUN = sum)
names(annual_total)[2] <- "annual_total_mm"
annual_counts <- aggregate(mean_mm ~ year, data = summary_results, FUN = length)
names(annual_counts)[2] <- "n_months"
annual_total <- merge(annual_total, annual_counts, by = "year")
annual_complete <- annual_total[annual_total$n_months == 12, ]

m_annual <- lm(annual_total_mm ~ year, data = annual_complete)
s_annual <- summary(m_annual)

summary_results$season_year_djfm <- ifelse(summary_results$month == 12,
                                            summary_results$year + 1, summary_results$year)
djfm <- summary_results[summary_results$month %in% c(12, 1, 2, 3), ]
djfm_annual <- aggregate(mean_mm ~ season_year_djfm, data = djfm, FUN = sum)
names(djfm_annual) <- c("year", "djfm_total_mm")
djfm_counts <- aggregate(mean_mm ~ season_year_djfm, data = djfm, FUN = length)
names(djfm_counts) <- c("year", "n_months")
djfm_annual <- merge(djfm_annual, djfm_counts, by = "year")
djfm_complete <- djfm_annual[djfm_annual$n_months == 4, ]

m_djfm <- lm(djfm_total_mm ~ year, data = djfm_complete)
s_djfm <- summary(m_djfm)

am <- summary_results[summary_results$month %in% c(4, 5), ]
am_annual <- aggregate(mean_mm ~ year, data = am, FUN = sum)
names(am_annual)[2] <- "am_total_mm"
am_counts <- aggregate(mean_mm ~ year, data = am, FUN = length)
names(am_counts)[2] <- "n_months"
am_annual <- merge(am_annual, am_counts, by = "year")
am_complete <- am_annual[am_annual$n_months == 2, ]

m_am <- lm(am_total_mm ~ year, data = am_complete)
s_am <- summary(m_am)

write.csv(annual_complete, file.path(DIR_CHIRPS, "chirps_annual_totals.csv"), row.names = FALSE)
write.csv(djfm_complete, file.path(DIR_CHIRPS, "chirps_djfm_totals.csv"), row.names = FALSE)
write.csv(am_complete, file.path(DIR_CHIRPS, "chirps_am_totals.csv"), row.names = FALSE)

trend_label <- function(slope, p) {
  if (p < 0.05) { if (slope < 0) "SIGNIFICANT DECLINE" else "SIGNIFICANT INCREASE" }
  else "NO SIGNIFICANT TREND"
}

cat(sprintf("Annual total : n=%d years | slope = %+.3f mm/yr | p = %.4f -> %s\n",
            nrow(annual_complete), coef(m_annual)[2],
            s_annual$coefficients[2,4], trend_label(coef(m_annual)[2], s_annual$coefficients[2,4])))
cat(sprintf("DJFM season  : n=%d years | slope = %+.3f mm/yr | p = %.4f -> %s\n",
            nrow(djfm_complete), coef(m_djfm)[2],
            s_djfm$coefficients[2,4], trend_label(coef(m_djfm)[2], s_djfm$coefficients[2,4])))
cat(sprintf("AM window    : n=%d years | slope = %+.3f mm/yr | p = %.4f -> %s\n",
            nrow(am_complete), coef(m_am)[2],
            s_am$coefficients[2,4], trend_label(coef(m_am)[2], s_am$coefficients[2,4])))

# ----------------------------------------------------------------------------
# Baseline vs Recent comparison (uses the SAME windows as the rest of the
# pipeline: CHIRPS_BASELINE_YEARS / CHIRPS_RECENT_YEARS from 00_config.R)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("PHASE 4: BASELINE vs RECENT COMPARISON\n")
cat("========================================\n\n")

compare_epochs <- function(df, value_col, label) {
  base_val <- mean(df[[value_col]][df$year %in% CHIRPS_BASELINE_YEARS])
  rec_val  <- mean(df[[value_col]][df$year %in% CHIRPS_RECENT_YEARS])
  diff_mm  <- rec_val - base_val
  diff_pct <- 100 * diff_mm / base_val
  cat(sprintf("%-8s: baseline = %7.1f mm | recent = %7.1f mm | diff = %+6.1f mm (%+5.1f%%)\n",
              label, base_val, rec_val, diff_mm, diff_pct))
  data.frame(variable = label, baseline_mm = base_val, recent_mm = rec_val,
             diff_mm = diff_mm, diff_pct = diff_pct)
}

epoch_comparison <- rbind(
  compare_epochs(annual_complete, "annual_total_mm", "Annual"),
  compare_epochs(djfm_complete, "djfm_total_mm", "DJFM"),
  compare_epochs(am_complete, "am_total_mm", "AM")
)

write.csv(epoch_comparison, file.path(DIR_CHIRPS, "chirps_epoch_comparison.csv"), row.names = FALSE)

cat("\nBaseline:", min(CHIRPS_BASELINE_YEARS), "-", max(CHIRPS_BASELINE_YEARS),
    "| Recent:", min(CHIRPS_RECENT_YEARS), "-", max(CHIRPS_RECENT_YEARS), "\n")
print(epoch_comparison, row.names = FALSE)
cat("\n05_chirps_download.R complete.\n")