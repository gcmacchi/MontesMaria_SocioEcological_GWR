# ============================================================================
# 06_chirps_seasonal.R
# Builds the DJFM (dry-season) and AM seasonal composites for the baseline
# and recent periods, and their percent-anomaly difference.
# ============================================================================
#
# REQUIRES: 00_config.R sourced first, and 05_chirps_download.R already run
#           (monthly cropped rasters present in DIR_CHIRPS_MONTHLY).
#
# DJFM: December (year-1) + January + February + March (year), matching the
#       seasonal-year convention used throughout this pipeline.
# AM:   April + May, same calendar year (no lag).
#
# DEPENDENCIES: terra
# ============================================================================

source("00_config.R")

library(terra)

load_month <- function(year, month) {
  f <- file.path(DIR_CHIRPS_MONTHLY, sprintf("CHIRPS_%d_%02d.tif", year, month))
  if (!file.exists(f)) stop("Missing monthly file: ", f)
  rast(f)
}

# ----------------------------------------------------------------------------
# DJFM: sum Dec(y-1)+Jan+Feb+Mar(y) per year, then average across the years
# of each epoch
# ----------------------------------------------------------------------------

build_djfm_mean <- function(years, label) {
  cat("Building DJFM composite for", label, "(years", min(years), "-", max(years), ")\n")
  yearly_rasters <- list()
  for (y in years) {
    r <- load_month(y - 1, 12) + load_month(y, 1) + load_month(y, 2) + load_month(y, 3)
    names(r) <- paste0("DJFM_", y)
    yearly_rasters[[as.character(y)]] <- r
    cat("  DJFM", y, ": mean =", round(global(r, "mean", na.rm=TRUE)[1,1], 2), "mm\n")
  }
  stack <- rast(yearly_rasters)
  mean_r <- mean(stack, na.rm = TRUE)
  names(mean_r) <- paste0("DJFM_mean_", label)
  mean_r
}

djfm_baseline <- build_djfm_mean(CHIRPS_BASELINE_YEARS, "baseline")
djfm_recent   <- build_djfm_mean(CHIRPS_RECENT_YEARS, "recent")

writeRaster(djfm_baseline, file.path(DIR_CHIRPS_SEASONAL, "DJFM_mean_baseline.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(djfm_recent, file.path(DIR_CHIRPS_SEASONAL, "DJFM_mean_recent.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("\nDJFM Baseline mean:", round(global(djfm_baseline, "mean", na.rm=TRUE)[1,1], 2), "mm\n")
cat("DJFM Recent mean:  ", round(global(djfm_recent, "mean", na.rm=TRUE)[1,1], 2), "mm\n\n")

# ----------------------------------------------------------------------------
# AM: sum Apr+May (same year), then average across the years of each epoch
# ----------------------------------------------------------------------------

build_am_mean <- function(years, label) {
  cat("Building AM composite for", label, "(years", min(years), "-", max(years), ")\n")
  yearly_rasters <- list()
  for (y in years) {
    r <- load_month(y, 4) + load_month(y, 5)
    names(r) <- paste0("AM_", y)
    yearly_rasters[[as.character(y)]] <- r
    cat("  AM", y, ": mean =", round(global(r, "mean", na.rm=TRUE)[1,1], 2), "mm\n")
  }
  stack <- rast(yearly_rasters)
  mean_r <- mean(stack, na.rm = TRUE)
  names(mean_r) <- paste0("AM_mean_", label)
  mean_r
}

am_baseline <- build_am_mean(CHIRPS_BASELINE_YEARS, "baseline")
am_recent   <- build_am_mean(CHIRPS_RECENT_YEARS, "recent")

writeRaster(am_baseline, file.path(DIR_CHIRPS_SEASONAL, "AM_mean_baseline.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(am_recent, file.path(DIR_CHIRPS_SEASONAL, "AM_mean_recent.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("\nAM Baseline mean:", round(global(am_baseline, "mean", na.rm=TRUE)[1,1], 2), "mm\n")
cat("AM Recent mean:  ", round(global(am_recent, "mean", na.rm=TRUE)[1,1], 2), "mm\n\n")

# ----------------------------------------------------------------------------
# Differences (Recent - Baseline), absolute and percent. The percent-
# difference rasters are the ones consumed downstream by 07_chirps_to_mau.R.
# ----------------------------------------------------------------------------

djfm_diff_mm   <- djfm_recent - djfm_baseline
djfm_diff_pct  <- (djfm_diff_mm / djfm_baseline) * 100

am_diff_mm  <- am_recent - am_baseline
am_diff_pct <- (am_diff_mm / am_baseline) * 100

writeRaster(djfm_diff_mm, file.path(DIR_CHIRPS_SEASONAL, "DJFM_diff_mm.tif"), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(djfm_diff_pct, file.path(DIR_CHIRPS_SEASONAL, "DJFM_diff_pct.tif"), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(am_diff_mm, file.path(DIR_CHIRPS_SEASONAL, "AM_diff_mm.tif"), overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(am_diff_pct, file.path(DIR_CHIRPS_SEASONAL, "AM_diff_pct.tif"), overwrite = TRUE, gdal = c("COMPRESS=LZW"))

cat("=== FINAL DIFFERENCES (Recent - Baseline) ===\n")
cat("DJFM: mean diff =", round(global(djfm_diff_mm, "mean", na.rm=TRUE)[1,1], 2), "mm (",
    round(global(djfm_diff_pct, "mean", na.rm=TRUE)[1,1], 1), "%)\n")
cat("AM:   mean diff =", round(global(am_diff_mm, "mean", na.rm=TRUE)[1,1], 2), "mm (",
    round(global(am_diff_pct, "mean", na.rm=TRUE)[1,1], 1), "%)\n")

cat("\n06_chirps_seasonal.R complete. Outputs saved to:", DIR_CHIRPS_SEASONAL, "\n")