# ============================================================================
# 07_chirps_to_mau.R
# Area-proportional extraction of CHIRPS precipitation anomalies onto the
# MAU grid, via fine-grid disaggregation.
# ============================================================================
#
# REQUIRES: 00_config.R sourced first, and 06_chirps_seasonal.R already run.
#
# Disaggregates CHIRPS anomaly rasters (0.05 deg native resolution) to a
# 50 m grid via nearest-neighbor resampling (no interpolation), then
# computes the mean of fine-grid cells falling within each MAU polygon.
# This approximates area-weighted extraction (proportional to the fraction
# of each native CHIRPS cell overlapping the MAU) without explicit
# polygon-pixel intersection geometry.
#
# DEPENDENCIES: terra
# ============================================================================

source("00_config.R")

library(terra)

if ("package:R.utils" %in% search()) detach("package:R.utils", unload = TRUE)

mau_grid <- vect(DEMOGRAPHY_SHP)

djfm_diff_pct <- rast(file.path(DIR_CHIRPS_SEASONAL, "DJFM_diff_pct.tif"))
am_diff_pct   <- rast(file.path(DIR_CHIRPS_SEASONAL, "AM_diff_pct.tif"))

# Crop to MAU extent (with buffer) before disaggregation, to avoid
# processing the full CHIRPS extent at fine resolution.
mau_proj   <- if (crs(mau_grid) != crs(djfm_diff_pct)) project(mau_grid, crs(djfm_diff_pct)) else mau_grid
mau_buffer <- buffer(mau_proj, width = 0.05)

djfm_cropped <- crop(djfm_diff_pct, mau_buffer)
am_cropped   <- crop(am_diff_pct, mau_buffer)

# Disaggregate: nearest-neighbor only, to avoid introducing interpolated
# values between native CHIRPS cells.
cat("Disaggregating DJFM and AM to fine grid (nearest-neighbor)...\n")
djfm_fine <- disagg(djfm_cropped, fact = CHIRPS_DISAGG_FACTOR, method = "near")
am_fine   <- disagg(am_cropped, fact = CHIRPS_DISAGG_FACTOR, method = "near")

cat("Fine grid resolution:", res(djfm_fine), "\n\n")

# Extract mean of fine-grid cells per MAU polygon.
mau_proj_fine <- if (crs(mau_grid) != crs(djfm_fine)) project(mau_grid, crs(djfm_fine)) else mau_grid

djfm_per_mau <- terra::extract(djfm_fine, mau_proj_fine, fun = mean, na.rm = TRUE)
am_per_mau   <- terra::extract(am_fine, mau_proj_fine, fun = mean, na.rm = TRUE)

mau_grid$CHIRPS_DJFM_pct <- djfm_per_mau[, 2]
mau_grid$CHIRPS_AM_pct   <- am_per_mau[, 2]

cat("CHIRPS_DJFM_pct summary:\n")
print(summary(mau_grid$CHIRPS_DJFM_pct))
cat("\nCHIRPS_AM_pct summary:\n")
print(summary(mau_grid$CHIRPS_AM_pct))

# NOTE: the .shp format truncates field names to 10 characters. The two
# fields above will be stored as CHIRPS_DJF and CHIRPS_AM_ in the output
# file -- printed explicitly below so this is never a silent surprise.

out_file <- file.path(DIR_BINDER, "mau_grid_with_chirps.shp")
writeVector(mau_grid, out_file, overwrite = TRUE)

cat("\nSaved:", out_file, "\n")
cat("Actual field names on disk:\n")
print(names(vect(out_file)))
cat("\n07_chirps_to_mau.R complete.\n")