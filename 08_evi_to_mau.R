# ============================================================================
# 08_evi_to_mau.R
# Adds the EVI difference (Recent minus Baseline, cross-sensor harmonized)
# to the MAU grid, completing the analysis dataset.
# ============================================================================
#
# REQUIRES: 00_config.R sourced first; 02_evi_epoch1.R, 03_evi_epoch2.R,
#           and 07_chirps_to_mau.R already run (mau_grid_with_chirps.shp
#           present in DIR_BINDER).
#
# No spatial disaggregation is required here, unlike the CHIRPS extraction
# in 07_chirps_to_mau.R: EVI is derived at native Landsat resolution
# (30 m), which is already finer than the MAU grid.
#
# Sign convention: Recent - Baseline (positive = greening, negative =
# browning).
#
# DEPENDENCIES: terra
# ============================================================================

source("00_config.R")

library(terra)

if ("package:R.utils" %in% search()) detach("package:R.utils", unload = TRUE)

mau_grid_path <- file.path(DIR_BINDER, "mau_grid_with_chirps.shp")
if (!file.exists(mau_grid_path)) {
  stop("Required input not found: ", mau_grid_path,
       "\nRun 07_chirps_to_mau.R first.")
}
mau_grid <- vect(mau_grid_path)

epoch1_file <- "I:/R_Workspace/LANDSAT_PROCESSOR/Landsat_pipeline_v2/02_evi/EPOCH1_EVI_MEAN_harmonized_nowater.tif"
epoch2_file <- "I:/R_Workspace/LANDSAT_PROCESSOR/Landsat_pipeline_v2/02_evi/EPOCH2_EVI_MEAN_nowater.tif"

if (!file.exists(epoch1_file)) stop("Missing: ", epoch1_file, " -- run 02_evi_epoch1.R first.")
if (!file.exists(epoch2_file)) stop("Missing: ", epoch2_file, " -- run 03_evi_epoch2.R first.")

epoch1 <- rast(epoch1_file)
epoch2 <- rast(epoch2_file)

if (!compareGeom(epoch1, epoch2, stopOnError = FALSE)) {
  cat("Geometry mismatch between Epoch 1 and Epoch 2 composites -- resampling.\n")
  epoch2 <- resample(epoch2, epoch1, method = "bilinear")
}

diff_evi <- epoch2 - epoch1
names(diff_evi) <- "DIFF_EVI"

writeRaster(diff_evi, file.path(DIR_EVI_OUT, "DIFF_EVI.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

# Zonal mean of the difference raster within each MAU polygon.
mau_proj <- if (crs(mau_grid) != crs(diff_evi)) project(mau_grid, crs(diff_evi)) else mau_grid
evi_per_mau <- terra::extract(diff_evi, mau_proj, fun = mean, na.rm = TRUE)
mau_grid$DIFF_EVI <- evi_per_mau[, 2]

cat("DIFF_EVI summary (Recent - Baseline):\n")
print(summary(mau_grid$DIFF_EVI))

out_file <- file.path(DIR_BINDER, "mau_grid_full.shp")
writeVector(mau_grid, out_file, overwrite = TRUE)

cat("\nSaved:", out_file, "\n")
cat("Field names on disk:\n")
print(names(vect(out_file)))

cat("\n08_evi_to_mau.R complete. This is the final analysis dataset,\n")
cat("consumed by 09_analysis.R.\n")