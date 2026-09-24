# =============================================================================
# 10_landcover_forest_mask.R
# Land-cover overlay, forest-mask construction, and forest-restricted
# OLS / Moran's I / GWR re-analysis (Montes de María, Colombia)
# =============================================================================
#
# Inputs (already cropped to study_area and reprojected to EPSG:32618):
#   ./data/mapBioMas/LC_2009_UTM.tif
#   ./data/mapBioMas/LC_2021_UTM.tif
#
# Requires objects/files produced by earlier scripts in the pipeline:
#   - DIFF_EVI_nowater.tif        (from 08_evi_to_mau.R / water-masking step)
#   - mau_grid_full.shp           (from 08_evi_to_mau.R)
#   - GWR_final_results_nowater.shp (from 09_analysis.R, full-area GWR)
#
# Outputs:
#   ./output/landcover/lc_composition_per_mau.csv
#   ./output/landcover/forest_pixel_stats_per_mau.csv
#   ./output/landcover/DIFF_EVI_forest_masked.tif
#   ./output/GWR/GWR_forest_masked_results.shp
#   ./output/GWR/gwr_forest_masked_local_coefficients.csv
#   ./output/GWR/ols_forest_masked_coefficients.csv
# =============================================================================

library(terra)
library(sf)
library(sp)
library(GWmodel)
library(car)
library(ape)

if ("package:R.utils" %in% search()) detach("package:R.utils", unload = TRUE)

# -----------------------------------------------------------------------
# 0. Paths and parameters
# -----------------------------------------------------------------------

PATH_LC_2009      <- "./data/mapBioMas/LC_2009_UTM.tif"
PATH_LC_2021      <- "./data/mapBioMas/LC_2021_UTM.tif"
PATH_DIFF_EVI     <- "./data/evi/DIFF_EVI_nowater.tif"
PATH_MAU_GRID     <- "./data/mau/mau_grid_full.shp"
PATH_GWR_FULL     <- "./output/GWR/GWR_final_results_nowater.shp"

DIR_OUT_LC        <- "./output/landcover"
DIR_OUT_GWR       <- "./output/GWR"
dir.create(DIR_OUT_LC,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT_GWR, recursive = TRUE, showWarnings = FALSE)

FOREST_CODE       <- 3      # MapBiomas "Forest" class used for the forest mask
FOREST_PCT_THRESH <- 10     # minimum % forest cover (2009 baseline) for MAU eligibility

# MapBiomas Colombia Collection 3.0 legend -> aggregated categories
LC_CATEGORY_MAP <- c(
  "3"  = "Forest",            "5"  = "Forest",            "6"  = "Forest",
  "9"  = "Agricultural",      "21" = "Agricultural",      "35" = "Agricultural",
  "11" = "Natural_nonforest",
  "23" = "Nonvegetated",      "24" = "Nonvegetated",
  "25" = "Nonvegetated",      "30" = "Nonvegetated",
  "33" = "Water"
)

# -----------------------------------------------------------------------
# 1. Load rasters and vector layers
# -----------------------------------------------------------------------

lc_2009  <- rast(PATH_LC_2009)
lc_2021  <- rast(PATH_LC_2021)
diff_evi <- rast(PATH_DIFF_EVI)
mau_grid <- vect(PATH_MAU_GRID)

stopifnot(same.crs(lc_2009, diff_evi))
stopifnot(same.crs(lc_2009, mau_grid))

# -----------------------------------------------------------------------
# 2. Land-cover composition per MAU (2009 and 2021)
# -----------------------------------------------------------------------

extract_lc_composition <- function(lc_raster, mau_vect, category_map, prefix) {
  ext_vals <- terra::extract(lc_raster, mau_vect)
  names(ext_vals) <- c("ID", "code")
  ext_vals$category <- category_map[as.character(ext_vals$code)]
  ext_vals$category[is.na(ext_vals$category)] <- "Unclassified"

  comp     <- table(ext_vals$ID, ext_vals$category)
  comp_pct <- prop.table(comp, margin = 1) * 100
  comp_df  <- as.data.frame.matrix(comp_pct)
  names(comp_df) <- paste0(prefix, "_", names(comp_df))
  comp_df$ID <- as.integer(rownames(comp_df))
  comp_df
}

lc09_comp <- extract_lc_composition(lc_2009, mau_grid, LC_CATEGORY_MAP, "LC09")
lc21_comp <- extract_lc_composition(lc_2021, mau_grid, LC_CATEGORY_MAP, "LC21")

lc_composition <- merge(lc09_comp, lc21_comp, by = "ID", all = TRUE)
write.csv(lc_composition,
          file.path(DIR_OUT_LC, "lc_composition_per_mau.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------
# 3. Forest mask (2009 baseline) and ΔEVI masking
# -----------------------------------------------------------------------

if (!compareGeom(lc_2009, diff_evi, stopOnError = FALSE)) {
  lc_2009_on_evi <- resample(lc_2009, diff_evi, method = "near")
} else {
  lc_2009_on_evi <- lc_2009
}

forest_mask      <- lc_2009_on_evi == FOREST_CODE
diff_evi_forest  <- ifel(forest_mask, diff_evi, NA)

writeRaster(diff_evi_forest,
            file.path(DIR_OUT_LC, "DIFF_EVI_forest_masked.tif"),
            overwrite = TRUE)

# -----------------------------------------------------------------------
# 4. Per-MAU forest-pixel statistics (2009 baseline)
# -----------------------------------------------------------------------

forest_mask_int <- as.numeric(forest_mask)  # 1 = forest, 0 = non-forest, NA = outside raster

n_forest_per_mau <- terra::extract(forest_mask_int, mau_grid, fun = sum,   na.rm = TRUE)
n_total_per_mau  <- terra::extract(forest_mask_int, mau_grid, fun = function(x, ...) sum(!is.na(x)))

forest_stats <- data.frame(
  ID        = 1:nrow(mau_grid),
  n_forest  = n_forest_per_mau[, 2],
  n_total   = n_total_per_mau[, 2]
)
forest_stats$pct_forest <- 100 * forest_stats$n_forest / forest_stats$n_total

write.csv(forest_stats,
          file.path(DIR_OUT_LC, "forest_pixel_stats_per_mau.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------
# 5. Attach forest statistics and forest-only ΔEVI to the GWR analysis grid
# -----------------------------------------------------------------------

gwr_full <- vect(PATH_GWR_FULL)

gwr_full$pct_forest <- forest_stats$pct_forest[match(1:nrow(gwr_full), forest_stats$ID)]

evi_forest_per_mau      <- terra::extract(diff_evi_forest, gwr_full, fun = mean, na.rm = TRUE)
gwr_full$DIFF_EVI_forest <- evi_forest_per_mau[, 2]

# -----------------------------------------------------------------------
# 6. Eligibility filter: MAUs with >= FOREST_PCT_THRESH% forest cover
#    and a valid forest-only ΔEVI value
# -----------------------------------------------------------------------

df <- as.data.frame(gwr_full, geom = "XY")

eligible <- df$pct_forest >= FOREST_PCT_THRESH & !is.na(df$DIFF_EVI_forest)
df_f     <- df[eligible, ]

cat(sprintf("Forest-eligible MAUs (>= %d%% forest cover): %d of %d\n",
            FOREST_PCT_THRESH, sum(eligible), nrow(df)))

# -----------------------------------------------------------------------
# 7. OLS on the forest-masked subset
# -----------------------------------------------------------------------

ols_forest <- lm(DIFF_EVI_forest ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_, data = df_f)
print(summary(ols_forest))

write.csv(as.data.frame(summary(ols_forest)$coefficients),
          file.path(DIR_OUT_GWR, "ols_forest_masked_coefficients.csv"))

# -----------------------------------------------------------------------
# 8. Moran's I on OLS residuals (forest-masked subset)
# -----------------------------------------------------------------------

coords_f  <- as.matrix(df_f[, c("x", "y")])
dist_mat  <- as.matrix(dist(coords_f))
inv_dist  <- 1 / dist_mat
diag(inv_dist) <- 0

moran_forest <- Moran.I(residuals(ols_forest), inv_dist, na.rm = TRUE)
print(moran_forest)

# -----------------------------------------------------------------------
# 9. GWR on the forest-masked subset
# -----------------------------------------------------------------------

sp_f <- SpatialPointsDataFrame(coords = coords_f, data = df_f)

bw_forest <- bw.gwr(
  DIFF_EVI_forest ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_,
  data     = sp_f,
  approach = "AICc",
  kernel   = "gaussian",
  adaptive = TRUE
)

gwr_forest <- gwr.basic(
  DIFF_EVI_forest ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_,
  data     = sp_f,
  bw       = bw_forest,
  kernel   = "gaussian",
  adaptive = TRUE
)

print(gwr_forest)

sdf_f <- gwr_forest$SDF
print(names(sdf_f))

# -----------------------------------------------------------------------
# 10. Local t-values, regime classification (sign reported separately)
# -----------------------------------------------------------------------

t_DEM  <- sdf_f$DENS_DIFF  / sdf_f$DENS_DIFF_SE
t_DJFM <- sdf_f$CHIRPS_DJF / sdf_f$CHIRPS_DJF_SE
t_AM   <- sdf_f$CHIRPS_AM_ / sdf_f$CHIRPS_AM__SE

sig_DEM  <- abs(t_DEM)  > 1.96
sig_DJFM <- abs(t_DJFM) > 1.96
sig_AM   <- abs(t_AM)   > 1.96
sig_CLIM <- sig_DJFM | sig_AM

regime <- ifelse(sig_DEM & sig_CLIM, "Dem x Clim",
           ifelse(sig_DEM & !sig_CLIM, "Dem",
            ifelse(!sig_DEM & sig_CLIM, "Clim", "N.s.")))

df_f$t_DEM      <- t_DEM
df_f$t_DJFM     <- t_DJFM
df_f$t_AM       <- t_AM
df_f$sign_DEM   <- sign(sdf_f$DENS_DIFF)
df_f$sign_DJFM  <- sign(sdf_f$CHIRPS_DJF)
df_f$sign_AM    <- sign(sdf_f$CHIRPS_AM_)
df_f$beta_DEM   <- sdf_f$DENS_DIFF
df_f$beta_DJFM  <- sdf_f$CHIRPS_DJF
df_f$beta_AM    <- sdf_f$CHIRPS_AM_
df_f$local_R2   <- sdf_f$Local_R2
df_f$REGIME     <- regime

# -----------------------------------------------------------------------
# 11. Export
# -----------------------------------------------------------------------

write.csv(df_f, file.path(DIR_OUT_GWR, "gwr_forest_masked_local_coefficients.csv"),
          row.names = FALSE)

out_vect <- vect(df_f, geom = c("x", "y"), crs = crs(gwr_full))
writeVector(out_vect,
            file.path(DIR_OUT_GWR, "GWR_forest_masked_results.shp"),
            overwrite = TRUE)

cat("\n10_landcover_forest_mask.R complete.\n")
cat(sprintf("Bandwidth (adaptive, AICc): %d\n", bw_forest))
cat(sprintf("Adjusted R2 (global OLS, forest-only): %.3f\n",
            summary(ols_forest)$adj.r.squared))