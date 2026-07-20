# 04_gwr_analysis.R
# Geographically Weighted Regression – Montes de María
# Author: Giancarlo Macchi Janica
# Date: 2025-11-23

library(terra)
library(ape)
library(GWmodel)
library(sp)

# --- CONFIGURATION ---
# Assumes:
# - diff_evi: raster (Epoch 2 - Epoch 1)
# - mau_grid: SpatVector (4 km dasymetric grid)
# - final_results: data.frame with climate anomalies
# - DENS_DIFF already in mau_grid

# --- STEP 1: Extract ΔEVI to MAU ---
mau_grid$mean_delta_evi <- extract(diff_evi, mau_grid, fun = mean, na.rm = TRUE)[,2]

# --- STEP 2: Merge Climate Data ---
analysis_df <- as.data.frame(mau_grid)
analysis_df <- merge(analysis_df, final_results, by = "id")

# --- STEP 3: Clean ---
analysis_clean <- na.omit(analysis_df[, c("mean_delta_evi", "DENS_DIFF",
                                          "DIFF_DGFM_perc", "DIFF_MGLA_perc")])

# --- STEP 4: Global OLS ---
model_ols <- lm(mean_delta_evi ~ DENS_DIFF + DIFF_DGFM_perc + DIFF_MGLA_perc,
                data = analysis_clean)
summary(model_ols)

# --- STEP 5: Moran's I ---
analysis_clean$residuals <- residuals(model_ols)
coords <- crds(centroids(mau_grid[as.numeric(rownames(analysis_clean)), ]))
dist_matrix <- as.matrix(dist(coords))
dist_inv <- 1 / dist_matrix
diag(dist_inv) <- 0
moran_test <- Moran.I(analysis_clean$residuals, dist_inv)
print(moran_test)

# --- STEP 6: GWR ---
sp_data <- SpatialPointsDataFrame(coords = coords,
                                  data = analysis_clean,
                                  proj4string = CRS(crs(mau_grid)))

bw_adaptive <- bw.gwr(mean_delta_evi ~ DENS_DIFF + DIFF_DGFM_perc + DIFF_MGLA_perc,
                       data = sp_data, approach = "AICc",
                       kernel = "gaussian", adaptive = TRUE)

gwr_model <- gwr.basic(mean_delta_evi ~ DENS_DIFF + DIFF_DGFM_perc + DIFF_MGLA_perc,
                       data = sp_data, bw = bw_adaptive,
                       kernel = "gaussian", adaptive = TRUE)

print(gwr_model)

# --- STEP 7: Regime Classification ---
results_sdf <- as.data.frame(gwr_model$SDF)
t_threshold <- 1.96
t_dem <- results_sdf$DENS_DIFF / results_sdf$DENS_DIFF_SE
t_am  <- results_sdf$DIFF_MGLA_perc / results_sdf$DIFF_MGLA_perc_SE

mau_grid_final <- mau_grid[as.numeric(rownames(analysis_clean)), ]
mau_grid_final$REGIME <- "H0"

mau_grid_final$REGIME[abs(t_dem) > t_threshold & abs(t_am) <= t_threshold] <- "H1"
mau_grid_final$REGIME[abs(t_am) > t_threshold & abs(t_dem) <= t_threshold] <- "H2"
mau_grid_final$REGIME[abs(t_dem) > t_threshold & abs(t_am) > t_threshold] <- "H1xH2"

writeVector(mau_grid_final, "outputs/figures/GWR_FINAL_REGIMES.shp", overwrite = TRUE)

# --- STEP 8: Export Summary ---
write.csv(summary(model_ols)$coefficients, "outputs/tables/OLS_Summary.csv")
write.csv(as.data.frame(gwr_model$SDF), "outputs/tables/GWR_Coefficients.csv")