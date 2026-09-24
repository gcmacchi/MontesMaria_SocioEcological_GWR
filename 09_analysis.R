# ============================================================================
# 09_analysis.R
# Final model: global OLS baseline and Geographically Weighted Regression
# (GWR). Predictors: demographic change, DJFM and AM precipitation
# anomalies. Response: cross-sensor-harmonized EVI change
# (Recent - Baseline).
# ============================================================================
#
# REQUIRES: 00_config.R sourced first, and 08_evi_to_mau.R already run
#           (mau_grid_full.shp present in DIR_BINDER).
#
# Field names reflect .shp 10-character truncation:
#   CHIRPS_DJFM_pct -> CHIRPS_DJF
#   CHIRPS_AM_pct   -> CHIRPS_AM_
#
# Regime classification: each MAU is classified purely on the statistical
# significance of its local coefficients (|t| > GWR_T_THRESHOLD), into
# four descriptive categories -- Dem (demographic-dominant), Clim
# (climate-dominant), Dem x Clim (co-dominant), N.s. (not significant).
# The sign (direction) of each significant coefficient is recorded
# separately and does NOT affect this classification, so that
# counter-intuitive or mixed-direction effects are reported rather than
# suppressed by the classification scheme itself.
#
# Model diagnostics: in addition to the OLS residual Moran's I test,
# this script also computes Moran's I on the GWR residuals (to quantify
# how much spatial autocorrelation the local model removes relative to
# the global model) and local variance inflation factors via
# GWmodel::gwr.collin.diagno() (to flag localities where the two
# precipitation predictors are not well separated at the chosen
# bandwidth).
#
# DEPENDENCIES: terra, sp, GWmodel, car, ape
# ============================================================================

source("00_config.R")

library(terra)
library(sp)
library(GWmodel)
library(car)
library(ape)

if ("package:R.utils" %in% search()) detach("package:R.utils", unload = TRUE)

mau_grid_path <- file.path(DIR_BINDER, "mau_grid_full.shp")
if (!file.exists(mau_grid_path)) {
  stop("Required input not found: ", mau_grid_path, "\nRun 08_evi_to_mau.R first.")
}

mau_grid <- vect(mau_grid_path)
df <- as.data.frame(mau_grid)

cat("Fields available:\n")
print(names(df))
cat("\n")

# ----------------------------------------------------------------------------
# Build the analysis dataset: response + predictors, NA-cleaned
# ----------------------------------------------------------------------------

analysis_df <- df[, c("DENS_DIFF", "CHIRPS_DJF", "CHIRPS_AM_", "DIFF_EVI")]
complete_idx <- complete.cases(analysis_df)
analysis_clean <- analysis_df[complete_idx, ]

cat("Total MAUs:", nrow(df), "  Complete cases:", nrow(analysis_clean), "\n\n")

# ============================================================================
# GLOBAL OLS MODEL
# ============================================================================

model_ols <- lm(DIFF_EVI ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_, data = analysis_clean)
s_ols <- summary(model_ols)
vifs <- vif(model_ols)

cat("=== GLOBAL OLS MODEL ===\n")
print(s_ols)
cat("\nVIF:\n")
print(vifs)

ols_coefs <- as.data.frame(s_ols$coefficients)
ols_coefs$term <- rownames(ols_coefs)
write.csv(ols_coefs, file.path(DIR_GWR, "ols_coefficients.csv"), row.names = FALSE)

# ============================================================================
# RESIDUAL SPATIAL AUTOCORRELATION CHECK (Moran's I) -- OLS
# ============================================================================

coords <- crds(centroids(mau_grid[complete_idx, ]))
dist_mat <- as.matrix(dist(coords))
dist_inv <- 1 / dist_mat
diag(dist_inv) <- 0
dist_inv[is.infinite(dist_inv)] <- 0

moran_ols <- Moran.I(residuals(model_ols), dist_inv)
cat("\n=== MORAN'S I ON OLS RESIDUALS ===\n")
cat("Observed:", moran_ols$observed, "  p-value:", moran_ols$p.value, "\n\n")

# ============================================================================
# GEOGRAPHICALLY WEIGHTED REGRESSION
# ============================================================================

sp_data <- SpatialPointsDataFrame(
  coords = coords,
  data = analysis_clean,
  proj4string = CRS(crs(mau_grid))
)

cat("Selecting optimal adaptive bandwidth (AICc)...\n")
bw_final <- bw.gwr(
  DIFF_EVI ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_,
  data = sp_data, approach = "AICc", kernel = "gaussian", adaptive = TRUE
)

gwr_model <- gwr.basic(
  DIFF_EVI ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_,
  data = sp_data, bw = bw_final, kernel = "gaussian", adaptive = TRUE
)

cat("\n=== GWR MODEL ===\n")
cat("Bandwidth (adaptive, n neighbors):", bw_final, "\n")
print(gwr_model)

sdf <- as.data.frame(gwr_model$SDF)
cat("\nSDF field names (verify SE column names before use):\n")
print(names(sdf))

t_dem  <- sdf$DENS_DIFF  / sdf$DENS_DIFF_SE
t_djfm <- sdf$CHIRPS_DJF / sdf$CHIRPS_DJF_SE
t_am   <- sdf$CHIRPS_AM_ / sdf$CHIRPS_AM__SE

n_mau <- nrow(sdf)

# ============================================================================
# GWR RESIDUAL DIAGNOSTICS (R3, Comment 7)
# ============================================================================

residuals_gwr <- sdf$residual

moran_gwr <- Moran.I(residuals_gwr, dist_inv, na.rm = TRUE)
cat("\n=== MORAN'S I ON GWR RESIDUALS ===\n")
cat("Observed:", moran_gwr$observed, "  p-value:", moran_gwr$p.value, "\n")
cat("Reduction relative to OLS residual Moran's I:",
    round(100 * (1 - moran_gwr$observed / moran_ols$observed), 1), "%\n\n")

# ============================================================================
# LOCAL COLLINEARITY DIAGNOSTIC (R3, Comment 8)
# ============================================================================

collin_diag <- gwr.collin.diagno(
  DIFF_EVI ~ DENS_DIFF + CHIRPS_DJF + CHIRPS_AM_,
  data = sp_data, bw = bw_final, kernel = "gaussian", adaptive = TRUE
)

cat("\n=== LOCAL VARIANCE INFLATION FACTORS ===\n")
print(summary(collin_diag$VIF))

vif_df <- as.data.frame(collin_diag$VIF)
names(vif_df) <- c("VIF_DENS_DIFF", "VIF_CHIRPS_DJF", "VIF_CHIRPS_AM")
write.csv(vif_df, file.path(DIR_GWR, "gwr_local_vif.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# Regime classification: significance-based, four descriptive categories.
# ----------------------------------------------------------------------------

dem_sig  <- abs(t_dem)  > GWR_T_THRESHOLD
clim_sig <- (abs(t_djfm) > GWR_T_THRESHOLD) | (abs(t_am) > GWR_T_THRESHOLD)

regime <- rep("N.s.", n_mau)
regime[dem_sig & !clim_sig] <- "Dem"
regime[!dem_sig & clim_sig] <- "Clim"
regime[dem_sig & clim_sig]  <- "Dem x Clim"

sign_dem  <- ifelse(sdf$DENS_DIFF  > 0, "positive", "negative")
sign_djfm <- ifelse(sdf$CHIRPS_DJF > 0, "positive", "negative")
sign_am   <- ifelse(sdf$CHIRPS_AM_ > 0, "positive", "negative")

# ----------------------------------------------------------------------------
# Export results to the MAU grid
# ----------------------------------------------------------------------------

mau_result <- mau_grid[complete_idx, ]
mau_result$beta_DEM  <- sdf$DENS_DIFF
mau_result$beta_DJFM <- sdf$CHIRPS_DJF
mau_result$beta_AM   <- sdf$CHIRPS_AM_
mau_result$local_R2  <- sdf$Local_R2
mau_result$t_DEM     <- t_dem
mau_result$t_DJFM    <- t_djfm
mau_result$t_AM      <- t_am
mau_result$REGIME    <- regime
mau_result$sign_DEM  <- sign_dem
mau_result$sign_DJFM <- sign_djfm
mau_result$sign_AM   <- sign_am
mau_result$resid_gwr <- residuals_gwr
mau_result$vif_DJFM  <- collin_diag$VIF[, 2]
mau_result$vif_AM    <- collin_diag$VIF[, 3]

out_file <- file.path(DIR_GWR, "GWR_final_results.shp")
writeVector(mau_result, out_file, overwrite = TRUE)
cat("\nShapefile written to:", out_file, "\n")

write.csv(as.data.frame(sdf), file.path(DIR_GWR, "gwr_local_coefficients.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------------

cat("\n=== REGIME DISTRIBUTION ===\n")
print(table(regime))

cat("\n=== SIGN BREAKDOWN AMONG SIGNIFICANT CELLS ===\n")
cat("DJFM significant: positive", sum(abs(t_djfm) > GWR_T_THRESHOLD & sign_djfm == "positive"),
    " negative", sum(abs(t_djfm) > GWR_T_THRESHOLD & sign_djfm == "negative"), "\n")
cat("AM significant:   positive", sum(abs(t_am) > GWR_T_THRESHOLD & sign_am == "positive"),
    " negative", sum(abs(t_am) > GWR_T_THRESHOLD & sign_am == "negative"), "\n")
cat("DEM significant:  positive", sum(abs(t_dem) > GWR_T_THRESHOLD & sign_dem == "positive"),
    " negative", sum(abs(t_dem) > GWR_T_THRESHOLD & sign_dem == "negative"), "\n")

cat("\n=== MODEL DIAGNOSTICS SUMMARY ===\n")
cat("OLS residual Moran's I:", moran_ols$observed, " (p =", moran_ols$p.value, ")\n")
cat("GWR residual Moran's I:", moran_gwr$observed, " (p =", moran_gwr$p.value, ")\n")
cat("Local VIF (DENS_DIFF):  min", min(collin_diag$VIF[,1]), " max", max(collin_diag$VIF[,1]), "\n")
cat("Local VIF (CHIRPS_DJF): min", min(collin_diag$VIF[,2]), " max", max(collin_diag$VIF[,2]), "\n")
cat("Local VIF (CHIRPS_AM_): min", min(collin_diag$VIF[,3]), " max", max(collin_diag$VIF[,3]), "\n")

cat("\nAll outputs saved to:", DIR_GWR, "\n")
cat("09_analysis.R complete. Pipeline finished.\n")