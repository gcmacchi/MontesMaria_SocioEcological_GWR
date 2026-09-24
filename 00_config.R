# ============================================================================
# 00_config.R
# Central configuration for the Montes de Maria reproducibility pipeline.
# Every other script in this repository sources this file first and reads
# its parameters from here -- no path, threshold, or date range is
# hard-coded a second time anywhere else.
#
# USAGE: set the working directory to the repository root (the folder
# containing this file and the "data" subfolder) before sourcing:
#   setwd("<repository root>")
#   source("00_config.R")
# ============================================================================

# ----------------------------------------------------------------------------
# 1. DIRECTORY STRUCTURE (all relative to the repository root)
# ----------------------------------------------------------------------------

DIR_ROOT       <- "."
DIR_DATA       <- file.path(DIR_ROOT, "data")
DIR_STUDY_AREA <- file.path(DIR_DATA, "study_area")
DIR_DEMOGRAPHY <- file.path(DIR_DATA, "demography")

STUDY_AREA_SHP <- file.path(DIR_STUDY_AREA, "study_area.shp")
DEMOGRAPHY_SHP <- file.path(DIR_DEMOGRAPHY, "demography.shp")

DIR_HARMONIZATION <- file.path(DIR_ROOT, "01_harmonization")
DIR_EVI            <- file.path(DIR_ROOT, "02_evi")
DIR_CHIRPS          <- file.path(DIR_ROOT, "03_chirps")
DIR_BINDER          <- file.path(DIR_ROOT, "04_binder")
DIR_GWR             <- file.path(DIR_ROOT, "05_gwr")

# Harmonization subfolders (mirrors the internal structure already used by
# 01_harmonization.R)
DIR_HARM_RAW    <- file.path(DIR_HARMONIZATION, "01_raw")
DIR_HARM_CROP   <- file.path(DIR_HARMONIZATION, "02_crop")
DIR_HARM_EVI    <- file.path(DIR_HARMONIZATION, "03_evi")
DIR_HARM_SYNTH  <- file.path(DIR_HARMONIZATION, "04_synthesis")

# EVI subfolders
DIR_EVI_RAW <- file.path(DIR_EVI, "01_raw")
DIR_EVI_OUT <- file.path(DIR_EVI, "02_composites")

# CHIRPS subfolders
DIR_CHIRPS_RAW      <- file.path(DIR_CHIRPS, "01_raw_gz")
DIR_CHIRPS_MONTHLY  <- file.path(DIR_CHIRPS, "02_monthly_cropped")
DIR_CHIRPS_SEASONAL <- file.path(DIR_CHIRPS, "03_seasonal")

# Create all output directories up front (input directories are checked,
# not created -- their absence is a fatal error, handled in Section 5).
for (d in c(DIR_HARM_RAW, DIR_HARM_CROP, DIR_HARM_EVI, DIR_HARM_SYNTH,
            DIR_EVI_RAW, DIR_EVI_OUT,
            DIR_CHIRPS_RAW, DIR_CHIRPS_MONTHLY, DIR_CHIRPS_SEASONAL,
            DIR_BINDER, DIR_GWR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------------------------------------------------------
# 2. STUDY AREA / DOWNLOAD BOUNDING BOX
# ----------------------------------------------------------------------------

# Generous WGS84 bounding box used for the initial STAC/CHIRPS download crop
# (~30 km buffer around the strict study area extent). The exact study-area
# boundary crop is applied afterward, per-dataset, using STUDY_AREA_SHP.
BBOX_WGS84 <- c(xmin = -75.75, ymin = 9.05, xmax = -74.45, ymax = 10.75)

# WRS-2 Path/Row for all Landsat scenes used in this study
LANDSAT_PATH <- 9L
LANDSAT_ROW  <- 53L

# ----------------------------------------------------------------------------
# 3. LANDSAT / EVI PARAMETERS
# ----------------------------------------------------------------------------

# Collection 2 Level-2 Surface Reflectance scale and offset
SR_SCALE  <- 0.0000275
SR_OFFSET <- -0.2

# Sensor-specific QA_PIXEL "clear land" codes
QA_GOOD_VALUE_ETM  <- 5440    # Landsat 5/7 (TM/ETM+)
QA_GOOD_VALUE_OLI  <- 21824   # Landsat 8/9 (OLI)

# Harmonization pixel-eligibility thresholds (used only in 01_harmonization.R)
L7_ATMOS_OPACITY_MAX   <- 250
L8_AEROSOL_GOOD_VALUES <- c(64, 66, 96)

# Local ETM+ -> OLI cross-sensor calibration (pixel-weighted pooled fit,
# 5 near-coincident scene pairs, n = 8,988,291 pixels, R2 = 0.9145).
# Derived by 01_harmonization.R; hard-coded here for use by downstream
# scripts (02_evi_epoch1.R) once derived and confirmed.
LOCAL_HARMONIZATION_INTERCEPT <- 0.020907
LOCAL_HARMONIZATION_SLOPE     <- 0.970029

# Reference external calibration (Roy et al. 2016, NIR band, CONUS),
# retained only for comparison in the harmonization synthesis report.
ROY_NIR_SLOPE     <- 0.8462
ROY_NIR_INTERCEPT <- 0.0412

# Epoch scene windows (Path/Row fixed above)
EPOCH1_YEARS <- 2008:2010   # baseline: Landsat 5/7
EPOCH2_YEARS <- 2021:2022   # recent: Landsat 8/9
CHECKPOINT_YEARS <- 2014:2015  # intermediate, trajectory visualization only

# ----------------------------------------------------------------------------
# 4. CHIRPS PARAMETERS
# ----------------------------------------------------------------------------

CHIRPS_YEARS  <- 1999:2024
CHIRPS_MONTHS <- 1:12

CHIRPS_BASELINE_YEARS <- 2003:2007   # centered on the 2005 census
CHIRPS_RECENT_YEARS   <- 2016:2020   # centered on the 2018 census

CHIRPS_URL_CAMER  <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0/camer-carib_monthly/tifs/"
CHIRPS_URL_GLOBAL <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/tifs/"

CHIRPS_NODATA <- -9999   # ocean / non-land sentinel value

# Crop criterion: a pixel is retained if its centroid falls within this
# buffer distance (meters) of the study area boundary.
CHIRPS_BUFFER_M <- 5000

# Physically plausible monthly precipitation bounds (mm); values outside
# this range are treated as data artefacts, not valid extreme rainfall.
CHIRPS_LOWER_BOUND <- 0
CHIRPS_UPPER_BOUND <- 2000

# Area-proportional extraction: disaggregation factor from native CHIRPS
# resolution (0.05 deg, ~5000 m) to the fine grid (50 m) used for zonal
# averaging onto the MAU grid.
CHIRPS_DISAGG_FACTOR <- 100

# ----------------------------------------------------------------------------
# 5. GWR / REGIME CLASSIFICATION PARAMETERS
# ----------------------------------------------------------------------------

GWR_T_THRESHOLD <- 1.96   # |t| > threshold defines local significance

# ----------------------------------------------------------------------------
# 6. INPUT FILE CHECKS (fail fast, with a clear message, if missing)
# ----------------------------------------------------------------------------

if (!file.exists(STUDY_AREA_SHP)) {
  stop("Required input not found: ", STUDY_AREA_SHP,
       "\nExpected location: ", DIR_STUDY_AREA, "/study_area.shp (+ .dbf/.shx/.prj)")
}

if (!file.exists(DEMOGRAPHY_SHP)) {
  stop("Required input not found: ", DEMOGRAPHY_SHP,
       "\nExpected location: ", DIR_DEMOGRAPHY, "/demography.shp (+ .dbf/.shx/.prj)")
}

cat("=== 00_config.R loaded ===\n")
cat("Repository root:", normalizePath(DIR_ROOT), "\n")
cat("Study area:", STUDY_AREA_SHP, "\n")
cat("Demography (MAU grid):", DEMOGRAPHY_SHP, "\n")
cat("Configuration OK.\n\n")