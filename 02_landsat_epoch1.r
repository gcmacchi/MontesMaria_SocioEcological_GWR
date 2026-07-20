# 02_landsat_epoch1.R
# Landsat 5/7 EVI Composite (Epoch 1: 1999–2003)
# Author: Giancarlo Macchi Janica
# Date: 2025-11-23

library(terra)
library(stringr)
library(fs)

# --- CONFIGURATION ---
ROOT_DIR <- "data/raw/landsat/epoch_1"
OUTPUT_DIR <- "data/processed"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# --- FUNCTIONS ---
create_cloud_mask <- function(qa_raster) {
  is_cloud <- (qa_raster %/% 8) %% 2
  is_shadow <- (qa_raster %/% 16) %% 2
  is_dilated <- (qa_raster %/% 2) %% 2
  mask <- ifel((is_cloud == 1 | is_shadow == 1 | is_dilated == 1), NA, 1)
  return(mask)
}

calc_evi <- function(blue, red, nir) {
  scale_factor <- 0.0000275
  offset <- -0.2
  b_blue <- (blue * scale_factor) + offset
  b_red <- (red * scale_factor) + offset
  b_nir <- (nir * scale_factor) + offset
  evi <- 2.5 * ((b_nir - b_red) / (b_nir + 6 * b_red - 7.5 * b_blue + 1))
  evi <- clamp(evi, -1, 1)
  return(evi)
}

# --- STEP 1: Common Extent ---
b1_files <- list.files(ROOT_DIR, pattern = "SR_B1.TIF$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if(length(b1_files) == 0) stop("No files found.")

common_ext <- ext(rast(b1_files[1]))
for(i in 2:length(b1_files)) {
  common_ext <- intersect(common_ext, ext(rast(b1_files[i])))
}
ref_rast <- rast(b1_files[1])
MASTER_GRID <- crop(ref_rast, common_ext)

# --- STEP 2: Process Scenes ---
qa_files <- list.files(ROOT_DIR, pattern = "QA_PIXEL.TIF$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
evi_stack <- list()

for (qa_path in qa_files) {
  b1_path <- str_replace(qa_path, regex("QA_PIXEL.TIF", ignore_case = TRUE), "SR_B1.TIF")
  b3_path <- str_replace(qa_path, regex("QA_PIXEL.TIF", ignore_case = TRUE), "SR_B3.TIF")
  b4_path <- str_replace(qa_path, regex("QA_PIXEL.TIF", ignore_case = TRUE), "SR_B4.TIF")

  if(!file.exists(b1_path)) next

  r_blue <- resample(crop(rast(b1_path), common_ext), MASTER_GRID, method = "bilinear")
  r_red  <- resample(crop(rast(b3_path), common_ext), MASTER_GRID, method = "bilinear")
  r_nir  <- resample(crop(rast(b4_path), common_ext), MASTER_GRID, method = "bilinear")
  r_qa   <- resample(crop(rast(qa_path), common_ext), MASTER_GRID, method = "near")

  cloud_mask <- create_cloud_mask(r_qa)
  r_blue <- mask(r_blue, cloud_mask)
  r_red  <- mask(r_red, cloud_mask)
  r_nir  <- mask(r_nir, cloud_mask)

  evi_scene <- calc_evi(r_blue, r_red, r_nir)
  evi_stack[[basename(qa_path)]] <- evi_scene
}

# --- STEP 3: Temporal Mean ---
stack_collection <- rast(evi_stack)
epoch_mean_evi <- mean(stack_collection, na.rm = TRUE)

writeRaster(epoch_mean_evi, file.path(OUTPUT_DIR, "EPOCH1_EVI_MEAN.tif"), overwrite = TRUE)