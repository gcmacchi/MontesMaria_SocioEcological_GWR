# 01_climate_analysis.R
# ERA5 Precipitation Anomalies – Montes de María
# Author: Giancarlo Macchi Janica
# Date: 2025-11-23

library(dplyr)
library(tidyr)
library(stringr)

# --- CONFIGURATION ---
# Assumes RAW_DATA_FINAL exists in the global environment
# (loaded from data/raw/RAW_DATA_FINAL.csv)

# --- STEP 1: Seasonal Aggregation ---
seasonal_totals <- RAW_DATA_FINAL %>%
  pivot_longer(cols = matches("^[A-Z]{3}_[0-9]{4}$"),
               names_to = "MonthYear", values_to = "Value") %>%
  separate(MonthYear, into = c("Month", "Year"), sep = "_") %>%
  mutate(Year = as.numeric(Year),
         Season = case_when(
           Month %in% c("MAY", "JUN", "JUL", "AUG") ~ "MGLA",
           Month %in% c("DEC", "JAN", "FEB", "MAR") ~ "DGFM",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(Season)) %>%
  mutate(SeasonYear = ifelse(Season == "DGFM" & Month == "DEC", Year + 1, Year)) %>%
  group_by(id, COMUN_, ROW_, Season, SeasonYear) %>%
  summarise(Seasonal_Total = sum(Value, na.rm = TRUE),
            n_months = n(), .groups = "drop") %>%
  filter(n_months == 4) %>%
  mutate(SeasonName = paste0(Season, "_", SeasonYear)) %>%
  pivot_wider(id_cols = c(id, COMUN_, ROW_),
              names_from = SeasonName,
              values_from = Seasonal_Total)

# --- STEP 2: 5-Year Moving Averages ---
period_2003 <- 2001:2005
period_2020 <- 2018:2022

analysis_data <- seasonal_totals %>%
  mutate(
    AVG_DGFM_2003 = rowMeans(dplyr::select(., any_of(paste0("DGFM_", period_2003))), na.rm = TRUE),
    AVG_MGLA_2003 = rowMeans(dplyr::select(., any_of(paste0("MGLA_", period_2003))), na.rm = TRUE),
    AVG_DGFM_2020 = rowMeans(dplyr::select(., any_of(paste0("DGFM_", period_2020))), na.rm = TRUE),
    AVG_MGLA_2020 = rowMeans(dplyr::select(., any_of(paste0("MGLA_", period_2020))), na.rm = TRUE)
  )

# --- STEP 3: Final Indices ---
final_results <- analysis_data %>%
  mutate(
    DIFF_DGFM_mm = AVG_DGFM_2020 - AVG_DGFM_2003,
    DIFF_MGLA_mm = AVG_MGLA_2020 - AVG_MGLA_2003,
    DIFF_DGFM_perc = ifelse(AVG_DGFM_2003 == 0, NA, (DIFF_DGFM_mm / AVG_DGFM_2003) * 100),
    DIFF_MGLA_perc = ifelse(AVG_MGLA_2003 == 0, NA, (DIFF_MGLA_mm / AVG_MGLA_2003) * 100)
  ) %>%
  dplyr::select(id, COMUN_, ROW_,
                AVG_DGFM_2003, AVG_DGFM_2020, DIFF_DGFM_mm, DIFF_DGFM_perc,
                AVG_MGLA_2003, AVG_MGLA_2020, DIFF_MGLA_mm, DIFF_MGLA_perc)

# --- EXPORT ---
write.csv(final_results, "data/processed/final_results.csv", row.names = FALSE)