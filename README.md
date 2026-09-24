# Montes de María — Reproducibility Repository

This repository contains the full processing pipeline supporting the demographic,
climatic, and vegetation-trajectory analysis of Montes de María (Colombian
Caribbean) reported in the manuscript. Every intermediate output (rasters, CSV
tables, and the final MAU-level dataset) can be regenerated from these scripts
and the two input files listed below.

## 1. Requirements

**R version:** 4.3 or later recommended.

**Required packages:**

| Package | Used for |
|---|---|
| `terra` | Raster/vector I/O, cropping, masking, zonal statistics |
| `sf` | Coordinate reprojection helpers |
| `rstac` | STAC search and asset retrieval (Microsoft Planetary Computer) |
| `httr` | CHIRPS file download |
| `R.utils` | Decompression of downloaded CHIRPS `.tif.gz` files |
| `sp` | Spatial data structures required by `GWmodel` |
| `GWmodel` | Geographically Weighted Regression (bandwidth selection and fit) |
| `car` | Variance Inflation Factor (multicollinearity check) |
| `ape` | Moran's I test on OLS residuals |

Install with:

```r
install.packages(c("terra", "sf", "rstac", "httr", "R.utils",
                    "sp", "GWmodel", "car", "ape"))
```

**Internet access** is required: scripts 02–05 download Landsat Collection 2
Level-2 imagery from the Microsoft Planetary Computer STAC API and CHIRPS v2.0
precipitation data from the Climate Hazards Center (UC Santa Barbara).

**Note on `R.utils` / `terra::extract()` conflict:** `R.utils` masks
`terra::extract()` when both packages are loaded. Scripts that need both call
`detach("package:R.utils", unload = TRUE)` before extraction, or qualify the
call as `terra::extract()`. This is already handled inside the scripts; it is
noted here only so the behavior is not mistaken for an error.

## 2. Required input data

Two files must be placed under `./data/` before running anything:

```
data/
  study_area/
    study_area.shp   (+ .dbf, .shx, .prj)
  demography/
    demography.shp   (+ .dbf, .shx, .prj)
```

- `study_area.shp` — polygon boundary of the Montes de María study area.
- `demography.shp` — the Minimum Analytical Unit (MAU) grid, with
  pre-computed dasymetric population fields (`DENS_2005`, `DENS_2018`,
  `DENS_DIFF`) already attached to each MAU polygon. Construction of this
  file (dasymetric population reconstruction) is outside the scope of this
  repository.

All other data — Landsat scenes and CHIRPS rasters — are downloaded
automatically by the scripts below.

## 3. Pipeline overview and execution order

Scripts must be run in the order listed. Each writes its outputs to a
dedicated subfolder; later scripts read from the outputs of earlier ones.

| # | Script | Purpose | Key output |
|---|---|---|---|
| 1 | `01_harmonization.R` | Derives the local ETM+ → OLI cross-sensor calibration from 5 near-coincident Landsat 7/8 scene pairs | Pooled calibration coefficients (intercept, slope), synthesis report |
| 2 | `02_evi_epoch1.R` | Downloads and processes the Epoch 1 (baseline) Landsat 5/7 scenes; applies the calibration from step 1 | `EPOCH1_EVI_MEAN_harmonized.tif`; `epoch1_pixel_log.csv` (defines the pixel-count target used in step 3) |
| 3 | `03_evi_epoch2.R` | Downloads and processes the Epoch 2 (recent) Landsat 8/9 scenes (no calibration needed — OLI is the native reference) | `EPOCH2_EVI_MEAN.tif` |
| 4 | `04_evi_checkpoint.R` *(optional)* | Builds an intermediate (2014–2015) EVI composite for trajectory visualization only; not used in the regression model | `CHECKPOINT_EVI_MEAN.tif` |
| 5 | `05_chirps_download.R` | Downloads the full monthly CHIRPS v2.0 record (1999–2024) and crops it to the study area | Monthly cropped rasters; `chirps_monthly_summary_table.csv` |
| 6 | `06_chirps_seasonal.R` | Builds the DJFM (dry-season) and AM seasonal composites for the baseline and recent periods, and their percent-anomaly difference | `DJFM_diff_pct.tif`, `AM_diff_pct.tif` |
| 7 | `07_chirps_to_mau.R` | Extracts the CHIRPS anomalies onto the MAU grid via area-proportional (fine-grid disaggregation) sampling | `mau_grid_with_chirps.shp` |
| 8 | `08_evi_to_mau.R` | Computes the EVI difference (Recent − Baseline) and extracts its zonal mean onto the MAU grid | `mau_grid_full.shp` — the complete analysis dataset |
| 9 | `09_analysis.R` | Fits the global OLS model, tests residual spatial autocorrelation (Moran's I), fits the Geographically Weighted Regression, and classifies each MAU into a descriptive regime | `GWR_final_results.shp`; `gwr_local_coefficients.csv`; `ols_coefficients.csv` |

Script 4 is independent of scripts 5–9 and can be skipped if only the
regression results (not the descriptive trajectory figure) are needed.

## 4. Key methodological parameters

These values are fixed in each script as named constants, not hard-coded
inline, so they can be located and audited directly in the source:

- **Cross-sensor calibration** (from step 1): `EVI_OLI = 0.020907 + 0.970029 × EVI_ETM+` (pooled, pixel-weighted fit; R² = 0.9145; n = 8,988,291 pixels from 5 scene pairs).
- **Cloud/shadow/fill masks:** `QA_PIXEL == 5440` (Landsat 5/7), `QA_PIXEL == 21824` (Landsat 8/9).
- **Surface reflectance scaling:** DN × 0.0000275 − 0.2 (Collection 2 Level-2 standard).
- **EVI formula:** 2.5 × (NIR − Red) / (NIR + 6×Red − 7.5×Blue + 1), clamped to [−1, 1].
- **Epoch 1 (baseline) window:** Landsat 5/7 scenes, January–March 2008–2010 (Path 9, Row 53).
- **Epoch 2 (recent) window:** Landsat 8/9 scenes, January–February 2021–2022 (Path 9, Row 53).
- **CHIRPS crop criterion:** a pixel is retained if its centroid falls within a 5 km buffer of the study area boundary (see comments in `05_chirps_download.R` for why this replaced earlier, less reproducible crop methods).
- **CHIRPS seasonal windows:** DJFM = December (year−1) + January–March (year); AM = April–May (same year, no lag). Baseline = 2003–2007; Recent = 2016–2020 (both centered on the respective census year).
- **Area-proportional CHIRPS extraction:** native 0.05° CHIRPS rasters are disaggregated to a 50 m grid (nearest-neighbor, no interpolation) before zonal averaging onto each MAU polygon.
- **GWR regime classification:** each MAU is classified by the statistical significance of its local coefficients (|t| > 1.96) into `Dem` (demographic-dominant), `Clim` (climate-dominant), `Dem x Clim` (co-dominant), or `N.s.` (not significant). The sign (direction) of each significant coefficient is recorded separately and does not affect this classification.

## 5. Known limitations, documented for transparency

- **Epoch depth asymmetry:** Landsat 7 (Epoch 1) carries a permanent Scan
  Line Corrector failure (since 31 May 2003), which removes roughly 22% of
  pixels per scene. Epoch 1 therefore uses 11 scenes over a 3-year window,
  while Epoch 2 (unaffected sensor) uses 11 scenes over a shorter window
  needed to reach a comparable cumulative valid-pixel count. This is a
  deliberate pixel-count-matching design, not an inconsistency.
- **Calibration window vs. analysis windows:** the 5 scene pairs used to
  derive the local cross-sensor calibration (step 1) fall in 2014–2016,
  outside both the Epoch 1 (2008–2010) and Epoch 2 (2021–2022) windows.
  The calibration is assumed stable across this gap; no calibration pairs
  closer to 2008–2010 were available in the archive.
- **Shapefile field-name truncation:** the `.shp` format truncates field
  names to 10 characters. `CHIRPS_DJFM_pct` and `CHIRPS_AM_pct` are stored
  as `CHIRPS_DJF` and `CHIRPS_AM_` in all output shapefiles. Field names
  are printed to the console at the point of creation in each script for
  verification.

## 6. Outputs relevant to the manuscript

The final analysis dataset (`mau_grid_full.shp`, step 8) and the GWR results
(`GWR_final_results.shp`, step 9) are the two files that reproduce, respectively,
Table [X] and Figure [X] of the manuscript. The regime distribution reported
in the manuscript is reproduced by `table(mau_result$REGIME)` at the end of
`09_analysis.R`.
