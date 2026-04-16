# Methods

## Analytical Framework

All three analyses follow the methodology of Rutterford et al. (2023) {cite}`rutterford2023sea`:

### 1. Data Preparation

- Fish abundance data (CPUE or counts) downloaded from source databases
- Environmental layers (SST, salinity) extracted from Bio-ORACLE {cite}`assis2018biooracle` at each sampling location using nearest-neighbour interpolation
- Hauls aggregated to 1×1° grid cells, computing mean CPUE per species and mean environmental values per cell
- Grid cells with fewer than 3 hauls excluded

### 2. Community Analysis (PCoA)

- Abundance data fourth-root transformed to reduce influence of dominant species
- Bray-Curtis dissimilarity matrix computed between all grid cells
- Principal Coordinates Analysis (PCoA) via classical multidimensional scaling
- Three axes retained, with variance explained reported for each

### 3. Regression Analysis

- Pearson correlation (univariate R²) between each PCoA axis and each environmental variable (SST, salinity, log₁₀ depth)
- Multiple linear regression (full model R²) using standardised environmental predictors
- The environmental variable with the highest univariate |r| on Axis 1 is identified as the "primary driver"

### 4. Reproduction Criterion

The claim is **confirmed** if SST has the highest univariate |r| on PCoA Axis 1 (p < 0.05).

## Environmental Data

All pipelines use Bio-ORACLE v2.0 {cite}`assis2018biooracle` surface layers:

| Variable | Bio-ORACLE dataset | Resolution |
|---|---|---|
| SST | `thetao_baseline_2000_2019_depthsurf` | 0.05° (~5.5 km) |
| Salinity | `so_baseline_2000_2019_depthsurf` | 0.05° (~5.5 km) |

Depth is taken from the survey data itself (haul depth), not from Bio-ORACLE.

## Software

| Package | Purpose |
|---|---|
| pandas | Data manipulation |
| numpy | Numerical computation |
| scipy | Bray-Curtis dissimilarity, Pearson correlation |
| scikit-learn | Linear regression, standardisation |
| xarray | NetCDF environmental data |
| matplotlib | Visualisation |
| requests | API access (DATRAS, SEANOE, Bio-ORACLE, WoRMS) |
