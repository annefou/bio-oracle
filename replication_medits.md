# Replication: MEDITS (Mediterranean Trawl)

This replication tests the SST claim in the **Mediterranean Sea** using **trawl survey data** from the MEDITS programme, restricted to continental shelf depths (≤200m).

## Nanopublications

The replication study and outcome nanopubs cover both Mediterranean datasets (ClimateFish + MEDITS):

::::{grid} 1 1 2 2
:gutter: 3

:::{card} Replication Study (Mediterranean)
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAVvRoVuJqBuWCK1XiDhVP8gS3QjO3k1k5CjPJ6Rdf_qQ
ClimateFish visual census + MEDITS trawl, different region and method
:::

:::{card} Replication Outcome — Confirmed
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAjEJL1PpNdE7lBYjHkOQidLG9iYmyRtW2Of5PCTAotMg
SST confirmed in Mediterranean (ClimateFish R²=0.615, MEDITS R²=0.367)
:::

::::

## Data

- **Source**: [JRC EU DCF MEDITS](https://data.jrc.ec.europa.eu/dataset/ef36af5d-eb4e-4a9f-9513-82a877fe71fe)
- **Method**: Standardised demersal trawl surveys
- **Species**: 230 (after filtering for haul frequency)
- **Grid cells**: 115 (1x1° aggregation, shelf only)
- **Period**: 2005-2024
- **Depth filter**: ≤ 200m (continental shelf)

## Why the Depth Filter Matters

Without the shelf depth filter, MEDITS hauls span 10-800m. At that range, **depth dominates** PCoA Axis 1 (R²=0.870) and SST is secondary (R²=0.165). This makes ecological sense — the transition from shelf to slope to bathyal communities is the strongest gradient in demersal fish assemblages.

By restricting to shelf depths (≤200m), we match Rutterford et al.'s scope and allow the SST signal to emerge:

| | Full depth range | Shelf only (≤200m) |
|---|---|---|
| Axis 1 primary | **Depth** (R²=0.870) | **SST** (R²=0.367) |
| SST on Axis 1 | R²=0.165 | R²=0.367 |
| Grid cells | 130 | 115 |

## Results

### PCoA Axis 1: SST is the primary driver (shelf)

```
                      Rutterford et al.    This replication
Region                NE Atlantic          Mediterranean
Method                Trawl                Trawl
Axis 1 variance       29%                  18.9%
SST R²                0.890                0.367
SST r                 —                    -0.605
Grid cells            193                  115
Species               198                  230
```

:::{figure} results/medits_pcoa.png
:name: medits-pcoa
PCoA ordination of Mediterranean shelf fish communities (MEDITS trawl data, depth ≤200m).
:::

### Variable Importance

:::{figure} results/medits_variable_importance.png
:name: medits-variable-importance
Comparison of variable importance: NE Atlantic (Rutterford) vs Mediterranean shelf (MEDITS).
:::

### Spatial Distribution

:::{figure} results/medits_pcoa_map.png
:name: medits-map
Spatial distribution of PCoA Axis 1 scores across the Mediterranean (shelf stations).
:::

## Regression Table

| Axis | Variable | r | R² | p-value |
|------|----------|---|----|----|
| PCoA1 | **SST** | **-0.605** | **0.367** | **<0.001** |
| PCoA1 | Salinity | -0.433 | 0.187 | <0.001 |
| PCoA1 | log(Depth) | -0.517 | 0.268 | <0.001 |
| PCoA2 | log(Depth) | -0.624 | 0.389 | <0.001 |

## Interpretation

The weaker SST effect (R²=0.367 vs 0.890) likely reflects the narrower SST range in the Mediterranean (5.6°C) compared to the NE Atlantic (10°C). With less temperature variation, the signal is weaker but still the dominant driver on Axis 1.

## Running

```bash
python replication_medits/01_download_medits.py
python replication_medits/02_spatial_join.py      # Includes shelf filter
python replication_medits/03_community_analysis.py
```
