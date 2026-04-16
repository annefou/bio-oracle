# Replication: ClimateFish (Mediterranean)

This replication tests the SST claim in a **different region** (Mediterranean Sea) using a **different data collection method** (visual census) with the ClimateFish database {cite}`azzurro2022climatefish`.

## Nanopublications

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

- **Source**: [ClimateFish](https://www.seanoe.org/data/00756/86784/) (CC-BY 4.0)
- **Method**: Visual census transects (SCUBA)
- **Species**: 15 climate indicator species (7 indigenous + 8 non-indigenous)
- **Sites**: 28 (aggregated from 2,222 transects)
- **Period**: 2009-2021
- **Countries**: 7 Mediterranean countries

## Results

### PCoA Axis 1: SST is the primary driver

```
                      Rutterford et al.    This replication
Region                NE Atlantic          Mediterranean
Axis 1 variance       29%                  45.2%
SST R²                0.890                0.615
SST r                 —                    -0.784
Grid cells            193                  28
Species               198                  15
```

:::{figure} results/pcoa_community.png
:name: climatefish-pcoa
PCoA ordination of Mediterranean fish communities (ClimateFish visual census data).
:::

### Variable Importance

:::{figure} results/variable_importance.png
:name: climatefish-variable-importance
Comparison of variable importance: NE Atlantic (Rutterford) vs Mediterranean (this study).
:::

## Regression Table

| Axis | Variable | r | R² | p-value |
|------|----------|---|----|----|
| PCoA1 | **SST** | **-0.784** | **0.615** | **<0.001** |
| PCoA1 | Salinity | -0.472 | 0.222 | 0.011 |
| PCoA1 | Depth | +0.087 | 0.008 | 0.659 |

## Interpretation

Despite using a different region, different species set (15 vs 198), and different sampling method (visual census vs trawl), SST emerges as the primary driver of community structure on Axis 1. The higher variance explained on Axis 1 (45.2% vs 29%) may reflect the smaller, more ecologically coherent species set.

## Running

```bash
python replication_climatefish/01_download_climatefish.py
python replication_climatefish/02_download_biooracle.py
python replication_climatefish/03_spatial_join.py
python replication_climatefish/04_community_analysis.py
python replication_climatefish/05_variable_importance.py
```
