# Reproduction: ICES DATRAS (NE Atlantic)

The reproduction uses the **same data source** as Rutterford et al. (2023) — ICES DATRAS trawl surveys across the NE Atlantic — to verify the original results.

## Nanopublications

::::{grid} 1 1 2 2
:gutter: 3

:::{card} Reproduction Study
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAzP6xzTxbXWC9hJDkdp2M5gq4di4cVggPvxRPnBh9-1k
Same data source (ICES DATRAS), same method as Rutterford et al.
:::

:::{card} Reproduction Outcome — Confirmed
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RA7f-x8JBwYZCid2GYcvFzkrGQWyJkVyDXa1pvMH2FgwM
SST confirmed as primary driver (R²=0.573, Axis 1: 29.2% variance)
:::

::::

## Data

- **Source**: [ICES DATRAS](https://datras.ices.dk/) (Database of Trawl Surveys)
- **Surveys**: 22 NE Atlantic bottom trawl surveys (2005-2018)
- **Hauls**: 38,012 (after filtering)
- **Grid cells**: 247 (1x1° aggregation, ≥3 hauls per cell)
- **Species**: 258 fish species (Actinopterygii + Elasmobranchii)

## Critical Filters

To match Rutterford et al.'s study design, three filters were applied:

1. **Exclude Baltic surveys** (BITS, SE-SOUND) — brackish salinity (~7 PSU) creates a gradient that overwhelms the SST signal
2. **Shelf depth ≤ 200m** — continental shelf only, matching Rutterford's scope
3. **Fish species only** — classified via [WoRMS](https://www.marinespecies.org/) API (Actinopterygii + Elasmobranchii), excluding invertebrates

## Results

### PCoA Axis 1: SST is the primary driver

```
                      Rutterford et al.    This reproduction
Axis 1 variance       29%                  29.2%
SST R²                0.890                0.573
Full model R²         —                    0.879
Grid cells            193                  247
Species               198                  258
```

:::{figure} results/datras_pcoa.png
:name: datras-pcoa
PCoA ordination of NE Atlantic fish communities coloured by SST, salinity, and depth.
:::

### Variable Importance

:::{figure} results/datras_variable_importance.png
:name: datras-variable-importance
Comparison of variable importance between Rutterford et al. (2023) and this reproduction.
:::

### Spatial Distribution

:::{figure} results/datras_pcoa_map.png
:name: datras-map
Spatial distribution of PCoA Axis 1 scores across the NE Atlantic.
:::

## Regression Table

| Axis | Variable | r | R² | p-value |
|------|----------|---|----|----|
| PCoA1 | **SST** | **-0.757** | **0.573** | **<0.001** |
| PCoA1 | Salinity | -0.399 | 0.159 | <0.001 |
| PCoA1 | log(Depth) | -0.642 | 0.412 | <0.001 |
| PCoA2 | SST | +0.645 | 0.416 | <0.001 |
| PCoA2 | **log(Depth)** | **+0.679** | **0.460** | **<0.001** |
| PCoA3 | Salinity | -0.223 | 0.050 | <0.001 |

## Running

```bash
python reproduction/01_download_datras.py          # ~3 hours, cached
python reproduction/02_download_biooracle.py
python reproduction/03_filter_and_join.py
python reproduction/04_community_analysis.py
```
