# SST as Primary Driver of Fish Community Structure

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19628179.svg)](https://doi.org/10.5281/zenodo.19628179)

Rutterford et al. (2023) {cite}`rutterford2023sea` demonstrated that **sea surface temperature (SST) is the primary driver of fish community structure** across the Northeast Atlantic continental shelf. This project tests that claim at three levels:

| Level | Dataset | Region | Method | Result |
|---|---|---|---|---|
| **Reproduction** | ICES DATRAS | NE Atlantic | Trawl surveys | **Confirmed** |
| **Replication** | ClimateFish | Mediterranean | Visual census | **Confirmed** |
| **Replication** | MEDITS | Mediterranean | Trawl surveys (shelf) | **Confirmed** |

## Original Claim

> "In conclusion, our study clearly identifies sea temperature as the primary driver of fish community composition across the Northeast Atlantic continental shelf."
>
> — Rutterford et al. (2023), *Global Change Biology*, 29(9), 2510-2521. [DOI: 10.1111/gcb.16633](https://doi.org/10.1111/gcb.16633)

## Approach

All three pipelines follow the same analytical framework:

1. Download fish abundance data (CPUE or counts)
2. Download Bio-ORACLE {cite}`assis2018biooracle` environmental layers (SST, salinity)
3. Spatial join and aggregation to 1x1° grid cells
4. PCoA on Bray-Curtis dissimilarity of fourth-root transformed abundances
5. Regression of PCoA axes against SST, salinity, and log(depth)

## Results at a Glance

:::{figure} results/summary_comparison.png
:name: summary-comparison
SST is the primary driver of PCoA Axis 1 across all three analyses. Left: SST R² on Axis 1. Centre: variance explained by Axis 1. Right: variable importance comparison.
:::

:::{figure} results/summary_maps.png
:name: summary-maps
Spatial distribution of PCoA Axis 1 scores across the NE Atlantic (DATRAS) and Mediterranean (ClimateFish, MEDITS). The colour gradient tracks SST in all three analyses.
:::

## FORRT Replication Chain on Science Live

The results are published as a chain of nanopublications — cryptographically signed, machine-readable, and citable semantic assertions:

::::{grid} 1 1 2 2
:gutter: 3

:::{card} Paper Quotation
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAzy03sZgyjiQONhr5-f0bgw_kPqdf1nddc9xbOGa-OxM
Rutterford et al. (2023) quotation on SST as primary driver
:::

:::{card} AIDA Sentence
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAtKboaVur5jc4_i0iSi4sk70v6ocJLAgdmN6MDxcnJtU
Sea surface temperature is the primary driver of fish community structure across Northeast Atlantic shelf seas.
:::

:::{card} FORRT Claim
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/np/RA-DG_PBY8ddwmtYQzkFuejXePhrmJgUZyP65mupugNNg
The testable claim registered in the FORRT replication framework
:::

:::{card} Reproduction Study (DATRAS)
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAzP6xzTxbXWC9hJDkdp2M5gq4di4cVggPvxRPnBh9-1k
Same data source (ICES DATRAS), same method — **Confirmed**
:::

:::{card} Reproduction Outcome
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RA7f-x8JBwYZCid2GYcvFzkrGQWyJkVyDXa1pvMH2FgwM
SST confirmed as primary driver (R²=0.573, Axis 1: 29.2%)
:::

:::{card} Replication Study (Mediterranean)
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAVvRoVuJqBuWCK1XiDhVP8gS3QjO3k1k5CjPJ6Rdf_qQ
ClimateFish + MEDITS shelf, different region — **Confirmed**
:::

:::{card} Replication Outcome
:link: https://platform.sciencelive4all.org/np/?uri=https://w3id.org/sciencelive/np/RAjEJL1PpNdE7lBYjHkOQidLG9iYmyRtW2Of5PCTAotMg
SST confirmed in Mediterranean (ClimateFish R²=0.615, MEDITS R²=0.367)
:::

::::

## Key Finding

SST is confirmed as the primary driver of fish community structure when analysis is restricted to **continental shelf communities** (depth ≤ 200m). This holds across two ocean regions (NE Atlantic, Mediterranean) and three independent datasets.

:::{note}
Without the shelf depth filter, depth overwhelms the SST signal in trawl surveys that sample across wide depth ranges. The claim is specific to shelf-depth fish communities.
:::

## Reproducibility

All analyses are containerised and can be reproduced:

```bash
docker build -t bio-oracle-replication .
docker run -v $(pwd)/results:/app/results bio-oracle-replication
```

The full pipeline code is available at [github.com/annefou/bio-oracle](https://github.com/annefou/bio-oracle) and archived on Zenodo at [doi.org/10.5281/zenodo.19628179](https://doi.org/10.5281/zenodo.19628179).

## Citation

> Fouilloux, A. (2026). *Bio-ORACLE Mediterranean Replication: SST as primary driver of fish community structure*. Zenodo. [https://doi.org/10.5281/zenodo.19628179](https://doi.org/10.5281/zenodo.19628179)
