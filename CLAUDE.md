# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Testing whether sea surface temperature (SST) is the primary driver of fish community structure, as claimed by Rutterford et al. (2023) for the NE Atlantic. Three levels of evidence, all confirmed:

1. **Reproduction** (`reproduction/`) — same DATRAS data, same method → **confirmed**
2. **Replication — ClimateFish** (`replication_climatefish/`) — Mediterranean visual census → **confirmed**
3. **Replication — MEDITS** (`replication_medits/`) — Mediterranean trawl, shelf ≤200m → **confirmed**

Results are published as a FORRT replication chain of nanopublications on Science Live.

## Commands

```bash
# Local setup
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# --- Reproduction (DATRAS, NE Atlantic) ---
python reproduction/01_download_datras.py          # ~3 hours first run, cached after
python reproduction/02_download_biooracle.py
python reproduction/03_filter_and_join.py           # Fish only, shelf depth, no Baltic
python reproduction/04_community_analysis.py

# --- Replication: ClimateFish (Mediterranean visual census) ---
python replication_climatefish/01_download_climatefish.py
python replication_climatefish/02_download_biooracle.py
python replication_climatefish/03_spatial_join.py
python replication_climatefish/04_community_analysis.py
python replication_climatefish/05_variable_importance.py

# --- Replication: MEDITS (Mediterranean trawl, shelf only) ---
python replication_medits/01_download_medits.py
python replication_medits/02_spatial_join.py         # Includes shelf depth filter ≤200m
python replication_medits/03_community_analysis.py

# Docker (recommended for reproducibility)
docker build -t bio-oracle-replication .
docker run -v $(pwd)/results:/app/results bio-oracle-replication
```

## Repository Structure

```
bio-oracle/
├── reproduction/              # DATRAS NE Atlantic (same data as Rutterford)
│   ├── 01_download_datras.py
│   ├── 02_download_biooracle.py
│   ├── 03_filter_and_join.py
│   └── 04_community_analysis.py
├── replication_climatefish/   # Mediterranean visual census
│   ├── 01_download_climatefish.py
│   ├── 02_download_biooracle.py
│   ├── 03_spatial_join.py
│   ├── 04_community_analysis.py
│   └── 05_variable_importance.py
├── replication_medits/        # Mediterranean trawl (shelf ≤200m)
│   ├── 01_download_medits.py
│   ├── 02_spatial_join.py
│   └── 03_community_analysis.py
├── data/                      # Shared data directory (cached downloads)
├── results/                   # All outputs (regression tables, plots)
├── Dockerfile
├── requirements.txt
├── CITATION.cff
├── codemeta.json
└── LICENSE
```

## Key Results

### Reproduction (DATRAS): SST confirmed as primary driver
- PCoA Axis 1: 29.2% variance, SST R²=0.573, full model R²=0.879
- Rutterford et al.: 29% variance, SST R²=0.890
- 247 grid cells, 258 fish species, 22 surveys, 2005-2018
- Critical filters: no Baltic (BITS, SE-SOUND), shelf ≤200m, fish only (WoRMS)

### Replication — ClimateFish: SST confirmed in Mediterranean
- PCoA Axis 1: 45.2% variance, SST r=-0.784, R²=0.615, p<0.001
- 15 species, 28 sites, visual census, 2009-2021

### Replication — MEDITS shelf: SST confirmed in Mediterranean trawl
- PCoA Axis 1: 18.9% variance, SST r=-0.605, R²=0.367, p<0.001
- 115 grid cells, 230 species, shelf depth ≤200m
- Without depth filter: depth dominates (R²=0.870) — claim only holds for shelf communities

## Parent Project

This repo lives under **ScienceLive/**. The main platform is at `../science-live-platform/`.

## Data Sources

- **ICES DATRAS**: https://datras.ices.dk/ (ICES Data Policy, open access)
- **ClimateFish**: https://www.seanoe.org/data/00756/86784/ (CC-BY 4.0)
- **MEDITS/JRC**: https://data.jrc.ec.europa.eu/dataset/ef36af5d-eb4e-4a9f-9513-82a877fe71fe
- **Bio-ORACLE ERDDAP**: https://erddap.bio-oracle.org/erddap/
- **Bio-ORACLE layers**: `thetao_baseline_2000_2019_depthsurf` (SST), `so_baseline_2000_2019_depthsurf` (salinity)
- **WoRMS**: https://www.marinespecies.org/ (species classification for fish filtering)

## Research Software Metadata

- `CITATION.cff` — citation metadata (CFF format)
- `codemeta.json` — CodeMeta software metadata
- `LICENSE` — MIT
