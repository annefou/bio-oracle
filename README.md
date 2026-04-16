# SST as Primary Driver of Fish Community Structure

Reproduction and replication of [Rutterford et al. (2023)](https://doi.org/10.1111/gcb.16633) testing whether sea surface temperature (SST) is the primary driver of fish community structure.

**Three independent analyses all confirm the claim for continental shelf communities.**

| Analysis | Dataset | Region | Method | SST R² (Axis 1) | Outcome |
|---|---|---|---|---|---|
| Reproduction | ICES DATRAS | NE Atlantic | Trawl | 0.573 | Confirmed |
| Replication | ClimateFish | Mediterranean | Visual census | 0.615 | Confirmed |
| Replication | MEDITS (shelf) | Mediterranean | Trawl | 0.367 | Confirmed |

## Quick Start

```bash
# Setup
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# Reproduction (DATRAS — NE Atlantic)
python reproduction/01_download_datras.py          # ~3 hours, cached
python reproduction/02_download_biooracle.py
python reproduction/03_filter_and_join.py
python reproduction/04_community_analysis.py

# Replication (ClimateFish — Mediterranean visual census)
python replication_climatefish/01_download_climatefish.py
python replication_climatefish/02_download_biooracle.py
python replication_climatefish/03_spatial_join.py
python replication_climatefish/04_community_analysis.py
python replication_climatefish/05_variable_importance.py

# Replication (MEDITS — Mediterranean trawl, shelf only)
python replication_medits/01_download_medits.py
python replication_medits/02_spatial_join.py
python replication_medits/03_community_analysis.py
```

## Docker

```bash
docker build -t bio-oracle-replication .
docker run -v $(pwd)/results:/app/results bio-oracle-replication
```

## Repository Structure

```
bio-oracle/
├── reproduction/              # DATRAS NE Atlantic (same data as Rutterford)
├── replication_climatefish/   # ClimateFish Mediterranean visual census
├── replication_medits/        # MEDITS Mediterranean trawl (shelf ≤200m)
├── data/                      # Shared data directory
├── results/                   # All outputs
├── myst.yml                   # Jupyter Book configuration
├── index.md                   # Book landing page
├── Dockerfile
└── requirements.txt
```

## Data Sources

- [ICES DATRAS](https://datras.ices.dk/) — NE Atlantic trawl surveys (ICES Data Policy)
- [ClimateFish](https://www.seanoe.org/data/00756/86784/) — Mediterranean visual census (CC-BY 4.0)
- [JRC MEDITS](https://data.jrc.ec.europa.eu/dataset/ef36af5d-eb4e-4a9f-9513-82a877fe71fe) — Mediterranean trawl surveys
- [Bio-ORACLE](https://bio-oracle.org/) — Marine environmental layers

## Citation

See [CITATION.cff](CITATION.cff) for citation metadata.

## License

Code: MIT | Content: CC-BY-4.0
