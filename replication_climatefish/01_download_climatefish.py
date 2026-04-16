"""
Step 1: Download ClimateFish abundance data from SEANOE.

ClimateFish (Azzurro et al., 2022) provides abundance counts for 15 fish species
proposed as candidate indicators of climate change in the Mediterranean Sea.
Data collected via visual census along 3,142 transects in 7 countries (2009-2021).

Source: https://www.seanoe.org/data/00756/86784/
License: CC-BY 4.0
"""

import os
from pathlib import Path

import pandas as pd
import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

CLIMATEFISH_URL = "https://www.seanoe.org/data/00756/86784/data/92824.csv"
FIELDS_URL = "https://www.seanoe.org/data/00756/86784/data/92825.csv"

SPECIES_COLUMNS = [
    "Coris_julis",
    "Epinephelus_marginatus",
    "Fistularia_commersonii",
    "Parupeneus_forskali",
    "Pempheris_rhomboidea",
    "Pterois_miles",
    "Sarpa_salpa",
    "Serranus_cabrilla",
    "Serranus_scriba",
    "Siganus_luridus",
    "Siganus_rivulatus",
    "Sparisoma_cretense",
    "Stephanolopis_diaspros",
    "Thalassoma_pavo",
    "Torquigener_flavimaculosus",
]


def download_file(url: str, dest: Path) -> Path:
    """Download a file if it doesn't already exist."""
    if dest.exists():
        print(f"  Already exists: {dest}")
        return dest
    print(f"  Downloading {url}")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"  Saved to {dest} ({len(resp.content):,} bytes)")
    return dest


def main():
    print("=== Step 1: Download ClimateFish data ===\n")

    # Download main dataset
    csv_path = download_file(CLIMATEFISH_URL, DATA_DIR / "climatefish_raw.csv")

    # Download field descriptions
    download_file(FIELDS_URL, DATA_DIR / "climatefish_fields.csv")

    # Load and inspect
    df = pd.read_csv(csv_path)
    print(f"\nDataset shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print(f"\nYear range: {df['Year'].min()} - {df['Year'].max()}")
    print(f"Countries: {sorted(df['Country'].unique())}")
    print(f"Number of transects: {len(df)}")

    # Check which species columns are present
    found_species = [c for c in df.columns if c in SPECIES_COLUMNS]
    missing_species = [c for c in SPECIES_COLUMNS if c not in df.columns]
    print(f"\nSpecies columns found: {len(found_species)}")
    if missing_species:
        print(f"  Missing (check column names): {missing_species}")

    # Basic quality summary
    print(f"\nCoordinate ranges:")
    print(f"  Latitude:  {df['Decimal_latitude'].min():.2f} to {df['Decimal_latitude'].max():.2f}")
    print(f"  Longitude: {df['Decimal_longitude'].min():.2f} to {df['Decimal_longitude'].max():.2f}")

    # Total individuals per species
    print(f"\nTotal individuals per species:")
    for sp in found_species:
        total = df[sp].sum()
        nonzero = (df[sp] > 0).sum()
        print(f"  {sp}: {total:,} individuals in {nonzero} transects")

    print(f"\nTotal individuals across all species: {df[found_species].sum().sum():,.0f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
