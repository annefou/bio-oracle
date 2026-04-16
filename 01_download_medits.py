"""
Step 1: Load and prepare MEDITS trawl survey data from JRC.

The JRC EU DCF dataset provides standardised demersal trawl survey data from
the Mediterranean (MEDITS programme, 1994-2024). We extract fish abundance
per haul with spatial and environmental metadata.

Source: https://data.jrc.ec.europa.eu/dataset/ef36af5d-eb4e-4a9f-9513-82a877fe71fe
License: Open access, no registration required
"""

import os
import zipfile
from pathlib import Path

import pandas as pd
import numpy as np
import requests

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

ZIP_URL = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/FAD/MEDBSsurvey/2025_MEDBSsurvey.zip"
ZIP_PATH = DATA_DIR / "2025_MEDBSsurvey.zip"

# Match Rutterford et al.'s modern study period
YEAR_MIN = 2005
YEAR_MAX = 2024

# Minimum number of hauls per species to include
MIN_HAULS = 30


def download_data():
    """Download JRC MEDITS data if not already present."""
    if ZIP_PATH.exists():
        print(f"  Already exists: {ZIP_PATH}")
        return

    print(f"  Downloading JRC MEDITS data ({ZIP_URL})...")
    resp = requests.get(ZIP_URL, timeout=600, stream=True)
    resp.raise_for_status()
    ZIP_PATH.write_bytes(resp.content)
    print(f"  Saved ({ZIP_PATH.stat().st_size / 1e6:.1f} MB)")


def extract_csv():
    """Extract demersal CSV files from zip."""
    ta_path = DATA_DIR / "Demersal" / "TA.csv"
    tb_path = DATA_DIR / "Demersal" / "TB.csv"
    if ta_path.exists() and tb_path.exists():
        print("  Already extracted")
        return ta_path, tb_path

    with zipfile.ZipFile(ZIP_PATH) as zf:
        zf.extract("Demersal/TA.csv", DATA_DIR)
        zf.extract("Demersal/TB.csv", DATA_DIR)
    print(f"  Extracted TA.csv and TB.csv")
    return ta_path, tb_path


def main():
    print("=== Step 1: Load MEDITS trawl survey data ===\n")

    download_data()
    ta_path, tb_path = extract_csv()

    # Load haul metadata (TA)
    print("\nLoading haul metadata (TA)...")
    ta = pd.read_csv(ta_path, low_memory=False)
    print(f"  Total hauls: {len(ta):,}")

    # Filter: MEDITS survey, valid hauls, study period
    ta = ta[
        (ta["name_of_survey"] == "MEDITS")
        & (ta["validity"] == "V")
        & (ta["year"] >= YEAR_MIN)
        & (ta["year"] <= YEAR_MAX)
    ].copy()
    print(f"  MEDITS valid hauls ({YEAR_MIN}-{YEAR_MAX}): {len(ta):,}")

    # Convert coordinates from DDMM.MM to decimal degrees
    def convert_coord(val):
        degrees = int(val / 100)
        minutes = val - degrees * 100
        return degrees + minutes / 60

    ta["lat"] = ta["shooting_latitude"].apply(convert_coord)
    ta["lon"] = ta["shooting_longitude"].apply(convert_coord)
    # Handle quadrant (1=NE, 2=NW, 3=SW, 4=SE)
    ta.loc[ta["shooting_quadrant"].isin([3, 4]), "lat"] *= -1
    ta.loc[ta["shooting_quadrant"].isin([2, 3]), "lon"] *= -1

    ta["depth"] = (ta["shooting_depth"] + ta["hauling_depth"]) / 2

    # Create unique haul ID
    ta["haul_id"] = (
        ta["country"] + "_" + ta["area"].astype(str) + "_"
        + ta["vessel"] + "_" + ta["year"].astype(str) + "_"
        + ta["haul_number"].astype(str)
    )

    print(f"  Countries: {sorted(ta['country'].unique())}")
    print(f"  Years: {ta['year'].min()} - {ta['year'].max()}")
    print(f"  Lat: {ta['lat'].min():.2f} to {ta['lat'].max():.2f}")
    print(f"  Lon: {ta['lon'].min():.2f} to {ta['lon'].max():.2f}")
    print(f"  Depth: {ta['depth'].min():.0f} to {ta['depth'].max():.0f} m")

    # Load catch data (TB)
    print("\nLoading catch data (TB)...")
    tb = pd.read_csv(tb_path, low_memory=False)

    # Filter: MEDITS, fish only (Ao = Osteichthyes, A = Actinopterygii)
    tb = tb[
        (tb["name_of_survey"] == "MEDITS")
        & (tb["catfau"].isin(["A", "Ao"]))
    ].copy()

    # Create matching haul ID
    tb["haul_id"] = (
        tb["country"] + "_" + tb["area"].astype(str) + "_"
        + tb["vessel"] + "_" + tb["year"].astype(str) + "_"
        + tb["haul_number"].astype(str)
    )

    # Create species name
    tb["species_name"] = tb["genus"] + "_" + tb["species"]

    # Keep only hauls in our filtered TA
    valid_hauls = set(ta["haul_id"])
    tb = tb[tb["haul_id"].isin(valid_hauls)]

    # Aggregate: total abundance per species per haul
    catch = tb.groupby(["haul_id", "species_name"])["nbtot"].sum().reset_index()
    catch.columns = ["haul_id", "species_name", "abundance"]

    print(f"  Fish catch records: {len(catch):,}")
    n_species_total = catch["species_name"].nunique()
    print(f"  Total fish taxa: {n_species_total}")

    # Filter species: require presence in at least MIN_HAULS hauls
    species_counts = catch.groupby("species_name")["haul_id"].nunique()
    common_species = species_counts[species_counts >= MIN_HAULS].index.tolist()
    catch = catch[catch["species_name"].isin(common_species)]
    print(f"  Species with >= {MIN_HAULS} hauls: {len(common_species)}")

    # Pivot to site x species matrix
    abundance_matrix = catch.pivot_table(
        index="haul_id", columns="species_name", values="abundance", fill_value=0
    )
    print(f"  Abundance matrix: {abundance_matrix.shape[0]} hauls x {abundance_matrix.shape[1]} species")

    # Merge with haul metadata
    haul_info = ta[["haul_id", "lat", "lon", "depth", "year", "country"]].drop_duplicates("haul_id")
    haul_info = haul_info.set_index("haul_id")

    result = abundance_matrix.join(haul_info, how="inner")
    print(f"  Matched hauls: {len(result):,}")

    # Save
    output_path = DATA_DIR / "medits_fish_hauls.csv"
    result.to_csv(output_path)
    print(f"\nSaved to {output_path}")

    # Save species list
    species_path = DATA_DIR / "medits_species_list.csv"
    species_summary = pd.DataFrame({
        "species": common_species,
        "n_hauls": [species_counts[s] for s in common_species],
    }).sort_values("n_hauls", ascending=False)
    species_summary.to_csv(species_path, index=False)
    print(f"Species list saved to {species_path}")

    # Summary
    print(f"\n=== Summary ===")
    print(f"  Hauls: {len(result):,}")
    print(f"  Species: {abundance_matrix.shape[1]}")
    print(f"  Countries: {sorted(result['country'].unique())}")
    print(f"  Period: {YEAR_MIN}-{YEAR_MAX}")
    print(f"  Depth range: {result['depth'].min():.0f}-{result['depth'].max():.0f} m")

    print("\nDone.")


if __name__ == "__main__":
    main()
