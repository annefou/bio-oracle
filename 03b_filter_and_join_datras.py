"""
Step 3b: Filter DATRAS data to match Rutterford et al. (2023) more closely.

Two key filters:
  1. Fish only — use WoRMS API to identify Actinopterygii (bony fish) and
     Elasmobranchii (sharks/rays), excluding invertebrates and other taxa.
  2. Continental shelf — restrict to hauls with depth <= 200m.

Then: spatial join with Bio-ORACLE + aggregate to 1x1° grid cells.
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

CACHE_FILE = DATA_DIR / "datras_fish_aphia_ids.json"

# Maximum shelf depth (Rutterford focused on continental shelf)
MAX_DEPTH = 200

# Exclude Baltic surveys (different salinity regime, not NE Atlantic shelf)
EXCLUDE_SURVEYS = {"BITS", "SE-SOUND"}

# Fish classes in WoRMS taxonomy
FISH_TAXA = {"Actinopterygii", "Actinopteri", "Teleostei",
             "Elasmobranchii", "Chondrichthyes", "Osteichthyes",
             "Myxini", "Petromyzonti"}


def classify_species(aphia_ids):
    """Use WoRMS API to identify which AphiaIDs are fish."""
    if CACHE_FILE.exists():
        with open(CACHE_FILE) as f:
            cached = json.load(f)
        # Convert string keys back to int
        cached = {int(k): v for k, v in cached.items()}
        missing = [a for a in aphia_ids if a not in cached]
        if not missing:
            print(f"  All {len(aphia_ids)} species classified (cached)")
            return cached
        print(f"  {len(cached)} cached, {len(missing)} to classify")
    else:
        cached = {}
        missing = list(aphia_ids)

    print(f"  Classifying {len(missing)} species via WoRMS API...")
    for i, aphia_id in enumerate(missing):
        if i > 0 and i % 50 == 0:
            print(f"    {i}/{len(missing)}...")

        try:
            url = f"https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{int(aphia_id)}"
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                # Walk the classification tree
                is_fish = False
                node = data
                while node:
                    name = node.get("scientificname", "")
                    if name in FISH_TAXA:
                        is_fish = True
                        break
                    node = node.get("child")
                cached[int(aphia_id)] = is_fish
            else:
                cached[int(aphia_id)] = False
        except Exception as e:
            print(f"    Warning: {aphia_id}: {e}")
            cached[int(aphia_id)] = False

        # Rate limit: WoRMS asks for max 1 req/sec
        time.sleep(0.3)

    # Save cache
    with open(CACHE_FILE, "w") as f:
        json.dump({str(k): v for k, v in cached.items()}, f)
    print(f"  Classification saved to {CACHE_FILE}")

    return cached


def extract_values(lats, lons, nc_path):
    """Extract Bio-ORACLE values at given coordinates using nearest neighbour."""
    ds = xr.open_dataset(nc_path)
    var_name = list(ds.data_vars)[0]
    var = ds[var_name]
    if "time" in var.dims:
        var = var.isel(time=0)

    values = []
    for lat, lon in zip(lats, lons):
        val = var.sel(latitude=lat, longitude=lon, method="nearest")
        values.append(float(val.values))
    ds.close()
    return values


def main():
    print("=== Step 3b: Filter + spatial join (DATRAS → Rutterford match) ===\n")

    # Load DATRAS haul data
    df = pd.read_csv(DATA_DIR / "datras_fish_hauls.csv", index_col=0)
    species_cols = [c for c in df.columns if c.startswith("sp_")]
    print(f"Loaded: {len(df):,} hauls, {len(species_cols)} species")

    # === Filter 0: Exclude Baltic surveys ===
    if EXCLUDE_SURVEYS:
        print(f"\n--- Filter 0: Exclude Baltic/non-Atlantic surveys ---")
        n_before = len(df)
        df = df[~df["Survey"].isin(EXCLUDE_SURVEYS)].copy()
        print(f"  Excluded: {EXCLUDE_SURVEYS}")
        print(f"  {n_before:,} → {len(df):,} hauls ({n_before - len(df):,} removed)")

    # === Filter 1: Depth ===
    print(f"\n--- Filter 1: Depth <= {MAX_DEPTH}m (continental shelf) ---")
    n_before = len(df)
    df = df[(df["depth"] > 0) & (df["depth"] <= MAX_DEPTH)].copy()
    print(f"  {n_before:,} → {len(df):,} hauls ({n_before - len(df):,} removed)")

    # === Filter 2: Fish species only ===
    print(f"\n--- Filter 2: Fish species only (WoRMS classification) ---")
    aphia_ids = [int(c.replace("sp_", "")) for c in species_cols]
    classification = classify_species(aphia_ids)

    fish_cols = [c for c in species_cols
                 if classification.get(int(c.replace("sp_", "")), False)]
    non_fish_cols = [c for c in species_cols if c not in fish_cols]
    print(f"  Fish species: {len(fish_cols)}")
    print(f"  Non-fish (removed): {len(non_fish_cols)}")

    # Drop non-fish columns
    df = df.drop(columns=non_fish_cols)
    species_cols = fish_cols

    # Remove species with no occurrences after depth filter
    species_sums = df[species_cols].sum()
    active = species_sums[species_sums > 0].index.tolist()
    dropped = [c for c in species_cols if c not in active]
    if dropped:
        df = df.drop(columns=dropped)
        species_cols = active
        print(f"  Species absent after depth filter: {len(dropped)} removed")
    print(f"  Final species: {len(species_cols)}")

    # === Spatial join with Bio-ORACLE ===
    print(f"\n--- Extracting Bio-ORACLE environment ---")

    sst_path = DATA_DIR / "biooracle_sst_neatl.nc"
    sal_path = DATA_DIR / "biooracle_salinity_neatl.nc"

    # Defragment DataFrame before adding columns
    df = df.copy()

    print("  Extracting SST...")
    df["sst"] = extract_values(df["lat"].values, df["lon"].values, sst_path)
    print(f"  SST range: {df['sst'].min():.1f} to {df['sst'].max():.1f} °C")

    print("  Extracting salinity...")
    df["salinity"] = extract_values(df["lat"].values, df["lon"].values, sal_path)
    print(f"  Salinity range: {df['salinity'].min():.1f} to {df['salinity'].max():.1f} PSU")

    # Drop invalid environment
    n_before = len(df)
    df = df.dropna(subset=["sst", "salinity"])
    df = df[(df["sst"] > 0) & (df["sst"] < 40)]
    df = df[(df["salinity"] > 0) & (df["salinity"] < 45)]
    print(f"  Dropped {n_before - len(df)} hauls with invalid environment")

    # === Aggregate to 1x1° grid cells ===
    print(f"\n--- Aggregating to 1x1° grid cells ---")
    df["grid_lat"] = np.floor(df["lat"]) + 0.5
    df["grid_lon"] = np.floor(df["lon"]) + 0.5
    df["grid_id"] = df["grid_lat"].astype(str) + "_" + df["grid_lon"].astype(str)

    grid_species = df.groupby("grid_id")[species_cols].mean()
    grid_env = df.groupby("grid_id")[["sst", "salinity", "depth"]].mean()
    grid_coords = df.groupby("grid_id")[["grid_lat", "grid_lon"]].first()
    grid_n = df.groupby("grid_id").size().rename("n_hauls")
    grid_n_surveys = df.groupby("grid_id")["Survey"].nunique().rename("n_surveys")

    grid = pd.concat([grid_species, grid_env, grid_coords, grid_n, grid_n_surveys], axis=1)

    # Require at least 3 hauls per grid cell
    grid = grid[grid["n_hauls"] >= 3]

    # Log-transform depth
    grid["log_depth"] = np.log10(grid["depth"])

    # Remove absent species
    sp_present = grid[species_cols].sum() > 0
    active_sp = sp_present[sp_present].index.tolist()
    removed_sp = [c for c in species_cols if c not in active_sp]
    if removed_sp:
        grid = grid.drop(columns=removed_sp)
    species_cols = active_sp

    # Save
    output_path = DATA_DIR / "datras_grid_cells_filtered.csv"
    grid.to_csv(output_path)

    print(f"\n=== Filtered dataset summary ===")
    print(f"  Grid cells: {len(grid)} (Rutterford: 193)")
    print(f"  Species: {len(species_cols)} (Rutterford: 198)")
    print(f"  SST: {grid['sst'].min():.1f} to {grid['sst'].max():.1f} °C "
          f"(range: {grid['sst'].max() - grid['sst'].min():.1f}°C)")
    print(f"  Salinity: {grid['salinity'].min():.1f} to {grid['salinity'].max():.1f} PSU")
    print(f"  Depth: {grid['depth'].min():.0f} to {grid['depth'].max():.0f} m "
          f"(max: {MAX_DEPTH}m)")
    print(f"  Hauls: {grid['n_hauls'].sum():,.0f}")
    print(f"\nSaved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
