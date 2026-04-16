"""
Step 3: Spatial join — match DATRAS hauls to Bio-ORACLE environmental values
and aggregate to 1x1 degree grid cells (following Rutterford et al.'s method).

Rutterford et al. aggregated CPUE to 1x1° grid cells, requiring at least
3 hauls from each survey operating in that area.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

LAYERS = {
    "sst": {
        "filename": "biooracle_sst_neatl.nc",
    },
    "salinity": {
        "filename": "biooracle_salinity_neatl.nc",
    },
}


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
    print("=== Step 3: Spatial join + grid aggregation (DATRAS) ===\n")

    # Load DATRAS haul data
    df = pd.read_csv(DATA_DIR / "datras_fish_hauls.csv", index_col=0)
    species_cols = [c for c in df.columns if c.startswith("sp_")]
    meta_cols = [c for c in df.columns if c not in species_cols]
    print(f"Loaded {len(df):,} hauls, {len(species_cols)} species")
    print(f"Metadata columns: {meta_cols}")

    # Extract SST
    print("\nExtracting SST at haul locations...")
    sst_path = DATA_DIR / "biooracle_sst_neatl.nc"
    df["sst"] = extract_values(df["lat"].values, df["lon"].values, sst_path)
    print(f"  SST range: {df['sst'].min():.1f} to {df['sst'].max():.1f} °C")

    # Extract salinity
    print("Extracting salinity...")
    sal_path = DATA_DIR / "biooracle_salinity_neatl.nc"
    df["salinity"] = extract_values(df["lat"].values, df["lon"].values, sal_path)
    print(f"  Salinity range: {df['salinity'].min():.1f} to {df['salinity'].max():.1f} PSU")

    # Drop hauls with missing environment (over land / outside Bio-ORACLE grid)
    n_before = len(df)
    df = df.dropna(subset=["sst", "salinity", "depth"])
    df = df[(df["sst"] > 0) & (df["sst"] < 40)]
    df = df[(df["salinity"] > 0) & (df["salinity"] < 45)]
    df = df[df["depth"] > 0]
    print(f"\nDropped {n_before - len(df)} hauls with invalid data, {len(df):,} remaining")

    # Aggregate to 1x1 degree grid cells (matching Rutterford et al.)
    print("\nAggregating to 1x1 degree grid cells...")
    df["grid_lat"] = np.floor(df["lat"]) + 0.5
    df["grid_lon"] = np.floor(df["lon"]) + 0.5
    df["grid_id"] = df["grid_lat"].astype(str) + "_" + df["grid_lon"].astype(str)

    # Mean CPUE per grid cell per species
    grid_species = df.groupby("grid_id")[species_cols].mean()

    # Mean environment per grid cell
    grid_env = df.groupby("grid_id")[["sst", "salinity", "depth"]].mean()

    # Grid cell coordinates + haul count
    grid_coords = df.groupby("grid_id")[["grid_lat", "grid_lon"]].first()
    grid_n = df.groupby("grid_id").size().rename("n_hauls")

    # Number of unique surveys per grid cell
    grid_n_surveys = df.groupby("grid_id")["Survey"].nunique().rename("n_surveys")

    grid = pd.concat([grid_species, grid_env, grid_coords, grid_n, grid_n_surveys], axis=1)

    # Filter: require at least 3 hauls per grid cell
    grid = grid[grid["n_hauls"] >= 3]
    print(f"  Grid cells (>= 3 hauls): {len(grid)}")
    print(f"  Total hauls represented: {grid['n_hauls'].sum():,.0f}")

    # Log-transform depth (following Rutterford)
    grid["log_depth"] = np.log10(grid["depth"])

    # Remove species that are absent from all remaining grid cells
    species_present = grid[species_cols].sum() > 0
    active_species = species_present[species_present].index.tolist()
    grid = grid.drop(columns=[c for c in species_cols if c not in active_species])
    print(f"  Active species in grid: {len(active_species)}")

    # Save
    output_path = DATA_DIR / "datras_grid_cells.csv"
    grid.to_csv(output_path)
    print(f"\nSaved to {output_path}")

    print(f"\n=== Grid cell summary ===")
    print(f"  Grid cells: {len(grid)}")
    print(f"  SST: {grid['sst'].min():.1f} to {grid['sst'].max():.1f} °C "
          f"(range: {grid['sst'].max() - grid['sst'].min():.1f}°C)")
    print(f"  Salinity: {grid['salinity'].min():.1f} to {grid['salinity'].max():.1f} PSU")
    print(f"  Depth: {grid['depth'].min():.0f} to {grid['depth'].max():.0f} m")
    print(f"  Species: {len(active_species)}")

    # Compare with Rutterford
    print(f"\n  Rutterford et al.: 193 grid cells, 198 species")
    print(f"  This study:        {len(grid)} grid cells, {len(active_species)} species")

    print("\nDone.")


if __name__ == "__main__":
    main()
