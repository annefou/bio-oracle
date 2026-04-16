"""
Step 2: Spatial join — match MEDITS hauls to Bio-ORACLE environmental values
and aggregate to 1x1 degree grid cells (following Rutterford et al.'s method).

Downloads Bio-ORACLE SST and salinity layers, then extracts values at each
haul location and aggregates to grid cells.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

ERDDAP_BASE = "https://erddap.bio-oracle.org/erddap/griddap"

# Mediterranean bounding box
LAT_MIN, LAT_MAX = 30.0, 46.0
LON_MIN, LON_MAX = -6.0, 36.5

LAYERS = {
    "sst": {
        "dataset_id": "thetao_baseline_2000_2019_depthsurf",
        "variable": "thetao_mean",
        "filename": "biooracle_sst_med.nc",
    },
    "salinity": {
        "dataset_id": "so_baseline_2000_2019_depthsurf",
        "variable": "so_mean",
        "filename": "biooracle_salinity_med.nc",
    },
}


def download_layer(key, info):
    """Download Bio-ORACLE layer if not present."""
    dest = DATA_DIR / info["filename"]
    if dest.exists():
        print(f"  Already exists: {dest}")
        return dest

    url = (
        f"{ERDDAP_BASE}/{info['dataset_id']}.nc?"
        f"{info['variable']}"
        f"[(last)]"
        f"[({LAT_MIN}):1:({LAT_MAX})]"
        f"[({LON_MIN}):1:({LON_MAX})]"
    )
    print(f"  Downloading {key}...")
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"  Saved ({len(resp.content):,} bytes)")
    return dest


def extract_values(lats, lons, nc_path, variable):
    """Extract Bio-ORACLE values at given coordinates."""
    ds = xr.open_dataset(nc_path)
    var = ds[variable]
    if "time" in var.dims:
        var = var.isel(time=0)

    values = []
    for lat, lon in zip(lats, lons):
        val = var.sel(latitude=lat, longitude=lon, method="nearest")
        values.append(float(val.values))
    ds.close()
    return values


def main():
    print("=== Step 2: Spatial join + grid aggregation ===\n")

    # Download Bio-ORACLE
    print("Downloading Bio-ORACLE layers...")
    for key, info in LAYERS.items():
        download_layer(key, info)

    # Load MEDITS data
    df = pd.read_csv(DATA_DIR / "medits_fish_hauls.csv", index_col=0)
    species_cols = [c for c in df.columns if c not in ["lat", "lon", "depth", "year", "country"]]
    print(f"\nLoaded {len(df):,} hauls, {len(species_cols)} species")

    # Extract SST
    print("\nExtracting SST...")
    sst_path = DATA_DIR / "biooracle_sst_med.nc"
    ds = xr.open_dataset(sst_path)
    sst_var = list(ds.data_vars)[0]
    df["sst"] = extract_values(df["lat"].values, df["lon"].values, sst_path, sst_var)
    print(f"  SST range: {df['sst'].min():.1f} to {df['sst'].max():.1f} °C")

    # Extract salinity
    print("Extracting salinity...")
    sal_path = DATA_DIR / "biooracle_salinity_med.nc"
    ds = xr.open_dataset(sal_path)
    sal_var = list(ds.data_vars)[0]
    df["salinity"] = extract_values(df["lat"].values, df["lon"].values, sal_path, sal_var)
    print(f"  Salinity range: {df['salinity'].min():.1f} to {df['salinity'].max():.1f} PSU")

    # Drop hauls with NaN environment (over land)
    n_before = len(df)
    df = df.dropna(subset=["sst", "salinity", "depth"])
    df = df[df["sst"] > 0]  # Remove any invalid SST
    print(f"\nDropped {n_before - len(df)} hauls with missing data")

    # Filter to continental shelf (following Rutterford et al.'s scope)
    MAX_DEPTH = 200
    print(f"\nFiltering to shelf depth (<= {MAX_DEPTH}m)...")
    n_before = len(df)
    df = df[(df["depth"] > 0) & (df["depth"] <= MAX_DEPTH)]
    print(f"  {n_before:,} → {len(df):,} hauls ({n_before - len(df):,} removed)")

    # Aggregate to 1x1 degree grid cells (matching Rutterford's method)
    print("\nAggregating to 1x1 degree grid cells...")
    df["grid_lat"] = df["lat"].round(0)
    df["grid_lon"] = df["lon"].round(0)
    df["grid_id"] = df["grid_lat"].astype(str) + "_" + df["grid_lon"].astype(str)

    # Mean abundance per grid cell
    grid_species = df.groupby("grid_id")[species_cols].mean()

    # Mean environment per grid cell
    grid_env = df.groupby("grid_id")[["sst", "salinity", "depth"]].mean()

    # Grid cell coordinates + metadata
    grid_coords = df.groupby("grid_id")[["grid_lat", "grid_lon"]].first()
    grid_n = df.groupby("grid_id").size().rename("n_hauls")

    grid = pd.concat([grid_species, grid_env, grid_coords, grid_n], axis=1)

    # Filter: require at least 3 hauls per grid cell
    grid = grid[grid["n_hauls"] >= 3]
    print(f"  Grid cells (>= 3 hauls): {len(grid)}")
    print(f"  Total hauls used: {grid['n_hauls'].sum():,.0f}")

    # Log-transform depth (following Rutterford)
    grid["log_depth"] = np.log10(grid["depth"])

    # Save
    output_path = DATA_DIR / "medits_grid_cells_shelf.csv"
    grid.to_csv(output_path)
    print(f"\nSaved to {output_path}")

    print(f"\n=== Grid cell summary ===")
    print(f"  Grid cells: {len(grid)}")
    print(f"  SST: {grid['sst'].min():.1f} to {grid['sst'].max():.1f} °C (range: {grid['sst'].max()-grid['sst'].min():.1f}°C)")
    print(f"  Salinity: {grid['salinity'].min():.1f} to {grid['salinity'].max():.1f} PSU")
    print(f"  Depth: {grid['depth'].min():.0f} to {grid['depth'].max():.0f} m")
    print(f"  Species: {len(species_cols)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
