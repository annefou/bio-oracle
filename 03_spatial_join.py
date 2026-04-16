"""
Step 3: Spatial join — match ClimateFish transects to Bio-ORACLE environmental values.

For each ClimateFish transect (lat/lon), extract the corresponding Bio-ORACLE
sea surface temperature and salinity values from the nearest grid cell.

Output: a joined CSV with species abundances + environmental variables per transect.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

DATA_DIR = Path("data")

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


def extract_biooracle_values(lats, lons, nc_path, variable):
    """Extract values from a Bio-ORACLE NetCDF at given lat/lon coordinates."""
    ds = xr.open_dataset(nc_path)
    var = ds[variable]

    # Squeeze out single-element time dimension if present
    if "time" in var.dims:
        var = var.isel(time=0)

    values = []
    for lat, lon in zip(lats, lons):
        val = var.sel(latitude=lat, longitude=lon, method="nearest")
        values.append(float(val.values))

    ds.close()
    return values


def main():
    print("=== Step 3: Spatial join (ClimateFish + Bio-ORACLE) ===\n")

    # Load ClimateFish data
    fish_path = DATA_DIR / "climatefish_raw.csv"
    df = pd.read_csv(fish_path)
    print(f"ClimateFish transects: {len(df)}")

    # Identify species columns present in the data
    species_cols = [c for c in SPECIES_COLUMNS if c in df.columns]
    if not species_cols:
        # Try case-insensitive match or check actual column names
        print(f"  Available columns: {list(df.columns)}")
        # Attempt flexible matching
        for expected in SPECIES_COLUMNS:
            matches = [c for c in df.columns if expected.lower() in c.lower()]
            if matches:
                species_cols.append(matches[0])
    print(f"Species columns matched: {len(species_cols)}")

    lats = df["Decimal_latitude"].values
    lons = df["Decimal_longitude"].values

    # Extract SST values
    print("\nExtracting SST from Bio-ORACLE...")
    sst_path = DATA_DIR / "biooracle_sst_med.nc"
    ds_sst = xr.open_dataset(sst_path)
    sst_var = list(ds_sst.data_vars)[0]
    print(f"  Variable name: {sst_var}")
    df["sst_mean"] = extract_biooracle_values(lats, lons, sst_path, sst_var)
    print(f"  SST range: {df['sst_mean'].min():.2f} to {df['sst_mean'].max():.2f} C")

    # Extract salinity values
    print("\nExtracting salinity from Bio-ORACLE...")
    sal_path = DATA_DIR / "biooracle_salinity_med.nc"
    ds_sal = xr.open_dataset(sal_path)
    sal_var = list(ds_sal.data_vars)[0]
    print(f"  Variable name: {sal_var}")
    df["salinity_mean"] = extract_biooracle_values(lats, lons, sal_path, sal_var)
    print(f"  Salinity range: {df['salinity_mean'].min():.2f} to {df['salinity_mean'].max():.2f} PSU")

    # Use the depth column from ClimateFish as a categorical proxy
    # (depth ranges: 0-3m, 5-10m, 11-20m, 21-30m)
    print(f"\nDepth categories in ClimateFish: {sorted(df['Depth'].unique())}")

    # Create a numeric depth midpoint for regression
    depth_map = {
        "0_3": 1.5,
        "1_3": 2.0,
        "5_10": 7.5,
        "11_20": 15.5,
        "21_30": 25.5,
        # Alternative formats
        "0-3": 1.5,
        "1-3": 2.0,
        "5-10": 7.5,
        "11-20": 15.5,
        "21-30": 25.5,
    }
    df["depth_midpoint"] = df["Depth"].map(depth_map)
    unmapped = df["depth_midpoint"].isna().sum()
    if unmapped > 0:
        print(f"  Warning: {unmapped} transects with unmapped depth values")
        print(f"  Unique depth values: {df['Depth'].unique()}")

    # Drop rows with missing environmental data
    env_cols = ["sst_mean", "salinity_mean", "depth_midpoint"]
    n_before = len(df)
    df = df.dropna(subset=env_cols + species_cols)
    print(f"\nDropped {n_before - len(df)} rows with missing data")
    print(f"Final dataset: {len(df)} transects")

    # Save joined dataset
    output_path = DATA_DIR / "fish_environment_joined.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved to {output_path}")

    # Summary statistics
    print(f"\nEnvironmental summary:")
    for col in env_cols:
        print(f"  {col}: mean={df[col].mean():.2f}, std={df[col].std():.2f}")

    print(f"\nSpecies prevalence (% transects with presence):")
    for sp in species_cols:
        prevalence = (df[sp] > 0).mean() * 100
        print(f"  {sp}: {prevalence:.1f}%")

    print("\nDone.")


if __name__ == "__main__":
    main()
