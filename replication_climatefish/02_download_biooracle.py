"""
Step 2: Download Bio-ORACLE environmental layers for the Mediterranean.

We download sea surface temperature (SST) and sea surface salinity (SSS)
from Bio-ORACLE v3 via the ERDDAP server. Depth (bathymetry) is derived
from the ETOPO dataset or Bio-ORACLE depth layers.

Bio-ORACLE: https://bio-oracle.org/
ERDDAP: https://erddap.bio-oracle.org/erddap/
Reference: Assis et al. (2018), DOI: 10.1111/geb.12693
"""

from pathlib import Path

import xarray as xr
import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

# Mediterranean bounding box (with some margin)
LAT_MIN, LAT_MAX = 30.0, 46.0
LON_MIN, LON_MAX = -6.0, 36.5

# Bio-ORACLE ERDDAP base
ERDDAP_BASE = "https://erddap.bio-oracle.org/erddap/griddap"

# Datasets: baseline 2000-2019 (overlaps with ClimateFish 2009-2021)
LAYERS = {
    "sst": {
        "dataset_id": "thetao_baseline_2000_2019_depthsurf",
        "variable": "thetao_mean",
        "filename": "biooracle_sst_med.nc",
        "description": "Sea surface temperature (mean, 2000-2019)",
    },
    "salinity": {
        "dataset_id": "so_baseline_2000_2019_depthsurf",
        "variable": "so_mean",
        "filename": "biooracle_salinity_med.nc",
        "description": "Sea surface salinity (mean, 2000-2019)",
    },
}


def download_layer(layer_key: str, layer_info: dict) -> Path:
    """Download a Bio-ORACLE layer for the Mediterranean region via ERDDAP."""
    dest = DATA_DIR / layer_info["filename"]
    if dest.exists():
        print(f"  Already exists: {dest}")
        return dest

    dataset_id = layer_info["dataset_id"]
    variable = layer_info["variable"]

    # ERDDAP griddap query for subset
    # Dataset has dimensions: time, latitude, longitude
    # We take the last available time step (climatological mean)
    url = (
        f"{ERDDAP_BASE}/{dataset_id}.nc?"
        f"{variable}"
        f"[(last)]"
        f"[({LAT_MIN}):1:({LAT_MAX})]"
        f"[({LON_MIN}):1:({LON_MAX})]"
    )

    print(f"  Downloading {layer_info['description']}...")
    print(f"  URL: {url[:120]}...")
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"  Saved to {dest} ({len(resp.content):,} bytes)")
    return dest


def main():
    print("=== Step 2: Download Bio-ORACLE environmental layers ===\n")
    print(f"Mediterranean bounding box:")
    print(f"  Lat: {LAT_MIN} to {LAT_MAX}")
    print(f"  Lon: {LON_MIN} to {LON_MAX}\n")

    for key, info in LAYERS.items():
        print(f"Layer: {key}")
        path = download_layer(key, info)

        # Inspect downloaded file
        ds = xr.open_dataset(path)
        print(f"  Variables: {list(ds.data_vars)}")
        print(f"  Dimensions: {dict(ds.dims)}")
        var = list(ds.data_vars)[0]
        print(f"  {var} range: {float(ds[var].min()):.2f} to {float(ds[var].max()):.2f}")
        ds.close()
        print()

    print("Done.")


if __name__ == "__main__":
    main()
