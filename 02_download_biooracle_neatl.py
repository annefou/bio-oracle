"""
Step 2: Download Bio-ORACLE environmental layers for the NE Atlantic.

Downloads sea surface temperature (SST) and salinity from Bio-ORACLE ERDDAP
for the NE Atlantic region matching Rutterford et al.'s study area
(southern Portugal to northern Norway).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

ERDDAP_BASE = "https://erddap.bio-oracle.org/erddap/griddap"

# NE Atlantic bounding box (Portugal to Norway)
LAT_MIN, LAT_MAX = 34.0, 72.0
LON_MIN, LON_MAX = -15.0, 30.0

LAYERS = {
    "sst": {
        "dataset_id": "thetao_baseline_2000_2019_depthsurf",
        "variable": "thetao_mean",
        "filename": "biooracle_sst_neatl.nc",
    },
    "salinity": {
        "dataset_id": "so_baseline_2000_2019_depthsurf",
        "variable": "so_mean",
        "filename": "biooracle_salinity_neatl.nc",
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
    print(f"  Downloading {key} from Bio-ORACLE ERDDAP...")
    print(f"  URL: {url}")
    resp = requests.get(url, timeout=600)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"  Saved ({len(resp.content) / 1e6:.1f} MB)")
    return dest


def main():
    print("=== Step 2: Download Bio-ORACLE layers (NE Atlantic) ===\n")
    print(f"Bounding box: {LAT_MIN}-{LAT_MAX}°N, {LON_MIN}-{LON_MAX}°E\n")

    for key, info in LAYERS.items():
        path = download_layer(key, info)
        ds = xr.open_dataset(path)
        var = list(ds.data_vars)[0]
        data = ds[var]
        if "time" in data.dims:
            data = data.isel(time=0)
        print(f"  {key}: shape={data.shape}, "
              f"range={float(data.min()):.2f} to {float(data.max()):.2f}")
        ds.close()
        print()

    print("Done.")


if __name__ == "__main__":
    main()
