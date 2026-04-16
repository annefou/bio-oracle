"""
Step 1: Download ICES DATRAS trawl survey data for the NE Atlantic.

Reproduces the dataset used by Rutterford et al. (2023): standardised CPUE
for fish species across NE Atlantic bottom trawl surveys, 2005-2018.

Data source: ICES DATRAS (Database of Trawl Surveys)
API: https://datras.ices.dk/WebServices/DATRASWebService.asmx
License: ICES Data Policy (open access, citation required)
"""

import os
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)
CACHE_DIR = DATA_DIR / "datras_cache"
CACHE_DIR.mkdir(exist_ok=True)

# DATRAS Web Service base URL
WS_BASE = "https://datras.ices.dk/WebServices/DATRASWebService.asmx"
NS = {"d": "ices.dk.local/DATRAS"}

# Study period (matching Rutterford et al. 2023)
YEAR_MIN = 2005
YEAR_MAX = 2018

# NE Atlantic bottom trawl surveys (23 surveys, Portugal to Norway)
# These cover the geographic range described by Rutterford et al.
SURVEYS = [
    "BITS",       # Baltic International Trawl Survey
    "BTS",        # Beam Trawl Survey (North Sea)
    "BTS-VIII",   # Beam Trawl Survey (Bay of Biscay)
    "DWS",        # Deep Water Survey
    "DYFS",       # Demersal Young Fish Survey
    "EVHOE",      # French Southern Atlantic BTS
    "FR-CGFS",    # French Channel Ground Fish Survey
    "FR-WCGFS",   # French Western Channel GFS
    "IE-IGFS",    # Irish Ground Fish Survey
    "NIGFS",      # Northern Ireland Ground Fish Survey
    "NL-BSAS",    # Netherlands Beam Trawl Survey
    "NS-IBTS",    # North Sea International BTS
    "NSSS",       # North Sea Sandeel Survey
    "PT-IBTS",    # Portuguese IBTS
    "ROCKALL",    # Scottish Rockall Survey (pre-2011)
    "SCOROC",     # Scottish Rockall Survey (post-2011)
    "SCOWCGFS",   # Scottish West Coast GFS (post-2011)
    "SE-SOUND",   # Sound Survey (Denmark/Sweden)
    "SNS",        # Sole Net Survey
    "SP-ARSA",    # Spanish Gulf of Cadiz
    "SP-NORTH",   # Spanish North Coast
    "SP-PORC",    # Spanish Porcupine Bank
    "SWC-IBTS",   # Scottish West Coast IBTS (pre-2011)
]

# Minimum hauls per species (across the entire dataset)
MIN_HAULS = 30


def fetch_xml(endpoint, params, retries=3, delay=5):
    """Fetch XML from DATRAS API with retry logic."""
    url = f"{WS_BASE}/{endpoint}"
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=300, verify=False)
            resp.raise_for_status()
            return ET.fromstring(resp.content)
        except (requests.RequestException, ET.ParseError) as e:
            if attempt < retries - 1:
                print(f"    Retry {attempt + 1}/{retries}: {e}")
                time.sleep(delay * (attempt + 1))
            else:
                print(f"    FAILED: {e}")
                return None


def get_survey_year_quarters(survey):
    """Get available year/quarter combinations for a survey."""
    root = fetch_xml("getSurveyYearList", {"survey": survey})
    if root is None:
        return []

    years = []
    for elem in root.findall(".//d:Year", NS):
        y = int(elem.text)
        if YEAR_MIN <= y <= YEAR_MAX:
            years.append(y)

    # For each year, get quarters
    combos = []
    for year in sorted(years):
        root = fetch_xml("getSurveyYearQuarterList",
                         {"survey": survey, "year": str(year)})
        if root is None:
            continue
        for elem in root.findall(".//d:Quarter", NS):
            combos.append((year, int(elem.text)))

    return combos


def parse_hh_xml(root):
    """Parse HH (haul) XML into a list of dicts."""
    rows = []
    for haul in root.findall(".//d:Cls_DatrasExchange_HH", NS):
        row = {}
        for child in haul:
            tag = child.tag.split("}")[-1]
            row[tag] = child.text.strip() if child.text else None
        rows.append(row)
    return rows


def parse_hl_xml(root):
    """Parse HL (catch) XML into a list of dicts."""
    rows = []
    for rec in root.findall(".//d:Cls_DatrasExchange_HL", NS):
        row = {}
        for child in rec:
            tag = child.tag.split("}")[-1]
            row[tag] = child.text.strip() if child.text else None
        rows.append(row)
    return rows


def download_survey_data(survey, year, quarter):
    """Download HH + HL data for a survey/year/quarter. Returns cached if available."""
    cache_hh = CACHE_DIR / f"{survey}_{year}_Q{quarter}_HH.csv"
    cache_hl = CACHE_DIR / f"{survey}_{year}_Q{quarter}_HL.csv"

    if cache_hh.exists() and cache_hl.exists():
        hh = pd.read_csv(cache_hh)
        hl = pd.read_csv(cache_hl)
        return hh, hl

    # Download HH
    root = fetch_xml("getHHdata",
                     {"survey": survey, "year": str(year), "quarter": str(quarter)})
    if root is None:
        return None, None

    hh_rows = parse_hh_xml(root)
    if not hh_rows:
        return None, None
    hh = pd.DataFrame(hh_rows)
    hh.to_csv(cache_hh, index=False)

    # Small delay between requests
    time.sleep(1)

    # Download HL
    root = fetch_xml("getHLdata",
                     {"survey": survey, "year": str(year), "quarter": str(quarter)})
    if root is None:
        return hh, None

    hl_rows = parse_hl_xml(root)
    if not hl_rows:
        hl = pd.DataFrame()
    else:
        hl = pd.DataFrame(hl_rows)
    hl.to_csv(cache_hl, index=False)

    time.sleep(1)
    return hh, hl


def main():
    print("=== Step 1: Download ICES DATRAS data ===\n")
    print(f"Period: {YEAR_MIN}-{YEAR_MAX}")
    print(f"Surveys: {len(SURVEYS)}\n")

    # Suppress SSL warnings (DATRAS certificate issue)
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    all_hauls = []
    all_catch = []
    survey_stats = []

    for s_idx, survey in enumerate(SURVEYS):
        print(f"\n[{s_idx + 1}/{len(SURVEYS)}] {survey}")

        # Get available year/quarter combinations
        combos = get_survey_year_quarters(survey)
        if not combos:
            print(f"  No data for {YEAR_MIN}-{YEAR_MAX}")
            continue

        print(f"  {len(combos)} year/quarter combinations")
        survey_haul_count = 0

        for year, quarter in combos:
            cache_tag = f"{survey}_{year}_Q{quarter}"
            hh, hl = download_survey_data(survey, year, quarter)

            if hh is None or hh.empty:
                continue
            if hl is None or hl.empty:
                continue

            # Filter valid hauls
            hh = hh[hh["HaulVal"].str.strip() == "V"].copy()
            if hh.empty:
                continue

            # Parse numeric columns
            for col in ["ShootLat", "ShootLong", "Depth", "HaulDur", "HaulNo", "Year"]:
                if col in hh.columns:
                    hh[col] = pd.to_numeric(hh[col], errors="coerce")

            # Create haul ID
            hh["haul_id"] = (
                survey + "_" + hh["Country"].str.strip() + "_"
                + hh["Ship"].str.strip() + "_" + hh["Year"].astype(int).astype(str) + "_"
                + hh["Quarter"].astype(str).str.strip() + "_"
                + hh["HaulNo"].astype(int).astype(str)
            )

            # Keep essential haul columns
            haul_cols = ["haul_id", "Survey", "Year", "Quarter", "Country",
                         "ShootLat", "ShootLong", "Depth", "HaulDur"]
            hh_clean = hh[[c for c in haul_cols if c in hh.columns]].copy()
            hh_clean = hh_clean.rename(columns={
                "ShootLat": "lat", "ShootLong": "lon",
                "Depth": "depth", "HaulDur": "haul_dur"
            })

            # Filter valid coordinates and depth
            hh_clean = hh_clean.dropna(subset=["lat", "lon"])
            hh_clean = hh_clean[
                (hh_clean["lat"] > 0) & (hh_clean["lat"] < 90)
                & (hh_clean["lon"] > -30) & (hh_clean["lon"] < 40)
            ]

            if hh_clean.empty:
                continue

            # Parse HL (catch) data
            for col in ["TotalNo", "HaulNo", "SubFactor"]:
                if col in hl.columns:
                    hl[col] = pd.to_numeric(hl[col], errors="coerce")

            # Valid_Aphia is the WoRMS species identifier
            if "Valid_Aphia" not in hl.columns:
                if "SpecCode" in hl.columns:
                    hl["Valid_Aphia"] = hl["SpecCode"]
                else:
                    continue

            hl["Valid_Aphia"] = pd.to_numeric(hl["Valid_Aphia"], errors="coerce")

            # Create matching haul ID in HL
            hl["Quarter_str"] = hl["Quarter"].astype(str).str.strip()
            hl["haul_id"] = (
                survey + "_" + hl["Country"].str.strip() + "_"
                + hl["Ship"].str.strip() + "_" + hl["Year"].astype(str).str.strip() + "_"
                + hl["Quarter_str"] + "_"
                + hl["HaulNo"].astype(int).astype(str)
            )

            # Filter to valid hauls
            valid_ids = set(hh_clean["haul_id"])
            hl = hl[hl["haul_id"].isin(valid_ids)]

            if hl.empty:
                continue

            # Aggregate: total abundance per species per haul
            # TotalNo is the raised number at each length class
            catch = (
                hl.groupby(["haul_id", "Valid_Aphia"])["TotalNo"]
                .sum()
                .reset_index()
            )
            catch.columns = ["haul_id", "aphia_id", "total_n"]

            # Calculate CPUE (number per hour)
            catch = catch.merge(
                hh_clean[["haul_id", "haul_dur"]].drop_duplicates("haul_id"),
                on="haul_id", how="left"
            )
            catch["cpue"] = catch["total_n"] / (catch["haul_dur"] / 60)
            catch = catch[catch["cpue"].notna() & (catch["cpue"] >= 0)]

            all_hauls.append(hh_clean)
            all_catch.append(catch[["haul_id", "aphia_id", "cpue"]])
            survey_haul_count += len(hh_clean)

            print(f"    {cache_tag}: {len(hh_clean)} hauls, "
                  f"{catch['aphia_id'].nunique()} species")

        if survey_haul_count > 0:
            survey_stats.append({"survey": survey, "hauls": survey_haul_count})
            print(f"  Total: {survey_haul_count} hauls")

    if not all_hauls:
        print("\nERROR: No data downloaded!")
        return

    # Combine all data
    print("\n\nCombining all surveys...")
    hauls = pd.concat(all_hauls, ignore_index=True)
    catch = pd.concat(all_catch, ignore_index=True)

    print(f"  Total hauls: {len(hauls):,}")
    print(f"  Total catch records: {len(catch):,}")
    print(f"  Total species (AphiaID): {catch['aphia_id'].nunique()}")

    # Remove duplicate haul IDs (in case of overlap)
    hauls = hauls.drop_duplicates("haul_id")

    # Mean CPUE per species per haul (in case of duplicates from subsampling)
    catch = catch.groupby(["haul_id", "aphia_id"])["cpue"].mean().reset_index()

    # Filter species: present in at least MIN_HAULS hauls
    species_haul_counts = catch[catch["cpue"] > 0].groupby("aphia_id")["haul_id"].nunique()
    common_species = species_haul_counts[species_haul_counts >= MIN_HAULS].index.tolist()
    catch = catch[catch["aphia_id"].isin(common_species)]
    print(f"  Species with >= {MIN_HAULS} hauls: {len(common_species)}")

    # Pivot to haul × species matrix (CPUE values)
    print("  Building abundance matrix...")
    abundance = catch.pivot_table(
        index="haul_id", columns="aphia_id", values="cpue", fill_value=0
    )
    # Rename columns to string for consistency
    abundance.columns = [f"sp_{int(c)}" for c in abundance.columns]
    print(f"  Abundance matrix: {abundance.shape[0]} hauls x {abundance.shape[1]} species")

    # Merge with haul metadata
    haul_info = hauls[["haul_id", "Survey", "Year", "Quarter", "Country",
                       "lat", "lon", "depth"]].drop_duplicates("haul_id")
    haul_info = haul_info.set_index("haul_id")

    result = abundance.join(haul_info, how="inner")
    print(f"  Matched hauls: {len(result):,}")

    # Save
    output_path = DATA_DIR / "datras_fish_hauls.csv"
    result.to_csv(output_path)
    print(f"\nSaved to {output_path}")

    # Save species list
    species_path = DATA_DIR / "datras_species_list.csv"
    species_df = pd.DataFrame({
        "aphia_id": common_species,
        "n_hauls": [species_haul_counts[s] for s in common_species],
    }).sort_values("n_hauls", ascending=False)
    species_df.to_csv(species_path, index=False)
    print(f"Species list saved to {species_path}")

    # Survey summary
    print(f"\n=== Summary ===")
    print(f"  Surveys with data: {len(survey_stats)}")
    for s in sorted(survey_stats, key=lambda x: x["hauls"], reverse=True):
        print(f"    {s['survey']:<12} {s['hauls']:>6} hauls")
    print(f"  Total hauls: {len(result):,}")
    print(f"  Species: {abundance.shape[1]}")
    print(f"  Period: {YEAR_MIN}-{YEAR_MAX}")
    print(f"  Lat range: {result['lat'].min():.1f} to {result['lat'].max():.1f}")
    print(f"  Lon range: {result['lon'].min():.1f} to {result['lon'].max():.1f}")
    print(f"  Depth range: {result['depth'].min():.0f}-{result['depth'].max():.0f} m")

    # Compare with Rutterford
    print(f"\n  Rutterford et al.: 23 surveys, 31,502 hauls, 198 species")
    print(f"  This download:     {len(survey_stats)} surveys, {len(result):,} hauls, "
          f"{abundance.shape[1]} species")

    print("\nDone.")


if __name__ == "__main__":
    main()
