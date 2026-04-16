"""
Step 4: Community analysis — PCoA on Bray-Curtis dissimilarity.

Following Rutterford et al. (2023), we perform Principal Coordinates Analysis (PCoA)
on the Bray-Curtis dissimilarity matrix of fish species abundance to identify the
primary axes of spatial variation in community structure.

We then test which environmental variables (SST, salinity, depth) are most strongly
associated with each PCoA axis using multiple linear regression.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

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

ENV_COLUMNS = ["sst_mean", "salinity_mean", "depth_midpoint"]


def pcoa(distance_matrix, n_components=3):
    """Classical multidimensional scaling (PCoA) from a distance matrix."""
    n = distance_matrix.shape[0]

    # Double centering
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (distance_matrix ** 2) @ H

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(B)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only positive eigenvalues
    positive = eigenvalues > 0
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]

    # Compute coordinates
    coords = eigenvectors[:, :n_components] * np.sqrt(eigenvalues[:n_components])

    # Proportion of variance explained
    total_var = eigenvalues.sum()
    var_explained = eigenvalues[:n_components] / total_var

    return coords, var_explained, eigenvalues


def aggregate_by_site(df, species_cols):
    """Aggregate transects by location to reduce noise.

    We group by rounded coordinates (0.5-degree grid cells) to create
    site-level abundance summaries, similar to Rutterford et al.'s
    1x1 degree grid cells.
    """
    df = df.copy()
    df["site_lat"] = (df["Decimal_latitude"] * 2).round() / 2
    df["site_lon"] = (df["Decimal_longitude"] * 2).round() / 2
    df["site_id"] = df["site_lat"].astype(str) + "_" + df["site_lon"].astype(str)

    # Mean abundance per site
    agg_species = df.groupby("site_id")[species_cols].mean()

    # Mean environment per site
    agg_env = df.groupby("site_id")[ENV_COLUMNS].mean()

    # Site coordinates
    agg_coords = df.groupby("site_id")[["site_lat", "site_lon"]].first()

    # Number of transects per site
    agg_n = df.groupby("site_id").size().rename("n_transects")

    result = pd.concat([agg_species, agg_env, agg_coords, agg_n], axis=1)
    return result


def main():
    print("=== Step 4: Community analysis (PCoA) ===\n")

    # Load joined data
    df = pd.read_csv(DATA_DIR / "fish_environment_joined.csv")
    species_cols = [c for c in SPECIES_COLUMNS if c in df.columns]
    print(f"Loaded {len(df)} transects, {len(species_cols)} species")

    # Aggregate by site
    sites = aggregate_by_site(df, species_cols)
    print(f"Aggregated to {len(sites)} sites (0.5-degree grid)")

    # Filter: require at least 3 transects per site
    sites = sites[sites["n_transects"] >= 3]
    print(f"After filtering (>=3 transects): {len(sites)} sites")

    if len(sites) < 10:
        print("WARNING: Very few sites. Results may not be robust.")

    # Compute Bray-Curtis dissimilarity
    print("\nComputing Bray-Curtis dissimilarity matrix...")
    abundance = sites[species_cols].values

    # Fourth-root transform to reduce influence of dominant species
    # (following Rutterford et al.'s approach)
    abundance_transformed = np.power(abundance + 1e-10, 0.25)

    bc_distances = pdist(abundance_transformed, metric="braycurtis")
    bc_matrix = squareform(bc_distances)
    print(f"  Distance matrix: {bc_matrix.shape}")

    # PCoA
    print("\nRunning PCoA...")
    coords, var_explained, eigenvalues = pcoa(bc_matrix, n_components=3)
    print(f"  Variance explained:")
    for i, ve in enumerate(var_explained):
        print(f"    PCoA axis {i+1}: {ve*100:.1f}%")
    print(f"  Total (3 axes): {var_explained.sum()*100:.1f}%")

    # Save PCoA results
    pcoa_df = pd.DataFrame(
        coords,
        columns=["PCoA1", "PCoA2", "PCoA3"],
        index=sites.index,
    )
    pcoa_df = pd.concat([pcoa_df, sites[ENV_COLUMNS + ["site_lat", "site_lon"]]], axis=1)
    pcoa_df.to_csv(RESULTS_DIR / "pcoa_axes.csv")

    # Correlations between PCoA axes and environmental variables
    print("\n=== Correlations: PCoA axes vs. environmental variables ===\n")
    print(f"{'Variable':<20} {'PCoA1':>10} {'PCoA2':>10} {'PCoA3':>10}")
    print("-" * 52)

    correlation_results = []
    for env_var in ENV_COLUMNS:
        row = {"variable": env_var}
        for i, axis in enumerate(["PCoA1", "PCoA2", "PCoA3"]):
            r, p = pearsonr(pcoa_df[axis], pcoa_df[env_var])
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            row[f"r_axis{i+1}"] = r
            row[f"p_axis{i+1}"] = p
            print(f"  {env_var if i == 0 else '':<20} r={r:+.3f} (p={p:.1e}) {sig}")
        correlation_results.append(row)
        print()

    corr_df = pd.DataFrame(correlation_results)
    corr_df.to_csv(RESULTS_DIR / "correlations.csv", index=False)

    # Identify primary driver for each axis
    print("=== Summary: Primary driver per PCoA axis ===\n")
    for i in range(3):
        axis_name = f"PCoA{i+1}"
        best_var = max(ENV_COLUMNS, key=lambda v: abs(
            pearsonr(pcoa_df[axis_name], pcoa_df[v])[0]
        ))
        r, p = pearsonr(pcoa_df[axis_name], pcoa_df[best_var])
        print(f"  {axis_name} ({var_explained[i]*100:.1f}% variance): "
              f"most associated with {best_var} (r={r:+.3f}, p={p:.1e})")

    # Plot: PCoA axes coloured by SST and salinity
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, (env_var, cmap, label) in zip(axes, [
        ("sst_mean", "RdYlBu_r", "SST (°C)"),
        ("salinity_mean", "YlGnBu", "Salinity (PSU)"),
        ("depth_midpoint", "viridis_r", "Depth midpoint (m)"),
    ]):
        sc = ax.scatter(
            pcoa_df["PCoA1"], pcoa_df["PCoA2"],
            c=pcoa_df[env_var], cmap=cmap, s=60, edgecolors="k", linewidth=0.3,
        )
        ax.set_xlabel(f"PCoA 1 ({var_explained[0]*100:.1f}%)")
        ax.set_ylabel(f"PCoA 2 ({var_explained[1]*100:.1f}%)")
        ax.set_title(f"Coloured by {label}")
        plt.colorbar(sc, ax=ax, label=label)

    fig.suptitle(
        "Mediterranean fish community structure (ClimateFish + Bio-ORACLE)\n"
        f"PCoA on Bray-Curtis dissimilarity — {len(sites)} sites, {len(species_cols)} species",
        fontsize=12,
    )
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "pcoa_community.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "pcoa_community.pdf", bbox_inches="tight")
    print(f"\nPlot saved to {RESULTS_DIR / 'pcoa_community.png'}")

    # Plot: Map of PCoA axis 1 scores
    fig, ax = plt.subplots(figsize=(10, 6))
    sc = ax.scatter(
        pcoa_df["site_lon"], pcoa_df["site_lat"],
        c=pcoa_df["PCoA1"], cmap="RdYlBu_r", s=80, edgecolors="k", linewidth=0.3,
    )
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title(f"PCoA Axis 1 scores ({var_explained[0]*100:.1f}% variance)")
    plt.colorbar(sc, ax=ax, label="PCoA1 score")
    plt.savefig(RESULTS_DIR / "pcoa_map.png", dpi=150, bbox_inches="tight")
    print(f"Map saved to {RESULTS_DIR / 'pcoa_map.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
