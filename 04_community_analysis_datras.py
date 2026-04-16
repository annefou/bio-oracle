"""
Step 4: Community analysis — PCoA + regression for DATRAS NE Atlantic data.

Reproduces Rutterford et al. (2023) analysis:
  - PCoA on Bray-Curtis dissimilarity of fourth-root transformed CPUE
  - Multiple linear regression of PCoA axes against SST, salinity, log(depth)
  - Univariate correlations to identify the primary driver per axis

Target results (Rutterford et al. 2023):
  PCoA Axis 1 (29% variance): SST       R² = 0.890
  PCoA Axis 2 (21% variance): Salinity   R² = 0.852
  PCoA Axis 3 (12% variance): Depth      R² = 0.460
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

ENV_COLUMNS = ["sst", "salinity", "log_depth"]
ENV_LABELS = {"sst": "SST (°C)", "salinity": "Salinity (PSU)", "log_depth": "log₁₀(Depth)"}


def pcoa(distance_matrix, n_components=3):
    """Classical multidimensional scaling (PCoA)."""
    n = distance_matrix.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (distance_matrix ** 2) @ H

    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    positive = eigenvalues > 0
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]

    coords = eigenvectors[:, :n_components] * np.sqrt(eigenvalues[:n_components])
    total_var = eigenvalues.sum()
    var_explained = eigenvalues[:n_components] / total_var

    return coords, var_explained


def main():
    print("=== Step 4: Community analysis (DATRAS — NE Atlantic) ===\n")

    # Use filtered data (fish only, shelf depth) if available
    filtered_path = DATA_DIR / "datras_grid_cells_filtered.csv"
    default_path = DATA_DIR / "datras_grid_cells.csv"
    data_path = filtered_path if filtered_path.exists() else default_path
    print(f"Using: {data_path}")

    grid = pd.read_csv(data_path, index_col=0)
    species_cols = [c for c in grid.columns if c.startswith("sp_")]
    print(f"Grid cells: {len(grid)}, Species: {len(species_cols)}")

    # Fourth-root transform (reduces influence of dominant species)
    abundance = grid[species_cols].values
    abundance_transformed = np.power(abundance + 1e-10, 0.25)

    # Bray-Curtis dissimilarity
    print("\nComputing Bray-Curtis dissimilarity...")
    bc = squareform(pdist(abundance_transformed, metric="braycurtis"))
    print(f"  Matrix: {bc.shape}")

    # PCoA
    print("Running PCoA...")
    coords, var_explained = pcoa(bc, n_components=3)
    print(f"  Variance explained:")
    for i, ve in enumerate(var_explained):
        print(f"    PCoA {i + 1}: {ve * 100:.1f}%")
    print(f"  Total (3 axes): {var_explained.sum() * 100:.1f}%")

    # Build results dataframe
    pcoa_df = pd.DataFrame(
        coords, columns=["PCoA1", "PCoA2", "PCoA3"], index=grid.index
    )
    pcoa_df = pd.concat([
        pcoa_df,
        grid[ENV_COLUMNS + ["grid_lat", "grid_lon", "n_hauls", "n_surveys"]]
    ], axis=1)
    pcoa_df.to_csv(RESULTS_DIR / "datras_pcoa_axes.csv")

    # === Regression analysis ===
    print("\n=== Regression: PCoA axes vs environment ===\n")

    all_results = []
    for axis in ["PCoA1", "PCoA2", "PCoA3"]:
        y = pcoa_df[axis].values
        X = pcoa_df[ENV_COLUMNS].values
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        model = LinearRegression().fit(X_scaled, y)
        r2_full = model.score(X_scaled, y)

        print(f"--- {axis} ---")
        print(f"  Full model R²: {r2_full:.3f}\n")
        print(f"  {'Variable':<20} {'r':>8} {'R²':>8} {'p-value':>12} {'Std.β':>10}")
        print(f"  {'-' * 58}")

        axis_result = {"axis": axis, "r2_full": r2_full, "vars": {}}
        for i, col in enumerate(ENV_COLUMNS):
            r, p = pearsonr(pcoa_df[axis], pcoa_df[col])
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {ENV_LABELS[col]:<20} {r:+.3f}   {r ** 2:.3f}  {p:>11.2e}  "
                  f"{model.coef_[i]:+.3f} {sig}")
            axis_result["vars"][col] = {"r": r, "r2": r ** 2, "p": p, "beta": model.coef_[i]}

        primary = max(ENV_COLUMNS, key=lambda v: abs(axis_result["vars"][v]["r"]))
        print(f"\n  -> Primary driver: {ENV_LABELS[primary]} "
              f"(r={axis_result['vars'][primary]['r']:+.3f})\n")
        all_results.append(axis_result)

    # === Comparison with Rutterford et al. (2023) ===
    print("=" * 60)
    print("=== COMPARISON: Rutterford et al. (2023) vs. this reproduction ===")
    print("=" * 60)

    print(f"\nRutterford et al. (NE Atlantic, 198 species, 193 grid cells, 23 surveys):")
    print(f"  PCoA 1 (29% var): SST        R² = 0.890")
    print(f"  PCoA 2 (21% var): Salinity    R² = 0.852")
    print(f"  PCoA 3 (12% var): Depth       R² = 0.460")

    print(f"\nThis reproduction (NE Atlantic, {len(species_cols)} species, "
          f"{len(grid)} grid cells):")
    for i, res in enumerate(all_results):
        primary = max(ENV_COLUMNS, key=lambda v: abs(res["vars"][v]["r"]))
        r2 = res["vars"][primary]["r2"]
        print(f"  PCoA {i + 1} ({var_explained[i] * 100:.0f}% var): "
              f"{ENV_LABELS[primary]:<15} R² = {r2:.3f}")

    # Reproduction outcome
    axis1_sst_r2 = all_results[0]["vars"]["sst"]["r2"]
    axis1_primary = max(ENV_COLUMNS, key=lambda v: abs(all_results[0]["vars"][v]["r"]))
    axis1_p = all_results[0]["vars"]["sst"]["p"]

    print(f"\n--- Reproduction assessment ---")
    print(f"  SST on Axis 1: R² = {axis1_sst_r2:.3f} (Rutterford: 0.890)")
    print(f"  Primary driver of Axis 1: {ENV_LABELS[axis1_primary]}")

    if axis1_primary == "sst" and axis1_p < 0.001:
        print(f"  REPRODUCTION: SUCCESSFUL — SST is the primary driver of Axis 1")
        if abs(axis1_sst_r2 - 0.890) < 0.1:
            print(f"  Effect size closely matches (within 0.1 R²)")
        else:
            print(f"  Effect size differs: {axis1_sst_r2:.3f} vs 0.890")
    elif axis1_primary == "sst":
        print(f"  REPRODUCTION: PARTIAL — SST is primary but weaker than reported")
    else:
        print(f"  REPRODUCTION: FAILED — {ENV_LABELS[axis1_primary]} dominates, not SST")

    # === Save regression table ===
    rows = []
    for res in all_results:
        for col in ENV_COLUMNS:
            v = res["vars"][col]
            rows.append({
                "axis": res["axis"], "variable": ENV_LABELS[col],
                "r": v["r"], "r2": v["r2"], "p_value": v["p"], "std_beta": v["beta"],
            })
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "datras_regression_table.csv", index=False)
    print(f"\nRegression table: {RESULTS_DIR / 'datras_regression_table.csv'}")

    # === Plots ===

    # 1. PCoA ordination coloured by environment
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (env, cmap, label) in zip(axes, [
        ("sst", "RdYlBu_r", "SST (°C)"),
        ("salinity", "YlGnBu", "Salinity (PSU)"),
        ("log_depth", "viridis_r", "log₁₀(Depth)"),
    ]):
        sc = ax.scatter(pcoa_df["PCoA1"], pcoa_df["PCoA2"],
                        c=pcoa_df[env], cmap=cmap, s=40, edgecolors="k", linewidth=0.3,
                        alpha=0.7)
        ax.set_xlabel(f"PCoA 1 ({var_explained[0] * 100:.1f}%)")
        ax.set_ylabel(f"PCoA 2 ({var_explained[1] * 100:.1f}%)")
        ax.set_title(f"Coloured by {label}")
        plt.colorbar(sc, ax=ax, label=label)
    fig.suptitle(
        f"NE Atlantic fish community structure (DATRAS + Bio-ORACLE)\n"
        f"PCoA on Bray-Curtis — {len(grid)} grid cells, {len(species_cols)} species\n"
        f"Reproducing Rutterford et al. (2023)",
        fontsize=12)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "datras_pcoa.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "datras_pcoa.pdf", bbox_inches="tight")
    print(f"\nPlot: {RESULTS_DIR / 'datras_pcoa.png'}")

    # 2. Variable importance comparison with Rutterford
    fig, axes_plot = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: Our R² per axis
    ax = axes_plot[0]
    x = np.arange(3)
    width = 0.25
    for i, col in enumerate(ENV_COLUMNS):
        r2_vals = [res["vars"][col]["r2"] for res in all_results]
        ax.bar(x + i * width, r2_vals, width, label=ENV_LABELS[col])
    ax.set_xticks(x + width)
    ax.set_xticklabels(["PCoA 1", "PCoA 2", "PCoA 3"])
    ax.set_ylabel("R² (univariate)")
    ax.set_title("This reproduction:\nVariable importance per axis")
    ax.legend()
    ax.set_ylim(0, 1)

    # Panel 2: Rutterford's R² per axis
    ax = axes_plot[1]
    rutterford_r2 = {
        "PCoA1": {"sst": 0.890, "salinity": 0.530, "log_depth": 0.170},
        "PCoA2": {"sst": 0.500, "salinity": 0.852, "log_depth": 0.090},
        "PCoA3": {"sst": 0.100, "salinity": 0.050, "log_depth": 0.460},
    }
    for i, col in enumerate(ENV_COLUMNS):
        r2_vals = [rutterford_r2[f"PCoA{j + 1}"][col] for j in range(3)]
        ax.bar(x + i * width, r2_vals, width, label=ENV_LABELS[col])
    ax.set_xticks(x + width)
    ax.set_xticklabels(["PCoA 1", "PCoA 2", "PCoA 3"])
    ax.set_ylabel("R² (univariate)")
    ax.set_title("Rutterford et al. (2023):\nVariable importance per axis")
    ax.legend()
    ax.set_ylim(0, 1)

    # Panel 3: Side-by-side Axis 1 comparison
    ax = axes_plot[2]
    rutterford_axis1 = [0.890, 0.530, 0.170]
    this_study_axis1 = [all_results[0]["vars"][c]["r2"] for c in ENV_COLUMNS]
    labels = [ENV_LABELS[c] for c in ENV_COLUMNS]
    x2 = np.arange(len(labels))
    ax.bar(x2 - 0.15, rutterford_axis1, 0.3,
           label="Rutterford et al.", color="steelblue")
    ax.bar(x2 + 0.15, this_study_axis1, 0.3,
           label="This reproduction", color="coral")
    ax.set_xticks(x2)
    ax.set_xticklabels(labels)
    ax.set_ylabel("R² (PCoA Axis 1)")
    ax.set_title("Axis 1 comparison:\nOriginal vs reproduction")
    ax.legend()
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "datras_variable_importance.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "datras_variable_importance.pdf", bbox_inches="tight")
    print(f"Plot: {RESULTS_DIR / 'datras_variable_importance.png'}")

    # 3. Map of PCoA1 scores
    fig, ax = plt.subplots(figsize=(14, 8))
    sc = ax.scatter(pcoa_df["grid_lon"], pcoa_df["grid_lat"],
                    c=pcoa_df["PCoA1"], cmap="RdYlBu_r", s=40, edgecolors="k",
                    linewidth=0.3, alpha=0.8)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title(f"PCoA Axis 1 ({var_explained[0] * 100:.1f}% variance) — "
                 f"NE Atlantic DATRAS surveys\n"
                 f"Reproducing Rutterford et al. (2023)")
    plt.colorbar(sc, ax=ax, label="PCoA1 score")
    ax.set_aspect("equal")
    plt.savefig(RESULTS_DIR / "datras_pcoa_map.png", dpi=150, bbox_inches="tight")
    print(f"Map: {RESULTS_DIR / 'datras_pcoa_map.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
