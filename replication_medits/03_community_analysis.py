"""
Step 3: Community analysis and variable importance for MEDITS data.

PCoA on Bray-Curtis dissimilarity + multiple linear regression to identify
the primary environmental drivers of fish community structure, following
Rutterford et al. (2023).
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

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
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
    print("=== Step 3: Community analysis (MEDITS shelf ≤200m) ===\n")

    grid = pd.read_csv(DATA_DIR / "medits_grid_cells_shelf.csv", index_col=0)
    species_cols = [c for c in grid.columns
                    if c not in ["sst", "salinity", "depth", "log_depth",
                                 "grid_lat", "grid_lon", "n_hauls"]]
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
        print(f"    PCoA {i+1}: {ve*100:.1f}%")
    print(f"  Total (3 axes): {var_explained.sum()*100:.1f}%")

    # Build results dataframe
    pcoa_df = pd.DataFrame(
        coords, columns=["PCoA1", "PCoA2", "PCoA3"], index=grid.index
    )
    pcoa_df = pd.concat([pcoa_df, grid[ENV_COLUMNS + ["grid_lat", "grid_lon", "n_hauls"]]], axis=1)
    pcoa_df.to_csv(RESULTS_DIR / "medits_pcoa_axes.csv")

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
        print(f"  {'-'*58}")

        axis_result = {"axis": axis, "r2_full": r2_full, "vars": {}}
        for i, col in enumerate(ENV_COLUMNS):
            r, p = pearsonr(pcoa_df[axis], pcoa_df[col])
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {ENV_LABELS[col]:<20} {r:+.3f}   {r**2:.3f}  {p:>11.2e}  {model.coef_[i]:+.3f} {sig}")
            axis_result["vars"][col] = {"r": r, "r2": r**2, "p": p, "beta": model.coef_[i]}

        primary = max(ENV_COLUMNS, key=lambda v: abs(axis_result["vars"][v]["r"]))
        print(f"\n  -> Primary: {ENV_LABELS[primary]} (r={axis_result['vars'][primary]['r']:+.3f})\n")
        all_results.append(axis_result)

    # === Comparison with Rutterford et al. ===
    print("=== Comparison with Rutterford et al. (2023) ===\n")
    print("Rutterford et al. (NE Atlantic, 198 species, 193 grid cells):")
    print("  PCoA 1 (29% var): SST        R²=0.890")
    print("  PCoA 2 (21% var): Salinity    R²=0.852")
    print("  PCoA 3 (12% var): Depth       R²=0.460\n")
    print(f"This study (Mediterranean, {len(species_cols)} species, {len(grid)} grid cells):")
    for i, res in enumerate(all_results):
        primary = max(ENV_COLUMNS, key=lambda v: abs(res["vars"][v]["r"]))
        r2 = res["vars"][primary]["r2"]
        print(f"  PCoA {i+1} ({var_explained[i]*100:.0f}% var): {ENV_LABELS[primary]:<15} R²={r2:.3f}")

    # Replication outcome
    axis1_primary = max(ENV_COLUMNS, key=lambda v: abs(all_results[0]["vars"][v]["r"]))
    r_sst = all_results[0]["vars"]["sst"]["r"]
    p_sst = all_results[0]["vars"]["sst"]["p"]
    print(f"\nPCoA Axis 1 — SST correlation: r={r_sst:+.3f}, p={p_sst:.2e}")
    if axis1_primary == "sst" and p_sst < 0.05:
        print("REPLICATION OUTCOME: CONFIRMED")
    elif axis1_primary != "sst":
        alt = ENV_LABELS[axis1_primary]
        print(f"REPLICATION OUTCOME: NOT CONFIRMED ({alt} is the primary driver, not SST)")
    else:
        print("REPLICATION OUTCOME: INCONCLUSIVE (SST not significant)")

    # === Save regression table ===
    rows = []
    for res in all_results:
        for col in ENV_COLUMNS:
            v = res["vars"][col]
            rows.append({
                "axis": res["axis"], "variable": ENV_LABELS[col],
                "r": v["r"], "r2": v["r2"], "p_value": v["p"], "std_beta": v["beta"],
            })
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "medits_regression_table.csv", index=False)

    # === Plots ===

    # 1. PCoA ordination coloured by environment
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (env, cmap, label) in zip(axes, [
        ("sst", "RdYlBu_r", "SST (°C)"),
        ("salinity", "YlGnBu", "Salinity (PSU)"),
        ("log_depth", "viridis_r", "log₁₀(Depth)"),
    ]):
        sc = ax.scatter(pcoa_df["PCoA1"], pcoa_df["PCoA2"],
                        c=pcoa_df[env], cmap=cmap, s=60, edgecolors="k", linewidth=0.3)
        ax.set_xlabel(f"PCoA 1 ({var_explained[0]*100:.1f}%)")
        ax.set_ylabel(f"PCoA 2 ({var_explained[1]*100:.1f}%)")
        ax.set_title(f"Coloured by {label}")
        plt.colorbar(sc, ax=ax, label=label)
    fig.suptitle(
        f"Mediterranean fish community structure (MEDITS + Bio-ORACLE)\n"
        f"PCoA on Bray-Curtis — {len(grid)} grid cells, {len(species_cols)} species",
        fontsize=12)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "medits_pcoa.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "medits_pcoa.pdf", bbox_inches="tight")
    print(f"\nPlot: {RESULTS_DIR / 'medits_pcoa.png'}")

    # 2. Variable importance comparison
    fig, axes_plot = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes_plot[0]
    x = np.arange(3)
    width = 0.25
    for i, col in enumerate(ENV_COLUMNS):
        r2_vals = [res["vars"][col]["r2"] for res in all_results]
        ax.bar(x + i * width, r2_vals, width, label=ENV_LABELS[col])
    ax.set_xticks(x + width)
    ax.set_xticklabels(["PCoA 1", "PCoA 2", "PCoA 3"])
    ax.set_ylabel("R² (univariate)")
    ax.set_title("Variable importance per PCoA axis\n(Mediterranean MEDITS)")
    ax.legend()
    ax.set_ylim(0, 1)

    ax = axes_plot[1]
    rutterford = [0.890, 0.852, 0.460]
    this_study = [all_results[0]["vars"][c]["r2"] for c in ENV_COLUMNS]
    labels = [ENV_LABELS[c] for c in ENV_COLUMNS]
    x = np.arange(len(labels))
    ax.bar(x - 0.15, rutterford, 0.3, label="Rutterford (NE Atlantic)", color="steelblue")
    ax.bar(x + 0.15, this_study, 0.3, label="This study (Mediterranean)", color="coral")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("R² (PCoA Axis 1)")
    ax.set_title("Axis 1: NE Atlantic vs Mediterranean")
    ax.legend()
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "medits_variable_importance.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "medits_variable_importance.pdf", bbox_inches="tight")
    print(f"Plot: {RESULTS_DIR / 'medits_variable_importance.png'}")

    # 3. Map of PCoA1 scores
    fig, ax = plt.subplots(figsize=(12, 6))
    sc = ax.scatter(pcoa_df["grid_lon"], pcoa_df["grid_lat"],
                    c=pcoa_df["PCoA1"], cmap="RdYlBu_r", s=80, edgecolors="k", linewidth=0.3)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title(f"PCoA Axis 1 ({var_explained[0]*100:.1f}% variance)")
    plt.colorbar(sc, ax=ax, label="PCoA1 score")
    plt.savefig(RESULTS_DIR / "medits_pcoa_map.png", dpi=150, bbox_inches="tight")
    print(f"Map: {RESULTS_DIR / 'medits_pcoa_map.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
