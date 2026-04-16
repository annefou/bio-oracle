"""
Step 5: Variable importance — regression analysis of PCoA axes vs environment.

Following Rutterford et al. (2023), we use multiple linear regression to test
which environmental variables (SST, salinity, depth) explain the most variation
in fish community structure (PCoA axes).

Key question: Does SST emerge as the primary driver of fish community structure
in the Mediterranean, as Rutterford et al. (2023) found for the NE Atlantic?
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

RESULTS_DIR = Path("results")

ENV_COLUMNS = ["sst_mean", "salinity_mean", "depth_midpoint"]
ENV_LABELS = {
    "sst_mean": "SST (°C)",
    "salinity_mean": "Salinity (PSU)",
    "depth_midpoint": "Depth (m)",
}


def regression_analysis(pcoa_df, axis_name, env_cols):
    """Run multiple linear regression: PCoA axis ~ environmental variables."""
    y = pcoa_df[axis_name].values
    X = pcoa_df[env_cols].values

    # Standardise predictors for comparable coefficients
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Full model
    model = LinearRegression()
    model.fit(X_scaled, y)
    r2 = model.score(X_scaled, y)

    # Per-variable R² (univariate)
    univariate_r2 = {}
    for i, col in enumerate(env_cols):
        r, p = pearsonr(pcoa_df[axis_name], pcoa_df[col])
        univariate_r2[col] = {"r": r, "r2": r ** 2, "p": p}

    # Standardised coefficients
    coefficients = {}
    for i, col in enumerate(env_cols):
        coefficients[col] = {
            "beta": model.coef_[i],
            "abs_beta": abs(model.coef_[i]),
        }

    return {
        "r2_full": r2,
        "univariate": univariate_r2,
        "coefficients": coefficients,
    }


def main():
    print("=== Step 5: Variable importance analysis ===\n")

    # Load PCoA results
    pcoa_df = pd.read_csv(RESULTS_DIR / "pcoa_axes.csv", index_col=0)
    print(f"Loaded {len(pcoa_df)} sites with PCoA scores\n")

    # Run regression for each axis
    all_results = []
    for axis in ["PCoA1", "PCoA2", "PCoA3"]:
        result = regression_analysis(pcoa_df, axis, ENV_COLUMNS)

        print(f"--- {axis} ---")
        print(f"  Full model R²: {result['r2_full']:.3f}")
        print()

        print(f"  {'Variable':<20} {'r':>8} {'R²':>8} {'p-value':>12} {'Std.beta':>10}")
        print(f"  {'-'*58}")
        for col in ENV_COLUMNS:
            uni = result["univariate"][col]
            coef = result["coefficients"][col]
            sig = "***" if uni["p"] < 0.001 else "**" if uni["p"] < 0.01 else "*" if uni["p"] < 0.05 else ""
            print(
                f"  {ENV_LABELS[col]:<20} {uni['r']:+.3f}   {uni['r2']:.3f}  "
                f"{uni['p']:>11.2e}  {coef['beta']:+.3f} {sig}"
            )

        # Primary driver
        primary = max(
            ENV_COLUMNS,
            key=lambda v: abs(result["univariate"][v]["r"]),
        )
        print(f"\n  -> Primary driver: {ENV_LABELS[primary]} "
              f"(r={result['univariate'][primary]['r']:+.3f})\n")

        all_results.append({"axis": axis, **result})

    # Comparison with Rutterford et al. (2023)
    print("\n=== Comparison with Rutterford et al. (2023) ===\n")
    print("Rutterford et al. (NE Atlantic, 198 species, trawl surveys):")
    print("  PCoA Axis 1 (29% variance): SST       (R² = 0.890)")
    print("  PCoA Axis 2 (21% variance): Salinity   (R² = 0.852)")
    print("  PCoA Axis 3 (12% variance): Depth      (R² = 0.460)")
    print()
    print("This study (Mediterranean, 15 species, visual census):")

    for res in all_results:
        axis = res["axis"]
        primary = max(
            ENV_COLUMNS,
            key=lambda v: abs(res["univariate"][v]["r"]),
        )
        r2 = res["univariate"][primary]["r2"]
        print(f"  {axis}: {ENV_LABELS[primary]:<15} (R² = {r2:.3f})")

    # Determine replication outcome
    axis1_primary = max(
        ENV_COLUMNS,
        key=lambda v: abs(all_results[0]["univariate"][v]["r"]),
    )
    if axis1_primary == "sst_mean":
        r2_sst = all_results[0]["univariate"]["sst_mean"]["r2"]
        print(f"\n  REPLICATION OUTCOME: CONFIRMED")
        print(f"  SST is the primary driver of community variation (Axis 1)")
        print(f"  in the Mediterranean, consistent with Rutterford et al. (2023).")
    else:
        print(f"\n  REPLICATION OUTCOME: NOT CONFIRMED")
        print(f"  {ENV_LABELS[axis1_primary]} (not SST) is the primary driver")
        print(f"  of community variation (Axis 1) in the Mediterranean.")

    # Save regression table
    rows = []
    for res in all_results:
        for col in ENV_COLUMNS:
            uni = res["univariate"][col]
            coef = res["coefficients"][col]
            rows.append({
                "axis": res["axis"],
                "variable": ENV_LABELS[col],
                "r": uni["r"],
                "r2": uni["r2"],
                "p_value": uni["p"],
                "std_beta": coef["beta"],
            })
    reg_df = pd.DataFrame(rows)
    reg_df.to_csv(RESULTS_DIR / "regression_table.csv", index=False)
    print(f"\nRegression table saved to {RESULTS_DIR / 'regression_table.csv'}")

    # Plot: variable importance comparison
    fig, axes_plot = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: Bar chart of R² per axis
    ax = axes_plot[0]
    x = np.arange(3)
    width = 0.25
    for i, col in enumerate(ENV_COLUMNS):
        r2_values = [res["univariate"][col]["r2"] for res in all_results]
        ax.bar(x + i * width, r2_values, width, label=ENV_LABELS[col])
    ax.set_xticks(x + width)
    ax.set_xticklabels(["PCoA 1", "PCoA 2", "PCoA 3"])
    ax.set_ylabel("R² (univariate)")
    ax.set_title("Variable importance per PCoA axis\n(Mediterranean — this study)")
    ax.legend()
    ax.set_ylim(0, 1)

    # Panel 2: Comparison with Rutterford et al.
    ax = axes_plot[1]
    rutterford_r2 = {
        "PCoA 1": {"SST (°C)": 0.890, "Salinity (PSU)": 0.852, "Depth (m)": 0.460},
    }
    # Just show Axis 1 comparison
    this_r2_axis1 = {
        ENV_LABELS[col]: all_results[0]["univariate"][col]["r2"] for col in ENV_COLUMNS
    }
    labels = list(ENV_LABELS.values())
    rutterford_vals = [0.890, 0.852, 0.460]
    this_vals = [this_r2_axis1[l] for l in labels]

    x = np.arange(len(labels))
    ax.bar(x - 0.15, rutterford_vals, 0.3, label="Rutterford et al. (NE Atlantic)", color="steelblue")
    ax.bar(x + 0.15, this_vals, 0.3, label="This study (Mediterranean)", color="coral")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("R² (univariate, PCoA Axis 1)")
    ax.set_title("PCoA Axis 1: NE Atlantic vs. Mediterranean")
    ax.legend()
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "variable_importance.png", dpi=150, bbox_inches="tight")
    plt.savefig(RESULTS_DIR / "variable_importance.pdf", bbox_inches="tight")
    print(f"Plot saved to {RESULTS_DIR / 'variable_importance.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
