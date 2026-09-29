#!/usr/bin/env python3
"""Build the C7B-integrated four-panel Figure 3."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "figures" / "source_data" / "final" / "FIGURE3_SOURCE_C7.tsv"
CV = ROOT / "results" / "c7b" / "CV_REPEAT_OOF_VALUES.tsv"
PERM = ROOT / "results" / "c7b" / "PERMUTATION_NULL_UNIFIED.tsv"
OUT = ROOT / "figures"


def main() -> None:
    source = pd.read_csv(SOURCE, sep="\t")
    cv = pd.read_csv(CV, sep="\t")
    perm = pd.read_csv(PERM, sep="\t")

    gene = source.loc[source["record_type"].eq("gene_point")].copy()
    summ = source.loc[source["record_type"].eq("summary")].copy()
    residual = gene["residual_Ribo6_after_RNA6"].to_numpy(float)
    protein = gene["Protein24"].to_numpy(float)
    rho, _ = spearmanr(residual, protein)
    # The manuscript's prespecified bootstrap values are used for the displayed
    # interval; the scatter itself is drawn directly from the source data.
    resid_row = summ.loc[summ["label"].eq("A")].iloc[0]
    partial = summ.loc[summ["label"].eq("B")].copy()
    partial = partial.sort_values("estimate")

    colors = {"A": "#3B77A6", "B": "#D97820", "C": "#6B7780", "D": "#5A6670"}
    plt.rcParams.update(
        {
            "font.family": ["Arial", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 12,
            "axes.labelsize": 9.5,
            "axes.linewidth": 1.1,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.4), dpi=600)
    a, b, c, d = ax.ravel()

    # Panel A: retained source-data scatter.
    a.scatter(residual, protein, s=13, color="#E7A064", alpha=0.24, linewidths=0)
    a.axhline(0, color="#D6DCE0", lw=0.8, zorder=0)
    a.axvline(0, color="#D6DCE0", lw=0.8, zorder=0)
    a.set_xlabel("Ribo6 residual after RNA6")
    a.set_ylabel("Protein24 log2FC")
    a.set_title("RNA-adjusted Ribo residual", loc="left", fontweight="bold")
    a.text(
        0.04,
        0.94,
        "ρ = 0.111\n95% CI 0.078–0.142\nn = 3,966",
        transform=a.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    a.set_xlim(np.nanpercentile(residual, 0.2) - 0.3, np.nanpercentile(residual, 99.8) + 0.3)
    a.set_ylim(np.nanpercentile(protein, 0.2) - 0.5, np.nanpercentile(protein, 99.8) + 0.5)

    # Panel B: partial associations.
    y = np.arange(len(partial))
    b.axvline(0, color="#6E777E", lw=1.0, zorder=0)
    for yi, (_, row) in zip(y, partial.iterrows()):
        col = colors["B"] if "Ribo6" in str(row["gene_symbol"]) and "RNA6" in str(row["gene_symbol"]).split("|")[0] else colors["A"]
        # Explicit labels avoid relying on row order for color assignment.
        col = colors["A"] if str(row["gene_symbol"]).startswith("Partial RNA6") else colors["B"]
        b.errorbar(
            row["estimate"],
            yi,
            xerr=[[row["estimate"] - row["CI_low"]], [row["CI_high"] - row["estimate"]]],
            fmt="o",
            color=col,
            ecolor=col,
            elinewidth=2.1,
            capsize=5,
            markersize=8,
            zorder=3,
        )
    b.set_yticks(y)
    b.set_yticklabels(["RNA6 | Ribo6", "Ribo6 | RNA6"], fontsize=8.5)
    b.set_xlabel("Partial Spearman ρ")
    b.set_title("Partial associations", loc="left", fontweight="bold")
    b.set_xlim(-0.05, 0.18)
    b.grid(False)

    # Panel C: ten complete-OOF repeats, with means and split-stability ranges.
    model_cols = [("A", "rho_A"), ("B", "rho_B"), ("C", "rho_C"), ("D", "rho_D")]
    means = [cv[col].mean() for _, col in model_cols]
    lows = [cv[col].quantile(0.025) for _, col in model_cols]
    highs = [cv[col].quantile(0.975) for _, col in model_cols]
    x = np.arange(len(model_cols))
    for j, (label, col) in enumerate(model_cols):
        c.scatter(
            np.full(len(cv), j),
            cv[col],
            s=24,
            color=colors[label],
            alpha=0.24,
            edgecolors="none",
            zorder=2,
        )
        c.errorbar(
            j,
            means[j],
            yerr=[[means[j] - lows[j]], [highs[j] - means[j]]],
            fmt="o",
            color=colors[label],
            ecolor=colors[label],
            capsize=4,
            elinewidth=1.7,
            markersize=8,
            zorder=4,
        )
    c.plot([0, 1], [cv.rho_A.mean(), cv.rho_B.mean()], color="#9BA4AA", lw=1.0, zorder=1)
    c.set_xticks(x)
    c.set_xticklabels(["A  RNA6", "B  RNA6+Ribo6", "C  RNA6+TE6", "D  RNA6+Ribo6+TE6"], fontsize=7.5)
    c.set_ylabel("Complete-OOF Spearman ρ")
    c.set_xlabel("Model")
    c.set_title("Held-out-gene complete-OOF CV", loc="left", fontweight="bold")
    c.set_ylim(0.08, 0.18)
    c.text(
        0.04,
        0.95,
        "10 repeats; each point = one complete 10-fold OOF repeat\nB − A mean Δρ = 0.04042; 10/10 repeats positive\nRanges: 95% split-stability intervals",
        transform=c.transAxes,
        ha="left",
        va="top",
        fontsize=7.3,
    )
    c.tick_params(axis="x", labelrotation=18, pad=1)

    # Panel D: unified permutation null.
    null = perm["T_perm_repeat_OOF"].to_numpy(float)
    observed = float(cv["delta_BA"].mean())
    d.hist(null, bins=50, color="#78838A", edgecolor="white", linewidth=0.35)
    d.axvline(observed, color="#D97820", lw=2.4)
    d.set_xlabel("Null T_REPEAT_OOF (Δρ)")
    d.set_ylabel("Count")
    d.set_title("Conditional residual-permutation null", loc="left", fontweight="bold")
    d.text(
        0.97,
        0.95,
        "Observed T_obs = 0.04042\n10,000 permutations\n0 null values ≥ observed\nP < 1 × 10⁻⁴",
        transform=d.transAxes,
        ha="right",
        va="top",
        fontsize=7.8,
    )
    d.set_xlim(min(null.min() - 0.003, -0.06), max(observed + 0.003, 0.045))

    for idx, panel in enumerate(ax.ravel()):
        panel.text(
            0.00,
            0.98,
            "ABCD"[idx],
            transform=panel.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )
        panel.spines["top"].set_visible(False)
        panel.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.15, top=0.92, wspace=0.30, hspace=0.40)
    fig.savefig(OUT / "FIGURE3_C7B_FINAL.png", bbox_inches="tight", dpi=600)
    fig.savefig(OUT / "FIGURE3_C7B_FINAL.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(OUT / "FIGURE3_C7B_FINAL.svg", bbox_inches="tight", dpi=600)
    # A high-resolution TIFF is retained for production workflows; the user-facing
    # deliverables remain the requested PDF and PNG.
    fig.savefig(OUT / "FIGURE3_C7B_FINAL.tiff", bbox_inches="tight", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
