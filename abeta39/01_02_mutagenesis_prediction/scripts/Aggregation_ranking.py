#!/usr/bin/env python3
"""
Aβ39 Variant Aggregation Propensity Analysis
============================================

Consensus analysis of Aβ39 peptide variants across four aggregation predictors
(TANGO, PASTA, AmyPred-FRL, CrossBeta): per-tool z-score normalisation, a
consensus score, identification of top aggregators/disruptors, and figures.

This is the loader-refactored version. The four per-tool readers previously
duplicated the same `for _, row in df.iterrows()` loop; they now share one
vectorised loader. At 65–465 variants this is not a speed fix (the loops cost
microseconds) — it removes ~40 lines of duplicated branching and the WT-skip
logic, which is the readability point the reviewers raised. The plotting code
below still uses iterrows for per-bar text annotation, which is the idiomatic
use and is left unchanged.

Author: Sergey Ilin
Date: May 2026
Version: 1.1
"""

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

# ========================== PLOT STYLE CONFIGURATION ==========================

matplotlib.rcParams["font.family"] = "monospace"

# Color palette (dark theme)
BG = "#0a0c14"
SURFACE = "#111520"
MUTED = "#64748b"
TEXT = "#e2e8f0"
WT_COL = "#4ff7c0"
AGG_COL = "#ef4444"  # Increased aggregation
DIS_COL = "#3b82f6"  # Decreased aggregation
ACC_COL = "#fbbf24"  # Highlight / consensus

TOOL_COLORS = {
    "TANGO": "#ef4444",
    "PASTA": "#3b82f6",
    "AmyPred-FRL": "#a78bfa",
    "CrossBeta": "#4ff7c0",
}

WT_NAME = "abeta39_WT"

# ========================== HELPER FUNCTIONS ==========================


def normalize_variant(name: str) -> str:
    """Standardise variant names (handles deletions and asterisks)."""
    if name in ("abeta39_E3", "abeta39_E3*", "abeta39_E3del"):
        return "abeta39_E3del"
    if name in ("abeta39_R5", "abeta39_R5*", "abeta39_R5del"):
        return "abeta39_R5del"
    return name


def short_name(full_name: str) -> str:
    """Remove 'abeta39_' prefix for cleaner labels."""
    return full_name.replace("abeta39_", "")


# ========================== DATA LOADING ==========================


def load_predictor(
    path: str,
    name_col: str,
    value_col: str,
    out_col: str,
    *,
    sep: str = ",",
    delta_from_wt: bool = True,
    coerce_numeric: bool = False,
) -> pd.DataFrame:
    """Load one predictor's table into a single-column, variant-indexed frame.

    Vectorised replacement for the four near-identical per-tool loops. The
    operation is the same for every tool:

      1. read the table,
      2. (optionally) coerce the value column to numeric,
      3. map raw names to canonical variant names,
      4. (optionally) express each value as a difference from the WT value,
      5. drop the WT row and index by variant.

    Parameters
    ----------
    name_col, value_col : the columns holding the variant name and the score.
    out_col             : the name to give the score column in the result.
    delta_from_wt       : subtract the WT score from every variant (TANGO,
                          PASTA, AmyPred-FRL). CrossBeta is used as an absolute
                          value, so pass False.
    coerce_numeric      : TANGO's column arrives as strings with whitespace.
    """
    df = pd.read_csv(path, sep=sep)

    if coerce_numeric:
        df[value_col] = pd.to_numeric(
            df[value_col].astype(str).str.strip(), errors="coerce"
        )

    df["variant"] = df[name_col].map(normalize_variant)

    values = df[value_col]
    if delta_from_wt:
        wt_value = df.loc[df[name_col] == WT_NAME, value_col].iloc[0]
        values = values - wt_value

    out = pd.DataFrame({"variant": df["variant"], out_col: values})
    out = out[df[name_col] != WT_NAME]  # drop WT after the delta is computed
    return out.set_index("variant")


def load_all_tools() -> pd.DataFrame:
    """Load and merge all four predictors on the shared variant index."""
    tango = load_predictor(
        "../tools_assessment/TANGO.txt",
        name_col="Sequence",
        value_col="Aggregation",
        out_col="TANGO",
        sep="\t",
        coerce_numeric=True,
    )
    pasta = load_predictor(
        "../tools_assessment/pasta.csv",
        name_col="Protein name",
        value_col="Best Energy",
        out_col="PASTA",
        sep=";",
    )
    amypred = load_predictor(
        "../tools_assessment/amypred.csv",
        name_col="Name",
        value_col="Probability",
        out_col="AmyPred-FRL",
    )
    crossbeta = load_predictor(
        "../tools_assessment/crossbeta.csv",
        name_col="Query_name",
        value_col="Average_protein_prediction",
        out_col="CrossBeta",
        sep=";",
        delta_from_wt=False,
    )

    return (
        tango.join(pasta, how="outer")
        .join(amypred, how="outer")
        .join(crossbeta, how="outer")
    )


# ========================== VISUALIZATION FUNCTIONS ==========================


def draw_top10(ax, idx_list, tool, direction, color, zscore):
    """Draw horizontal bar plot for top-10 aggregators or disruptors."""
    ax.set_facecolor(SURFACE)
    vals = zscore.loc[idx_list, tool].values
    names = [short_name(n) for n in idx_list]

    order = np.argsort(vals)[::-1] if direction == "agg" else np.argsort(vals)
    vals = vals[order]
    names = [names[i] for i in order]

    y = np.arange(len(vals))
    ax.barh(y, vals, color=color, alpha=0.25, height=0.55, zorder=1)
    ax.hlines(y, 0, vals, color=color, lw=1.2, alpha=0.5, zorder=2)
    ax.scatter(vals, y, color=color, s=45, zorder=3, linewidths=0)

    for i, v in enumerate(vals):
        ha = "left" if v >= 0 else "right"
        off = 0.05 if v >= 0 else -0.05
        ax.text(
            v + off, i, f"{v:+.2f}", ha=ha, va="center",
            fontsize=6.5, fontfamily="monospace", color=TEXT, alpha=0.85,
        )

    ax.axvline(0, color=WT_COL, lw=0.8, linestyle="--", alpha=0.4, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8, fontfamily="monospace", color=TEXT)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7)

    for sp in ax.spines.values():
        sp.set_edgecolor("#1e2540")

    ax.set_title(
        tool, color=TOOL_COLORS[tool], fontsize=10,
        fontfamily="monospace", fontweight="bold", pad=8,
    )
    ax.set_xlabel(
        "z-score (within tool)", color=MUTED, fontsize=7.5,
        fontfamily="monospace", labelpad=4,
    )
    ax.xaxis.set_major_locator(MaxNLocator(5))


def plot_selected_variants(variants_list, title, filename, zscore, tools):
    """Plot grouped bar chart for manually selected variants."""
    full_names = ["abeta39_" + v for v in variants_list]
    existing = [name for name in full_names if name in zscore.index]
    missing = set(full_names) - set(existing)

    if missing:
        print(f"Warning: Missing variants: {missing}")
    if not existing:
        print(f"No available variants for '{title}'")
        return

    df_plot = zscore.loc[existing, tools].copy()
    short_names = [short_name(name) for name in existing]

    n_vars = len(existing)
    fig, ax = plt.subplots(figsize=(max(9, n_vars * 1.15), 6.2), facecolor=BG)
    ax.set_facecolor(SURFACE)

    x = np.arange(n_vars)
    width = 0.2
    offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

    for i, tool in enumerate(tools):
        ax.bar(
            x + offsets[i], df_plot[tool], width, label=tool,
            color=TOOL_COLORS[tool], edgecolor="none", alpha=0.9,
        )

    ax.axhline(0, color=WT_COL, lw=1, linestyle="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, fontsize=9.5, fontfamily="monospace", color=TEXT)
    ax.set_ylabel("z-score", color=MUTED, fontsize=10, fontfamily="monospace")
    ax.tick_params(axis="y", colors=MUTED, labelsize=9)

    for sp in ax.spines.values():
        sp.set_edgecolor("#1e2540")

    ax.legend(
        loc="upper right", frameon=True, framealpha=0.2,
        edgecolor=MUTED, facecolor=BG, fontsize=9, labelcolor=TEXT,
    )
    ax.set_title(
        title, color=TEXT, fontsize=12,
        fontfamily="monospace", fontweight="bold", pad=12,
    )

    fig.tight_layout()
    fig.savefig("../../images/" + filename, dpi=170, bbox_inches="tight", facecolor=BG)
    print(f"Saved: {filename}")
    plt.close(fig)


# ========================== MAIN ANALYSIS ==========================


def compute_consensus(raw, tools):
    """Z-score each tool (PASTA inverted: lower energy = stronger aggregation)
    and return (zscore_frame, raw_with_consensus_column)."""
    raw_norm = raw.copy()
    raw_norm["PASTA"] = -raw_norm["PASTA"]

    zscore = pd.DataFrame(index=raw_norm.index)
    for tool in tools:
        col = raw_norm[tool]
        zscore[tool] = (col - col.mean()) / col.std()

    raw = raw.copy()
    raw["consensus_z"] = zscore.mean(axis=1)
    return zscore, raw


def print_text_summary(zscore, raw, tools, top10_agg, top10_dis):
    for tool in tools:
        print(f"\n{'=' * 65}")
        print(f" {tool} - TOP 10 AGGREGATORS")
        print(f"{'=' * 65}")
        for i, idx in enumerate(top10_agg[tool], 1):
            print(
                f" {i:2d}. {short_name(idx):<20}  "
                f"z={zscore.loc[idx, tool]:+.3f}   raw={raw.loc[idx, tool]:+.4f}"
            )
        print(f"\n {tool} - TOP 10 DISRUPTORS")
        print(f"{'-' * 65}")
        for i, idx in enumerate(top10_dis[tool], 1):
            print(
                f" {i:2d}. {short_name(idx):<20}  "
                f"z={zscore.loc[idx, tool]:+.3f}   raw={raw.loc[idx, tool]:+.4f}"
            )


def plot_consensus_ranking(raw, n):
    """Full-width consensus ranking bar chart."""
    fig1, ax0 = plt.subplots(figsize=(20, 7.5), facecolor=BG)
    ax0.set_facecolor(SURFACE)

    cons_sorted = raw["consensus_z"].sort_values(ascending=False)
    labels = [short_name(name) for name in cons_sorted.index]
    colors_bar = [AGG_COL if v >= 0 else DIS_COL for v in cons_sorted.values]

    bars = ax0.bar(range(len(cons_sorted)), cons_sorted.values,
                   color=colors_bar, width=0.75, zorder=2)

    for i, v in enumerate(cons_sorted.values):
        if i < 10 or i >= len(cons_sorted) - 10:
            bars[i].set_edgecolor(ACC_COL)
            bars[i].set_linewidth(1.3)
        else:
            bars[i].set_alpha(0.65)

    ax0.axhline(0, color=WT_COL, lw=0.8, linestyle="--", alpha=0.5)

    for i, (label, v) in enumerate(zip(labels, cons_sorted.values)):
        if i < 10 or i >= len(cons_sorted) - 10:
            va = "bottom" if v >= 0 else "top"
            off = 0.04 if v >= 0 else -0.04
            ax0.text(i, v + off, label, ha="center", va=va, fontsize=6.2,
                     fontfamily="monospace", color=ACC_COL, rotation=90,
                     fontweight="bold")

    ax0.set_xlim(-0.8, len(cons_sorted) - 0.2)
    ax0.set_xticks([])
    ax0.set_ylabel("Consensus z-score", color=MUTED, fontsize=9.5,
                   fontfamily="monospace", labelpad=6)
    ax0.set_title(
        f"Aβ39 · Consensus Ranking of All {n} Variants\n"
        "Consensus = average z-score across four tools",
        color=TEXT, fontsize=12, fontfamily="monospace", fontweight="bold", pad=15,
    )

    agg_p = mpatches.Patch(color=AGG_COL, label="Increased aggregation")
    dis_p = mpatches.Patch(color=DIS_COL, label="Decreased aggregation")
    top_p = mpatches.Patch(facecolor="none", edgecolor=ACC_COL, lw=1.5,
                           label="Top-10 (highlighted)")
    ax0.legend(handles=[agg_p, dis_p, top_p], loc="upper right",
               frameon=True, framealpha=0.25, edgecolor=MUTED,
               facecolor=BG, fontsize=9, labelcolor=TEXT)

    fig1.tight_layout(pad=1.8)
    fig1.savefig("../../images/Aβ39_ranking_consensus.png", dpi=170,
                 bbox_inches="tight", facecolor=BG)
    print("Saved: Aβ39_ranking_consensus.png")
    plt.close(fig1)


def plot_top10_panels(zscore, tools, top10, direction, suptitle, color, filename):
    """Four-panel top-10 figure (aggregators or disruptors)."""
    fig, axes = plt.subplots(1, 4, figsize=(20, 7.5), facecolor=BG)
    fig.subplots_adjust(wspace=0.45, left=0.06, right=0.96, top=0.85, bottom=0.12)
    fig.suptitle(suptitle, color=color, fontsize=13,
                 fontfamily="monospace", fontweight="bold", y=0.96)
    for ax, tool in zip(axes, tools):
        draw_top10(ax, top10[tool], tool, direction, TOOL_COLORS[tool], zscore)
    fig.savefig(f"../../images/{filename}", dpi=170, bbox_inches="tight", facecolor=BG)
    print(f"Saved: {filename}")
    plt.close(fig)


def main():
    raw = load_all_tools().dropna()
    n = len(raw)
    print(f"Loaded {n} variants common to all four prediction tools.\n")

    tools = ["TANGO", "PASTA", "AmyPred-FRL", "CrossBeta"]
    zscore, raw = compute_consensus(raw, tools)

    top10_agg = {t: zscore[t].nlargest(10).index for t in tools}
    top10_dis = {t: zscore[t].nsmallest(10).index for t in tools}

    print_text_summary(zscore, raw, tools, top10_agg, top10_dis)

    plot_consensus_ranking(raw, n)
    plot_top10_panels(
        zscore, tools, top10_agg, "agg",
        "Aβ39 · Top-10 Aggregators by Individual Tool\n(z-score within each predictor)",
        AGG_COL, "Aβ39_top10_aggregators.png",
    )
    plot_top10_panels(
        zscore, tools, top10_dis, "dis",
        "Aβ39 · Top-10 Disruptors by Individual Tool\n(z-score within each predictor)",
        DIS_COL, "Aβ39_top10_disruptors.png",
    )

    plot_selected_variants(
        ["E3G", "E22Q", "R5L", "E3Q", "R5Q", "V12I"],
        "Aβ39 · Selected Aggregators", "Aβ39_selected_aggregators.png", zscore, tools,
    )
    plot_selected_variants(
        ["G25D", "H13P", "G25S", "R5P", "H6R", "H6Q"],
        "Aβ39 · Selected Disruptors", "Aβ39_selected_disruptors.png", zscore, tools,
    )

    print("\nAnalysis completed successfully.")


if __name__ == "__main__":
    main()
