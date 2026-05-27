#!/usr/bin/env python3
"""
Aβ39 Variant Aggregation Propensity Analysis
============================================

Comprehensive consensus analysis of Aβ39 peptide variants based on four
aggregation prediction tools: TANGO, PASTA, AmyPred-FRL, and CrossBeta.

The script performs z-score normalization per tool, calculates a consensus
score, identifies top aggregators and disruptors, and generates publication-ready
visualizations.

Author: [Your Name]
Date: May 2026
Version: 1.0
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from typing import List, Dict

# ========================== PLOT STYLE CONFIGURATION ==========================

matplotlib.rcParams['font.family'] = 'monospace'

# Color palette (dark theme)
BG = '#0a0c14'
SURFACE = '#111520'
MUTED = '#64748b'
TEXT = '#e2e8f0'
WT_COL = '#4ff7c0'
AGG_COL = '#ef4444'      # Increased aggregation
DIS_COL = '#3b82f6'      # Decreased aggregation
ACC_COL = '#fbbf24'      # Highlight / consensus

TOOL_COLORS = {
    'TANGO': '#ef4444',
    'PASTA': '#3b82f6',
    'AmyPred-FRL': '#a78bfa',
    'CrossBeta': '#4ff7c0',
}

# ========================== HELPER FUNCTIONS ==========================

def normalize_variant(name: str) -> str:
    """Standardize variant names (handles deletions and asterisks)."""
    if name in ('abeta39_E3', 'abeta39_E3*', 'abeta39_E3del'):
        return 'abeta39_E3del'
    if name in ('abeta39_R5', 'abeta39_R5*', 'abeta39_R5del'):
        return 'abeta39_R5del'
    return name


def short_name(full_name: str) -> str:
    """Remove 'abeta39_' prefix for cleaner labels."""
    return full_name.replace('abeta39_', '')


# ========================== DATA LOADING FUNCTIONS ==========================

def load_tango(path: str) -> pd.DataFrame:
    """Load and process TANGO results."""
    df = pd.read_csv(path, sep='\t')
    df['Aggregation'] = pd.to_numeric(df['Aggregation'].astype(str).str.strip(), errors='coerce')
    
    wt_value = df.loc[df['Sequence'] == 'abeta39_WT', 'Aggregation'].iloc[0]
    
    data = []
    for _, row in df.iterrows():
        if row['Sequence'] == 'abeta39_WT':
            continue
        variant = normalize_variant(row['Sequence'])
        data.append({'variant': variant, 'TANGO': row['Aggregation'] - wt_value})
    
    return pd.DataFrame(data).set_index('variant')


def load_pasta(path: str) -> pd.DataFrame:
    """Load and process PASTA results."""
    df = pd.read_csv(path, sep=';')
    wt_value = df.loc[df['Protein name'] == 'abeta39_WT', 'Best Energy'].iloc[0]
    
    data = []
    for _, row in df.iterrows():
        if row['Protein name'] == 'abeta39_WT':
            continue
        variant = normalize_variant(row['Protein name'])
        data.append({'variant': variant, 'PASTA': row['Best Energy'] - wt_value})
    
    return pd.DataFrame(data).set_index('variant')


def load_amypred(path: str) -> pd.DataFrame:
    """Load and process AmyPred-FRL results."""
    df = pd.read_csv(path)
    wt_value = df.loc[df['Name'] == 'abeta39_WT', 'Probability'].iloc[0]
    
    data = []
    for _, row in df.iterrows():
        if row['Name'] == 'abeta39_WT':
            continue
        variant = normalize_variant(row['Name'])
        data.append({'variant': variant, 'AmyPred-FRL': row['Probability'] - wt_value})
    
    return pd.DataFrame(data).set_index('variant')


def load_crossbeta(path: str) -> pd.DataFrame:
    """Load and process CrossBeta results."""
    df = pd.read_csv(path, sep=';')
    wt_value = df.loc[df['Query_name'] == 'abeta39_WT', 'Average_protein_prediction'].iloc[0]
    
    data = []
    for _, row in df.iterrows():
        if row['Query_name'] == 'abeta39_WT':
            continue
        variant = normalize_variant(row['Query_name'])
        data.append({'variant': variant, 'CrossBeta': row['Average_protein_prediction']})
    
    return pd.DataFrame(data).set_index('variant')


# ========================== VISUALIZATION FUNCTIONS ==========================

def draw_top10(ax, idx_list: List[str], tool: str, direction: str, color: str):
    """Draw horizontal bar plot for top-10 aggregators or disruptors."""
    ax.set_facecolor(SURFACE)
    vals = zscore.loc[idx_list, tool].values
    names = [short_name(n) for n in idx_list]
    
    order = np.argsort(vals)[::-1] if direction == 'agg' else np.argsort(vals)
    vals = vals[order]
    names = [names[i] for i in order]
    
    y = np.arange(len(vals))
    ax.barh(y, vals, color=color, alpha=0.25, height=0.55, zorder=1)
    ax.hlines(y, 0, vals, color=color, lw=1.2, alpha=0.5, zorder=2)
    ax.scatter(vals, y, color=color, s=45, zorder=3, linewidths=0)
    
    for i, v in enumerate(vals):
        ha = 'left' if v >= 0 else 'right'
        off = 0.05 if v >= 0 else -0.05
        ax.text(v + off, i, f'{v:+.2f}', ha=ha, va='center',
                fontsize=6.5, fontfamily='monospace', color=TEXT, alpha=0.85)
    
    ax.axvline(0, color=WT_COL, lw=0.8, linestyle='--', alpha=0.4, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8, fontfamily='monospace', color=TEXT)
    ax.tick_params(axis='y', length=0, pad=4)
    ax.tick_params(axis='x', colors=MUTED, labelsize=7)
    
    for sp in ax.spines.values():
        sp.set_edgecolor('#1e2540')
    
    ax.set_title(tool, color=TOOL_COLORS[tool], fontsize=10, 
                 fontfamily='monospace', fontweight='bold', pad=8)
    ax.set_xlabel('z-score (within tool)', color=MUTED, fontsize=7.5,
                  fontfamily='monospace', labelpad=4)
    ax.xaxis.set_major_locator(MaxNLocator(5))


def plot_selected_variants(variants_list: List[str], title: str, filename: str, 
                          bar_color: str):
    """Plot grouped bar chart for manually selected variants."""
    full_names = ['abeta39_' + v for v in variants_list]
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
    offsets = [-1.5*width, -0.5*width, 0.5*width, 1.5*width]
    
    for i, tool in enumerate(tools):
        offset = offsets[i]
        ax.bar(x + offset, df_plot[tool], width, label=tool,
               color=TOOL_COLORS[tool], edgecolor='none', alpha=0.9)
    
    ax.axhline(0, color=WT_COL, lw=1, linestyle='--', alpha=0.6)
    
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, fontsize=9.5, fontfamily='monospace', color=TEXT)
    ax.set_ylabel('z-score', color=MUTED, fontsize=10, fontfamily='monospace')
    ax.tick_params(axis='y', colors=MUTED, labelsize=9)
    
    for sp in ax.spines.values():
        sp.set_edgecolor('#1e2540')
    
    ax.legend(loc='upper right', frameon=True, framealpha=0.2, 
              edgecolor=MUTED, facecolor=BG, fontsize=9, labelcolor=TEXT)
    
    ax.set_title(title, color=TEXT, fontsize=12, 
                 fontfamily='monospace', fontweight='bold', pad=12)
    
    fig.tight_layout()
    fig.savefig('../../images/'+filename, dpi=170, bbox_inches='tight', facecolor=BG)
    print(f"✓ Saved: {filename}")
    plt.close(fig)


# ========================== MAIN ANALYSIS ==========================

def main():
    global zscore, tools, raw   # for use in plotting functions
    
    # Load data
    tango = load_tango('../tools_assessment/TANGO.txt')
    pasta = load_pasta('../tools_assessment/pasta.csv')
    amypred = load_amypred('../tools_assessment/amypred.csv')
    crossbeta = load_crossbeta('../tools_assessment/crossbeta.csv')
    
    # Merge all tools
    raw = tango.join(pasta, how='outer') \
               .join(amypred, how='outer') \
               .join(crossbeta, how='outer')
    
    raw = raw.dropna()
    n = len(raw)
    
    print(f"Loaded {n} variants common to all four prediction tools.\n")
    
    # Z-score normalization
    tools = ['TANGO', 'PASTA', 'AmyPred-FRL', 'CrossBeta']
    raw_norm = raw.copy()
    raw_norm['PASTA'] = -raw_norm['PASTA']   # Invert PASTA (lower energy = stronger aggregation)
    
    zscore = pd.DataFrame(index=raw_norm.index)
    for tool in tools:
        zscore[tool] = (raw_norm[tool] - raw_norm[tool].mean()) / raw_norm[tool].std()
    
    raw['consensus_z'] = zscore.mean(axis=1)
    
    # Top-10 per tool
    top10_agg = {t: zscore[t].nlargest(10).index for t in tools}
    top10_dis = {t: zscore[t].nsmallest(10).index for t in tools}
    
    # ========================== TEXT SUMMARY ==========================
    
    for tool in tools:
        print(f"\n{'═' * 65}")
        print(f" {tool} — TOP 10 AGGREGATORS")
        print(f"{'═' * 65}")
        for i, idx in enumerate(top10_agg[tool], 1):
            print(f" {i:2d}. {short_name(idx):<20}  z={zscore.loc[idx, tool]:+.3f}   "
                  f"raw={raw.loc[idx, tool]:+.4f}")
        
        print(f"\n {tool} — TOP 10 DISRUPTORS")
        print(f"{'─' * 65}")
        for i, idx in enumerate(top10_dis[tool], 1):
            print(f" {i:2d}. {short_name(idx):<20}  z={zscore.loc[idx, tool]:+.3f}   "
                  f"raw={raw.loc[idx, tool]:+.4f}")
    
    # ========================== PLOTS ==========================
    
    # 1. Consensus Ranking
    fig1, ax0 = plt.subplots(figsize=(20, 7.5), facecolor=BG)
    ax0.set_facecolor(SURFACE)
    
    cons_sorted = raw['consensus_z'].sort_values(ascending=False)
    labels = [short_name(n) for n in cons_sorted.index]
    colors_bar = [AGG_COL if v >= 0 else DIS_COL for v in cons_sorted.values]
    
    bars = ax0.bar(range(len(cons_sorted)), cons_sorted.values,
                   color=colors_bar, width=0.75, zorder=2)
    
    for i, v in enumerate(cons_sorted.values):
        if i < 10 or i >= len(cons_sorted) - 10:
            bars[i].set_edgecolor(ACC_COL)
            bars[i].set_linewidth(1.3)
        else:
            bars[i].set_alpha(0.65)
    
    ax0.axhline(0, color=WT_COL, lw=0.8, linestyle='--', alpha=0.5)
    
    for i, (label, v) in enumerate(zip(labels, cons_sorted.values)):
        if i < 10 or i >= len(cons_sorted) - 10:
            va = 'bottom' if v >= 0 else 'top'
            off = 0.04 if v >= 0 else -0.04
            ax0.text(i, v + off, label, ha='center', va=va, fontsize=6.2,
                     fontfamily='monospace', color=ACC_COL, rotation=90, fontweight='bold')
    
    ax0.set_xlim(-0.8, len(cons_sorted) - 0.2)
    ax0.set_xticks([])
    ax0.set_ylabel('Consensus z-score', color=MUTED, fontsize=9.5, 
                   fontfamily='monospace', labelpad=6)
    
    ax0.set_title(f'Aβ39 · Consensus Ranking of All {n} Variants\n'
                  'Consensus = average z-score across four tools', 
                  color=TEXT, fontsize=12, fontfamily='monospace', 
                  fontweight='bold', pad=15)
    
    # Legend
    agg_p = mpatches.Patch(color=AGG_COL, label='Increased aggregation')
    dis_p = mpatches.Patch(color=DIS_COL, label='Decreased aggregation')
    top_p = mpatches.Patch(facecolor='none', edgecolor=ACC_COL, lw=1.5, 
                           label='Top-10 (highlighted)')
    ax0.legend(handles=[agg_p, dis_p, top_p], loc='upper right', 
               frameon=True, framealpha=0.25, edgecolor=MUTED, 
               facecolor=BG, fontsize=9, labelcolor=TEXT)
    
    fig1.tight_layout(pad=1.8)
    fig1.savefig('../../images/Aβ39_ranking_consensus.png', dpi=170, bbox_inches='tight', facecolor=BG)
    print("✓ Saved: Aβ39_ranking_consensus.png")
    plt.close(fig1)
    
    # 2. Top-10 Aggregators
    fig2, axes2 = plt.subplots(1, 4, figsize=(20, 7.5), facecolor=BG)
    fig2.subplots_adjust(wspace=0.45, left=0.06, right=0.96, top=0.85, bottom=0.12)
    fig2.suptitle('Aβ39 · Top-10 Aggregators by Individual Tool\n'
                  '(z-score within each predictor)', color=AGG_COL, 
                  fontsize=13, fontfamily='monospace', fontweight='bold', y=0.96)
    
    for ax, tool in zip(axes2, tools):
        draw_top10(ax, top10_agg[tool], tool, 'agg', TOOL_COLORS[tool])
    
    fig2.savefig('../../images/Aβ39_top10_aggregators.png', dpi=170, bbox_inches='tight', facecolor=BG)
    print("✓ Saved: Aβ39_top10_aggregators.png")
    plt.close(fig2)
    
    # 3. Top-10 Disruptors
    fig3, axes3 = plt.subplots(1, 4, figsize=(20, 7.5), facecolor=BG)
    fig3.subplots_adjust(wspace=0.45, left=0.06, right=0.96, top=0.85, bottom=0.12)
    fig3.suptitle('Aβ39 · Top-10 Disruptors by Individual Tool\n'
                  '(z-score within each predictor)', color=DIS_COL, 
                  fontsize=13, fontfamily='monospace', fontweight='bold', y=0.96)
    
    for ax, tool in zip(axes3, tools):
        draw_top10(ax, top10_dis[tool], tool, 'dis', TOOL_COLORS[tool])
    
    fig3.savefig('../../images/Aβ39_top10_disruptors.png', dpi=170, bbox_inches='tight', facecolor=BG)
    print("✓ Saved: Aβ39_top10_disruptors.png")
    plt.close(fig3)
    
    # 4. Selected variants
    agg_selected = ['E3G', 'E22Q', 'R5L', 'E3Q', 'R5Q', 'V12I']
    dis_selected = ['G25D', 'H13P', 'G25S', 'R5P', 'H6R', 'H6Q']
    
    plot_selected_variants(agg_selected,
                           'Aβ39 · Selected Aggregators',
                           'Aβ39_selected_aggregators.png',
                           AGG_COL)
    
    plot_selected_variants(dis_selected,
                           'Aβ39 · Selected Disruptors',
                           'Aβ39_selected_disruptors.png',
                           DIS_COL)
    
    print("\nAnalysis completed successfully.")


if __name__ == "__main__":
    main()
