"""
Visualization of Aβ42 mutations — sequence alignment style.

Input file: mutations_list.txt  (one mutation per line, format p.Ala673Gly)
Output file: abeta42_mutations_alignment.png

Usage: python3 abeta42_alignment_plotting.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
import re
import os

INPUT_FILE  = "../additional/mutations_list.txt"
OUTPUT_FILE = "abeta42_mutations_alignment_dark.png"
DPI         = 150

COL_VARIANT    = '#4ff7c0'
COL_WT_BG      = '#111520'
COL_WT_TEXT    = '#4ff7c0'
COL_WT_STRIPE  = '#161b2e'
COL_DASH       = '#334166'
COL_NUM        = '#64748b'
COL_HOTSPOT    = '#1a2240'
COL_PATHOGENIC = '#ef4444'
COL_PROTECTIVE = '#22c55e'

HOTSPOTS = [
    (15, 23, 'CHC / CAA region'),
    (28, 42, 'Hydrophobic core'),
]

NAMED_MUTATIONS = {
    (22, 'G'): ('Arctic',    'pathogenic'),
    (22, 'Q'): ('Dutch',     'pathogenic'),
    (23, 'N'): ('Iowa',      'pathogenic'),
    (2,  'T'): ('Icelandic', 'protective'),
}

WT      = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
OFFSET  = 671
SEQ_LEN = len(WT)

AA3TO1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
    'Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
    'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
    'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
    'Ter':'*'
}

def conv(s):
    """Convert three-letter amino acid codes to one-letter."""
    for three, one in AA3TO1.items():
        s = s.replace(three, one)
    return s


def parse_mutations(filepath):
    """Parse mutation list from file."""
    mutations = []
    seen = set()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)', line)
            if not m:
                continue
            wt_aa   = conv(m.group(1))
            app_pos = int(m.group(2))
            mut_aa  = conv(m.group(3))
            pos     = app_pos - OFFSET
            if pos < 1 or pos > SEQ_LEN or mut_aa == '*':
                continue
            key = (pos, mut_aa)
            if key in seen:
                continue
            seen.add(key)
            mutations.append({'pos': pos, 'wt_aa': wt_aa, 'mut_aa': mut_aa})
    mutations.sort(key=lambda x: (x['pos'], x['mut_aa']))
    print(f"Loaded mutations: {len(mutations)}")
    return mutations


def draw(mutations, output_file):
    """Draw the mutation alignment plot."""
    n_rows = len(mutations)
    CHAR_W = 0.50; CHAR_H = 0.22; LEFT_PAD = 1.8; RIGHT_PAD = 0.3; TOP_PAD = 0.95
    FIG_W = LEFT_PAD + SEQ_LEN * CHAR_W + RIGHT_PAD
    FIG_H = TOP_PAD + n_rows * CHAR_H + 0.4

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W); ax.set_ylim(0, FIG_H); ax.axis('off')
    fig.patch.set_facecolor('#0a0c14')

    def cx(p): return LEFT_PAD + (p - 0.5) * CHAR_W
    mut_top = FIG_H - TOP_PAD
    def ry(i): return mut_top - (i + 0.5) * CHAR_H

    mono      = FontProperties(family='monospace', size=7.5)
    mono_bold = FontProperties(family='monospace', size=7.5, weight='bold')

    ax.text(0.15, FIG_H - 0.22, "Aβ42 Mutation Map",
            fontsize=10, fontweight='bold', color='#e2e8f0', va='top')

    mut_p = mpatches.Patch(color=COL_VARIANT,    label='Mutation')
    pa_p  = mpatches.Patch(color=COL_PATHOGENIC, label='Pathogenic (named)')
    pr_p  = mpatches.Patch(color=COL_PROTECTIVE, label='Protective (named)')
    hs_p  = mpatches.Patch(facecolor=COL_HOTSPOT, edgecolor='#fbbf24', lw=0.8, label='Aggregation hot-spot')
    ax.legend(handles=[mut_p, pa_p, pr_p, hs_p], loc='upper right',
              bbox_to_anchor=(FIG_W - 0.05, FIG_H - 0.05),
              fontsize=7.5, framealpha=0.95, edgecolor='#1e2540', facecolor='#111520', labelcolor='#e2e8f0',
              handlelength=1.2, handleheight=0.9)

    num_y = FIG_H - 0.55
    for i in range(SEQ_LEN):
        p = i + 1
        if p == 1 or p % 5 == 0 or p == SEQ_LEN:
            ax.text(cx(p), num_y, str(p), fontsize=5.8,
                    ha='center', va='center', color=COL_NUM, fontproperties=mono)

    wt_y = FIG_H - TOP_PAD + CHAR_H * 0.5
    ax.add_patch(plt.Rectangle((LEFT_PAD, wt_y - CHAR_H * 0.48),
        SEQ_LEN * CHAR_W, CHAR_H * 0.96, facecolor=COL_WT_BG, edgecolor='#1e2540', lw=0.7, zorder=0))
    for i in range(SEQ_LEN):
        if (i // 5) % 2 == 0:
            ax.add_patch(plt.Rectangle((LEFT_PAD + i * CHAR_W, wt_y - CHAR_H * 0.48),
                CHAR_W, CHAR_H * 0.96, facecolor=COL_WT_STRIPE, edgecolor='none', zorder=1))

    hs_full_bottom = mut_top - n_rows * CHAR_H - CHAR_H * 0.1
    for hs_start, hs_end, hs_name in HOTSPOTS:
        hx0 = LEFT_PAD + (hs_start - 1) * CHAR_W
        hx1 = LEFT_PAD + hs_end * CHAR_W
        ax.add_patch(plt.Rectangle((hx0, wt_y - CHAR_H * 0.48), hx1 - hx0, CHAR_H * 0.96,
            facecolor=COL_HOTSPOT, edgecolor='none', alpha=0.9, zorder=1))
        ax.add_patch(plt.Rectangle((hx0, hs_full_bottom), hx1 - hx0,
            wt_y + CHAR_H * 0.48 - hs_full_bottom,
            facecolor=COL_HOTSPOT, edgecolor='none', alpha=0.5, zorder=0))
        ax.text((hx0 + hx1) / 2, wt_y + CHAR_H * 0.62, hs_name,
                fontsize=5.5, ha='center', va='bottom', color='#fbbf24', style='italic')

    ax.text(LEFT_PAD - 0.12, wt_y, "WT", fontsize=7, fontweight='bold',
            color=COL_WT_TEXT, ha='right', va='center')
    for i, aa in enumerate(WT):
        ax.text(cx(i + 1), wt_y, aa, fontsize=7.5, ha='center', va='center',
                color=COL_WT_TEXT, fontproperties=mono_bold, zorder=2)

    sep_y = mut_top - CHAR_H * 0.05
    ax.plot([LEFT_PAD - 0.15, LEFT_PAD + SEQ_LEN * CHAR_W + 0.05],
            [sep_y, sep_y], color='#1e2540', lw=0.6, ls='--')

    for row_i, mut in enumerate(mutations):
        y = ry(row_i); pos = mut['pos']; mut_aa = mut['mut_aa']
        named_key = (pos, mut_aa)
        if named_key in NAMED_MUTATIONS:
            _, effect = NAMED_MUTATIONS[named_key]
            color = COL_PATHOGENIC if effect == 'pathogenic' else COL_PROTECTIVE
        else:
            color = COL_VARIANT
        for i in range(SEQ_LEN):
            ax.text(cx(i + 1), y, '–', fontsize=7, ha='center', va='center',
                    color=COL_DASH, fontproperties=mono)
        ax.text(cx(pos), y, mut_aa, fontsize=7.5, ha='center', va='center',
                color=color, fontproperties=mono_bold, zorder=3)
        if named_key in NAMED_MUTATIONS:
            name, _ = NAMED_MUTATIONS[named_key]
            ax.text(LEFT_PAD - 1.65, y, name, fontsize=5.8, ha='left', va='center',
                    color=color, style='italic')

    plt.savefig(output_file, dpi=DPI, facecolor='#0a0c14')
    print(f"Saved: {os.path.abspath(output_file)}")
    plt.close()


if __name__ == '__main__':
    if not os.path.exists(INPUT_FILE):
        print(f"Error: file '{INPUT_FILE}' not found.")
        exit(1)
    mutations = parse_mutations(INPUT_FILE)
    draw(mutations, OUTPUT_FILE)