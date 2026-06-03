# Step 1: Amyloid Prediction Tools Analysis

Analysis of Aβ42 variants using six external amyloidogenicity predictors, followed by consensus ranking and multi-tool visualisations.

---

## Repository Structure

```
1_step_prediction_tools_analysis/
├── predictions/              # Raw output files downloaded from web tools
│   ├── tango/                # tango_input_aggregation.txt
│   ├── pasta/                # pasta.csv
│   ├── waltz/                # Waltz_bestoverall/highsens/highspecif.txt
│   ├── cross-beta/           # cross-beta_result.csv + per-variant PNGs
│   ├── amypred/              # amypred-frl.csv
│   └── amylogram/            # AmyloGram_results_ed.csv
├── notebooks/                # One Jupyter notebook per tool + integration steps
└── results/                  # Output PNGs, one subfolder per tool
    ├── tango/
    ├── pasta/
    ├── waltz/
    ├── cross-beta/
    ├── amypred/
    ├── amylogram/
    └── rank/
```

---

## Requirements

```
numpy>=2.3
matplotlib>=3.10
scipy>=1.16
pandas>=2.3
```

Install:

```bash
pip install -r requirements.txt
```

Each notebook is run **from its own tool subfolder within `notebooks/`**, so all relative paths (`../predictions/...`) resolve correctly. Notebooks are independent and can be run in any order except `amyloid_ranking` and `consensus_heatmap`, which aggregate results from the other tools and must be run last.

---

## Execution

```bash
bash run_analysis.sh
```

Where `run_analysis.sh` contains:

```bash
#!/bin/bash
set -e
cd notebooks

for nb in tango_aggregation_analysis waltz_region_map_analysis \
           cross-beta_analysis amypred_analysis amylogram_analysis \
           pasta_analysis vem amyloid_ranking consensus_heatmap; do
    echo "Running $nb..."
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=600 \
        ${nb}.ipynb \
        --output ${nb}_executed.ipynb
done
echo "Done."
```

> Run from the repository root. Executed notebooks with outputs are saved alongside the originals as `*_executed.ipynb` and do not overwrite the source.

---

## Substeps

### 1. TANGO — Aggregation Score

**Notebook:** `notebooks/tango_aggregation_analysis.ipynb`

**Description:** Reads the TANGO tab-separated output (`tango_input_aggregation.txt`) containing a total aggregation score for each Aβ42 variant. Computes Δ relative to the Wildtype score (WT = 1526.48) and classifies each mutation into five groups:

| Group | Δ score |
|---|---|
| Strongly destabilizing | Δ < −50 |
| Weakly destabilizing | −50 ≤ Δ < −5 |
| Neutral | −5 ≤ Δ ≤ +5 |
| Enhancing | +5 < Δ ≤ +30 |
| Strongly enhancing | Δ > +30 |

Produces two plots: a horizontal lollipop chart of Δ scores and a five-column classification table.

**Cautions:**

| Item | Detail |
|---|---|
| `WT_SCORE` | Derived at runtime from the file (`Wildtype_Abeta_42` row). If the wildtype entry is renamed or absent, the script will crash. |
| Thresholds (`THRESH`) | Hardcoded as `{strong_down: -50, weak_down: -5, neutral_hi: +5, up: +30}`. Adjust at the top of the data-loading cell if different classification boundaries are needed. |
| Input path | `../predictions/tango/tango_input_aggregation.txt` — tab-separated, must have columns `Sequence` and `Aggregation`. |
| Output directory | Figures are saved to the **current working directory** of the notebook (i.e. `notebooks/`). Move to `results/tango/` manually or via the script. |

**Results** (saved to `results/tango/`):

| File | Content |
|---|---|
| `tango_lollipop.png` | Horizontal lollipop — Δ aggregation per variant, coloured by group |
| `tango_groups.png` | Five-column classification table with Δ values |

---

### 2. WALTZ — Amyloid Region Map

**Notebook:** `notebooks/waltz_region_map_analysis.ipynb`

**Description:** Parses three WALTZ output files (Best Overall, High Sensitivity, High Specificity) and draws a per-variant amyloid region map for each mode. Each row represents a variant; coloured rectangles mark predicted amyloid-forming segments along the 42-residue sequence. Segment opacity encodes the WALTZ score (higher score → brighter). Variants lacking any predicted region are highlighted in red.

**Cautions:**

| Item | Detail |
|---|---|
| `SEQ_LEN = 42` | Hardcoded — change if analysing a peptide of different length. |
| Parser | Custom `parse_waltz()` function; expects WALTZ's native `>name / Positions header / start-end\tseq\tscore` format. Any structural change in the tool's output will break parsing. |
| Three separate files | All three (`Waltz_bestoverall.txt`, `Waltz_highsens.txt`, `Waltz_highspecif.txt`) must be present in `predictions/waltz/`. |
| `score_to_alpha` range | Mapped from 83–100; if WALTZ scores outside this range appear, all blocks will be clipped to min/max opacity. |

**Results** (saved to `results/waltz/`):

| File | Content |
|---|---|
| `waltz_region_map.png` | All three modes combined in one figure |
| `waltz_region_map_best_overall.png` | Best Overall mode only |
| `waltz_region_map_high_sensitivity.png` | High Sensitivity mode only |
| `waltz_region_map_high_specificity.png` | High Specificity mode only |

---

### 3. Cross-beta — Per-residue Amyloid Score

**Notebook:** `notebooks/cross-beta_analysis.ipynb`

**Description:** Reads a CSV containing per-residue Cross-beta amyloid scores for all 66 Aβ42 variants. Constructs a (66 × 42) score matrix and computes Δ relative to Wildtype. Produces four plots: an absolute-score heatmap, a Δ-score heatmap, a lollipop of average Δ per variant, and overlaid line profiles for the top 5 stabilizing and top 5 destabilizing mutants.

**Cautions:**

| Item | Detail |
|---|---|
| `THRESH = 0.005` | Neutrality threshold for colour-coding lollipop points. Adjust if the score distribution is very narrow or wide. |
| CSV separator | Semicolon-delimited (`sep=';'`). The `Amino_acids_score` column contains Python-literal lists of `{aa: score}` dicts, parsed with `ast.literal_eval`. |
| Input path | `../predictions/cross-beta/cross-beta_result.csv` |
| Per-variant PNGs | Raw Cross-beta result graphs (one PNG per variant) are stored in `predictions/cross-beta/imgs cross-beta/` — these are the tool's own output, not generated by the notebook. |

**Results** (saved to `results/cross-beta/`):

| File | Content |
|---|---|
| `cross-beta_heatmap_abs.png` | Absolute per-residue score heatmap (plasma colormap) |
| `cross-beta_heatmap_delta.png` | Δ score vs Wildtype heatmap (diverging colormap) |
| `cross-beta_lollipop.png` | Average Δ score per variant |
| `cross-beta_profiles.png` | Line profiles for top 5 destabilizing / stabilizing mutants |

---

### 4. AmyPred-FRL — Amyloid Probability

**Notebook:** `notebooks/amypred_analysis.ipynb`

**Description:** Reads the AmyPred-FRL CSV (`amypred-frl.csv`) with one amyloid probability per variant. Computes ΔP relative to the Wildtype probability and classifies variants as decreasing, neutral, or increasing amyloidogenicity. Produces a lollipop chart and a grouped bar chart.

**Cautions:**

| Item | Detail |
|---|---|
| `THRESH = 0.003` | Neutrality threshold; smaller than for other tools because AmyPred-FRL probabilities are discretised and differences are small. |
| CSV format | Semicolon-delimited, UTF-8-BOM (`encoding='utf-8-sig'`); columns `Name` and `Probability`. |
| Input path | `../predictions/amypred/amypred-frl.csv` |

**Results** (saved to `results/amypred/`):

| File | Content |
|---|---|
| `amypred_lollipop.png` | ΔProbability lollipop (blue = decreasing, red = increasing amyloidogenicity) |
| `amypred_bars.png` | Grouped bar chart of absolute probabilities |

---

### 5. AmyloGram — Amyloid Probability

**Notebook:** `notebooks/amylogram_analysis.ipynb`

**Description:** Reads the edited AmyloGram CSV (`AmyloGram_results_ed.csv`). The tool returns a highly discretised probability (few unique values), so the notebook first inventories distinct values and prints their counts before plotting. Produces a lollipop and a grouped bar chart analogous to AmyPred-FRL.

**Cautions:**

| Item | Detail |
|---|---|
| `THRESH = 0.001` | Neutrality threshold set very low due to output discreteness. |
| Edited CSV | The file `AmyloGram_results_ed.csv` is a manually consolidated version of `AmyloGram_results1.csv` and `AmyloGram_results2.csv` (two separate AmyloGram batch submissions). Do not swap it for the raw files without re-merging. |
| Input path | `../predictions/amylogram/AmyloGram_results_ed.csv`; columns `Input name` and `Amyloid probability`. |

**Results** (saved to `results/amylogram/`):

| File | Content |
|---|---|
| `amylogram_lollipop.png` | ΔProbability lollipop |
| `amylogram_groups.png` | Grouped bar chart |

---

### 6. PASTA 2.0 — Best Energy & β-Strand Content

**Notebook:** `notebooks/pasta_analysis.ipynb`

**Description:** Reads the PASTA CSV and extracts two metrics per variant: `Best Energy` (kcal/mol; more negative = more stable amyloid) and `% β-Strand`. Computes Δ for each relative to Wildtype. Sign convention: ΔEnergy > 0 means the mutation destabilises the amyloid (less favourable energy); Δβ < 0 means fewer β-strands (lower aggregation propensity). Produces three plots: ΔEnergy lollipop, Δβ-strand lollipop, and a combined scatter.

**Cautions:**

| Item | Detail |
|---|---|
| Decimal separator | The CSV uses commas as decimal separators; the notebook replaces them with periods before parsing (`content.replace(',', '.')`). Do not pre-process the file externally. |
| CSV separator | Semicolon-delimited. |
| `THRESH = 0.1` | Applied independently to ΔEnergy and Δβ-strand for colour-coding. |
| Input path | `../predictions/pasta/pasta.csv` |

**Results** (saved to `results/pasta/`):

| File | Content |
|---|---|
| `pasta_energy.png` | ΔBest Energy lollipop |
| `pasta_beta.png` | Δ% β-strand lollipop |
| `pasta_absolute.png` | Scatter: absolute Best Energy vs % β-strand per variant |
| `vem_pasta.png` | Variant Effect Map for PASTA (generated by `vem.ipynb`, stored here) |

---

### 7. VEM — Variant Effect Map

**Notebook:** `notebooks/vem.ipynb`

**Description:** Constructs two position × amino-acid matrices (one for TANGO Δ Aggregation, one for PASTA Δ Best Energy) using only single-point substitution variants (parsed by regex `^([A-Z])(\d+)([A-Z])`). Each cell shows the predicted effect of substituting the wildtype amino acid at a given position with a given alternative. Named clinical variants (Arctic, Dutch, Iowa, Icelandic) that match the pattern are included automatically.

**Cautions:**

| Item | Detail |
|---|---|
| `WT_SEQ` | Hardcoded as `'DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'`. Must match the sequence used for predictions. |
| `AA_ORDER` | 20 canonical amino acids in alphabetical order by one-letter code: `ACDEFGHIKLMNPQRSTVWY`. Rows with no data remain NaN and are shown as grey cells. |
| Multi-substitution variants | Variants whose names do not match `^([A-Z])(\d+)([A-Z])` are silently skipped (e.g. compound mutants or non-standard naming). |
| Output location | Figures saved directly in `notebooks/`; moved to `results/tango/` and `results/pasta/` respectively. |

**Results:**

| File | Destination |
|---|---|
| `vem_tango.png` | `results/tango/` |
| `vem_pasta.png` | `results/pasta/` |

---

### 8. Amyloid Ranking — Consensus Score

**Notebook:** `notebooks/amyloid_ranking.ipynb`

**Description:** Integrates four tools (TANGO, PASTA, AmyPred-FRL, CrossBeta). For each tool, variants are ranked by their Δ metric; ranks are then averaged into a consensus score. Variants at the top of the consensus ranking are predicted to most consistently increase aggregation across tools; those at the bottom most consistently decrease it. Outputs a consensus lollipop, a top-10 aggregators chart, and a top-10 disruptors chart.

**Cautions:**

| Item | Detail |
|---|---|
| Must run after tools 1–4 | Reads raw prediction files directly — does not depend on figures, but all four prediction CSVs must be present. |
| Tool colours | Hardcoded in `TOOL_COLORS` dict (`TANGO: red`, `PASTA: blue`, `AmyPred-FRL: purple`, `CrossBeta: teal`). |
| Rank direction | For TANGO and AmyPred-FRL, higher Δ = more aggregating (rank ascending). For PASTA, higher ΔEnergy = more destabilising amyloid (rank ascending). For CrossBeta, lower Δavg = more destabilising (rank descending). This sign convention is hardcoded in the loading functions. |
| Waltz and AmyloGram excluded | These tools do not produce a single numerical Δ per variant suitable for ranking (Waltz: region-based; AmyloGram: too discretised). |

**Results** (saved to `results/rank/`):

| File | Content |
|---|---|
| `ranking_consensus.png` | Full consensus lollipop for all variants |
| `ranking_top10_aggregators.png` | Top 10 variants most consistently predicted to increase aggregation |
| `ranking_top10_disruptors.png` | Top 10 variants most consistently predicted to decrease aggregation |

---

### 9. Consensus Heatmap

**Notebook:** `notebooks/consensus_heatmap.ipynb`

**Description:** Builds a (variants × tools) heatmap of normalised Δ values for the same four tools as the ranking step. Variants and tools are hierarchically clustered (Ward linkage on Euclidean distances) so that mutations with similar cross-tool profiles are grouped together. Clinically relevant variants (Arctic E22G, Dutch E22Q, Iowa D23N, Icelandic A2T) are labelled distinctly.

**Cautions:**

| Item | Detail |
|---|---|
| Must run last | Depends on all four prediction files. |
| `CLINICAL` set | Hardcoded set of full variant names used for special label colour. Update if clinical variants are added or renamed. |
| Normalisation | Each tool's Δ column is z-score normalised before clustering and display. Raw values are not shown in the heatmap. |
| Clustering | `scipy.cluster.hierarchy.linkage` with `method='ward'` and `metric='euclidean'`. Changing the method may reorder rows/columns significantly. |

**Results** (saved to `results/rank/`):

| File | Content |
|---|---|
| `consensus_heatmap.png` | Clustered heatmap of normalised Δ values across four tools |