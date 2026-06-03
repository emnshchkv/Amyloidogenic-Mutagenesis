# MD Analysis Pipeline — Step Documentation

Full data is avaible at:
https://drive.google.com/drive/folders/1QIL2TLIcyQhrr-qpx6V0-1W6LD2xOE7r?usp=sharing

---

## Step 1: Aβ42 Mutant Pentamer Analysis

### Requirements

Install dependencies:

```bash
pip install MDAnalysis matplotlib numpy seaborn scipy
```

### Description of the Process

This step performs a comprehensive molecular dynamics (MD) analysis of an Aβ42 mutant pentamer. The notebook loads a GROMACS trajectory and runs the following analyses sequentially:

1. **RMSD** — Root mean square deviation of Cα, backbone, and all-atom selections relative to the starting structure; tracks global structural stability over time.
2. **RMSF** — Per-residue root mean square fluctuation of Cα atoms; identifies flexible and rigid regions. Mutation positions in all five monomers are highlighted.
3. **Radius of gyration (Rg)** — Compactness metric calculated for every frame; displayed as a time series and distribution.
4. **DSSP** — Secondary structure (helix/sheet/coil) propensity per residue, averaged over all trajectory frames.
5. **Cα–Cα distance map** — Average inter-residue distance matrix sampled every 10th frame, with monomer boundaries annotated.
6. **Hydrogen bonds** — Backbone N–H···O hydrogen bonds counted per frame (cutoff 3.5 Å, angle ≥ 120°).
7. **PCA** — Principal component analysis on Cα coordinates; cumulated variance plot and PC1 vs PC2 projection.
8. **Ramachandran plot** — φ/ψ dihedral angle distribution for the full ensemble.
9. **Convergence assessment** — Running average of Rg, block averaging, and first-half vs. second-half comparison.
10. **Summary** — Printed statistics for all metrics; mutation characterisation (F20L context hardcoded in print statements).

### Execution

```bash
bash run_mutant_analysis.sh
```

Where `run_mutant_analysis.sh` contains:

```bash
#!/bin/bash
# Run from the mutant simulation directory (where md.gro and md_clean.xtc reside)
echo "20" | jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=3600 \
    abeta42_mutant_analysis.ipynb \
    --output abeta42_mutant_analysis_executed.ipynb
```

> The `echo "20"` pipes the mutation position (residue 20, F20L) to the `input()` call inside the notebook. Adjust the value if analysing a different mutation site.

### Cautions

| Item | Detail |
|---|---|
| `GRO_FILE = "./md.gro"` | Hardcoded — must be run from the directory containing the trajectory files, or the path must be updated. |
| `XTC_FILE = "./md_clean.xtc"` | Expects a **solvent-stripped** trajectory. Running on a raw (solvent-included) XTC will produce incorrect RMSF/distance maps and is much slower. |
| `MUT_POS = int(input())` | Interactive prompt — must be piped or entered manually. For the F20L mutant use `20`. Changing this value shifts all highlighted mutation positions in the plots. |
| `monomer_size = 42` | Hardcoded in the distance-map cell. Change if the peptide length differs. |
| Distance map sampling `[::10]` | Every 10th frame is used for the Cα–Cα map to reduce compute time. Increase sampling density (e.g. `[::5]`) for higher accuracy at the cost of runtime. |
| `d_a_cutoff=3.5`, `d_h_a_angle_cutoff=120` | H-bond geometric criteria. Standard values; modify only if using a non-standard force field. |

### What Is Gained

Nine publication-quality PNG figures (300 dpi) saved in the **current working directory**:

| File | Content |
|---|---|
| `01_rmsd.png` | RMSD time series (Cα, backbone, all-atom) |
| `02_rmsf.png` | Per-residue RMSF with mutation site markers |
| `03_radius_of_gyration.png` | Rg time series + distribution |
| `04_secondary_structure.png` | DSSP propensity stacked bar chart |
| `05_ca_distance_map.png` | Average Cα–Cα distance map |
| `06_hbonds.png` | Backbone H-bond count time series + distribution |
| `07_pca_analysis.png` | PCA cumulated variance + PC1 vs PC2 scatter |
| `08_ramachandran.png` | Ramachandran plot for full ensemble |
| `09_convergence_assessment.png` | Running average, block averaging, half-trajectory comparison |

A printed summary with mean ± std for RMSD, Rg, RMSF, secondary structure fractions, and convergence status is also produced in the notebook output.

---

## Step 2: WT vs Mutant Comparison

### Requirements

- Same Python packages as Step 1
- **Two** sets of GROMACS trajectory files:
  - Wild-type: `../md_wt/md.gro` and `../md_wt/md_clean.xtc`
  - Mutant: `./md.gro` and `./md_clean.xtc`
- Two interactive inputs at runtime:
  1. Mutation position (`MUT_POS`, integer)
  2. Mutation label (`MUT_LABEL`, string, e.g. `F20L`)

### Description of the Process

This step loads both the wild-type (WT) and mutant trajectories and performs a **side-by-side comparative analysis** across seven metrics. For each metric, statistical significance is assessed with an independent-samples t-test:

1. **RMSD** — Time series and distribution overlaid for WT vs mutant; Δ mean and t-test p-value reported.
2. **RMSF** — Per-residue profiles averaged over all five monomers (mean ± SEM). A ΔRMSF bar chart (mutant − WT) highlights regions of altered flexibility around the mutation site.
3. **Radius of gyration** — Overlaid time series and distributions; t-test on the full Rg vectors.
4. **DSSP** — Stacked bar plots for each system separately, followed by ΔHelix and ΔSheet difference charts.
5. **Contact maps** — Per-monomer average Cα–Cα distance maps for WT and mutant on a shared colour scale, plus a difference map (Δ distance, RdBu palette).
6. **Hydrogen bonds** — Backbone H-bond time series and distributions overlaid; t-test statistics.
7. **PCA** — Cumulated variance curves and PC1 vs PC2 scatter plots coloured by time (Blues for WT, Purples for mutant).

### Execution

```bash
bash run_comparison.sh
```

Where `run_comparison.sh` contains:

```bash
#!/bin/bash
# Run from the mutant simulation directory
printf "20\nF20L\n" | jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=7200 \
    comparison_wt_vs_mutant.ipynb \
    --output comparison_wt_vs_mutant_executed.ipynb
```

> `printf "20\nF20L\n"` supplies both `input()` calls in order: mutation residue number, then the label string. Adjust for other mutations (e.g. `"3\nE3K\n"`).

### Cautions

| Item | Detail |
|---|---|
| `GRO_WT = "../md_wt/md.gro"` | Hardcoded relative path — the notebook **must** be run from the mutant directory, with the WT simulation one level up in `md_wt/`. |
| `GRO_MUT = "./md.gro"` | Expects mutant files in the current directory. |
| `N_MON = 5`, `N_RES = 42` | Hardcoded pentamer dimensions. Update both constants if the oligomeric state or peptide length changes. |
| `MUT_LABEL` (second `input()`) | Used as a filename suffix and plot label. Avoid spaces or special characters (e.g. use `F20L` not `F→L20`). |
| `C_WT = "#2E86AB"`, `C_MUT = "#6A0572"` | Colour scheme hardcoded at the top of the notebook. Change here if different colours are needed for publication. |
| Distance map sampling `[::10]` | Same caveat as Step 1 — every 10th frame sampled for performance. |
| t-test assumption | `scipy.stats.ttest_ind` assumes normally distributed, independent samples. For very short trajectories this may not hold; consider Mann–Whitney U as an alternative. |
| Timeout `7200 s` | The comparison notebook processes two full trajectories; 2 hours is a safe default for ~500 ns simulations. Adjust `--ExecutePreprocessor.timeout` for longer runs. |

### What Is Gained

Fourteen publication-quality PNG figures (300 dpi) saved in the **mutant working directory**, all named with the pattern `NN_metric_wt_vs_<MUT_LABEL>.png`:

| File | Content |
|---|---|
| `01_rmsd_wt_vs_f20l.png` | RMSD time series + distribution, WT vs mutant |
| `02_rmsf_wt_vs_f20l.png` | Per-residue RMSF profiles (mean ± SEM) |
| `02b_delta_rmsf_wt_vs_f20l.png` | ΔRMSF bar chart (mutant − WT) |
| `03_rg_wt_vs_f20l.png` | Rg time series + distribution |
| `04_dssp_wt_vs_f20l.png` | DSSP stacked bars for WT and mutant |
| `04b_delta_dssp_wt_vs_f20l.png` | ΔHelix and ΔSheet difference charts |
| `05_contact_maps_wt_vs_f20l.png` | Cα–Cα distance maps + difference map |
| `06_hbonds_wt_vs_f20l.png` | H-bond time series + distribution |
| `07_pca_wt_vs_f20l.png` | PCA cumulated variance + PC1/PC2 scatter |

> File names use the lowercased `MUT_LABEL` as a suffix; for `F20L` they will be `..._f20l.png`.

Printed statistics for each metric include mean ± std for both systems, Δ mean, and t-test results (with significance stars: `***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `ns`).
