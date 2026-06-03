# Prediction and Prioritisation of Mutations Influencing the Amyloidogenic Properties of Aβ and PrP

A single computational pipeline — *in silico* saturation mutagenesis → multi-tool consensus amyloidogenicity ranking → all-atom molecular dynamics of selected variants → trajectory analysis — applied across three amyloidogenic substrates to test whether sequence-based aggregation predictors translate into measurable conformational change.

**Authors:** _Elizaveta Menshikova_ (Aβ42), _Sergey Ilin_ (Aβ39), _Daniil Spirin_ (PrP)
**Supervisors:** _Sukhanova Xenia_
**Institution:** Bioinformatics Institute, 2025–2026

---

## Contents

- [Background](#background)
  - [Goal](#goal)
  - [Objectives](#objectives)
- [Methods](#methods)
- [Workflow steps](#workflow-steps)
- [The three subprojects](#the-three-subprojects)
- [System requirements](#system-requirements)
- [Getting the data](#getting-the-data)
- [Repository structure](#repository-structure)
- [Key results](#key-results)
- [Conclusions](#conclusions)
- [References](#references)

---

## Background

Amyloid fibrils form when a polypeptide populates a cross-β architecture stabilised by a steric-zipper backbone hydrogen-bond network. Point substitutions shift the free-energy balance between the soluble ensemble and the aggregation-prone state, which is why a handful of single mutations in Aβ (e.g. the Arctic E22G, Dutch E22Q, Iowa D23N) and in PrP convert a normally clearing peptide into a clinically aggressive one. Sequence-based predictors (TANGO, PASTA 2.0, WALTZ, AmyloGram, AGGRESCAN, AmyPred-FRL, Cross-β) each encode a different physical proxy for that balance — β-sheet propensity, pairing free energy, hydrophobic-pattern matching — but no single tool is reliable alone, and a prediction of "more aggregation-prone" is not the same as an observed structural change.

This project asks whether the **consensus** of several predictors, used to *prioritise* mutations, picks out substitutions that actually perturb the fold under explicit-solvent molecular dynamics.

### Goal

To establish and validate a reusable computational pipeline that prioritises amyloidogenicity-modulating point mutations and tests the top candidates against all-atom MD, demonstrated across structurally distinct amyloidogenic systems.

### Objectives

1. Generate the complete (or region-restricted) single-substitution variant library for each target.
2. Score every variant with multiple predictors and combine them into a per-tool z-score consensus ranking.
3. Select top aggregators / disruptors spanning different sequence regions and effect magnitudes.
4. Simulate wild-type and selected variants under a common CHARMM36m protocol.
5. Quantify the structural consequences (RMSD, RMSF, Rg, DSSP, H-bond/salt-bridge networks, contact maps, PCA) and relate them back to the predicted ranking.

## Methods

![Graphical abstract of the workflow](docs/images/graphical_abstract.png)

> _Replace `docs/images/graphical_abstract.png` with a single figure showing the four steps feeding three substrates. This is the figure the reviewers asked for; it is the one thing that makes the "one project, three systems" structure legible at a glance._

The same four-step workflow is applied to every substrate; only the starting structure, the mutated region, and the oligomeric state differ.

## Workflow steps

| Step | What happens | Output |
|---|---|---|
| **01 — Mutagenesis** | Enumerate single substitutions (saturation, or restricted to the amyloidogenic core) and write FASTA. | variant FASTA |
| **02 — Prediction** | Submit variants to the predictor panel; normalise each tool to a within-tool z-score; average into a consensus; rank. | consensus table + ranking figures |
| **03 — MD run** | Build WT + selected variants, energy-minimise, equilibrate (NVT/NPT), production MD in explicit solvent (CHARMM36m, TIP3P). | trajectories |
| **04 — MD analysis** | RMSD, RMSF, Rg, DSSP, H-bonds, salt bridges, contact maps, PCA, convergence; WT-vs-mutant comparison. | analysis figures |

## The three subprojects

The substrates were chosen to stress the pipeline on deliberately different structural regimes:

| Subproject | System | Mutated region | Oligomeric state simulated | Selected variants |
|---|---|---|---|---|
| [`abeta42/`](abeta42/) | Aβ42 | full sequence (saturation, filtered to known variants) | pentameric β-sheet assembly | F20L, G33R, H14R |
| [`abeta39/`](abeta39/) | Aβ39 | amyloidogenic regions | monomer **and** tetramer | G25D (+ WT) |
| [`prp/`](prp/) | PrP | C-terminal core, residues 170–195 | monomer (C-terminal core 23–231) | V180R |

Each subproject directory contains its own README with the substrate-specific rationale and per-step documentation.

## System requirements

**Software**

| Component | Version | Used in |
|---|---|---|
| GROMACS | 2024.3 (CUDA 12.8) | steps 03–04 |
| CHARMM36m / CGenFF | charmm36-feb2026, cgenff-5.0 | step 03 |
| Python | ≥ 3.10 | steps 01, 02, 04 |
| R | ≥ 4.2 | step 02 (PrP) |
| Predictors | TANGO, PASTA 2.0, WALTZ, AmyloGram, AGGRESCAN, AmyPred-FRL, Cross-β | step 02 |

Exact Python dependencies are pinned in [`requirements.txt`](requirements.txt).

**Hardware**

GPU strongly recommended for production MD (simulations were run with CUDA 12.8). Prediction and analysis steps run comfortably on a CPU-only workstation.

## Getting the data

To keep the repository lightweight, **large binary assets are not stored in git** (see [Repository structure](#repository-structure) for what this excludes and why). They live in a shared Google Drive folder and are fetched with a script:

```bash
pip install gdown
python scripts/download_data.py            # everything
python scripts/download_data.py --only prp # one subproject
```

The CHARMM36m force field is a standard external asset; the script can also fetch it, but you may already have it shipped with your GROMACS install.

## Repository structure

```
.
├── README.md                  ← you are here (the unified story)
├── requirements.txt
├── pyproject.toml             ← ruff config (lint + quote/import normalisation)
├── .pre-commit-config.yaml    ← ruff + nbstripout, runs on commit
├── .gitignore                 ← excludes heavy data (see below)
├── scripts/
│   └── download_data.py       ← pulls heavy assets from Google Drive
├── docs/
│   └── images/                ← small figures for this README only
├── abeta42/
│   ├── README.md
│   ├── 01_mutagenesis/   { README.md, scripts/, data/ }
│   ├── 02_prediction/    { README.md, scripts/, notebooks/, data/, figures/ }
│   ├── 03_md_run/        { README.md, scripts/(mdp), structures/ }
│   └── 04_md_analysis/   { README.md, notebooks/, figures/ }
├── abeta39/   (same skeleton)
└── prp/       (same skeleton)
```

**Not tracked in git** (fetched via `download_data.py`):

| Asset | Why excluded |
|---|---|
| `charmm36-feb2026_cgenff-5.0.ff/` | ~15 MB standard external force field, was duplicated in every branch |
| `*.xtc`, `*.tpr`, `*.gro`, `*.pdb` (trajectories/structures) | large, regenerable from the MD protocol |
| `*.gif` trajectory animations | tens of MB each |
| notebook cell outputs | stripped on commit; figures live in `figures/` instead |

## Key results

_2–3 sentences per substrate with the headline finding and a link to the subproject. Example:_

- **Aβ42** — Of the three simulated variants, H14R showed progressive Cα RMSD drift after ~140 ns (→ ~13 Å) consistent with continuous remodelling of the pentamer, while F20L was structurally near-neutral. See [`abeta42/`](abeta42/).
- **Aβ39** — _<headline finding for G25D, monomer vs tetramer>_. See [`abeta39/`](abeta39/).
- **PrP** — V180R, the top predicted disruptor, _<observed structural effect on the H2–H3 region>_. See [`prp/`](prp/).

## Conclusions

_2–4 sentences: did the consensus ranking pick out variants that measurably perturbed the structure under MD? Where did prediction and simulation agree, and where did they diverge? What does that say about using sequence-based consensus to prioritise candidates for the (far more expensive) MD step?_

## References

_Numbered list — predictor papers (TANGO, PASTA 2.0, WALTZ, AmyloGram, AGGRESCAN, AmyPred-FRL, Cross-β), CHARMM36m, GROMACS, MDAnalysis, and the clinical-mutation sources (UniProt variant records)._
