# Mutagenesis and Predictions

## Overview

The `Mutagenesis_and_Predictions/` directory is dedicated to **in silico mutagenesis** and **computational prediction** of amyloidogenic properties for Aβ (Amyloid-β) peptides. It integrates sequence-based prediction tools and custom scripts to generate mutant variants and rank them according to their predicted aggregation propensity.

This module supports the identification of mutations that may enhance or suppress the amyloidogenic potential of the peptide, complementing the all-atom molecular dynamics simulations performed in the `MD/` directory.

## Directory Structure

- **`scripts/`** — Core Python scripts for mutagenesis and ranking:
  - `make_mutants.py` — Generates FASTA-file with WT sequence and all possible single point mutants in amyloidogenic regions.
  - `Aggregation_ranking.py` — Aggregates predictions from multiple amyloidogenicity tools and produces a final ranking of mutants.

- **`tools_assessment/`** — Contains input and output files from external amyloid prediction tools:
  - `mutants.fasta` — FASTA sequences of generated mutants.
  - Prediction results from **TANGO**, **PASTA**, **AmyPred**, and **CrossBeta** (`.csv` and `.txt` files).
  - These files serve as intermediate data for consensus scoring.

## Workflow and Usage

The pipeline follows a clear two-step process:

### 1. Mutant Generation
```bash
cd Mutagenesis_and_Predictions/scripts
python make_mutants.py
```

### 2. Aggregation Propensity Ranking
```bash
cd Mutagenesis_and_Predictions/scripts
python Aggregation_ranking.py
```

This script performs the following steps:

- Parses output files from multiple amyloid aggregation prediction algorithms (TANGO, PASTA, AmyPred, CrossBeta).
- Normalizes and integrates the scores into a unified consensus metric.
- Produces a ranked list of mutants according to their predicted amyloidogenic potential.

The resulting images can be seen in the ../images folder.
