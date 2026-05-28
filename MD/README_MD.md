# Molecular Dynamics Simulations

## Overview

The `MD/` directory contains all necessary files and scripts to perform all-atom molecular dynamics (MD) simulations of **wild-type (WT)** and **G25D mutant** Aβ39 peptides in both monomeric and tetrameric forms using **GROMACS 2024.3**.

## Directory Structure

- **`mdps/`** — Contains standardized GROMACS parameter files (`.mdp`) for all simulation stages:
  - `minim.mdp` — Energy minimization
  - `ions.mdp` — Ion addition
  - `nvt.mdp` — NVT equilibration
  - `npt.mdp` — NPT equilibration
  - `md.mdp` — Production MD run

- **`MD_monomers/`** and **`MD_tetramers/`** — Simulation directories for monomers and tetramers respectively. Each contains subfolders for WT and G25D variants.

- **`charmm36-feb2026_cgenff-5.0.ff`** — CHARMM36 force field (February 2026 release) with CGenFF parameters.

## Running Simulations via Makefile

Each simulation variant (WT/G25D, monomer/tetramer) includes an identical **`Makefile`** that automates the entire GROMACS workflow.

### Main Targets

- `make setup` — Performs full system preparation (PDB → topology → solvation → ionization → energy minimization → NVT → NPT equilibration) up to the production `.tpr` file.
- `make run_md` — Launches the production molecular dynamics simulation in a `screen` session.
- `make process` — Post-processing of the trajectory (PBC correction, centering on protein, compact unit cell).
- `make protein_only` — Creates protein-only trajectory and topology (removes solvent and ions).
- `make all` — Runs `setup` + `run_md` (recommended for initial launch).
- `make clean` — Removes temporary and intermediate files.

### Usage Example

```bash
cd MD/MD_monomers/WT
make -j1 setup          # Prepare the system
make run_md             # Start production MD
make process            # Process trajectory
```

# Requirements
- GROMACS 2024.3 compiled with CUDA 12.8 support (recommended for GPU acceleration)
- CHARMM36 force field (included)

All simulations use consistent parameters defined in the mdps/ directory, ensuring reproducibility across all systems.
