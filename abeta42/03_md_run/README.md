# Molecular dynamics with GROMACS and CHARMM36m: a complete protocol

**This guide provides a step-by-step protocol for simulating the amyloid-beta 42 (Aβ42) peptide pentamer using GROMACS with the CHARMM36m force field and CHARMM-modified TIP3P water on CPU-only hardware.**

---

## 0. GROMACS Installation

Set up a conda environment with GROMACS 2026.0 (CPU-only, no MPI) using `mamba` for faster dependency resolution.

```bash
conda install -n base -c conda-forge mamba
```

Create a new environment and install GROMACS:

```bash
mamba create -n gromacs_cpu -c conda-forge python=3.12 "gromacs=2026.0=nompi*"
```

Activate the environment:

```bash
conda activate gromacs_cpu
```

---

## 1. System preparation: from PDB to topology

### Generating the topology with pdb2gmx

The clean PDB files are ready for topology generation and are stored at `structures`. The CHARMM36m force field is available at the `charmm36-feb2026_cgenff-5.0.ff` directory. Place the `.ff` directory in your working directory.

```bash
gmx pdb2gmx -f protein.pdb -o processed.pdb -water tip3p -ignh -his
```

The command produces three outputs: `processed.gro` (coordinates), `topol.top` (system topology), and `posre.itp` (position restraint parameters)

---

## 2. Box setup and solvation for a disordered peptide

**Recommended approach:** Use at least **1.2 nm padding** from any protein atom to the box edge, but critically, consider the fully extended conformation, not just the initial compact structure.Use a **rhombic dodecahedron**, which saves ~29% solvent volume compared to a cubic box with the same minimum image distance:

```bash
gmx editconf -f processed.pdb -o boxed.pdb -c -d 1.2 -bt dodecahedron
```

### Solvation

```bash
gmx solvate -cp boxed.pdb -cs spc216.gro -o solvated.pdb -p topol.top
```

### Adding ions for charge neutralization and physiological ionic strength

Files with minimal parameters for genion are stored at `mdp`.

```bash
gmx grompp -f ions.mdp -c solvated.pdb -p topol.top -o ions.tpr -maxwarn 2
gmx genion -s ions.tpr -o solvated_ions.pdb -p topol.top \
  -pname NA -nname CL -neutral -conc 0.15
```

When prompted, **select the SOL group** (solvent molecules that will be replaced by ions).

**Always verify** that the `[molecules]` section of `topol.top` has been updated correctly with the reduced SOL count and added NA/CL entries.

---

## 3. Energy minimization: removing bad contacts

Energy minimization relaxes steric clashes and unfavorable geometry introduced during system construction. Steepest descent is the standard algorithm — it is robust, guaranteed to move downhill, and does not require second derivatives.

**Run the minimization:**

```bash
gmx grompp -f minim.mdp -c solvated_ions.pdb -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em -nt 2
```

---

## 4. NVT equilibration: stabilizing temperature

The NVT (canonical ensemble) step equilibrates the system temperature while holding volume constant. Position restraints on protein heavy atoms allow the solvent to relax around the peptide without disturbing its initial conformation.

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt -nt 2
```

---

## 5. NPT equilibration: stabilizing pressure and density

After temperature is stable, pressure coupling is activated to equilibrate the system density.

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt -nt 2
```

---

## 6. Production MD: sampling the conformational ensemble

The production run removes position restraints and lets the peptide evolve freely.

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2
gmx mdrun -deffnm md -nt 2
```

### How long to simulate Aβ42

Aβ42 is an IDP that samples a vast conformational landscape. **Conventional MD of 100 ns is insufficient** for meaningful conformational ensemble characterization. Literature recommendations for the Aβ42 monomer:

- **Minimum practical**: 500 ns for basic structural characterization
- **Recommended**: Multiple independent trajectories of ≥500 ns each (e.g., 3–5 × 500 ns), totaling several microseconds of aggregate sampling
- **Gold standard**: Replica exchange MD (REMD) with ≥100 ns per replica across 24–64 temperature replicas
- **Convergence check**: Compare calculated observables (Rg distribution, secondary structure content, chemical shifts) with experimental data

For a master's project on CPU-only hardware, **3 × 500 ns independent trajectories** is a practical and publishable target. This produces 1.5 μs of aggregate sampling and allows estimation of statistical uncertainty.

---

## 8. Preparation for trajectory analysis

Before any analysis, correct for periodic boundary conditions:

```bash
# Make molecules whole, remove jumps, center protein
gmx trjconv -f md.xtc -s md.tpr -pbc whole -o step1.xtc
gmx trjconv -f step1.xtc -s md.tpr -pbc nojump -o step2.xtc
gmx trjconv -f step2.xtc -s md.tpr -pbc mol -center -ur compact -o md_clean.xtc
# Select "Protein" for centering, "System" for output
```

---