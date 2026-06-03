# Step 1: Mutagenesis and Sequence Preparation

## Requirements
**System Requirements:**
* OS: Windows / WSL (Ubuntu)
* Python version: 3.14

**Software & Packages:**
* **Python Packages:** Built-in standard libraries only (`argparse`, `os`, `sys`). No third-party dependencies are required for this step.

## Description of the process
This step focuses on generating a library of mutant sequences and preparing them for downstream analysis. The workflow consists of two parts:
1. **Mutagenesis:** Generating single point mutations (substitutions) within a specified target consensus region of the Wild-Type (WT) Prion Protein (PrP). The algorithm introduces specific "breaker" amino acids at each position of the target region.
2. **Trimming:** Removing the unstructured and signal peptide regions from both the N- and C-termini of all generated sequences to isolate the structurally relevant core of the protein.

## Command example for execution
All scripts must be executed from the `results/01_mutagenesis/scripts/` subdirectory to ensure relative paths resolve correctly.

### 1. Generating Mutants
Run `01_01_generate_mutants.py` to create the mutation library.
By default, the script uses the breaker amino acids P, D, E, K, and R, and saves the output in the same directory as the input file with the name `prp_sequences.fasta`.

**Standard execution:**
```bash
python 01_01_generate_mutants.py --input "../structures/wt_prion.fasta" --start 170 --end 193
```

**Custom execution (specifying custom breakers and output name):**
```bash
python 01_01_generate_mutants.py --input "../structures/wt_prion.fasta" --start 170 --end 193 --breakers P E --output "custom_mutants.fasta"
```

### 2. Sequence Trimming
Run `01_02_trim_sequences.py` to cut off the N- and C-terminal regions.

Execution:
```bash
python 01_02_trim_sequences.py --input "../structures/prp_sequences.fasta" --start 23 --end 231
```

## Cautions needed to be highlighted
* **Coordinate System:** Both scripts utilize 1-based biological coordinates, meaning you input the exact residue numbers as they appear in the biological sequence (starting at 1, not 0).

* **Trimming Logic:** In `01_02_trim_sequences.py`, the `--start` and `--end` arguments define the region you want to KEEP, not the number of amino acids to remove. For example, passing `--start 23 --end 231` removes residues 1-22 and everything after 231.

* **Input Integrity:** The WT FASTA file for `generate_mutants.py` must contain only a single sequence. The header name inside the file will be ignored, as the script automatically labels the baseline sequence as `>WT`.

## What is gained
Execution of this pipeline populates the `../structures/` subfolder with the following artifacts:

* `prp_sequences.fasta`: A multi-FASTA file containing the full-length WT sequence and all generated single point mutants.

* `prp_sequences_23-to-231.fasta`: The trimmed multi-FASTA file containing the isolated structural core of the proteins. The filename is dynamically generated based on the input coordinates to maintain a clear audit trail. This file is directly used as the input for the downstream Amyloid Prediction step.
