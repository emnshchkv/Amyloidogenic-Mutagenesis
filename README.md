Repository for the prediction and prioritization of mutations that influence the amyloidogenic properties of Aβ and PrP. Project conducted at the Bioinformatics Institute, 2026.

# Amyloidogenic Mutagenesis Script

This script performs targeted mutagenesis to reduce amyloidogenic potential by identifying enhancer amino acids and replacing them with beta-breaker sequences.

## Features

- **Enhancer Detection**: Automatically finds enhancer amino acids (F, Y, W, V, L, I, Q, N, G) in specified regions
- **Multiple Beta-Breaker Types**: 
  - Single amino acids: R, P
  - Dipeptides: WY, WM
  - Pentapeptides: KLVFF, LPFFD
- **Mutation Strategies**:
  - Single mutations (one enhancer at a time)
  - Combinatorial mutations (multiple enhancers simultaneously)
  - Fixed combinations (all to P, all to R, alternating P/R)
- **Flexible Input**: Accepts both FASTA files and sequence strings
- **Region-Specific**: Target specific regions of interest

## Installation

The script requires Python 3.6 or higher. No additional packages needed - uses only standard library.

```bash
chmod +x amyloid_mutagenesis.py
```

## Usage

### Basic Usage

```bash
# Using a sequence string
python amyloid_mutagenesis.py --sequence "MKVLIVLLIPLASAPTVIGVK" --region "5:10,15:20" --output mutations.fasta

# Using a FASTA file
python amyloid_mutagenesis.py --fasta protein.fasta --region "1:50" --output mutations.fasta
```

### Advanced Usage

```bash
# Multiple regions with custom settings
python amyloid_mutagenesis.py \
  --fasta protein.fasta \
  --region "1:20,35:50,70:85" \
  --output mutations.fasta \
  --no-dipeptides \
  --max-combinations 2

# Only single mutations and fixed combinations
python amyloid_mutagenesis.py \
  --sequence "MKVLIVLLIPLASAPTVIGVK" \
  --region "1:21" \
  --output mutations.fasta \
  --no-dipeptides \
  --no-pentapeptides \
  --no-combinatorial
```

## Command Line Arguments

### Required Arguments

- `--sequence` or `--fasta`: Input sequence (mutually exclusive)
  - `--sequence`: Protein sequence as string
  - `--fasta`: Path to FASTA file
- `--region`: Regions of interest (format: "start:end,start:end") - 1-indexed
- `--output`: Output FASTA file path

### Optional Arguments

- `--no-single`: Exclude single amino acid mutations
- `--no-dipeptides`: Exclude dipeptide insertions  
- `--no-pentapeptides`: Exclude pentapeptide insertions
- `--no-combinatorial`: Exclude combinatorial mutations
- `--no-fixed`: Exclude fixed combination mutations
- `--max-combinations`: Maximum simultaneous mutations (default: 3)

## Example: Basic Amyloid-β Analysis

```bash
# Analyze the amyloid-β peptide
python amyloid_mutagenesis.py \
  --sequence "MKVLIVLLIPLASAPTVIGVK" \
  --region "17:21,30:42" \
  --output Abeta42_mutations.fasta
```

## Output Format

The script generates a FASTA file containing:

1. **Original Sequence**: The input sequence unchanged
2. **Single Mutations**: Each enhancer replaced individually with R or P
3. **Dipeptide Insertions**: Each enhancer replaced with WY or WM
4. **Pentapeptide Insertions**: Each enhancer replaced with KLVFF or LPFFD
5. **Combinatorial Mutations**: Multiple enhancers mutated simultaneously
6. **Fixed Combinations**: All enhancers to P, all to R, or alternating P/R

### FASTA Header Format

- `>Original_Sequence`
- `>Mutant_XXXX_Description`
  - Single: `V4R` (V at position 4 to R)
  - Dipeptide: `F5_WY` (F at position 5 to WY)
  - Combinatorial: `V4R_F5P` (Multiple mutations)
  - Fixed: `ALL_TO_P_V4P_F5P` (All enhancers to P)

## Performance Notes

- Single mutations: 2 × number of enhancers
- Dipeptide insertions: 2 × number of enhancers  
- Pentapeptide insertions: 2 × number of enhancers
- Combinatorial: Exponential with max_combinations setting
- Total mutations can be large for proteins with many enhancers

Use `--max-combinations 2` for large proteins to limit output size.

## Logging

The script generates a detailed log file (`amyloid_mutagenesis.log` by default) containing:
- Sequence and region validation details
- Enhancer detection results
- Mutation generation statistics
- Performance metrics
- Any warnings or errors encountered

Use `--log` to specify a custom log file path and `--verbose` for debug-level console output.

## Error Handling

The script includes comprehensive error handling with specific exit codes:
- Exit code 2: Sequence validation error
- Exit code 3: Region validation error
- Exit code 4: FASTA file reading error
- Exit code 5: Output file writing error
- Exit code 6: Argument validation error
- Exit code 7: General mutagenesis error
- Exit code 130: Script interrupted by user (Ctrl+C)
- Exit code 1: Unexpected error