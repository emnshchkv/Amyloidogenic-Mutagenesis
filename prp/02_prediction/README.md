# Step 2: Amyloidogenicity Prediction and Consensus Scoring

## Requirements
**System Requirements:**
* OS: Windows / WSL (Ubuntu)
* R version: 4.5.1
* Python version: 3.14

**Software & Packages:**
* **TANGO Core Executable:** Must be downloaded independently (e.g., `tango2_3_1.winXP_`).
* **R Packages:** `tidyverse`, `seqinr`, `appnn`, `AmyloGram`, `readxl`.
* **Python Packages:** Built-in standard libraries (`argparse`, `os`, `shutil`, `subprocess`, `sys`).

## Description of the process
During this step, mutant Prion Protein (PrP) sequences generated in the previous mutagenesis step are analyzed to evaluate their amyloidogenic potential compared to the Wild-Type (WT). The pipeline utilizes five different prediction algorithms:
1. **PASTA 2.0** (Energy-based profile)
2. **TANGO** (Statistical mechanics algorithm for cross-beta aggregation)
3. **APPNN** (Neural network-based amyloid propensity predictor)
4. **AmyloGram** (Machine learning predictor based on n-gram analysis)
5. **AGGRESCAN** (In vivo aggregation propensity scale)

Data from these independent predictors are parsed, filtered for the target structural region, and compared against the WT baseline. Finally, a consensus analysis calculates an integrative Z-score for each mutation, allowing for the prioritization of the most potent amyloid breakers.

## Command example for execution
All scripts must be executed from the `scripts/` subdirectory to ensure relative paths resolve correctly.

### 1. TANGO Automation (Python)
The `02_02_Tango_run.py` script requires the path to the TANGO executable. You can run it in one of three ways:

**Option A (Explicit Path via CLI argument):**
```bash
python 02_02_Tango_run.py --tango "C:\Path\To\Tango.exe" --fasta "../data/prp_sequences.fasta" --output "../data/prp_tango"
```

**Option B (Using Windows Environment Variable):**
```DOS
set TANGO_PATH="C:\Path\To\Tango.exe"
python 02_02_Tango_run.py
```

**Option C (If TANGO is in the system PATH):**
```bash
python 02_02_Tango_run.py
```

### R Data Analysis Scripts
Run the R scripts sequentially to process predictor outputs and generate visualizations. Note that due to Windows OS limitations, parallel processing (e.g., mclapply) is not supported, and scripts are executed sequentially.

```bash
Rscript 02_01_PASTA_analysis.R
Rscript 02_02_Tango_analysis.R
Rscript 02_03_APPNN_analysis.R
Rscript 02_04_AmyloGram_analysis.R
Rscript 02_05_AGGRESCAN_analysis.R
Rscript 02_06_Consensus_analysis.R
```

## Cautions needed to be highlighted
To ensure reproducibility and pipeline integrity, several variables are hardcoded within the R scripts. If adapting this pipeline for a different protein or region, these must be manually modified in the configuration block of each script:

* **Relative Paths:** All data extraction and saving operations rely on relative paths (`../data/`, `../figures/`). The scripts will automatically create these directories if they do not exist, but the scripts must be launched from inside the `scripts/` folder.

* `wt_name` **("WT"):** Used as the baseline reference. The naming convention in your input fasta/excel files must strictly match this string.

* `start_pos = 148` and `end_pos = 173`: These represent the structural region of interest (originally PrP 170-195 aa). They are hardcoded as 148-173 because the 22-amino-acid N-terminal signal peptide was removed during the previous structural preparation step, shifting the index.

* `prp_pasta20`: PASTA2.0 is a web tool, before running the respective script ensure that PASTA2.0 output files are unzipped into the `prp_pasta20` folder, no other modifications are needed. Data uploaded to the repository is exemplary.

## What is gained
Execution of this pipeline populates the repository subfolders with the following artifacts:

* `../data/` directory:

    * Raw TANGO calculation outputs (text files).

    * 5 intermediate `.csv` files (`pasta2_results.csv`, `tango_results.csv, etc`. containing the $\Delta$ metrics relative to WT

    * `.consensus_results.csv`: The final integrative dataframe containing aligned Z-scores and the overall Integrative Score for every sequence

* `.../figures/` directory:

    * High-resolution png waterfall plots for each of the 5 individual tools showing mutation effects.

    * `Consensus_Waterfall.png`: Visual ranking of the most effective amyloid breakers based on the mean Z-score.

    * `Consensus_Heatmap.png`: A heatmap of the top 30 candidates, illustrating the contribution and agreement of each specific tool to the final score.
