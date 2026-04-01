# Early Parkinson’s Revealed by Unlocking Longitudinal Omics at Population Scale

https://doi.org/10.64898/2026.03.12.26348299

# Abstract
Many diseases begin developing years before symptoms appear1–3, yet biospecimens from these early stages are rarely available. We developed Chronos, a framework that uses privacy-preserving tokenization4 to link archived plasma samples with longitudinal clinical records, enabling the modeling of molecular trajectories across time. Starting with >100 million archived, routine-donation samples from 3 million plasma donors, we assembled a longitudinal Parkinson’s disease cohort and profiled 2,609 samples from 348 cases and 348 matched controls using four proteomics platforms, covering more than 25,000 proteoforms. We reproduced proteomic signatures from clinically-phenotyped cohorts and revealed early, coordinated alterations in a CXCL12, cell adhesion, and integrin signaling network years before the estimated onset of PD. We used protein ratios to predict future diagnosis, achieving a maximum cross-validated area under the curve of 0.76 and replicated the findings in up to 5 independent cohorts. Chronos enables disease detection before clinical manifestation by prioritizing longitudinal molecular changes over symptoms, and provides a general framework to reconstruct chronic and acute disease trajectories from large plasma collections.

# Background

This repository will contain the source code supporting the analysis and findings from the manuscript _Early Parkinson’s Revealed by Unlocking Longitudinal Omics at Population Scale_ currently available as a medRxiv [preprint](https://www.medrxiv.org/content/10.64898/2026.03.12.26348299v1).

# Status

All code artifacts are currently undergoing legal review and will be made available upon completion of the review and publication of the manuscript.

# Access

The Chronos-PD datasets generated and analyzed during the Chronos-PD Phase 1 study will be made available to the research community upon publication. 

A web portal providing interactive access to the data is available at [https://chronos-pd-proteomics-explorer.share.connect.posit.cloud/](https://chronos-pd-proteomics-explorer.share.connect.posit.cloud/)

For further data access inquiries please contact: [chronospd\_data\_access@alkahest.com](mailto:chronospd_data_access@alkahest.com)

# Chronos PD Phase 1 Manuscript - Code Repository

This repository contains the code to generate all figures and tables for the Chronos Parkinson's Disease Phase 1 manuscript.

## Overview

This is the **visualization and analysis code** repository for our manuscript on Parkinson's Disease biomarkers in the Chronos study. Each figure and table in the manuscript has a corresponding R script that can be run to reproduce the output.

**Key characteristics:**

- Self-contained R scripts for each figure/table panel
- File-based data access (no database dependencies)
- Modular design with shared utility functions
- Reproducible outputs in publication-ready formats (PDF, TSV)

## Repository Structure

```

src/
├── 20_analysis/                       # Analysis workflows/pipelines
│   ├── joint_modelling/               # Joint modeling of longitudinal protein and survival outcomes
│   ├── lmm/                           # Linear mixed models screening for differential protein expression
│   ├── non_linear_mixed/              # Non-linear trajectory modeling (spline-based mixed effects)
│   ├── protein_ratio_screen/          # Systematic screening of all pairwise protein ratios
│   ├── sepa/                          # Stepwise Enrichment Pathway Analysis (SEPA)
│   └── target_validation/             # External cohort validation (PPMI, GNPC, UK Biobank)
└── 40_manuscript/                     # Individual figure/table generation scripts
    ├── main_figures/                  # Main figures scripts
    │   ├── figure_1/
    │   ├── figure_2/
    │   ├── figure_3/
    │   └── figure_4/
    ├── main_tables/                   # Main tables scripts
    ├── supplementary_figures/         # Supplementary figure scripts
    │   ├── figure_s1/
    │   ├── figure_s3/
    │   ├── figure_s4/
    │   ├── figure_s5/
    │   ├── figure_s6/
    │   ├── figure_s7/
    │   └── figure_s8/
    └── supplementary_tables/          # Supplementary tables scripts

packages/manuscript/                   # Shared R package with utility functions
└── R/                                 # Function definitions
    ├── create_forest_plot.R
    ├── create_volcano_plot.R
    ├── plot_disease_map.R
    └── ...

_input_data/                           # Input data directory
checkpoint_datasets/                   # Intermediate analysis outputs (git-ignored)
_output/                               # Output directory
config.yaml                            # Configuration (input/output paths, settings)
.Rprofile                              # Auto-loads manuscript package
renv.lock                              # R package dependencies (managed by renv)
```

## Code Organization

This repository contains two types of scripts:

### 1. Analysis Workflows/Pipelines (`src/20_analysis/`)

Complex, multi-step analyses that generate intermediate datasets and multiple outputs. These scripts **have dependencies** and must be run in order.

**Key pipelines:**

- **`lmm/`** - Linear mixed model screening across all proteins and platforms (preDx, postDx, combined)
  - Generates: Differential expression results used by downstream figures
  
- **`non_linear_mixed/`** - Tests whether spline-based models improve fit over linear trajectories
  - Generates: Model selection results (BIC comparisons)
  
- **`joint_modeling/`** - Joint models combining longitudinal protein data with survival outcomes
  - Subdirectories: `PD/` (Parkinson's diagnosis), `comorbidity/` (other conditions)
  
- **`sepa/`** - Pathway enrichment analysis (GO, KEGG, Reactome)
  - Run in order: `SEPA_data_preparation.R` → `SEPA_classic.R` → `PA_heatmap_network.R`
  - Generates: Figure 3c, 3d, Figure S5a, Table S6
  
- **`protein_ratio_screen/`** - Snakemake workflow for systematic screening of all protein pairs
  - See `src/20_analysis/protein_ratio_screen/README.md` for detailed documentation
  
- **`target_validation/`** - Validation in external cohorts (PPMI, GNPC, UK Biobank)

### 2. Manuscript Figure/Table Scripts (`src/40_manuscript/`)

Self-contained scripts that each generate a single figure or table panel. These scripts are **independent** and can be run in any order.

**Characteristics:**
- Each script corresponds to one manuscript output (e.g., `figure_1b.R` → Figure 1, panel b)
- Scripts use `here::i_am()` for location-independent execution
- All input data paths constructed via `get_config()$input_dir`
- Outputs saved automatically to `_out/` with organized subdirectories

## Quick Start

### Prerequisites

**System Requirements:**

  - R version 4.4.3 or higher
  - RStudio (recommended) or command-line R
  - Unix-like environment (Linux/macOS) or WSL on Windows

**R Package Management:**

This project uses `renv` for reproducible package management. All dependencies are specified in `renv.lock`.

### Installation

  1. **Clone the repository:**

   ```bash
   git clone https://github.com/AlkahestInc/chronos_pd_phase_1_manuscript.git
   cd chronos_pd_phase_1_manuscript
   ```

  2. **Restore R package environment:**

   ```r
   # Open R/RStudio in the project directory
   renv::restore()
   ```

   This will install all required packages with the exact versions used for the manuscript.

  3. **Configure data paths:**

   The scripts expect input data in the location specified in `config.yaml`:
   
   ```yaml
   input_dir: "_input_data"
   ```

   **For reviewers:** The full input datasets are not included in this repository due to size and privacy constraints. To reproduce figures, you would need access to the processed data files. Contact the corresponding author for data access requests.

### Running a Figure Script

Each figure script is self-contained and can be run independently:

**Example: Generate Figure 1b**

```bash
# From repository root
Rscript src/40_manuscript/main_figures/figure_1/figure_1b.R
```

Or in RStudio:

```r
# Open and run the script
source("src/40_manuscript/main_figures/figure_1/figure_1b.R")
```

**Output location:**

By default, outputs are saved locally to:
  - Figures: `_out/main_figures/`
  - Tables: `_out/main_tables/` or `_out/supplementary_tables/`

### Running Analysis Pipelines

Analysis pipelines generate intermediate datasets that are used by figure scripts. Run pipeline scripts in the specified order.

**Example: SEPA Pathway Enrichment (generates Figure 3c, 3d)**

```bash
# Must run in this order:
Rscript src/20_analysis/sepa/SEPA_data_preparation.R
Rscript src/20_analysis/sepa/SEPA_classic.R
Rscript src/20_analysis/sepa/PA_heatmap_network.R
```

**Example: Linear Mixed Model Screening**

```bash
# Run for each platform (can be run in parallel)
Rscript src/20_analysis/lmm/alamar.preDx.R
Rscript src/20_analysis/lmm/alamar.postDx.R
Rscript src/20_analysis/lmm/alamar.R
# etc ...
```

**Example: Protein Ratio Screen (Snakemake workflow)**

```bash
cd src/20_analysis/protein_ratio_screen
conda activate ./envs/snakemake

# Run small test
snakemake --profile profiles/local train_models_fixed_tiny

# Run full analysis (requires HPC and HPC profile to be setup)
snakemake --profile profiles/local all
```

**Pipeline outputs:**
- Intermediate results saved to `checkpoint_datasets/` (git-ignored)
- These checkpoints are used as inputs for figure generation scripts

## Code Conventions

This codebase tries follows strict coding standards for reproducibility:

  - **Tidyverse style guide** with base R pipes (`|>`)
  - **Explicit namespacing:** All functions use `package::function()` notation (no `library()` calls)
  - **Reproducible paths:** All paths via `here::here()` or `get_config()` (no `setwd()`)
  - **Idempotent code:** Select columns by name, never by index
  - **Publication quality:** Figures saved at 600 DPI with consistent styling


## Shared Utility Functions

Common visualizations and data processing are abstracted into the `manuscript` package.
The package is automatically loaded by `.Rprofile` when you start R in this project.

## Data Access and Privacy

**Input Data:**

  - Processed datasets are stored in Apache Parquet format
  - Data files are not included in this repository
  - Contact the corresponding author for data access requests

**Output Data:**

  - Generated figures (PDF format, 600 DPI)
  - Generated tables (TSV format, optionally gzipped)
  - Outputs are saved to `_out/` directory (git-ignored)

**Privacy Note:**

  - All patient identifiers have been removed from analysis datasets
  - Only aggregate statistics and de-identified visualizations are generated

## Configuration

The repository is configured via `config.yaml`:

```yaml
input_dir: "_input_data"
output_dir: "_out"
log_level: "INFO"
```

Access configuration in scripts via `get_config()$input_dir`, `get_config()$output_dir`, etc.

## Pre-commit Hooks

The repository uses pre-commit hooks to enforce code quality:

  - **Code formatting:** Automatic styling via `styler`
  - **Code linting:** Static analysis via `lintr`
  - **File size checks:** Blocks large files (>200KB)
  - **Data artifact prevention:** Blocks accidental commits of data files

To install hooks:

```bash
pre-commit install
```

## Troubleshooting

**"Error: No such file or directory"**

  - Check that `config.yaml` points to correct data location
  - Verify input data files exist

**"Package not found"**

  - Run `renv::restore()` to install all dependencies

**"Cannot open connection"**

  - Ensure you have read permissions for input data directory
  - Check that output directory exists and is writable

## Citation

If you use this code, please cite:

**Early Parkinson’s Revealed by Unlocking Longitudinal Omics at Population Scale**

*Chunmiao Feng, Idit Kosti, Yuelong Guo, Ying Wang, Nathan S. Watson-Haigh, Balint File, Nhi Hin, Tibor Nanasi, Jingyu Guo, Rad Suchecki, Rick Tearle, Kathleen Koborsi, Kim Dang, Rashi Saxena, Arnaud Teichert, Shalini Padmanabhan, Brit Mollenhauer, Samuel. M. Goldman, Tony Wyss-Coray, Karoly Nikolich, Scott Lohr, Benoit Lehallier*

medRxiv 2026.03.12.26348299; doi: https://doi.org/10.64898/2026.03.12.26348299
