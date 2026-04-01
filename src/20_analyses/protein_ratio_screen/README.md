# Protein Ratio Screening Analysis

## Overview

**Purpose**: Systematic screening of all pairwise protein ratios to identify biomarkers predictive of Parkinson's disease progression.

**Approach**: For each protein pair (A/B), fit logistic regression models using cross-validation to evaluate predictive performance across different time windows relative to disease onset.

**Output**: AUC performance metrics for all protein ratio combinations across multiple temporal windows.

---

## Table of Contents

1. [Software Dependencies](#software-dependencies)
2. [Quick Start](#quick-start)
3. [Workflow Overview](#workflow-overview)
4. [Configuration](#configuration)
5. [Running the Workflow](#running-the-workflow)
6. [Output Structure](#output-structure)
7. [Post-Analysis](#post-analysis)

---

## Software Dependencies

- Conda/Mamba package manager
- Snakemake workflow management system

---

## Quick Start

### Installation

```bash
git clone https://github.com/AlkahestInc/chronos_pd_phase_1_manuscript.git

# Navigate to workflow directory
cd chronos_pd_phase_1_manuscript/src/20_analyses/protein_ratio_screen

# Install conda (if needed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
conda config --set channel_priority strict

# Create Snakemake environment
conda env create --prefix ./envs/snakemake --file envs/snakemake.yaml
conda activate ./envs/snakemake
```

### Validation Run

```bash
# Test with 3 protein pairs (~5 minutes)
snakemake --profile profiles/local train_models_fixed_tiny --dry-run
snakemake --profile profiles/local train_models_fixed_tiny
```

### Full Production Run

```bash
# Edit config.yaml to select parameter sets (see Configuration section)

# Dry run to verify DAG
snakemake --profile profiles/local all --dry-run

# Execute workflow
snakemake --profile profiles/local all
```

---

## Workflow Overview

### Pipeline Stages

**Stage 1: Data Preparation** → Extract and merge proteomics datasets
- Rules: `get_olink_data`, `get_somalogic_data`, `get_alamar_data`, `extract_assay_ids`
- Outputs: 3 parquet files + assay ID list
- Runtime: Typically ~1-2 minutes

**Stage 2: Train/Test Splitting** → Create temporal stratified splits
- Rule: `initial_split`
- Strategy: Split donors by time window relative to disease onset
- Outputs: `training.parquet` and `testing.parquet` per parameter set
- Runtime: Typically ~1-2 minutes

**Stage 3: Model Training** → Fit all pairwise protein ratio models
- Rule: `fit_model`
- For each denominator protein:
  - Calculate ratios with all numerator proteins
  - Fit logistic regression with V-fold cross-validation
  - Evaluate on held-out test set
- Outputs: Performance metrics per protein pair
- Runtime: weeks to months of CPU time and days to weeks of walltime depending on parallelization

### DAG Visualization

```bash
snakemake --dag all | dot -Tpng > dag.png
```

---

## Configuration

### config.yaml Reference

#### Output Directory

```yaml
out_parent_dir: "output"  # Where all results are stored
```

#### Parameter Sets

Define temporal windows and data splitting strategies:

```yaml
parameter_sets:
  - window_start: -8          # Years before disease onset (start)
    window_end: -6            # Years before disease onset (end)
    strategy: middle          # Splitting strategy: "middle", "earliest", "latest"
    training_proportion: 1.0  # Fraction of data for training (0-1)
```

**Common Configurations:**

- **Pre-diagnostic window**: `window_start: -8, window_end: -6`
- **Recent window**: `window_start: -3, window_end: -1`
- **All pre-diagnostic**: `window_start: -Inf, window_end: -1`
- **Post-diagnostic**: `window_start: 0, window_end: Inf`

**Note**: Negative values indicate years before disease onset; positive values indicate years after.

#### Cross-Validation Settings

```yaml
vfold_cv:
  partitions: 3  # Number of folds (typically 3-10)
  repeats: 3     # Number of repeat iterations
```

**Note**: Total CV iterations = partitions × repeats (e.g., 3×3 = 9 models per protein pair)

#### Job Control

```yaml
jobs:
  max_numerators_per_job: 99999999  # Limit numerators per job (for memory control)
  debug_level: info                 # Log verbosity: debug, info, warn, error

training_data:
  permute: false                    # Set true for negative control (permute labels)

n_test_denominators: 50             # For equally_spaced test target only
```

### Profile Configuration

Profiles control execution environment (e.g. local vs. cluster).

**Local profile** (`profiles/local/config.yaml`):

```yaml
use-conda: true
cores: 2
jobs: 2
resources:
  mem_mb: 12000
show-failed-logs: true
latency-wait: 60
```

---

## Running the Workflow

### Test Targets (Development)

```bash
# 1. Tiny test: 3 small jobs (~1-2 min)
snakemake --profile profiles/local train_models_fixed_tiny

# 2. Small test: 1 small + 1 medium + 1 large job (~5-10 min)
snakemake --profile profiles/local train_models_fixed_small

# 3. Equally spaced: n jobs spanning runtime distribution
snakemake --profile profiles/local train_models_equally_spaced
```

### Production Targets

```bash
# Run Stage 1 only (data preparation)
snakemake --profile profiles/local prepare_data

# Run Stages 1+2 (data + splits)
snakemake --profile profiles/local prepare_splits

# Run all stages (full workflow)
snakemake --profile profiles/local all
```

### Resuming Failed Runs

Snakemake automatically resumes from where it left off:

```bash
# Rerun with same command - only missing outputs will be generated
snakemake --profile profiles/local all
```

---

## Output Structure

### Directory Layout

```
output/
├── assets/                                      # Stage 1 outputs
│   ├── olink.parquet
│   ├── somalogic.parquet
│   ├── alamar.parquet
│   ├── assay_ids.txt                            # List of all protein assays
│   ├── sessionInfo.txt                          # R package versions
│   └── window_start=-8/                         # Stage 2 output: per-parameter set
│       └── window_end=-6/
│           └── strategy=middle/
│               └── train_prop=1.0/              # Proportion of data in training data for initial split
│                   ├── training.parquet         # Initial split for training data
│                   └── testing.parquet          # Initial split for testing data
├── logs/                                        # Logs (STDERR/STDOUT) for model fitting jobs 
│   └── window_start=-8/
│       └── window_end=-6/
│           └── strategy=middle/
│               └── train_prop=1.0/
│                   └── denominator/
│                       ├── OID40001.log         # Log file for denominator=OID40001
│                       ├── seq.9997.12.log      # Log file for denominator=seq.9997.12
│                       └── ...                  # One log file per protein denominator
└── results/                                     # Stage 3 outputs: HIVE-like key=value directory structure
    └── window_start=-8/
        └── window_end=-6/
            └── strategy=middle/
                └── train_prop=1.0/
                    └── denominator/
                        ├── OID40001.parquet     # Results for denominator=OID40001
                        ├── seq.9997.12.parquet  # Results for denominator=seq.9997.12
                        └── ...                  # One result file per protein denominator
```

### Output File Schema

Each `denominator/*.parquet` file contains model performance metrics for all protein ratios using that denominator.

**Key columns**:

| Column                | Type    | Description                                              |
|-----------------------|---------|----------------------------------------------------------|
| `numerator`           | string  | Numerator protein assay ID                               |
| `denominator`         | string  | Denominator protein assay ID (same for all rows in file) |
| `pair_id`             | string  | Unique identifier: "numerator_o_denominator"             |
| `auc_cv_mean`         | double  | Mean AUC across CV folds                                 |
| `auc_cv_sd`           | double  | Standard deviation of AUC across CV folds                |
| `auc_assessment_mean` | double  | Mean AUC on test set                                     |
| `auc_assessment_sd`   | double  | SD of AUC on test set                                    |
| `n_folds`             | integer | Total CV folds attempted                                 |
| `n_folds_ok`          | integer | CV folds with successful model fit                       |

**Quality Filtering**: Consider filtering to `n_folds_ok == n_folds` for reliable results. Rows where `n_folds_ok < n_folds` indicate model fitting issues (e.g., convergence failure, perfect separation).

---

## Post-Analysis

### Loading Results in R

The workflow produces parquet files organized in a hive-like directory structure, allowing efficient querying with Apache Arrow.

#### Basic Loading

```r
# Load a specific parameter set
results <- arrow::open_dataset(
  "output/results/window_start=-8/window_end=-6/strategy=middle/train_prop=1.0/denominator/",
  format = "parquet"
)

results |>
  dplyr::filter(
    auc_assessment_mean >= 0.7
  ) |>
  dplyr::collect()
```

#### Loading Multiple Parameter Sets

```r
# Load all parameter sets with automatic partitioning
results <- arrow::open_dataset(
  "output/results/",
  format = "parquet",
  hive_style = TRUE,
  partitioning = arrow::hive_partition(
    window_start = arrow::int32(),
    window_end = arrow::int32(),
    strategy = arrow::string(),
    train_prop = arrow::float64()
  )
)

results |>
  dplyr::filter(
    window_start == -8,
    auc_assessment_mean >= 0.7
  ) |>
  dplyr::collect()
```
