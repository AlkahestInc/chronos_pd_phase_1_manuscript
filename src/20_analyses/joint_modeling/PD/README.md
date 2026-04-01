# Joint Modeling (JM) Analysis for Time-to-PD Prediction

This folder contains R scripts for performing Joint Modeling analysis to identify protein biomarkers whose longitudinal trajectories predict future Parkinson's Disease (PD) diagnosis.

> **Note**: Scripts process one protein at a time and are designed for array job execution (SLURM/HPC) or sequential local processing. Each script uses file-based data access and centralized functions from `packages/manuscript/`.

## Overview

The Joint Modeling (JM) approach combines:

- **Longitudinal data**: Protein measurements from plasma samples collected over time
- **Survival data**: Time-to-event data for PD diagnosis

The analysis uses Bayesian joint models implemented with JMbayes2 to identify proteins whose changing levels are associated with time to PD diagnosis.

**Summary by Platform:**
- **SomaScan**: ~400+ proteins (filtered by q_val < 0.05 from screen model)
- **Olink**: ~100+ proteins (filtered by q_val < 0.05 from screen model)
- **Alamar**: 14 proteins (filtered by q_val_all < 0.05 from screen model)

## Analysis Workflow

The analysis workflow for time-to-PD prediction:

```
1. Data preparation (completed in prior analysis steps)
   → Screen model results with pre-diagnosis samples
   ↓
2. [Platform].R scripts (somalogic.R, olink.R, alamar.R)
   → Process one protein at a time (specified by protein_index)
   → Filter proteins by screen model significance (q_val < 0.05)
   → Fit Joint Models for each protein
   → Outputs to get_config()$output_dir (somalogic) or get_config()$checkpoint_dir (olink, alamar)
```

---

## Folder Structure

```
PD/
├── README.md          # This documentation file
├── somalogic.R        # JM analysis for SomaScan platform (~400+ proteins)
├── olink.R            # JM analysis for Olink platform (~100+ proteins)
└── alamar.R           # JM analysis for Alamar platform (14 proteins)
```

---

## Main Analysis Scripts

### 1. somalogic.R

**Purpose**: Fit Joint Models for SomaScan proteomics data to predict time-to-PD diagnosis.

**Analysis Logic**:

- Joint Modeling functions auto-loaded from `packages/manuscript/` via `.Rprofile`
- Load critical variables and RWD summary data
- Create sample-level data with PD onset times
- Filter to include only pre-diagnosis samples:
  ```r
  dplyr::filter(sample_collection_date < Event_date)
  ```
- Load SomaScan proteomics data and screen model results
- Filter proteins by screen model significance (q_val < 0.05)
- For the specified protein index:
  - Create `joint_data` object with longitudinal and survival components
  - Set up Joint Model specification:
    - **Longitudinal submodel**: Mixed-effects model (LME) with random intercept and slope
      - Formula: `protein ~ gender + cluster + recent_plasma_donations + age_protein`
      - Random effects: `~ 1 + age_protein | patient_id`
    - **Survival submodel**: Cox proportional hazards model
      - Formula: `Surv(age_entry, Time_age, Status) ~ gender`
      - Entry time: `age_entry` (age at first medical record)
      - Event time: `Time_age` (age at PD diagnosis or censoring)
    - **Association**: Value-based (current protein level predicts PD diagnosis risk)
  - Configure MCMC settings (4 chains, 30,000 iterations, 3,000 warmup, thin=2)
  - Fit joint model using JMbayes2 package
  - Save full results and summary to output directory

**Input Files** (from `get_config()$input_dir` in versioned parquet/RDS format):

- `name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet` - Critical variables with PD onset dates
- `name=phase1_claim_summary/version=none/dataset.parquet` - RWD summary with first/last medical record dates
- `name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet` - Screen model results (filtered by q_val < 0.05)
- `name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds` - SomaScan proteomics data

**Output Files** (to `get_config()$output_dir/joint_modeling/PD/`):

- `somalogic.screening.protein.[protein_index].rds` - Complete model object
- `somalogic.screening.protein.[protein_index].summary.rds` - Model summary statistics

**How to Run:**

```bash
# Modify protein_index in script (line 6), then:
cd ~/chronos-pd_phase1_ms
Rscript src/20_analyses/joint_modeling/PD/somalogic.R

# For SLURM array job (recommended for multiple proteins):
# sbatch --array=1-400 run_somalogic.slurm
```

**Features**:
- One protein per run (set via `protein_index` variable)
- Designed for parallelization via SLURM array jobs
- MCMC parallelization: 4 chains across available cores
- Outputs: Full model object + summary statistics

**Number of proteins**: ~400+ proteins (filtered by q_val < 0.05 from screen model)

---

### 2. olink.R

**Purpose**: Fit Joint Models for Olink proteomics data to predict time-to-PD diagnosis.

**Analysis Logic**: Same as somalogic.R, but for Olink platform.

**Input Files** (from `get_config()$input_dir` in versioned parquet/RDS format):

- `name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet` - Critical variables with PD onset dates
- `name=phase1_claim_summary/version=none/dataset.parquet` - RWD summary with first/last medical record dates
- `name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet` - Screen model results (filtered by q_val < 0.05)
- `name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds` - Olink proteomics data

**Output Files** (to `get_config()$checkpoint_dir/joint_modeling/PD/`):

- `olink.screening.protein.[protein_index].rds` - Complete model object
- `olink.screening.protein.[protein_index].summary.rds` - Model summary statistics

**How to Run:**

```bash
# Modify protein_index in script (line 6), then:
cd ~/chronos-pd_phase1_ms
Rscript src/20_analyses/joint_modeling/PD/olink.R

# For SLURM array job (recommended for multiple proteins):
# sbatch --array=1-100 run_olink.slurm
```

**Features**:
- One protein per run (set via `protein_index` variable)
- Designed for parallelization via SLURM array jobs
- MCMC parallelization: 4 chains across available cores
- Outputs: Full model object + summary statistics

**Number of proteins**: ~100+ proteins (filtered by q_val < 0.05 from screen model)

---

### 3. alamar.R

**Purpose**: Fit Joint Models for Alamar proteomics data to predict time-to-PD diagnosis.

**Analysis Logic**: Same as somalogic.R, but for Alamar platform.

**Input Files** (from `get_config()$input_dir` in versioned parquet/RDS format):

- `name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet` - Critical variables with PD onset dates
- `name=phase1_claim_summary/version=none/dataset.parquet` - RWD summary with first/last medical record dates
- `name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet` - Screen model results (filtered by q_val_all < 0.05)
- `name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29464/dataset.rds` - Alamar proteomics data

**Output Files** (to `get_config()$checkpoint_dir/joint_modeling/PD/`):

- `alamar.screening.protein.[protein_index].rds` - Complete model object (no separate summary file)

**How to Run:**

```bash
# Modify protein_index in script (line 6), then:
cd ~/chronos-pd_phase1_ms
Rscript src/20_analyses/joint_modeling/PD/alamar.R

# For SLURM array job or sequential processing:
# sbatch --array=1-14 run_alamar.slurm
# (or run sequentially with protein_index = 1, 2, ..., 14)
```

**Features**:
- One protein per run (set via `protein_index` variable)
- Small protein set suitable for sequential processing
- MCMC parallelization: 4 chains across available cores
- Output: Full model object only (no separate summary file)

**Number of proteins**: 14 proteins (filtered by q_val_all < 0.05 from screen model)

---

## File Organization

**Output Directory Structure**:
```
output_dir/joint_modeling/PD/              # Somalogic outputs
├── somalogic.screening.protein.1.rds
├── somalogic.screening.protein.1.summary.rds
├── somalogic.screening.protein.2.rds
├── somalogic.screening.protein.2.summary.rds
└── ...

checkpoint_dir/joint_modeling/PD/          # Olink & Alamar outputs
├── olink.screening.protein.1.rds
├── olink.screening.protein.1.summary.rds
├── olink.screening.protein.2.rds
├── olink.screening.protein.2.summary.rds
├── alamar.screening.protein.1.rds
├── alamar.screening.protein.2.rds
└── ...
```

**Note**:
- Somalogic uses `get_config()$output_dir`
- Olink and Alamar use `get_config()$checkpoint_dir`
- Alamar only saves full model objects (no separate summary files)

---

## Sample Inclusion Criteria

**PD patients**: All plasma samples collected BEFORE estimated PD onset date
**Controls**: All plasma samples (no PD diagnosis)

This approach maximizes:
- Statistical power (more observations per subject)
- Model convergence (better MCMC stability)
- Biomarker discovery (greater sensitivity)

### Data Filtering Logic

```r
# Create event indicator
Event_date = ifelse(group == "pd", estimated_pd_onset_date, last_date_id)

# Include only pre-diagnosis samples
dplyr::filter(sample_collection_date < Event_date)

# Survival outcome
Status = ifelse(group == "pd", 1, 0)  # 1 = PD event, 0 = censored
```

---

## Key Statistical Concepts

### Joint Modeling Approach

Joint models simultaneously analyze:

- **Longitudinal submodel**: How protein levels change over time
  - Linear mixed-effects model (LME) with random intercept and slope
  - Accounts for within-subject correlation over multiple time points
- **Survival submodel**: Time-to-PD diagnosis
  - Cox proportional hazards model with left truncation (delayed entry)
  - Entry time: Age at first medical record (`age_entry`)
  - Event time: Age at PD diagnosis or last follow-up (`Time_age`)
  - Event indicator: `Status` (1 = PD, 0 = censored)
- **Association parameter**: Links protein trajectory to PD risk
  - "Value" association: Current protein level predicts hazard of PD diagnosis
  - Interpretation: Higher/lower protein → increased/decreased PD risk

### Model Convergence

- **Rhat**: Gelman-Rubin convergence diagnostic (target: < 1.05)
- **MCMC**: Markov Chain Monte Carlo (30,000 iterations, 4 chains, warmup = 3,000, thin = 2)
- **Stability**: More samples per subject → better convergence

### Interpretation

- **Estimate**: Log hazard ratio (association between protein and PD risk)
- **Hazard Ratio (HR)**: exp(estimate) - multiplicative effect on hazard
  - HR > 1: Higher protein → increased PD risk
  - HR < 1: Higher protein → decreased PD risk
- **P-value**: Statistical significance of association
- **Threshold**: Typically p < 0.01 for biomarker identification

---

## Implementation Notes

### Code Architecture

- **Functions**: Joint modeling functions (`joint_data`, `JointModelSpec`, `fit_joint_models`, etc.) are loaded automatically from `packages/manuscript/` via `.Rprofile`
- **No duplication**: All scripts use the same centralized function library
- **Configuration**: Uses `get_config()` from `config.yaml` for input/output paths
  - `get_config()$input_dir` - Published input data (default: `_input_data`)
  - `get_config()$output_dir` - Analysis outputs (default: `_out`)
  - `get_config()$checkpoint_dir` - Intermediate checkpoints (default: `_out`)

### Execution

- **One protein per run**: Each script processes a single protein specified by `protein_index` (line 6)
- **Array jobs**: Designed for HPC SLURM array job submission to parallelize across proteins
- **Sequential processing**: Can also run locally by modifying `protein_index` and running sequentially

### Parallelization

Scripts use `future::multisession` to parallelize MCMC chains:

- 4 MCMC chains run in parallel
- Uses `availableCores()` to automatically detect available cores
- Each protein is processed independently (embarrassingly parallel across proteins)

### Output Naming Convention

```
[platform].screening.protein.[protein_index].rds
[platform].screening.protein.[protein_index].summary.rds
```

Example: `somalogic.screening.protein.1.rds`

- `somalogic` - Platform (somalogic/olink/alamar)
- `screening` - Analysis type (screening proteins from screen model)
- `protein` - Target type
- `1` - Protein index (1-based indexing)

### Required R Packages

- `JMbayes2` - Joint modeling framework
- `furrr`, `future` - Parallel processing
- `nlme` - Linear mixed-effects models
- `survival` - Cox proportional hazards models
- Tidyverse packages (`dplyr`, `tidyr`, `arrow`, etc.)
- `here` - Path management

All packages managed via `renv` for reproducibility.

---

## Data Dependencies

### Configuration

File paths are managed via `config.yaml`:

```yaml
default:
  input_dir: "_input_data"     # Published input data
  output_dir: "_out"           # Analysis outputs
  checkpoint_dir: "_out"       # Intermediate checkpoints
```

Access in R scripts:

```r
input_dir <- get_config()$input_dir
output_dir <- get_config()$output_dir
checkpoint_dir <- get_config()$checkpoint_dir
```

### Input Data Sources

All input files are loaded from `get_config()$input_dir` in versioned parquet/RDS format:

- **Critical variables table**: PD onset dates, demographics, sample metadata
- **RWD summary**: First/last medical record dates for survival analysis
- **Proteomics data**: Platform-specific three-tibble format (meta_data, prot_data, annot_data)
- **Screen model results**: Master summary datasets with significance filtering (q_val or q_val_all < 0.05)

---

## Downstream Usage

Results from these scripts are used in:

1. **Comorbidity Effects Analysis**: Located in `../comorbidity/`

   - Uses PD-associated proteins identified here (filtered by p < 0.005, rhat < 1.05)
   - Analyzes relationship between comorbidities and protein trajectories
   - See `../comorbidity/model_processing/joint_modeling_preparation.R` for protein selection logic

2. **Master Summary Compilation**:

   - Individual protein results need to be combined
   - Calculate adjusted p-values (Benjamini-Hochberg correction)
   - Generate joint model master summary datasets for publication

---

## Implementation Summary

### Key Features

- **File-based data access**: All inputs loaded from `get_config()$input_dir` in versioned format
- **Centralized functions**: Joint modeling functions auto-loaded from `packages/manuscript/` via `.Rprofile`
- **One protein per run**: Each script processes a single protein specified by `protein_index`
- **MCMC parallelization**: 4 chains run in parallel using `future::multisession`
- **HPC-ready**: Designed for SLURM array job submission but also runs locally

### Platform Differences

- **Somalogic**: Uses `get_config()$output_dir`, saves full model + summary
- **Olink**: Uses `get_config()$checkpoint_dir`, saves full model + summary
- **Alamar**: Uses `get_config()$checkpoint_dir`, saves full model only
- **Filtering**: Somalogic/Olink use `q_val < 0.05`, Alamar uses `q_val_all < 0.05`

### Configuration

All scripts use `config.yaml` for path management:

```yaml
default:
  input_dir: "_input_data"       # Published input data
  output_dir: "_out"             # Analysis outputs (somalogic)
  checkpoint_dir: "_out"         # Intermediate checkpoints (olink, alamar)
```

---
