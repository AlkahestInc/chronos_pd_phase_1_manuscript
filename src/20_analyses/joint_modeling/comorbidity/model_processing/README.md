# Joint Modeling (JM) Analysis for Comorbidity Effects

This folder contains R scripts for performing Joint Modeling analysis to assess the relationship between comorbidities and PD-associated protein biomarkers.

> **Note**: Scripts have been updated to be **click-and-run** without requiring SLURM or command-line arguments. Each script automatically processes all diagnosis × protein combinations with built-in progress tracking and error handling.

## Overview

The Joint Modeling (JM) approach combines:

- **Longitudinal data**: Protein measurements from plasma samples collected over time
- **Survival data**: Time-to-event data for chronic condition diagnoses

The analysis uses Bayesian joint models implemented with JMbayes2.

## Analysis Workflow

The analysis follows this sequential workflow:

```
1. joint_modeling_preparation.R
   → Outputs to: get_config()$output_dir/joint_modeling/comorbidity/model_processing/
   ↓
   [Copy/publish datasets to get_config()$input_dir]
   ↓
2. [Platform]_PD_comorbidity.R OR [Platform]_other_comorbidity.R
   → Outputs to: get_config()$checkpoint_dir/joint_modeling/comorbidity/pd_comorbidity/ (or other_comorbidity/)
   (each script automatically processes all diagnosis × protein combinations)
```

**Two parallel analyses:**
- **PD Comorbidities**: 50 conditions with adjusted p < 0.05 (G20 excluded from modeling as it represents PD itself)
- **Other Comorbidities**: 52 non-PD-associated chronic conditions (adjusted p ≥ 0.05, n ≥ 30)

---

## Main Analysis Scripts

### 1. joint_modeling_preparation.R

**Purpose**: Data preparation and feature engineering for Joint Modeling analysis.

**Analysis Logic**:

- Load critical variables and proteomics data from three platforms (SomaScan, Olink, Alamar) via file-based access
- Load pre-prepared chronic condition ICD-10 codes (3-digit format with descriptions) from `ccir_icd10_3digits`
  - This file should contain: `patient_id`, `min_diagnosis_date`, `diagnosis_code_3digits`, `icd10_description_3digits`
- Calculate odds ratios (OR) for chronic conditions between PD vs control groups using Fisher's exact test
- Identify significantly enriched PD-related comorbidities (adjusted p < 0.05) - **50 conditions**
  - Note: G20 (Parkinson's disease code) is included in preparation data but excluded from comorbidity modeling
- Identify non-PD-related conditions (adjusted p ≥ 0.05, n ≥ 30) - **52 conditions**
- Prepare survival tables with time-to-event data using claim summary data
- Merge sample collection information with clinical records
- Filter protein biomarkers combining:
  - **Source 1**: Proteins from joint model analysis (rhat < 1.05, p < 0.005)
    - Group by uniprot_id, select best platform (minimum p_val)
  - **Source 2**: Proteins from screen model with pre-onset significance (q_val_pre < 0.05)
    - Group by uniprot_id, select best platform (minimum p_val_pre)
  - Combine both sources and remove duplicates
- Create analysis-ready datasets for Joint Modeling (longitudinal + survival data)
- Save all outputs to `get_config()$output_dir` in versioned parquet format with summary messages

**Input Files** (from `get_config()$input_dir` in versioned parquet/RDS format):

- `name=Chronos_phase1_critical_variables_table/version=29432/dataset.parquet`
- `name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds`
- `name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds`
- `name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29649/dataset.rds`
- `name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet`
- `name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet`
- `name=chronos_alamar_joint_model-pd_master_summary-pd_master_summary_2/version=29821/dataset.parquet`
- `name=chronos_somalogic_merged-pd_master_summary_2/version=29828/dataset.parquet`
- `name=chronos_olink_merged-pd_master_summary_2/version=29829/dataset.parquet`
- `name=chronos_alamar_merged-pd_master_summary_2/version=29830/dataset.parquet`
- `name=ccir_icd10_3digits/version=none/dataset.parquet` - Pre-prepared chronic condition ICD-10 codes
  - Required columns: `patient_id`, `min_diagnosis_date`, `diagnosis_code_3digits`, `icd10_description_3digits`
  - Contains 3-digit ICD-10 codes with first diagnosis dates and descriptions for each patient
- `name=patient_claims_converted/version=none/dataset.parquet` - ICD-9/ICD-10 converted claim records
  - Required columns: `patient_id`, `diagnosis_date`
  - Used to create claim summary with first/last diagnosis dates for survival analysis

**Output Files** (saved to `get_config()$output_dir/joint_modeling/comorbidity/model_processing/` in versioned parquet format):

- `name=OR_phase1_chronic_updated0120/version=none/dataset.parquet` - Odds ratios for all chronic conditions with adjusted p-values
- `name=OR_phase1_chronic_ref/version=none/dataset.parquet` - ICD-10 reference descriptions (merged from ccir_icd10_3digits)
- `name=protein_list/version=none/dataset.parquet` - List of proteins to test (probe_id, uniprot_id, gene_symbol, platform)
  - Combines proteins from joint model (rhat < 1.05, p < 0.005) and screen model (q_val_pre < 0.05)
  - For each source, selects best platform per UniProt ID (minimum p-value)
  - Final list contains distinct proteins across both sources
- `name=JM_lme_somalogic/version=none/dataset.parquet` - SomaScan longitudinal data
- `name=JM_surv_somalogic/version=none/dataset.parquet` - SomaScan survival data
- `name=JM_lme_olink/version=none/dataset.parquet` - Olink longitudinal data
- `name=JM_surv_olink/version=none/dataset.parquet` - Olink survival data
- `name=JM_lme_alamar/version=none/dataset.parquet` - Alamar longitudinal data
- `name=JM_surv_alamar/version=none/dataset.parquet` - Alamar survival data
- `name=JMcondition_lme_[platform]/version=none/dataset.parquet` - Longitudinal data for non-PD conditions
- `name=JMcondition_surv_[platform]/version=none/dataset.parquet` - Survival data for non-PD conditions

**How to Run:**

```bash
# In RStudio: Open the script and click "Source"
# Or from command line:
cd ~/chronos-pd_phase1_ms
Rscript src/20_analyses/joint_modeling/comorbidity/model_processing/joint_modeling_preparation.R
```

**Output Summary**: The script prints a summary at the end showing:
- Output directory location
- File organization structure (versioned format)
- Reminder to publish files to data publishing repository

**Note**: These files are saved locally to `get_config()$output_dir/joint_modeling/comorbidity/model_processing/`. For the modeling scripts to use them, they need to be copied/published to `get_config()$input_dir` (default: `_input_data`). A separate data publishing step or script should handle this transfer.

---

### 2. PD Comorbidity Analysis Scripts

#### Alamar_PD_comorbidity.R / Olink_PD_comorbidity.R / Somalogic_PD_comorbidity.R

**Purpose**: Fit Joint Models for PD-associated comorbidities across all diagnosis × protein combinations.

**Analysis Logic**:

- Joint Modeling functions auto-loaded from `packages/manuscript/` via `.Rprofile`
- Load prepared JM data for the respective platform from `get_config()$input_dir`
- Loop through all diagnosis codes × proteins automatically
- For each combination:
  - Filter data for specific chronic diagnosis code (3-digit ICD-10)
  - Create `joint_data` object with longitudinal and survival components
  - Set up Joint Model specification:
    - **Longitudinal submodel**: Mixed-effects model (LME) with random intercept and slope
      - Formula: `protein ~ age_protein + gender + cluster + recent_plasma_donations`
      - Random effects: `~ 1 + age_protein | PATIENT_ID`
    - **Survival submodel**: Cox proportional hazards model
      - Formula: `Surv(age_entry, Time_age, Status) ~ gender`
    - **Association**: Value-based (current protein level predicts event risk)
  - Configure MCMC settings (4 chains, 30,000 iterations, 3,000 warmup, thin=2)
  - Fit joint model using JMbayes2 package
  - Save results to `get_config()$checkpoint_dir`
- Skip already-completed combinations
- Handle errors gracefully and continue processing

**Input Files** (from `get_config()$input_dir` in versioned parquet format):

- `name=protein_list/version=none/dataset.parquet` - List of proteins to test (dynamically generated from filtering criteria)
- `name=JM_lme_[platform]/version=none/dataset.parquet` - Platform-specific longitudinal data
- `name=JM_surv_[platform]/version=none/dataset.parquet` - Platform-specific survival data

**Output Files** (to `get_config()$checkpoint_dir/joint_modeling/comorbidity/pd_comorbidity/`):

- `full_results/[platform]_JM.[diagnosis_code].[protein_index].[protein_id].RDS` - Complete model object
- `summary_only/[platform]_JM.[diagnosis_code].[protein_index].[protein_id].RDS` - Model summary statistics

**Number of combinations per platform:**
- Alamar: 50 diagnoses × 9 proteins = **450 combinations**
- Olink: 50 diagnoses × 11 proteins = **550 combinations**
- Somalogic: 50 diagnoses × 16 proteins = **800 combinations**

**How to Run:**

```bash
# In RStudio: Open the script and click "Source"
# Or from command line:
cd ~/chronos-pd_phase1_ms
Rscript src/20_analyses/joint_modeling/comorbidity/model_processing/Somalogic_PD_comorbidity.R
```

**Features**:
- Automatic progress tracking with `[5/800]` style counters
- Skips already-completed runs (checks for existing output files)
- Error handling - continues if one combination fails
- Final summary showing completed/skipped/failed counts

---

### 3. Other Comorbidity Analysis Scripts

#### Alamar_other_comorbidity.R / Olink_other_comorbidity.R / Somalogic_other_comorbidity.R

**Purpose**: Fit Joint Models for non-PD-associated chronic conditions (control analysis).

**Analysis Logic**: Same as PD scripts, but uses `JMcondition_lme_*` and `JMcondition_surv_*` datasets (52 non-PD conditions).

**Input Files** (from `get_config()$input_dir` in versioned parquet format):

- `name=protein_list/version=none/dataset.parquet`
- `name=JMcondition_lme_[platform]/version=none/dataset.parquet` - Non-PD condition longitudinal data
- `name=JMcondition_surv_[platform]/version=none/dataset.parquet` - Non-PD condition survival data

**Output Files** (to `get_config()$checkpoint_dir/joint_modeling/comorbidity/other_comorbidity/`):

- `full_results/[platform]_JM.[diagnosis_code].[protein_index].[protein_id].RDS`
- `summary_only/[platform]_JM.[diagnosis_code].[protein_index].[protein_id].RDS`

**Number of combinations per platform:**
- Alamar: 52 diagnoses × 9 proteins = **468 combinations**
- Olink: 52 diagnoses × 11 proteins = **572 combinations**
- Somalogic: 52 diagnoses × 16 proteins = **832 combinations**

**How to Run:**

```bash
# In RStudio: Open the script and click "Source"
# Or from command line:
Rscript src/20_analyses/joint_modeling/comorbidity/model_processing/Somalogic_other_comorbidity.R
```

---

## File Organization

```
model_processing/
├── joint_modeling_preparation.R       # Step 1: Data preparation
├── Somalogic_PD_comorbidity.R         # Step 2: JM for PD comorbidities
├── Olink_PD_comorbidity.R
├── Alamar_PD_comorbidity.R
├── Somalogic_other_comorbidity.R      # Step 2: JM for other conditions
├── Olink_other_comorbidity.R
├── Alamar_other_comorbidity.R
└── README.md                          # This file
```

**Output Directory Structure**:
```
output_dir/joint_modeling/comorbidity/
└── model_processing/                     # From joint_modeling_preparation.R
    ├── name=protein_list/version=none/dataset.parquet
    ├── name=JM_lme_[platform]/version=none/dataset.parquet
    ├── name=JM_surv_[platform]/version=none/dataset.parquet
    ├── name=JMcondition_lme_[platform]/version=none/dataset.parquet
    ├── name=JMcondition_surv_[platform]/version=none/dataset.parquet
    ├── name=OR_phase1_chronic_updated0120/version=none/dataset.parquet
    └── name=OR_phase1_chronic_ref/version=none/dataset.parquet

checkpoint_dir/joint_modeling/comorbidity/
├── pd_comorbidity/                       # From *_PD_comorbidity.R scripts
│   ├── full_results/
│   │   ├── somalogic_JM.[code].[idx].[protein].RDS
│   │   ├── olink_JM.[code].[idx].[protein].RDS
│   │   └── alamar_JM.[code].[idx].[protein].RDS
│   └── summary_only/
│       └── [same structure as full_results]
└── other_comorbidity/                    # From *_other_comorbidity.R scripts
    ├── full_results/
    └── summary_only/
```

---

## Key Statistical Concepts

### Joint Modeling Approach
Joint models simultaneously analyze:

- **Longitudinal submodel**: How protein levels change over time
  - Linear mixed-effects model (LME) with random intercept and slope
  - Accounts for within-subject correlation
- **Survival submodel**: Time-to-event (diagnosis) analysis
  - Cox proportional hazards model
  - Allows for censored observations
- **Association parameter**: Links protein trajectory to event risk
  - "Value" association: Current protein level predicts hazard

### Model Convergence

- **Rhat**: Gelman-Rubin convergence diagnostic (target: < 1.05)
- **MCMC**: Markov Chain Monte Carlo (30,000 iterations, 4 chains)
- **DIC/WAIC/LPML**: Model fit statistics

### Interpretation

- **Estimate**: Log hazard ratio (association between protein and event risk)
- **Hazard Ratio (HR)**: exp(estimate) - multiplicative effect on hazard
- **Risk change**: Percent change in risk per 0.1 unit protein change

---

## Data Dependencies

### Input Data Sources

**For joint_modeling_preparation.R:**
- All input files are loaded from `get_config()$input_dir` in versioned parquet/RDS format

**For modeling scripts:**
- Published input files via `get_config()$input_dir` in versioned parquet format
- All input files should be prepared by data publishing scripts

### Configuration

File paths are managed via `config.yaml`:
```yaml
default:
  input_dir: "_input_data"  # Published input data
  output_dir: "_out"                                   # Analysis outputs
  checkpoint_dir: "_out"                               # Intermediate checkpoints
```

Access in R scripts:
```r
input_dir <- get_config()$input_dir
output_dir <- get_config()$output_dir
checkpoint_dir <- get_config()$checkpoint_dir
```

---

## Implementation Notes

### Code Architecture

- **Functions**: Joint modeling functions (`joint_data`, `JointModelSpec`, `fit_joint_models`, etc.) are loaded automatically from `packages/manuscript/` via `.Rprofile`
- **No duplication**: All scripts use the same centralized function library
- **Click-and-run**: Scripts automatically process all combinations without command-line arguments
- **Configuration**: Uses `get_config()` from `config.yaml` for input/output paths
  - `get_config()$input_dir` - Published input data (default: `_input_data`)
  - `get_config()$output_dir` - Analysis outputs (used by joint_modeling_preparation.R, default: `_out`)
  - `get_config()$checkpoint_dir` - Intermediate checkpoints (used by modeling scripts, default: `_out`)

### Execution Features

- **Automatic looping**: Each script processes all diagnosis × protein combinations
- **Resume capability**: Skips already-completed combinations (checks for existing output files)
- **Error handling**: Catches errors and continues processing remaining combinations
- **Progress tracking**: Real-time progress indicators (`[5/800]` for Somalogic PD comorbidity)
- **Summary reporting**: Final counts of completed/skipped/failed runs

### Parallelization

Scripts use `future::multisession` to parallelize MCMC chains across available cores. Each diagnosis × protein combination runs sequentially, but within each model:
- 4 MCMC chains run in parallel
- Uses `availableCores() - 1` to leave one core free

### Output Naming Convention

```
[platform]_JM.[diagnosis_code].[protein_index].[protein_id].RDS
```

Example: `alamar_JM.E10.1.10000168.RDS`
- `alamar` - Platform (somalogic/olink/alamar)
- `E10` - 3-digit ICD-10 diagnosis code
- `1` - Protein index (for sorting)
- `10000168` - Protein ID (probe_id)

### Required R Packages

- `JMbayes2` - Joint modeling framework
- `furrr`, `future` - Parallel processing
- `nlme` - Linear mixed-effects models
- `survival` - Cox proportional hazards models
- Tidyverse packages (`dplyr`, `tidyr`, `ggplot2`, etc.)

All packages managed via `renv` for reproducibility.

---

