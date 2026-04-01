# GNPC Target Validation

This directory contains R scripts for validating Chronos protein ratio and single protein biomarkers in selected **GNPC Harmony V1 dataset** cohorts.

## Overview

**Analysis environment**: AD Workbench RStudio (requires approved access to GNPC Harmony V1 dataset)

**Cohorts**: Five GNPC cohorts are analyzed:

- **C**, **F**, **J**, **L**, **Q** (all use plasma samples; J and Q noted as "plasma" in filenames)

**Purpose**: Validate protein biomarkers discovered in Chronos datasets (Alamar/Somalogic/Olink platforms) using independent GNPC cohorts.

## Directory Structure

```
README.md                                                        # This file
README.html                                                      # Rendered documentation

data_preparation/                                                # Data preprocessing and sample matching
├── gnpc_cfjlq_3_tibble_prep.rmd                                 # Main pipeline (knit in RStudio)
├── cohort_C_matching_sample_prep.R                              # MatchIt matching for Cohort C
├── cohort_F_matching_sample_prep.R                              # MatchIt matching for Cohort F
├── cohort_J_plasma_sample_matching_sample_prep.R                # MatchIt matching for Cohort J (plasma)
├── cohort_L_matching_sample_prep.R                              # MatchIt matching for Cohort L
└── cohort_Q_plasma_sample_matching_sample_prep.R                # MatchIt matching for Cohort Q (plasma)

data_qc/                                                         # Quality control reports (knit in RStudio)
├── Cohort_C_QC_matching_sample_summary_vfin.rmd                 # QC report for Cohort C
├── Cohort_F_QC_matching_sample_summary_vfin.rmd                 # QC report for Cohort F
├── Cohort_J_plasma_samples_QC_matching_sample_summary_vfin.rmd  # QC report for Cohort J
├── Cohort_L_QC_matching_sample_summary_vfin.rmd                 # QC report for Cohort L
└── Cohort_Q_plasma_samples_QC_matching_sample_summary_vfin.rmd  # QC report for Cohort Q

protein_ratio_validation/                                        # Protein ratio biomarker validation
├── GNPC_cohort_CFJLQ_protein_ratio_validation.R                 # Validates ratio biomarkers (A/B)
└── GNPC_cohort_CFJLQ_protein_ratio_ind_prot_validation.R        # Validates individual proteins (A, B separately)

single_protein_validation/                                       # Single protein biomarker validation
└── GNPC_cohort_CFJLQ_single_protein_validation.R                # Validates individual protein biomarkers
```


## Workflow Overview

### 1. Data Preparation

The data preparation pipeline processes raw GNPC Harmony V1 data (requires approved access in AD Workbench) and creates standardized analysis-ready datasets.

**Step 1a: Core preprocessing** (`data_preparation/gnpc_cfjlq_3_tibble_prep.rmd`)

Transforms raw GNPC proteomics + clinical data into three tidy tibbles:

1. **Analytes tibble** - assay-level information
2. **Proteomics tibble** - quantified protein/peptide measurements
3. **Clinical/metadata tibble** - sample-level patient and clinical information

**How to run**: Open in RStudio and click "Knit" button

**Outputs**: Three RDS files per cohort saved to `/home/rstudio/data_prep/{cohort}/` (by default)

**Step 1b: Sample matching** (cohort-specific R scripts)

Uses `MatchIt::matchit` for propensity score matching to create balanced PD case/control subcohorts.

**How to run**: Source the R scripts in RStudio or run interactively

**Matching strategy**: 

- **Matching ratio**: 1:1 or 1:2 (PD cases : controls) based on cohort sample sizes
- **Predictors**: Typically age_at_visit, sex, and AD status
  - Note: One cohort lacks AD status completely and uses only age/sex
  - Samples with missing predictor values are removed before matching
- **Study design handling**:
  - Cross-sectional (CS) cohorts: Use all samples
  - Longitudinal (LD) cohorts: Use most recent visit per participant
- **Sample type**: All cohorts use plasma samples only (CSF excluded)

**Outlier detection**: All cohorts use both PCA and median plate-based outlier detection

**Outputs**: Matched sample RDS files saved to `/home/rstudio/data_prep/gnpc_cohort_cleanup_matched_data/`

### 2. Quality Control

**Purpose**: Assess quality of sample matching before validation analysis

**Reports** (`data_qc/*.rmd`):

- Matching quality assessment for each cohort
- Sample distribution summaries
- Covariate balance checks

**How to run**: Open each .rmd file in RStudio and click "Knit" to generate HTML reports

### 3. Validation Analysis

**Protein ratio validation** (`protein_ratio_validation/`)

Validates protein ratio biomarkers (e.g., ProteinA/ProteinB) discovered in Chronos datasets.

- **Ratio validation**: Tests whether the ratio (A/B) predicts disease status in GNPC cohorts
  - Script: `GNPC_cohort_CFJLQ_protein_ratio_validation.R`
- **Component validation**: Tests whether individual proteins (A, B separately) predict disease status
  - Script: `GNPC_cohort_CFJLQ_protein_ratio_ind_prot_validation.R`

**How to run**: Source scripts in RStudio or run interactively

**Single protein validation** (`single_protein_validation/`)

Validates individual protein biomarkers from Chronos Alamar/Somalogic/Olink platforms using linear regression models.

- Script: `GNPC_cohort_CFJLQ_single_protein_validation.R`

**How to run**: Source script in RStudio or run interactively

## Required Inputs

### External data (must be available before running)

**GNPC Harmony V1 raw data**:

- Requires approved access in AD Workbench environment
- Location: `/home/rstudio/GNPC Harmonised Dataset v1/`
- **Important**: GNPC data cannot be downloaded and must remain in AD Workbench

**Protein lists from Chronos discovery** (external inputs from upstream Chronos linear mixed model analyses):

- Single protein lists (LMM results by platform, RDS format in `/home/rstudio/input_data/single_protein_validation/`):
  - `chronos_somalogic_protein_list_all.rds` - Somalogic LMM results
  - `chronos_alamar_protein_list_all.rds` - Alamar LMM results
  - `chronos_olink_protein_list_all.rds` - Olink LMM results
- Protein ratio list (Parquet format, accessed via `get_config()$input_dir`):
  - `Chronos_phase1_protein_ratio_list_1` (version 29452)
  - Path: `{get_config()$input_dir}/name=Chronos_phase1_protein_ratio_list_1/version=29452/dataset.parquet`
  - Contains UniProt IDs and gene symbols for mapping (no separate mapping file needed)
  - **Note**: Only list_1 is used; lists 2 and 3 mentioned in older documentation are deprecated

**Somalogic annotation** (generated from GNPC raw data):

- `annotation_tbl_v6.rds` (location: `/home/rstudio/data_prep/annotation/`)
  - Human protein-filtered Somalogic vendor annotation table
  - Generated from: `/home/rstudio/GNPC Harmonised Dataset v1/Clinical/SomalogicAnalyteInfoV1_anonymized.csv`

### Generated during workflow

**From Step 1a** (`gnpc_cfjlq_3_tibble_prep.rmd`):

- Three tibbles per cohort (analytes, proteomics, clinical/metadata) saved as RDS files

**From Step 1b** (matching scripts):

- Matched sample datasets per cohort saved as RDS files

## Usage

### Prerequisites

1. Access to AD Workbench RStudio with approved GNPC Harmony V1 dataset access
2. Required R packages: `MatchIt`, `tidyverse`, `arrow`, `glue`, `knitr`, `ggrepel` (see Required R Packages section below)
3. External protein lists and metadata available (see Required Inputs above)

### Step-by-step execution

**Step 1: Data Preparation**

```r
# In RStudio, open and knit:
data_preparation/gnpc_cfjlq_3_tibble_prep.rmd

# Then source each matching script:
source("data_preparation/cohort_C_matching_sample_prep.R")
source("data_preparation/cohort_F_matching_sample_prep.R")
source("data_preparation/cohort_J_plasma_sample_matching_sample_prep.R")
source("data_preparation/cohort_L_matching_sample_prep.R")
source("data_preparation/cohort_Q_plasma_sample_matching_sample_prep.R")
```

**Step 2: Quality Control**

```r
# In RStudio, open and knit each QC report:
data_qc/Cohort_C_QC_matching_sample_summary_vfin.rmd
# Repeat for F, J, L, Q
```

**Step 3: Validation Analysis**

```r
# Protein ratio validation
source("protein_ratio_validation/GNPC_cohort_CFJLQ_protein_ratio_validation.R")
source("protein_ratio_validation/GNPC_cohort_CFJLQ_protein_ratio_ind_prot_validation.R")

# Single protein validation
source("single_protein_validation/GNPC_cohort_CFJLQ_single_protein_validation.R")
```

### Output locations

Validation results are automatically written to subdirectories:

- `protein_ratio_validation/protein_ratio_validation_result/` - Ratio validation results
- `protein_ratio_validation/ind_protein_validation_result/` - Individual protein component results
- `single_protein_validation/result/` - Single protein validation results

Output directories are created automatically by the scripts.

## Required R Packages

Core requirements:

- `MatchIt` - Propensity score matching
- `tidyverse` - Data manipulation (includes `dplyr`, `tidyr`, `ggplot2`, `readr`, `purrr`)
- `arrow` - Parquet file reading
- `glue` - Path construction
- `knitr` - R Markdown rendering
- `ggrepel` - Text label positioning in plots

Additional packages may be required depending on specific analyses.

## Technical Notes

**Cohort selection process**:

1. Initial QC identifies cohorts with: missingness, outliers, plate effects, hidden factors
2. Cohorts with obvious hidden factors or very different demographic profiles are excluded
3. Remaining cohorts undergo matching and validation analysis

**Cohort variability explained**:

- **Study design**: Mix of cross-sectional (CS) and longitudinal (LD) cohorts
- **Visit selection**: LD cohorts use most recent visit per participant
- **Sample type**: All validation uses plasma samples only (Cohorts J and Q are specifically noted as plasma; all others also use plasma)
- **Predictor availability**: One cohort completely lacks AD status, so matching uses only age and sex
- **Matching ratios**: 1:1 or 1:2 (cases:controls) based on available sample sizes

**Outlier detection**:

- All cohorts use both PCA and median plate-based outlier detection
- In 4 of 5 cohorts, median plate method found no additional outliers beyond PCA
- Only Cohort F scripts show both methods explicitly; others show PCA only as median plate was redundant

**Data format notes**:

- Original analysis used RDS files for some protein lists
- Current scripts primarily use Parquet format accessed through versioned datasets at `get_config()$input_dir`
- Only protein_ratio_list_1 is used (lists 2 and 3 are deprecated)

**Data security**: All GNPC Harmony V1 data must remain in AD Workbench and cannot be downloaded or transferred outside the secure environment.

**File-based access**: This is a visualization/validation repository. All data inputs must come from files - no direct database queries or internal data extraction tools.

## Related Documentation

- UK Biobank validation: `src/20_analyses/target_validation/uk_biobank/`
- Manuscript tables using these results: `src/40_manuscript/supplementary_tables/table_s11.R`
