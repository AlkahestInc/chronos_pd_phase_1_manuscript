## Script Details

### ppmi_data_preprocessing.R

Comprehensive cleaning and curation of PPMI Project 9000 plasma data

- Loads raw PPMI proteomics data from Excel file
- Updates UniProt IDs and gene symbols using master annotations
- Processes donor metadata: demographics, race, education, smoking, alcohol, medical conditions, falls
- Processes sample metadata: visit dates, vitals, age calculations, PD-related timepoints
- Creates structured output: `dataProt`, `pheno.samp`, `pheno.don`, `nom` tables

### ppmi_data_qc.R

Quality control analysis for PPMI proteomics data

- Missingness analysis for metadata and proteomics
- PCA analysis for outlier detection
- Sample median-based outlier detection
- Generates summary tables by cohort and visit
- Creates baseline (BL) and longitudinal (LD) datasets

### ppmi_chronos_pd_validation.R

Quality control analysis for PPMI proteomics data

- Loads `PPMI_p9000_Plasma_MS` dataset and generation LD dataset
- Loads somalogic and olink, alamarbio summary tables
- PD hits validation PPMI LD samples based on linear mixed model
- Plot protein trajectory
- Creates Chronos hit validation results
