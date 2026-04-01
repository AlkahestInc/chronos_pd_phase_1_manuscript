# Non-Linear Mixed Effects Model Analysis

This directory contains scripts for evaluating non-linear trajectory models of protein biomarkers in the Chronos Phase 1 Parkinson's Disease study.

## Purpose

The analysis tests whether non-linear mixed-effects models (using spline functions) provide better fit than linear models when modeling protein trajectories over time. This addresses the hypothesis that protein changes may follow non-linear patterns during disease progression.

## Contents

- `somalogic.R` - Non-linear modeling for SomaLogic platform proteins
- `alamar.R` - Non-linear modeling for Alamar platform proteins
- `olink.R` - Non-linear modeling for Olink platform proteins

## Analysis Workflow

Each script follows the same workflow:

1. **Load Data**: Reads analysis-ready proteomics data for the respective platform
2. **Identify Candidate Proteins**: Loads previously identified significant proteins (q-value < 0.05) from linear mixed model screening
3. **Fit Non-Linear Models**: Tests two types of non-linear spline models:
   - **Linear splines** with varying knot positions (range: -7 to 5 years)
   - **Natural cubic splines** with varying degrees of freedom (df range: 2-6)
4. **Model Selection**: Compares non-linear model fit to baseline linear model using BIC (Bayesian Information Criterion)
5. **Save Results**: Writes model results to checkpoint directory

## Model Comparison Strategy

- Models are compared using BIC difference (diff_lsBIC or diff_lmeBIC)
- Proteins with BIC improvement > 2 are considered to show meaningful non-linear patterns
- Selected proteins are further filtered for significant group differences

## Outputs

Results are saved to the checkpoint directory:

```
checkpoint_datasets/non_linear_mixed/
├── somalogic/
│   ├── somalogic_all_ls.rds  # Linear spline results
│   └── somalogic_all_nc.rds  # Natural cubic spline results
├── alamar/
│   ├── alamar_all_ls.rds
│   └── alamar_all_nc.rds
└── olink/
    ├── olink_all_ls.rds
    └── olink_all_nc.rds
```

## Execution

Each script can be run independently:

```r
source("src/20_analyses/non_linear_mixed/somalogic.R")
source("src/20_analyses/non_linear_mixed/alamar.R")
source("src/20_analyses/non_linear_mixed/olink.R")
```

Scripts use parallel processing with all available cores (`mc.cores = parallelly::availableCores()`) for computational efficiency.

## Interpretation

- **No proteins with BIC > 2**: Linear trajectory is sufficient
- **Proteins with BIC > 2**: Non-linear trajectory provides meaningful improvement
- Proteins showing both improved fit AND significant group differences are candidates for further investigation
