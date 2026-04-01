## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/target_validation/ppmi/ppmi_chronos_pd_validation.R")

# Input: Pre-processed PPMI dataset
p9k <- glue::glue(
  "{get_config()$checkpoint_dir}/name={name}/version={version}/dataset.rds",
  name = "PPMI_p9000_Plasma_MS",
  version = 29985
) |>
  readr::read_rds()

# Output
output_datafile <- glue::glue(
  "{get_config()$checkpoint_dir}/name={name}/version={version}/dataset.parquet",
  name = "PPMI_pd_hits_lme_results",
  version = "Feb8_2026"
)
dir.create(dirname(output_datafile), showWarnings = FALSE, recursive = TRUE)


# patch: add ID column to dataProt
p9k$dataProt$ID <- rownames(p9k$dataProt)
# patch: add protein column to nom and separate uniprot into individual uniprot IDs when there are complex
p9k$nom <- p9k$nom |>
  dplyr::mutate(protein = paste(LABEL, "NPX", sep = "_")) |>
  tidyr::separate_rows(UNIPROT)

# patch: add age_at_first_visit variable to donor table
p9k$pheno.don$age_at_first_visit <- NA
for (k in 1:nrow(p9k$pheno.don)) {
  i <- which(p9k$pheno.samp$PATNO == p9k$pheno.don$PATNO[k])
  p9k$pheno.don$age_at_first_visit[k] <- min(p9k$pheno.samp$age_at_visit[i], na.rm = T)
}


# generate LD dataset, all visits ####

# baseline meta using earliest date per subject ####
PPMI_pd_meta_BL <- p9k$pheno.samp |>
  dplyr::group_by(PATNO) |>
  dplyr::filter(cohort %in% c("control", "PD") & visit_date == min(visit_date)) |> # BL defined earliest visit for each subject
  dplyr::left_join(p9k$pheno.don |>
    dplyr::select(
      PATNO, ENROLL_STATUS, sex, race, education_years, smoking, alcohol, handedness,
      hypertension, diabetes, lipid_cholesterol, depression_anxiety, cognitive_deficit,
      pddiag, age_at_pd_diag, early_onset_pd, age_at_first_visit
    ), by = "PATNO") |>
  dplyr::rename(ID = SampleId) |>
  dplyr::ungroup()

# LD meta
PPMI_pd_meta_LD <- p9k$pheno.samp |>
  dplyr::filter(cohort %in% c("control", "PD")) |>
  dplyr::group_by(PATNO) |>
  dplyr::mutate(bl_subjects = ifelse(any(EVENT_ID == "BL"), "with_bl", "no_bl")) |>
  dplyr::ungroup() |>
  dplyr::left_join(p9k$pheno.don |>
    dplyr::select(
      PATNO, ENROLL_STATUS, sex, race, education_years, smoking, alcohol, handedness,
      hypertension, diabetes, lipid_cholesterol, depression_anxiety, cognitive_deficit,
      pddiag, age_at_pd_diag, early_onset_pd, age_at_first_visit
    ), by = "PATNO") |>
  dplyr::rename(ID = SampleId) |>
  dplyr::ungroup()


# LD proteomics
PPMI_pd_protein_LD <- p9k$dataProt |>
  dplyr::filter(ID %in% PPMI_pd_meta_LD$ID)
dim(PPMI_pd_protein_LD) # 651 1464
rownames(PPMI_pd_protein_LD) <- NULL

# BL proteomics
PPMI_pd_protein_BL <- PPMI_pd_protein_LD |>
  dplyr::filter(ID %in% PPMI_pd_meta_BL$ID)

# merged LD dataset ####
merged_data_LD <- PPMI_pd_protein_LD |>
  dplyr::inner_join(PPMI_pd_meta_LD, by = "ID")
dim(merged_data_LD) # 651 1506


# protein names
protein_names <- unique(paste0(p9k$nom$LABEL, "_NPX"))
length(unique(protein_names))
# analysis ready table generation after proteomics QC ####
analysis_ready_table <- merged_data_LD |>
  dplyr::filter(EVENT_ID %in% c("BL", "V02", "V04", "V06", "V10")) |> # number of samples>=20
  dplyr::mutate(
    cohort = factor(cohort, levels = c("control", "PD"))
  )
table(analysis_ready_table$cohort)
length(unique(analysis_ready_table$PATNO[which(analysis_ready_table$cohort == "PD")]))
length(unique(analysis_ready_table$PATNO[which(analysis_ready_table$cohort == "control")]))


# master summary table, one per protein, when there are multiple assys, most significant one was selected
master_summary <- glue::glue(
  "{get_config()$input_dir}/name={name}/version={version}/dataset.parquet",
  name = "chronos_all_tables_merged_by_uniprot-pd_master_summary_2",
  version = 29975
) |>
  arrow::open_dataset() |>
  dplyr::collect()

Chronos_pd_hits <- master_summary$uniprot_id[which(master_summary$q_val_all_somalogic < 0.05 |
  master_summary$q_val_all_olink < 0.05 |
  master_summary$q_val_all_alamar < 0.05)]

## linear mixed modeling mixed model formula ####
formula1.1 <- paste0("protein~", paste0("cohort", "*", "EVENT_ID", "+", "sex", "+", "age_at_first_visit"))
formula2 <- paste0("~1|PATNO")

# difference at baseline --BL as the reference level ####
lme.fit1 <- LME_results_P1(
  dataset = analysis_ready_table,
  protein_names = protein_names,
  ind_cov = c("cohort", "age_at_first_visit", "sex", "EVENT_ID"),
  ID = "PATNO",
  fixedFormula = formula1.1,
  randomFormula = formula2,
  method = "REML"
)

# plot results
# baseline difference
LME_results_P1_plot(
  dataset = lme.fit1,
  coef_names = "cohortPD",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.2, xlim_max = 1.2, ylim_min = 0, ylim_max = 16, plot = TRUE, topN = 10, pvalue = 0.05
)


PPMI_pd_hits_lme_BL <- lme.fit1 |>
  dplyr::filter(coef.names == "cohortPD") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

PPMI_pd_hits_lme_BL |>
  dplyr::arrange(p.value) |>
  dplyr::select(Value, p.value, adjP, SYMBOL) |>
  dplyr::filter(p.value < 0.05)

# compare these proteins in PPMI and Chronos
ggplot2::ggplot(lme.fit1 |> dplyr::filter(coef.names == "cohortPD"), ggplot2::aes(x = Value, y = -log10(p.value))) +
  ggplot2::geom_point() +
  ggplot2::geom_point(data = lme.fit1 |> dplyr::filter(coef.names == "cohortPD" & p.value < 0.05), col = "red") +
  ggplot2::geom_point(data = lme.fit1 |> dplyr::filter(coef.names == "cohortPD" & protein %in% PPMI_pd_hits_lme_BL$protein), col = "blue") +
  ggplot2::theme_minimal() +
  ggplot2::xlab("log2FC between PD and control at bl") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 12),
    axis.text.y = ggplot2::element_text(size = 12)
  )

# difference in change between V10 and baseline between PD and control
LME_results_P1_plot(
  dataset = lme.fit1,
  coef_names = "cohortPD:EVENT_IDV10",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.3, xlim_max = 1.3, ylim_min = 0, ylim_max = 16, plot = "TRUE", topN = 15, pvalue = 0.05
)

PPMI_pd_hits_lme_V10_BL <- lme.fit1 |>
  dplyr::filter(coef.names == "cohortPD:EVENT_IDV10") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

# compare these proteins in PPMI and Chronos
ggplot2::ggplot(lme.fit1 |> dplyr::filter(coef.names == "cohortPD:EVENT_IDV10"), ggplot2::aes(x = Value, y = -log10(p.value))) +
  ggplot2::geom_point() +
  ggplot2::geom_point(data = lme.fit1 |> dplyr::filter(coef.names == "cohortPD:EVENT_IDV10" & p.value < 0.05), col = "red") +
  ggplot2::geom_point(data = lme.fit1 |> dplyr::filter(coef.names == "cohortPD:EVENT_IDV10" &
    protein %in% PPMI_pd_hits_lme_V10_BL$protein), col = "blue") +
  ggplot2::theme_minimal() +
  ggplot2::xlab("log2FC difference in changes between PD and control from V10 to bl") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_text(size = 12), # X-axis text size
    axis.text.y = ggplot2::element_text(size = 12) # Y-axis text size
  )

# difference at V02 ####
# reset reference level of EVENT ID to have V02 as the reference level
analysis_ready_table$EVENT_ID <- factor(analysis_ready_table$EVENT_ID, levels = c("V02", "BL", "V04", "V06", "V10"))

lme.fit2 <- LME_results_P1(
  dataset = analysis_ready_table,
  protein_names = protein_names,
  ind_cov = c("cohort", "age_at_first_visit", "sex", "EVENT_ID"),
  ID = "PATNO",
  fixedFormula = formula1.1,
  randomFormula = formula2,
  method = "REML"
)

# plot results
# V02 difference
LME_results_P1_plot(
  dataset = lme.fit2,
  coef_names = "cohortPD",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.2, xlim_max = 1.2, ylim_min = 0, ylim_max = 16, plot = "TRUE", topN = 10, pvalue = 0.05
)

# save the pd hits measured and validated in PPMI p9k
PPMI_pd_hits_lme_V02 <- lme.fit2 |>
  dplyr::filter(coef.names == "cohortPD") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

# difference at V04 ####
analysis_ready_table$EVENT_ID <- factor(analysis_ready_table$EVENT_ID, levels = c("V04", "BL", "V02", "V06", "V10"))
lme.fit3 <- LME_results_P1(
  dataset = analysis_ready_table,
  protein_names = protein_names,
  ind_cov = c("cohort", "age_at_first_visit", "sex", "EVENT_ID"),
  ID = "PATNO",
  fixedFormula = formula1.1,
  randomFormula = formula2,
  method = "REML"
)

# plot results
# V04 difference
LME_results_P1_plot(
  dataset = lme.fit3,
  coef_names = "cohortPD",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.2, xlim_max = 1.2, ylim_min = 0, ylim_max = 16, plot = "TRUE", topN = 10, pvalue = 0.05
)

# save the pd hits measured and validated in PPMI p9k
PPMI_pd_hits_lme_V04 <- lme.fit3 |>
  dplyr::filter(coef.names == "cohortPD") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

PPMI_pd_hits_lme_V04 |>
  dplyr::arrange(p.value) |>
  dplyr::select(Value, p.value, adjP, SYMBOL) |>
  dplyr::filter(p.value < 0.05)


# difference at V06 ####
analysis_ready_table$EVENT_ID <- factor(analysis_ready_table$EVENT_ID, levels = c("V06", "BL", "V02", "V04", "V10"))

lme.fit4 <- LME_results_P1(
  dataset = analysis_ready_table,
  protein_names = protein_names,
  ind_cov = c("cohort", "age_at_first_visit", "sex", "EVENT_ID"),
  ID = "PATNO",
  fixedFormula = formula1.1,
  randomFormula = formula2,
  method = "REML"
)

# plot results
# V06 difference
LME_results_P1_plot(
  dataset = lme.fit4,
  coef_names = "cohortPD",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.2, xlim_max = 1.2, ylim_min = 0, ylim_max = 16, plot = "TRUE", topN = 10, pvalue = 0.05
)

# save the pd hits measured and validated in PPMI p9k
PPMI_pd_hits_lme_V06 <- lme.fit4 |>
  dplyr::filter(coef.names == "cohortPD") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

PPMI_pd_hits_lme_V06 |>
  dplyr::arrange(p.value) |>
  dplyr::select(Value, p.value, adjP, SYMBOL) |>
  dplyr::filter(p.value < 0.05)

# difference at V10 ####
analysis_ready_table$EVENT_ID <- factor(analysis_ready_table$EVENT_ID, levels = c("V10", "BL", "V02", "V04", "V06"))

lme.fit5 <- LME_results_P1(
  dataset = analysis_ready_table,
  protein_names = protein_names,
  ind_cov = c("cohort", "age_at_first_visit", "sex", "EVENT_ID"),
  ID = "PATNO",
  fixedFormula = formula1.1,
  randomFormula = formula2,
  method = "REML"
)

# plot results
# V10 difference
LME_results_P1_plot(
  dataset = lme.fit5,
  coef_names = "cohortPD",
  variable_names = c("Value", "Std.Error", "p.value", "protein"),
  xlim_min = -1.2, xlim_max = 1.2, ylim_min = 0, ylim_max = 16, plot = "TRUE", topN = 10, pvalue = 0.05
)

# save the pd hits measured and validated in PPMI p9k
PPMI_pd_hits_lme_V10 <- lme.fit5 |>
  dplyr::filter(coef.names == "cohortPD") |>
  dplyr::inner_join(p9k$nom, by = "protein") |>
  dplyr::filter(UNIPROT %in% Chronos_pd_hits) |>
  dplyr::mutate(adjP = p.adjust(p.value, method = "BH"))

# save the pd hits measured and validated in PPMI p9k
PPMI_lme_pd_hits <- PPMI_pd_hits_lme_BL |>
  dplyr::mutate(visit = "BL") |>
  dplyr::bind_rows(PPMI_pd_hits_lme_V02 |> dplyr::mutate(visit = "V02")) |>
  dplyr::bind_rows(PPMI_pd_hits_lme_V04 |> dplyr::mutate(visit = "V04")) |>
  dplyr::bind_rows(PPMI_pd_hits_lme_V06 |> dplyr::mutate(visit = "V06")) |>
  dplyr::bind_rows(PPMI_pd_hits_lme_V10 |> dplyr::mutate(visit = "V10"))

PPMI_lme_pd_hits |>
  arrow::write_parquet(output_datafile)

#####
# Generate plots
#####
analysis_ready_table$EVENT_ID <- factor(analysis_ready_table$EVENT_ID, levels = c("BL", "V02", "V04", "V06", "V10"))

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_ITGAV_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "loess", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()


ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_DDC_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "loess", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_CRIM1_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_CD276_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_NEFL_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_SMOC1_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_ITGA11_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()

ggplot2::ggplot(
  analysis_ready_table,
  ggplot2::aes(x = EVENT_ID, y = P9000_Plasma_SCARF2_NPX, group = PATNO, col = cohort)
) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::geom_line(alpha = 0.3, linewidth = 0.5) +
  ggplot2::geom_smooth(ggplot2::aes(group = cohort), method = "lm", linewidth = 1.2, linetype = "solid") +
  ggplot2::theme_minimal()
