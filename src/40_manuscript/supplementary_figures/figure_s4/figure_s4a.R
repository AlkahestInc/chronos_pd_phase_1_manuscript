## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s4/figure_s4a.R")

GNPC_C_result_file <- glue::glue("{get_config()$input_dir}/name=GNPC_C_single_protein_validation_result/version=none/dataset.parquet")
GNPC_F_result_file <- glue::glue("{get_config()$input_dir}/name=GNPC_F_single_protein_validation_result/version=none/dataset.parquet")
GNPC_J_result_file <- glue::glue("{get_config()$input_dir}/name=GNPC_J_single_protein_validation_result/version=none/dataset.parquet")
GNPC_L_result_file <- glue::glue("{get_config()$input_dir}/name=GNPC_L_single_protein_validation_result/version=none/dataset.parquet")
GNPC_Q_result_file <- glue::glue("{get_config()$input_dir}/name=GNPC_Q_single_protein_validation_result/version=none/dataset.parquet")
PPMI_validated_results_lme_file <- glue::glue("{get_config()$input_dir}/name=PPMI_pd_hits_lme_results/version=Feb8_2026/dataset.parquet")

soma_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet")
olink_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet")
alamar_all_file <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")

assay_id2_gene_symbol_file <- glue::glue("{get_config()$input_dir}/name=protein_platform_uniprot_ids/version=29404/dataset.parquet")
ukb_single_protein_adjusted_linear_regression_file <- glue::glue("{get_config()$input_dir}/name=ukb_single_protein_adjusted_linear_regression/version=20251216_060400_+0000/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s4/figure_s4a.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

#####
# Load Chronos data
#####
soma_all <- soma_all_file |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  tidyr::separate_rows(uniprot_id)
olink_all <- olink_all_file |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  tidyr::separate_rows(uniprot_id)
alamar_all <- alamar_all_file |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  tidyr::separate_rows(uniprot_id)

#####
# Load GNPC data
#####
GNPC_C_result <- GNPC_C_result_file |>
  arrow::open_dataset() |>
  dplyr::collect()
GNPC_F_result <- GNPC_F_result_file |>
  arrow::open_dataset() |>
  dplyr::collect()
GNPC_J_result <- GNPC_J_result_file |>
  arrow::open_dataset() |>
  dplyr::collect()
GNPC_L_result <- GNPC_L_result_file |>
  arrow::open_dataset() |>
  dplyr::collect()
GNPC_Q_result <- GNPC_Q_result_file |>
  arrow::open_dataset() |>
  dplyr::collect()

#####
# Load UKB data
#####
ukb_single_protein_adjusted_linear_regression <- ukb_single_protein_adjusted_linear_regression_file |>
  arrow::open_dataset() |>
  dplyr::collect()

#####
# Load PPMI data
#####
PPMI_validated_results_lme <- PPMI_validated_results_lme_file |>
  arrow::open_dataset() |>
  dplyr::collect()

assay_id2_gene_symbol <- assay_id2_gene_symbol_file |>
  arrow::open_dataset() |>
  dplyr::collect()

all_combined <- dplyr::bind_rows(
  soma_all |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, std_err, cohens_d, p_val, q_val) |>
    dplyr::mutate(
      platform = "SomaScanV5",
      phase_sample = "all"
    ) |>
    dplyr::distinct(),
  olink_all |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, std_err, cohens_d, p_val, q_val) |>
    dplyr::mutate(
      platform = "Olink_Explore_HT",
      phase_sample = "all"
    ) |>
    dplyr::distinct(),
  alamar_all |>
    dplyr::select(probe_id, uniprot_id, gene_symbol, estimate, std_err, cohens_d, p_val, q_val) |>
    dplyr::mutate(
      platform = "Alamar_NULISAseq_CNS",
      phase_sample = "all"
    ) |>
    dplyr::distinct()
)

Chronos_hits <- all_combined |>
  dplyr::filter(
    q_val < 0.05
  )


GNPC_measured_proteins <- GNPC_C_result |>
  dplyr::filter(prot_uniprot %in% Chronos_hits$uniprot_id) |>
  dplyr::select(prot_uniprot) |>
  dplyr::distinct()
PPMI_measured_proteins <- PPMI_validated_results_lme |>
  dplyr::filter(UNIPROT %in% Chronos_hits$uniprot_id) |>
  dplyr::select(UNIPROT) |>
  dplyr::distinct()

### adding UKBB model results ####
tbl <- ukb_single_protein_adjusted_linear_regression |>
  # Rename the Protein column as it;s really the ID for the purposes of the analysis
  dplyr::rename(
    id = Protein
  ) |>
  # Start mapping the ID to gene symbols so we can then map to uniprot
  #   the id is actually a concatenation of gene symbols, so lets use that as the basis for our conversion
  dplyr::mutate(
    gene_symbol = id,
  ) |>
  # Lets simplify things by only taking the column forward that we need, you could adjust/drop this as needed
  # The gene symbols are delimited by "_" in the UKB, so split them out so we have a gene symbol per row
  tidyr::separate_longer_delim(gene_symbol, delim = "_") |>
  # Lets mark these as coming from the results table
  dplyr::mutate(IN_UKB_RESULTS = TRUE) |>
  # Lets map the gene symbols to uniprot ids using Rick's mapping table
  dplyr::left_join(
    assay_id2_gene_symbol |>
      # Lets only use those we know are from olink...we might not need this?
      dplyr::filter(platform == "olink") |>
      # We are only interested in gene_symbol to uniprot_id mappings
      dplyr::distinct(
        gene_symbol,
        uniprot_id,
      ) |>
      # There may be 1:M mappings between uniprot_id:gene_symbol. Rick collapses these with a ":" separator, lets separate them out
      tidyr::separate_longer_delim(gene_symbol, delim = ":") |>
      # Lets mark these as coming from the mapping table
      dplyr::mutate(IN_MAPPING = TRUE),
    by = dplyr::join_by(gene_symbol)
  ) |>
  # Converte any NA's to FALSE
  dplyr::mutate(
    IN_MAPPING = tidyr::replace_na(IN_MAPPING, FALSE),
    IN_UKB_RESULTS = tidyr::replace_na(IN_UKB_RESULTS, FALSE),
  ) |>
  dplyr::relocate(
    id, gene_symbol, uniprot_id, IN_UKB_RESULTS, IN_MAPPING
  )

ukb_single_protein_adjusted_linear_regression_updated <- tbl |>
  dplyr::filter(!is.na(uniprot_id)) |>
  dplyr::distinct(id, gene_symbol, uniprot_id, beta, std.error, statistic, p.value, conf.low, conf.high, neg_log_p)
ukb_measured_proteins <- ukb_single_protein_adjusted_linear_regression_updated |>
  dplyr::filter(uniprot_id %in% Chronos_hits$uniprot_id) |>
  dplyr::select(uniprot_id) |>
  dplyr::distinct()

.prep_dataset <- function(data, estimate, pval, id, prefix) {
  data |>
    dplyr::select({{estimate}}, {{pval}}, {{id}}) |>
    dplyr::rename(
      !!paste0(prefix, "_estimate") := {{estimate}},
      !!paste0(prefix, "_p_val") := {{pval}},
      uniprot_id = {{id}}
    ) |>
    dplyr::group_by(uniprot_id) |>
    dplyr::slice_min(!!rlang::sym(paste0(prefix, "_p_val"))) |>
    dplyr::ungroup() |>
    dplyr::distinct()
}

master_summary <- all_combined |>
  dplyr::filter(q_val < 0.05) |>
  tidyr::separate_longer_delim(uniprot_id, delim = "[|;:]") |>
  dplyr::mutate(uniprot_id = trimws(uniprot_id)) |>
  dplyr::select(estimate, p_val, uniprot_id, gene_symbol) |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::rename(
    Chronos_estimate = estimate,
    Chronos_p_val = p_val
  ) |>
  dplyr::left_join(
    .prep_dataset(PPMI_validated_results_lme, Value, p.value, UNIPROT, "PPMI"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(GNPC_C_result, pdYes_estimate, pdYes_p_value, prot_uniprot, "GNPC_C"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(GNPC_F_result, pdYes_estimate, pdYes_p_value, prot_uniprot, "GNPC_F"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(GNPC_J_result, pdYes_estimate, pdYes_p_value, prot_uniprot, "GNPC_J"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(GNPC_L_result, pdYes_estimate, pdYes_p_value, prot_uniprot, "GNPC_L"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(GNPC_Q_result, pdYes_estimate, pdYes_p_value, prot_uniprot, "GNPC_Q"),
    by = "uniprot_id"
  ) |>
  dplyr::left_join(
    .prep_dataset(
      ukb_single_protein_adjusted_linear_regression_updated,
      beta,
      p.value,
      uniprot_id,
      "UKB"
    ),
    by = "uniprot_id"
  ) |>
  dplyr::filter(
    PPMI_p_val < 0.05 | GNPC_C_p_val < 0.05 | GNPC_J_p_val < 0.05 | GNPC_F_p_val < 0.05 | GNPC_L_p_val < 0.05 | GNPC_Q_p_val < 0.05 | UKB_p_val < 0.05
  )

# create a draframe for heatmap 1
mat <- tibble::tibble(
  "Chronos all" = -log10(master_summary$Chronos_p_val) * sign(master_summary$Chronos_estimate),
  "PPMI p9k" = -log10(master_summary$PPMI_p_val) * sign(master_summary$PPMI_estimate),
  "GNPC C" = -log10(master_summary$GNPC_C_p_val) * sign(master_summary$GNPC_C_estimate),
  "GNPC F" = -log10(master_summary$GNPC_F_p_val) * sign(master_summary$GNPC_F_estimate),
  "GNPC J" = -log10(master_summary$GNPC_J_p_val) * sign(master_summary$GNPC_J_estimate),
  "GNPC L" = -log10(master_summary$GNPC_L_p_val) * sign(master_summary$GNPC_L_estimate),
  "GNPC Q" = -log10(master_summary$GNPC_Q_p_val) * sign(master_summary$GNPC_Q_estimate),
  "UKB-PD" = -log10(master_summary$UKB_p_val) * sign(master_summary$UKB_estimate)
) |>
  as.matrix()

color_fun <- circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
cohort <- rep(c("Chronos","PPMI","GNPC","UKB"), c(1,1,5,1))
cohort_cols <- c(
  Chronos = "#222222",
  PPMI    = "#555555",
  GNPC    = "#888888",
  UKB     = "#BBBBBB"
)
top_anno <- ComplexHeatmap::HeatmapAnnotation(
  Cohort = cohort,
  col = list(Cohort = cohort_cols),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = grid::unit(3, "mm")
)
hmap <- mat |>
  ComplexHeatmap::Heatmap(
    col = color_fun,
    column_split = cohort,
    top_annotation = top_anno,
    cluster_columns = FALSE,
    column_names_rot = 45,
    show_heatmap_legend = FALSE
  )

# Create a manual legend object
lgd = ComplexHeatmap::Legend(
  title = "-log10(p-value)",
  direction = "vertical",
  title_position = "topcenter",
  legend_height = grid::unit(0.20, "npc"),
  col_fun = color_fun,
  at = c(-5, 0, 5)
)

pdf(
  output_figure_file,
  width = 8.0,
  height = 6.0
)
ComplexHeatmap::draw(
  hmap,
  padding = grid::unit(c(2, 2, 2, 30), "mm")
)
ComplexHeatmap::draw(
  lgd,
  x = grid::unit(0.95, "npc"),
  y = grid::unit(0.4, "npc"),
  just = c("right", "bottom")
)
dev.off()
