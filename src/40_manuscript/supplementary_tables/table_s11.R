## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_tables/table_s11.R")

# Define input files
input_chronos_data <- glue::glue("{get_config()$input_dir}/name=chronos_protein_ratio_list1_statistics/version=202510/dataset.parquet")
input_ukb_data <- glue::glue("{get_config()$input_dir}/name=ukb_validation_protein_ratio_lists_1_3/version=202510/dataset.parquet")
input_ppmi_data <- glue::glue("{get_config()$input_dir}/name=ppmi_protein_ratio_list1_statistics/version=202510/dataset.parquet")
input_gnpc_cohortC_data <- glue::glue("{get_config()$input_dir}/name=gnpc_cohort_c_protein_ratio_validation_results/version=202510/dataset.parquet")
input_gnpc_cohortJ_data <- glue::glue("{get_config()$input_dir}/name=gnpc_cohort_j_protein_ratio_validation_results/version=202510/dataset.parquet")
input_gnpc_cohortL_data <- glue::glue("{get_config()$input_dir}/name=gnpc_cohort_l_protein_ratio_validation_results/version=202510/dataset.parquet")
input_annotation_data <- glue::glue("{get_config()$input_dir}/name=protein_platform_uniprot_ids/version=29805/dataset.parquet")
protein_ratio_list_1_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_protein_ratio_list_1/version=29452/dataset.parquet")

# Define output file
output_table_file <- glue::glue("{get_config()$output_dir}/supplementary_tables/table_s11.tsv")
dir.create(dirname(output_table_file), showWarnings = FALSE, recursive = TRUE)

.keep_most_significant <- function(data, protein_pair_col = "protein_pair", p_value_col = "p_value") {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(protein_pair_col))) |>
    dplyr::slice_min(order_by = .data[[p_value_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

# ============================================================================
# Load and Process Data
# ============================================================================

# Load Chronos data
chronos_full <- input_chronos_data |>
  arrow::open_dataset() |>
  dplyr::rename(
    protein_pair = ratio
  ) |>
  dplyr::select(
    protein_pair,
    window_start,
    window_end,
    estimate,
    odds_ratio,
    or_lower_ci,
    or_upper_ci,
    p_value
  ) |>
  dplyr::collect() |>
  dplyr::mutate(cohort_label = glue::glue("Chronos_{window_start}_{window_end}"))


# Get window order for organizing columns (ascending by window_start to match R3)
window_order <- chronos_full |>
  dplyr::distinct(window_start, window_end, cohort_label) |>
  dplyr::arrange(window_start) |> # Ascending order like R3
  dplyr::pull(cohort_label)

# Pivot Chronos estimates (log odds) for heatmap using consistent labels
chronos_wide_est <- chronos_full |>
  dplyr::select(protein_pair, cohort_label, estimate, p_value) |>
  tidyr::pivot_wider(
    names_from = cohort_label,
    values_from = c(estimate, p_value),
    names_glue = "{.value}_{cohort_label}"
  )

# Separate estimates and p-values
estimate_cols <- paste0("estimate_", window_order)
pvalue_cols <- paste0("p_value_", window_order)

estimates_wide <- chronos_wide_est |>
  dplyr::select(protein_pair, dplyr::all_of(estimate_cols)) |>
  dplyr::rename_with(~ stringr::str_remove(.x, "^estimate_"), .cols = -protein_pair)

pvalues_wide <- chronos_wide_est |>
  dplyr::select(protein_pair, dplyr::all_of(pvalue_cols)) |>
  dplyr::rename_with(~ stringr::str_remove(.x, "^p_value_"), .cols = -protein_pair)

# Rename ratio column for consistency
estimates_wide <- estimates_wide |> dplyr::rename(ratio = protein_pair)
pvalues_wide <- pvalues_wide |> dplyr::rename(ratio = protein_pair)

# Load UKB data
ukb_full <- input_ukb_data |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::filter(stringr::str_detect(File, "list_1_"))

# Filter UKB for windows used in final plots (n_samples > 150)
ukb_filtered <- ukb_full |>
  dplyr::filter(window_id %in% c("[-9,-5)", "[-5,-1)")) |>
  dplyr::filter(n_samples > 150) |>
  dplyr::rename(
    protein_pair = pair_id,
    estimate = log_or,
    or_lower_ci = conf.low,
    or_upper_ci = conf.high,
    p_value = p.value
  ) |>
  dplyr::select(
    protein_pair, window_id, window_start, window_end, n_samples,
    estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value
  )

# Create consistent cohort labels for UKB (matching forest plot style)
ukb_filtered <- ukb_filtered |>
  dplyr::mutate(cohort_label = glue::glue("UKB_{window_start}_{window_end}"))

# UKB window order (ascending by window_start to match R3: UKB_-9_-5, then UKB_-5_-1)
ukb_window_order <- ukb_filtered |>
  dplyr::distinct(window_start, window_end, cohort_label) |>
  dplyr::arrange(window_start) |> # Ascending order like R3
  dplyr::pull(cohort_label)

# Pivot UKB estimates using consistent labels
ukb_wide_est <- ukb_filtered |>
  dplyr::select(protein_pair, cohort_label, estimate, p_value) |>
  tidyr::pivot_wider(
    names_from = cohort_label,
    values_from = c(estimate, p_value),
    names_glue = "{.value}_{cohort_label}"
  )

# Separate estimates and p-values
estimate_cols <- paste0("estimate_", ukb_window_order)
pvalue_cols <- paste0("p_value_", ukb_window_order)

ukb_estimates_wide <- ukb_wide_est |>
  dplyr::select(protein_pair, dplyr::all_of(estimate_cols)) |>
  dplyr::rename_with(~ stringr::str_remove(.x, "^estimate_"), .cols = -protein_pair)

ukb_pvalues_wide <- ukb_wide_est |>
  dplyr::select(protein_pair, dplyr::all_of(pvalue_cols)) |>
  dplyr::rename_with(~ stringr::str_remove(.x, "^p_value_"), .cols = -protein_pair)

# Join UKB data
estimates_wide <- estimates_wide |>
  dplyr::left_join(ukb_estimates_wide, by = c("ratio" = "protein_pair"))
pvalues_wide <- pvalues_wide |>
  dplyr::left_join(ukb_pvalues_wide, by = c("ratio" = "protein_pair"))

# Load validation cohorts
ppmi_full <- input_ppmi_data |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    odds_ratio = or,
    or_lower_ci = exp(or_lower_ci),
    or_upper_ci = exp(or_upper_ci)
  ) |>
  dplyr::select(protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value)

gnpc_cohortC_full <- input_gnpc_cohortC_data |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::filter(protein_pair_tbl_name == "chronos_PD_ratio_list1") |>
  dplyr::mutate(
    odds_ratio = or,
    or_lower_ci = exp(or_lower_ci),
    or_upper_ci = exp(or_upper_ci)
  ) |>
  dplyr::select(protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value)

gnpc_cohortJ_full <- input_gnpc_cohortJ_data |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::filter(protein_pair_tbl_name == "chronos_PD_ratio_list1") |>
  dplyr::mutate(
    odds_ratio = or,
    or_lower_ci = exp(or_lower_ci),
    or_upper_ci = exp(or_upper_ci)
  ) |>
  dplyr::select(protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value)

gnpc_cohortL_full <- input_gnpc_cohortL_data |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::filter(protein_pair_tbl_name == "chronos_PD_ratio_list1") |>
  dplyr::mutate(
    odds_ratio = or,
    or_lower_ci = exp(or_lower_ci),
    or_upper_ci = exp(or_upper_ci)
  ) |>
  dplyr::select(protein_pair, estimate, odds_ratio, or_lower_ci, or_upper_ci, p_value)

# ============================================================================
# Process Data for Table S11
# ============================================================================

# Define cohort configurations for joining estimates
cohort_configs <- list(
  "PPMI" = list(data = ppmi_full),
  "GNPC C" = list(data = gnpc_cohortC_full |> .keep_most_significant()),
  "GNPC J" = list(data = gnpc_cohortJ_full |> .keep_most_significant()),
  "GNPC L" = list(data = gnpc_cohortL_full |> .keep_most_significant())
)

# Loop through cohorts and join data
for (cohort_name in names(cohort_configs)) {
  cohort_data <- cohort_configs[[cohort_name]]$data

  # Join estimates
  estimates_wide <- estimates_wide |>
    dplyr::left_join(
      cohort_data |>
        dplyr::select(protein_pair, !!cohort_name := estimate),
      by = c("ratio" = "protein_pair")
    )

  # Join p-values (using backticks for names with spaces)
  pvalues_wide <- pvalues_wide |>
    dplyr::left_join(
      cohort_data |>
        dplyr::select(protein_pair, !!paste0(cohort_name, "_p_value") := p_value),
      by = c("ratio" = "protein_pair")
    )
}

# Verify alignment
stopifnot(all(estimates_wide$ratio == pvalues_wide$ratio))

# Manually fix ITGAV complex - external validation does not target this complex
pvalues_wide <- pvalues_wide |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("UKB"),
      ~ dplyr::if_else(ratio == "seq.15573.110_o_seq.20187.10", NA_real_, .x)
    )
  )
estimates_wide <- estimates_wide |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("UKB"),
      ~ dplyr::if_else(ratio == "seq.15573.110_o_seq.20187.10", NA_real_, .x)
    )
  )

# Get all protein pairs
all_rows <- pvalues_wide |> dplyr::pull(ratio)

# Function to create long format data for forest plots
.create_forest_data <- function(
    protein_pairs,
    chronos_full,
    ukb_filtered,
    ppmi_full,
    gnpc_cohortC_full,
    gnpc_cohortJ_full,
    gnpc_cohortL_full
) {
  forest_data_list <- list()

  for (pair in protein_pairs) {
    # Chronos data for this pair
    chronos_forest <- chronos_full |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::arrange(desc(window_start)) |>
      dplyr::mutate(
        cohort = glue::glue("Chronos_{window_start}_{window_end}"),
        cohort_group = "Chronos"
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # UKB data for this pair (filtered windows, ordered most recent first for forest plot)
    ukb_forest <- ukb_filtered |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::mutate(
        cohort = glue::glue("UKB_{window_start}_{window_end}"),
        cohort_group = "UKB"
      ) |>
      dplyr::arrange(desc(window_start)) |> # Most recent first
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # Ensure both UKB windows are present in correct order
    expected_ukb <- data.frame(
      protein_pair = pair,
      cohort = c("UKB_-5_-1", "UKB_-9_-5"),
      cohort_group = "UKB",
      window_start = c(-5, -9),
      window_end = c(-1, -5)
    )
    ukb_forest <- expected_ukb |>
      dplyr::left_join(ukb_forest |> dplyr::select(-protein_pair, -cohort_group, -window_start, -window_end),
                       by = "cohort"
      )

    # PPMI data for this pair
    ppmi_forest <- ppmi_full |>
      dplyr::filter(protein_pair == pair) |>
      dplyr::mutate(
        cohort = "PPMI",
        cohort_group = "PPMI",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no PPMI data, create empty row (like R4)
    if (nrow(ppmi_forest) == 0) {
      ppmi_forest <- data.frame(
        protein_pair = pair,
        cohort = "PPMI",
        cohort_group = "PPMI",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    # GNPC cohorts for this pair (keep most significant)
    gnpc_c_forest <- gnpc_cohortC_full |>
      dplyr::filter(protein_pair == pair) |>
      .keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC C",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC C data, create empty row (like R4)
    if (nrow(gnpc_c_forest) == 0) {
      gnpc_c_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC C",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    gnpc_j_forest <- gnpc_cohortJ_full |>
      dplyr::filter(protein_pair == pair) |>
      .keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC J",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC J data, create empty row (like R4)
    if (nrow(gnpc_j_forest) == 0) {
      gnpc_j_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC J",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    gnpc_l_forest <- gnpc_cohortL_full |>
      dplyr::filter(protein_pair == pair) |>
      .keep_most_significant() |>
      dplyr::mutate(
        cohort = "GNPC L",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA
      ) |>
      dplyr::select(
        protein_pair, cohort, cohort_group, window_start, window_end,
        odds_ratio, or_lower_ci, or_upper_ci, p_value
      )

    # If no GNPC L data, create empty row (like R4)
    if (nrow(gnpc_l_forest) == 0) {
      gnpc_l_forest <- data.frame(
        protein_pair = pair,
        cohort = "GNPC L",
        cohort_group = "GNPC",
        window_start = NA,
        window_end = NA,
        odds_ratio = NA,
        or_lower_ci = NA,
        or_upper_ci = NA,
        p_value = NA
      )
    }

    # Combine all cohorts for this pair (reverse order for coord_flip display, like R4)
    pair_data <- dplyr::bind_rows(
      gnpc_l_forest,
      gnpc_j_forest,
      gnpc_c_forest,
      ppmi_forest,
      ukb_forest,
      chronos_forest
    )

    forest_data_list[[pair]] <- pair_data
  }

  # Combine all pairs
  dplyr::bind_rows(forest_data_list)
}

# Create forest plot data for all pairs
forest_plot_all <- .create_forest_data(
  protein_pairs = all_rows,
  chronos_full,
  ukb_filtered,
  ppmi_full,
  gnpc_cohortC_full,
  gnpc_cohortJ_full,
  gnpc_cohortL_full
)

# ============================================================================
# Add Annotations and Create Table S11
# ============================================================================

# Load annotation data
annotation_table <- arrow::open_dataset(input_annotation_data) |>
  dplyr::collect()

somalogic2annotations <- annotation_table |>
  dplyr::filter(platform == "somalogic" & platform_version == "11k")
olink2annotations <- annotation_table |>
  dplyr::filter(platform == "olink" & platform_version == "5k")
alamar2annotations <- annotation_table |>
  dplyr::filter(platform == "alamar")

ids2annotations <- dplyr::bind_rows(
  somalogic2annotations,
  olink2annotations,
  alamar2annotations
)
# Manually fix ITGAV/ITGB3 as a complex
ids2annotations <- ids2annotations |>
  dplyr::filter(probe_id != "seq.20187.10") |>
  dplyr::bind_rows(
    tibble::tibble(
      platform = "somalogic",
      platform_version = "11k",
      probe_id = "seq.20187.10",
      uniprot_id = "P05106_P06756",
      gene_symbol = "ITGAV_ITGB3"
    )
  )

input_metadata <- protein_ratio_list_1_file |>
  arrow::open_dataset() |>
  dplyr::select(protein_pair) |>
  dplyr::collect() |>
  tidyr::separate(protein_pair, into = c("protein1", "protein2"), sep = "_o_", remove = FALSE) |>
  dplyr::mutate(
    platform1 = dplyr::case_when(
      stringr::str_starts(protein1, "OID") ~ "Olink",
      stringr::str_starts(protein1, "seq") ~ "Somalogic",
      .default = "AlamarBio"
    ),
    platform2 = dplyr::case_when(
      stringr::str_starts(protein2, "OID") ~ "Olink",
      stringr::str_starts(protein2, "seq") ~ "Somalogic",
      .default = "AlamarBio"
    )
  ) |>
  dplyr::left_join(
    ids2annotations |>
      dplyr::select(probe_id, gene_symbol, uniprot_id) |>
      dplyr::rename(protein1 = probe_id, gene_symbol1 = gene_symbol, uniprot1 = uniprot_id),
    by = "protein1"
  ) |>
  dplyr::left_join(
    ids2annotations |>
      dplyr::select(probe_id, gene_symbol, uniprot_id) |>
      dplyr::rename(protein2 = probe_id, gene_symbol2 = gene_symbol, uniprot2 = uniprot_id),
    by = "protein2"
  )

# Add gene symbols to estimates and pvalues
estimates_with_genes_unfiltered <- estimates_wide |>
  dplyr::left_join(
    input_metadata |> dplyr::select(protein_pair, gene_symbol1, gene_symbol2),
    by = c("ratio" = "protein_pair")
  ) |>
  dplyr::select(ratio, gene_symbol1, gene_symbol2, dplyr::everything())

pvalues_with_genes_unfiltered <- pvalues_wide |>
  dplyr::left_join(
    input_metadata |> dplyr::select(protein_pair, gene_symbol1, gene_symbol2),
    by = c("ratio" = "protein_pair")
  ) |>
  dplyr::select(ratio, gene_symbol1, gene_symbol2, dplyr::everything())

# Add gene symbols to forest plot data
forest_plot_all_with_genes <- forest_plot_all |>
  dplyr::left_join(
    input_metadata |> dplyr::select(protein_pair, gene_symbol1, gene_symbol2),
    by = "protein_pair"
  ) |>
  dplyr::select(protein_pair, gene_symbol1, gene_symbol2, dplyr::everything())

# Convert estimates and pvalues to long format
estimates_long_unfiltered <- estimates_with_genes_unfiltered |>
  tidyr::pivot_longer(
    cols = -c(ratio, gene_symbol1, gene_symbol2),
    names_to = "measure",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    study = dplyr::case_when(
      stringr::str_detect(measure, "^Chronos_") ~ "Chronos",
      stringr::str_detect(measure, "^UKB_") ~ "UKB",
      stringr::str_detect(measure, "^PPMI$") ~ "PPMI",
      stringr::str_detect(measure, "^GNPC ") ~ "GNPC",
      .default = "Other"
    ),
    window_info = dplyr::case_when(
      stringr::str_detect(measure, "^Chronos_") ~ stringr::str_extract(measure, "[-\\d]+_[-\\d]+$"),
      stringr::str_detect(measure, "^UKB_") ~ stringr::str_extract(measure, "[-\\d]+_[-\\d]+$"),
      .default = NA_character_
    )
  )

pvalues_long_unfiltered <- pvalues_with_genes_unfiltered |>
  tidyr::pivot_longer(
    cols = -c(ratio, gene_symbol1, gene_symbol2),
    names_to = "measure",
    values_to = "p_value"
  ) |>
  dplyr::mutate(
    measure = stringr::str_remove(measure, "_p_value$"),
    study = dplyr::case_when(
      stringr::str_detect(measure, "^Chronos_") ~ "Chronos",
      stringr::str_detect(measure, "^UKB_") ~ "UKB",
      stringr::str_detect(measure, "^PPMI$") ~ "PPMI",
      stringr::str_detect(measure, "^GNPC ") ~ "GNPC",
      .default = "Other"
    ),
    window_info = dplyr::case_when(
      stringr::str_detect(measure, "^Chronos_") ~ stringr::str_extract(measure, "[-\\d]+_[-\\d]+$"),
      stringr::str_detect(measure, "^UKB_") ~ stringr::str_extract(measure, "[-\\d]+_[-\\d]+$"),
      .default = NA_character_
    )
  )

# Combine estimates and p-values
heatmap_combined_unfiltered <- estimates_long_unfiltered |>
  dplyr::full_join(
    pvalues_long_unfiltered |> dplyr::select(ratio, measure, p_value),
    by = c("ratio", "measure")
  )

# Prepare forest data for export
forest_export_all <- forest_plot_all_with_genes |>
  dplyr::mutate(
    measure = cohort,
    # Extract window_info from measure name to match heatmap format
    window_info = dplyr::case_when(
      stringr::str_detect(cohort, "^Chronos_") ~ stringr::str_extract(cohort, "[-\\d]+_[-\\d]+$"),
      stringr::str_detect(cohort, "^UKB_") ~ stringr::str_extract(cohort, "[-\\d]+_[-\\d]+$"),
      .default = NA_character_
    )
  ) |>
  dplyr::select(
    ratio = protein_pair,
    gene_symbol1,
    gene_symbol2,
    measure,
    study = cohort_group,
    window_info,
    odds_ratio,
    or_lower_ci,
    or_upper_ci,
    p_value
  )

# Merge heatmap and forest data to create Table S11
comprehensive_export_all <- forest_export_all |>
  dplyr::full_join(
    heatmap_combined_unfiltered |> dplyr::select(ratio, measure, estimate),
    by = c("ratio", "measure")
  ) |>
  dplyr::select(
    ratio,
    study,
    measure,
    window_info,
    estimate,
    odds_ratio,
    or_lower_ci,
    or_upper_ci,
    p_value
  )

# Add annotation to create final Table S11
comprehensive_export_all <- input_metadata |>
  dplyr::select(protein_pair, platform1, platform2, uniprot1, uniprot2, gene_symbol1, gene_symbol2) |>
  dplyr::inner_join(
    comprehensive_export_all,
    dplyr::join_by(protein_pair == ratio)
  ) |>
  dplyr::rename(ratio = protein_pair) |>
  dplyr::select(
    ratio,
    `platform proteoform 1` = platform1,
    `platform proteoform 2` = platform2,
    `uniprot ID proteoform 1` = uniprot1,
    `uniprot ID proteoform 2` = uniprot2,
    `Gene symbol proteoform 1` = gene_symbol1,
    `Gene symbol proteoform 2` = gene_symbol2,
    Study = study,
    Cohort = measure,
    `Analysis Window` = window_info,
    Estimate = estimate,
    odds_ratio,
    or_lower_ci,
    or_upper_ci,
    Pvalue = p_value
  ) |>
  dplyr::mutate(
    `platform proteoform 1` = dplyr::recode_values(
      `platform proteoform 1`,
      "AlamarBio" ~ "AlamarBio-CNS-120",
      "Olink" ~ "Olink-HT-5k",
      "Somalogic" ~ "Somalogic-11k"
    ),
    `platform proteoform 2` = dplyr::recode_values(
      `platform proteoform 2`,
      "AlamarBio" ~ "AlamarBio-CNS-120",
      "Olink" ~ "Olink-HT-5k",
      "Somalogic" ~ "Somalogic-11k"
    )
  )

comprehensive_export_all |>
  readr::write_tsv(file = output_table_file)
