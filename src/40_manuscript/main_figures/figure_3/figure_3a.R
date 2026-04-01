## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_3/figure_3a.R")

input_lmm_result_files <- list(
  alamar = glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet"),
  biognosys = glue::glue("{get_config()$input_dir}/name=chronos_biognosys_all-pd_master_summary/version=29864/dataset.parquet"),
  olink = glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet"),
  somalogic = glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet")
)

input_proteomics_files <- list(
  alamar = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_alamar_analysis_ready_three_tibble/version=29649/dataset.rds"),
  biognosys = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_biognosys_analysis_ready_three_tibble/version=29639/dataset.rds"),
  olink = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds"),
  somalogic = glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_somalogic_analysis_ready_three_tibble/version=29650/dataset.rds")
)

output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_3a.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

###########################################
#### Load Phase 1 Significant Proteins ####
###########################################
signifproteins_phase1 <- purrr::map(
  input_lmm_result_files,
  ~ {
    arrow::open_dataset(.x) |>
      dplyr::select(
        probe_id,
        estimate,
        cohens_d,
        p_val,
        q_val,
        uniprot_id,
        gene_symbol
      ) |>
      dplyr::filter(
        q_val < 0.05
      ) |>
      dplyr::collect()
  }
) |>
  dplyr::bind_rows() |>
  dplyr::distinct()

############################################
#### Load/Process Three-Tibble Datasets ####
############################################
message("Loading proteomics datasets")
proteomics_data <- purrr::map(
  input_proteomics_files,
  ~ {
    three_tbl <- readr::read_rds(.x)
    tbl <- three_tbl$meta_data |>
      # Get the minimum metadata required for this script
      dplyr::select(
        ID,
        group,
        years_since_pd_onset
      ) |>
      dplyr::inner_join(
        three_tbl$prot_data |>
          # Only get a subset of data for the proteins we need
          dplyr::select(
            ID,
            dplyr::any_of(signifproteins_phase1[["probe_id"]])
          ),
        by = "ID"
      )

    return(tbl)
  },
  progress = TRUE
)

######################################################
#### Generate Loess-Smoothed Protein Trajectories ####
######################################################
predictions <- proteomics_data |>
  purrr::imap(
    ~ {
      message(glue::glue("Processing proteomics dataset: { .y }"))

      predict_all_proteins(
        .x,
        protein_cols = signifproteins_phase1[["probe_id"]],
        years_pred = seq(-7, 7, 0.1)
      )
    }
  )

##################################
#### Prepare Data for Heatmap ####
##################################
mat <- predictions |>
  purrr::compact() |>
  dplyr::bind_rows(.id = "platform") |>
  tidyr::unnest_longer(c(years_pred, difference)) |>
  dplyr::select(probe_id, years_pred, difference) |>
  dplyr::filter(
    !is.na(difference)
  ) |>
  tidyr::pivot_wider(
    names_from  = years_pred,
    values_from = difference
  ) |>
  tibble::column_to_rownames("probe_id") |>
  as.matrix()

##########################
#### Generate Heatmap ####
##########################
#  Recreate the row ordering that gplots::heatmap.2() generates
.heatmap.2.clustering <- function(mat, fun = rowMeans) {
  Rowv <- fun(mat, na.rm = TRUE)
  distr <- dist(mat)
  hcr <- hclust(distr)
  reorderfun <- function(d, w) reorder(d, w)
  ddr <- as.dendrogram(hcr)
  ddr <- reorderfun(ddr, Rowv) |> rev()

  return(ddr)
}
.col_labels <- function(mat, breaks = NULL, pos_prefix = "", suffix = "") {
  if (is.null(breaks)) {
    breaks <- -20:20
  }

  curr_breaks <- as.numeric(colnames(mat))
  new_labels <- rep("", ncol(mat))

  for (b in breaks) {
    idx <- which(abs(curr_breaks - b) == 0)
    new_labels[idx] <- if (b == 0) {
      "0"
    } else if (b > 0) {
      paste0(pos_prefix, b, suffix)
    } else {
      paste0(b, suffix)
    }
  }

  return(new_labels)
}

row_title_text <- glue::glue(
  "<span style='font-size:12pt'>{n_probe} Probe IDs</span><br>
   <span style='font-size:8pt'>({n_uniprot} UniProt IDs)</span>",
  n_probe = nrow(mat),
  n_uniprot = signifproteins_phase1 |>
    dplyr::filter(probe_id %in% rownames(mat)) |>
    dplyr::distinct(uniprot_id) |>
    nrow()
)

p <- mat |>
  ComplexHeatmap::Heatmap(
    name = "difference",  # Legend title
    col = circlize::colorRamp2(
      breaks = c(-0.5, 0, 0.5),
      colors = c("deepskyblue", "black", "yellow")
    ),

    # Clustering
    cluster_rows = .heatmap.2.clustering(mat),
    cluster_columns = FALSE,
    column_names_rot = 0,

    # Row settings
    row_title = ComplexHeatmap::gt_render(row_title_text),
    show_row_names = FALSE,

    # Column settings
    column_title = "Years to estimated PD onset",
    column_title_side = "bottom",
    column_labels = .col_labels(mat, breaks = c(-7, 0, 7), pos_prefix = "+"),

    # Legend
    heatmap_legend_param = list(
      title = "",
      direction = "horizontal"
    )
  )

pdf(output_figure_file, width = 5.0, height = 4.5)
p |>
  ComplexHeatmap::draw(
    heatmap_legend_side = "bottom"
  )
dev.off()
