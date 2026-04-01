## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/main_figures/figure_3/figure_3b.R")

# Define input files
alamar_all <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet")
alamar_pre <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_pre_onset-pd_master_summary/version=29813/dataset.parquet")
alamar_post <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_post_onset-pd_master_summary/version=29814/dataset.parquet")
alamar_jm <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_joint_model-pd_master_summary/version=29821/dataset.parquet")
soma_all <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet")
soma_pre <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_pre_onset-pd_master_summary/version=29807/dataset.parquet")
soma_post <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_post_onset-pd_master_summary/version=29808/dataset.parquet")
soma_jm <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet")
olink_all <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet")
olink_pre <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_pre_onset-pd_master_summary/version=29810/dataset.parquet")
olink_post <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_post_onset-pd_master_summary/version=29811/dataset.parquet")
olink_jm <- glue::glue("{get_config()$input_dir}/name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet")

# Define output file
output_figure_file <- glue::glue("{get_config()$output_dir}/main_figures/figure_3b.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


## merge results from multiple platforms ####
### mixed model results from all samples ####
all_combined <- soma_all |>
  arrow::read_parquet() |>
  dplyr::select(
    estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
  ) |>
  dplyr::distinct() |>
  dplyr::bind_rows(
    olink_all |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  ) |>
  dplyr::bind_rows(
    alamar_all |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  )

all_combined <- all_combined |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::mutate(adjP_combined = p.adjust(p_val, method = "BH"))

### mixed model from preDx samples  ####
preDx_combined <- soma_pre |>
  arrow::read_parquet() |>
  dplyr::select(
    estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
  ) |>
  dplyr::distinct() |>
  dplyr::bind_rows(
    olink_pre |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  ) |>
  dplyr::bind_rows(
    alamar_pre |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  )

preDx_combined <- preDx_combined |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::mutate(adjP_combined = p.adjust(p_val, method = "BH"))

### mixed model from postDx samples  ####
postDx_combined <- soma_post |>
  arrow::read_parquet() |>
  dplyr::select(
    estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
  ) |>
  dplyr::distinct() |>
  dplyr::bind_rows(
    olink_post |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  ) |>
  dplyr::bind_rows(
    alamar_post |>
      arrow::read_parquet() |>
      dplyr::select(
        estimate, cohens_d, p_val, q_val, uniprot_id, gene_symbol
      ) |>
      dplyr::distinct()
  )

postDx_combined <- postDx_combined |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::mutate(adjP_combined = p.adjust(p_val, method = "BH"))

#### table for scatter plot of effect size in commonly measured proteins
#### between whole sample dataset and preDx sample dataset ####
commone_combined_all_preDx <- all_combined |>
  dplyr::mutate(sign_sig_all = -log10(p_val) * sign(estimate)) |>
  dplyr::rename(
    estimate_all      = estimate,
    cohens_d_all      = cohens_d,
    p_val_all         = p_val,
    adjP_combined_all = q_val          # use original q_value
  ) |>
  dplyr::select(
    sign_sig_all, estimate_all, cohens_d_all, p_val_all,
    adjP_combined_all, uniprot_id, gene_symbol
  ) |>
  dplyr::inner_join(
    preDx_combined |>
      dplyr::mutate(sign_sig_pre = -log10(p_val) * sign(estimate)) |>
      dplyr::rename(
        estimate_pre      = estimate,
        cohens_d_pre      = cohens_d,
        p_val_pre         = p_val,
        adjP_combined_pre = q_val
      ) |>
      dplyr::select(
        sign_sig_pre, estimate_pre, cohens_d_pre,
        p_val_pre, adjP_combined_pre, uniprot_id
      ),
    by = "uniprot_id"
  ) |>
  dplyr::mutate(
    significance = ifelse(
      adjP_combined_all < 0.05 & adjP_combined_pre < 0.05,
      "both_significant",
      ifelse(
        adjP_combined_all < 0.05 & adjP_combined_pre >= 0.05,
        "whole_significant",
        ifelse(
          adjP_combined_all >= 0.05 & adjP_combined_pre < 0.05,
          "preDx_significant",
          "not_significant"
        )
      )
    )
  )

#### table for risk biomarker forest plot ####
risk_combined <- soma_jm |>
  arrow::read_parquet() |>
  dplyr::filter(rhat < 1.05) |>
  dplyr::select(
    probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val,
    rhat, risk, risk_lower, risk_upper
  ) |>
  dplyr::mutate(platform = "SomaScanV5") |>
  dplyr::bind_rows(
    olink_jm |>
      arrow::read_parquet() |>
      dplyr::filter(rhat < 1.05) |>
      dplyr::select(
        probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val,
        rhat, risk, risk_lower, risk_upper
      ) |>
      dplyr::mutate(platform = "Olink_HT")
  ) |>
  dplyr::bind_rows(
    alamar_jm |>
      arrow::read_parquet() |>
      dplyr::filter(rhat < 1.05) |>
      dplyr::select(
        probe_id, uniprot_id, gene_symbol, estimate, p_val, q_val,
        rhat, risk, risk_lower, risk_upper
      ) |>
      dplyr::mutate(platform = "Alamarbio_CNS")
  )

risk_combined <- risk_combined |>
  dplyr::group_by(uniprot_id) |>
  dplyr::slice_min(p_val) |>
  dplyr::ungroup() |>
  dplyr::mutate(gene_symbol = reorder(gene_symbol, risk))

## volcano plots and scatter plot ####
p1 <- preDx_combined |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = estimate,
      y = -log10(p_val)
    )
  ) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::xlim(-0.9, 0.9) +
  ggplot2::ylim(0, 12) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text  = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 15, face = "bold"),
    aspect.ratio = 1
  ) +
  ggplot2::labs(
    title = "Pre-PD Onset Samples\n(n=1,748 from 481 individuals)",
    x = "log2FC between PD and control"
  )

p2 <- commone_combined_all_preDx |>
  dplyr::mutate(
    plot_order = dplyr::case_when(
      significance == "whole_significant"  ~ 1,
      significance == "preDx_significant"  ~ 2,
      significance == "both_significant"   ~ 3,
      TRUE                                 ~ 0
    )
  ) |>
  dplyr::arrange(plot_order) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = -log10(p_val_all),
      y = -log10(p_val_pre),
      color = significance
    )
  ) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::scale_color_manual(
    values = c(
      "whole_significant"  = "purple",
      "preDx_significant"  = "blue",
      "both_significant"   = "red"
    ),
    na.value = "grey",
    guide = "none",
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    lty = 2
  ) +
  ggrepel::geom_label_repel(
    data = commone_combined_all_preDx |>
      dplyr::filter(
        significance %in% c("both_significant", "preDx_significant") | -log10(p_val_all) > 7
      ),
    ggplot2::aes(
      label = gene_symbol
    ),
    fill = scales::alpha("white", 0.6),
    color = "black",
    size = 3,
    label.size = NA,
    box.padding = 0.2,
    min.segment.length = 0.2,
    max.overlaps = Inf,
    seed = 1234
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text  = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 15, face = "bold"),
    aspect.ratio = 1
  ) +
  ggplot2::labs(
    x = "significance based on whole dataset",
    y = "significance based on preDx dataset"
  )

# volcano plot
p3 <- all_combined |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = estimate,
      y = -log10(p_val)
    )
  ) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::xlim(-0.9, 0.9) +
  ggplot2::ylim(0, 12) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text  = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 15, face = "bold"),
    aspect.ratio = 1
  ) +
  ggplot2::labs(
    title = "Pre- and post-PD Onset Samples\n(n=2,609 from 696 individuals)",
    x = "log2FC between PD and control"
  )

# Create blank placeholder plot for adjusting spacing between plot panels
blank <- cowplot::ggdraw()

# Plot all plot panels
p <- cowplot::plot_grid(
  p1,
  blank,
  p2,
  blank,
  p3,
  blank,
  nrow       = 1,
  rel_widths = c(0.8, 0.7, 1, 0.5, 0.8, 0.3)
)

ggplot2::ggsave(
  plot     = p,
  filename = output_figure_file,
  width    = 16,
  height   = 6,
  units    = "in"
)
