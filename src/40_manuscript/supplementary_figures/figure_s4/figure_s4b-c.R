## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s4/figure_s4b-c.R")

ppmi_input_file <- glue::glue("{get_config()$input_dir}/name=PPMI_p9000_Plasma_MS/version=29985/dataset.rds")
olink_input_file <- glue::glue("{get_config()$input_dir}/name=Chronos_phase1_proteomics_olink_analysis_ready_three_tibble/version=29465/dataset.rds")
ukb_input_file <- glue::glue("{get_config()$input_dir}/name=ukb_pd-proteomics-ddc_itgav/version=none/dataset.parquet")

output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s4/figure_s4b-c.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)


#####
# Process the PPMI Data
#####
p9k <- ppmi_input_file |>
  readr::read_rds()

# Extract person-level metadata
p9k[["pheno.don"]] <- p9k[["pheno.don"]] |>
  tibble::as_tibble() |>
  dplyr::select(
    PATNO,
    group = cohort
  ) |>
  dplyr::mutate(
    group = dplyr::recode_values(
      group,
      "PD" ~ "Case",
      "control" ~ "Control"
    )
  ) |>
  # Retain only those people designated as PD or control
  dplyr::filter(
    group %in% c("Case", "Control")
  )

# Extract sample-level metadata
p9k[["pheno.samp"]] <- p9k[["pheno.samp"]] |>
  tibble::as_tibble() |>
  dplyr::select(
    SampleId,
    PATNO,
    EVENT_ID
  ) |>
  dplyr::filter(
    EVENT_ID %in% c("BL", "V02", "V04", "V06", "V10")
  ) |>
  # Restrict samples to just those we retained
  dplyr::semi_join(
    p9k[["pheno.don"]],
    by = "PATNO"
  )

# Extract proteomics data for the samples of interest
p9k[["dataProt"]] <- p9k[["dataProt"]] |>
  tibble::as_tibble(rownames = "ID") |>
  dplyr::select(
    SampleId = ID,
    P9000_Plasma_ITGAV_NPX,
    P9000_Plasma_DDC_NPX
  ) |>
  dplyr::semi_join(
    p9k[["pheno.samp"]],
    by = "SampleId"
  )

ppmi_plotting_data <- p9k[["pheno.samp"]] |>
  dplyr::left_join(
    p9k[["pheno.don"]],
    by = "PATNO"
  ) |>
  dplyr::left_join(
    p9k[["dataProt"]],
    by = "SampleId"
  )

.make_summary_string <- function(
  data,
  id_col,
  group_col,
  sample_col,
  case_label = "Case",
  control_label = "Control"
) {
  id_col <- rlang::enquo(id_col)
  group_col <- rlang::enquo(group_col)
  sample_col <- rlang::enquo(sample_col)

  # Count unique patients per group
  patient_counts <- data |>
    dplyr::distinct(!!id_col, !!group_col) |>
    dplyr::count(!!group_col) |>
    dplyr::mutate(
      group = dplyr::case_when(
        stringr::str_detect(!!group_col, stringr::regex("case|pd", ignore_case = TRUE)) ~ case_label,
        stringr::str_detect(!!group_col, stringr::regex("control", ignore_case = TRUE)) ~ control_label,
        .default = as.character(!!group_col)
      )
    ) |>
    dplyr::select(group, n) |>
    tibble::deframe()

  total_samples <- data |>
    dplyr::distinct(!!sample_col) |>
    nrow()

  glue::glue(
    "{patient_counts[case_label]} PDs and {patient_counts[control_label]} controls\n",
    "{scales::label_comma()(total_samples)} samples"
  )
}

ppmi_point_size <- 0.25
ppmi_point_alpha <- 0.3
ppmi_line_width <- 0.1
ppmi_line_alpha <- 0.25
ppmi_lowess_width <- 1.0

p_ppmi_ddc <- ppmi_plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = EVENT_ID,
      y = P9000_Plasma_DDC_NPX,
      group = PATNO,
      colour = group
    )
  ) +
  ggplot2::geom_point(
    size = ppmi_point_size,
    alpha = ppmi_point_alpha
  ) +
  ggplot2::geom_line(
    linewidth = ppmi_line_width,
    alpha = ppmi_line_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = group
    ),
    method = "loess",
    linewidth = ppmi_lowess_width,
    linetype = "solid",
    formula = y ~ x
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Event ID",
    y = "NPX"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

p_ppmi_itgav <- ppmi_plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = EVENT_ID,
      y = P9000_Plasma_ITGAV_NPX,
      group = PATNO,
      colour = group
    )
  ) +
  ggplot2::geom_point(
    size = ppmi_point_size,
    alpha = ppmi_point_alpha
  ) +
  ggplot2::geom_line(
    linewidth = ppmi_line_width,
    alpha = ppmi_line_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = group
    ),
    method = "loess",
    linewidth = ppmi_lowess_width,
    linetype = "solid",
    formula = y ~ x
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Event ID",
    y = "NPX"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom"
  )

#####
# Process the UKB Data
#####
ukb_plotting_data <- ukb_input_file |>
  arrow::open_dataset() |>
  dplyr::collect()

ukb_point_size <- 0.25
ukb_point_alpha <- 0.3
ukb_line_width <- 0.1
ukb_line_alpha <- 0.25
ukb_lowess_width <- 1.0

p_ukb_ddc <- ukb_plotting_data |>
  dplyr::filter(
    gene_symbol == "DDC"
  ) |>
  dplyr::mutate(
    group = dplyr::recode_values(
      group,
      "case.eid" ~ "Case",
      "control.eid" ~ "Control",
      unmatched = "error"
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = YEARS_SINCE_ONSET,
      y = NPX,
      colour = group
    )
  ) +
  ggplot2::coord_cartesian(
    xlim = c(NA, 5)
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      group = Participant_ID_hash
    ),
    alpha = ukb_line_alpha,
    linewidth = ukb_line_width
  ) +
  ggplot2::geom_point(
    size = ukb_point_size,
    alpha = ukb_point_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = group),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = ukb_lowess_width
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Years Since PD Onset",
    y = "NPX"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

p_ukb_itgav <- ukb_plotting_data |>
  dplyr::filter(
    gene_symbol == "ITGAV"
  ) |>
  dplyr::mutate(
    group = dplyr::recode_values(
      group,
      "case.eid" ~ "Case",
      "control.eid" ~ "Control",
      unmatched = "error"
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = YEARS_SINCE_ONSET,
      y = NPX,
      colour = group
    )
  ) +
  ggplot2::coord_cartesian(
    xlim = c(NA, 5)
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      group = Participant_ID_hash
    ),
    alpha = ukb_line_alpha,
    linewidth = ukb_line_width
  ) +
  ggplot2::geom_point(
    size = ukb_point_size,
    alpha = ukb_point_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = group),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = ukb_lowess_width
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Years Since PD Onset",
    y = "NPX"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

#####
# Process the Chronos-PD Data
#####
chronos_olink_three_tibble <- olink_input_file |>
  readr::read_rds()

chronos_olink_three_tibble[["prot_data"]] <- chronos_olink_three_tibble[["prot_data"]] |>
  dplyr::select(
    ID,
    OID44598, # DDC
    OID43878 # ITGAV
  )

chronos_olink_three_tibble[["meta_data"]] <- chronos_olink_three_tibble[["meta_data"]] |>
  dplyr::select(
    ID,
    group,
    years_since_pd_onset,
    patient_id
  )

chronos_plotting_data <- chronos_olink_three_tibble[["prot_data"]] |>
  dplyr::inner_join(
    chronos_olink_three_tibble[["meta_data"]],
    by = "ID"
  ) |>
  dplyr::mutate(
    group = dplyr::recode_values(
      group,
      "pd" ~ "Case",
      "control" ~ "Control",
      unmatched = "error"
    )
  )

chronos_point_size <- 0.25
chronos_point_alpha <- 0.3
chronos_line_width <- 0.1
chronos_line_alpha <- 0.25
chronos_lowess_width <- 1.0

p_chronos_ddc <- chronos_plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = OID44598,
      colour = group
    )
  ) +
  ggplot2::geom_line(
    ggplot2::aes(group = patient_id),
    alpha = chronos_line_alpha,
    linewidth = chronos_line_width
  ) +
  ggplot2::geom_point(
    size = chronos_point_size,
    alpha = chronos_point_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = group),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = chronos_lowess_width
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Years Since PD Onset",
    y = "NPX"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )


p_chronos_itgav <- chronos_plotting_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = years_since_pd_onset,
      y = OID43878,
      colour = group
    )
  ) +
  ggplot2::geom_line(
    ggplot2::aes(group = patient_id),
    alpha = chronos_line_alpha,
    linewidth = chronos_line_width
  ) +
  ggplot2::geom_point(
    size = chronos_point_size,
    alpha = chronos_point_alpha
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = group),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = chronos_lowess_width
  ) +
  ggplot2::scale_color_manual(
    name = "Group",
    values = c(
      Case = "darkorange",
      Control = "dodgerblue"
    )
  ) +
  ggplot2::labs(
    x = "Years Since PD Onset",
    y = "NPX"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

label_top <- cowplot::ggdraw() +
  cowplot::draw_label(
    "DDC",
    x = 0,
    y = 0.65,
    hjust = 0,
    fontface = "bold",
    size = 12
  ) +
  cowplot::draw_label(
    "Aromatic L-amino\nacid decarboxylase",
    x = 0,
    y = 0.53,
    hjust = 0,
    size = 10
  )

label_bottom <- cowplot::ggdraw() +
  cowplot::draw_label(
    "ITGAV",
    x = 0,
    y = 0.62,
    hjust = 0,
    fontface = "bold",
    size = 12
  ) +
  cowplot::draw_label(
    "Integrin Alpha V",
    x = 0,
    y = 0.53,
    hjust = 0,
    size = 10
  )

ppmi_summary <- .make_summary_string(
  ppmi_plotting_data,
  PATNO,
  group,
  SampleId,
  case_label = "PD",
  control_label = "Control"
)

ukb_summary <- .make_summary_string(
  ukb_plotting_data,
  Participant_ID_hash,
  group,
  Sample_ID_hash,
  case_label = "Case",
  control_label = "Control"
)

chronos_summary <- .make_summary_string(
  chronos_plotting_data,
  patient_id,
  group,
  ID,
  case_label = "Case",
  control_label = "Control"
)

ppmi_label <- cowplot::ggdraw() +
  cowplot::draw_label(
    ppmi_summary,
    hjust = 0.5,
    size = 10
  )
ukb_label <- cowplot::ggdraw() +
  cowplot::draw_label(
    ukb_summary,
    hjust = 0.5,
    size = 10
  )
chronos_label <- cowplot::ggdraw() +
  cowplot::draw_label(
    chronos_summary,
    hjust = 0.5,
    size = 10
  )

row0 <- cowplot::plot_grid(
  ppmi_label + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = ""),
  ukb_label + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = "", y = ""),
  chronos_label + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = "", y = ""),
  ncol = 3,
  align = "h"
)
row1 <- cowplot::plot_grid(
  p_ppmi_ddc + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = ""),
  p_ukb_ddc + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = "", y = ""),
  p_chronos_ddc + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = "", y = ""),
  ncol = 3,
  align = "h"
)
row2 <- cowplot::plot_grid(
  p_ppmi_itgav + ggplot2::theme(legend.position = "none"),
  p_ukb_itgav + ggplot2::theme(legend.position = "none") + ggplot2::labs(y = ""),
  p_chronos_itgav + ggplot2::theme(legend.position = "none") + ggplot2::labs(y = ""),
  ncol = 3,
  align = "h"
)
row3 <- cowplot::plot_grid(
  NULL,
  cowplot::get_legend(p_ppmi_ddc),
  NULL,
  ncol = 3,
  align = "h"
)

row0_labeled <- cowplot::plot_grid(
  NULL,
  row0,
  ncol = 2,
  rel_widths = c(0.13, 1)
)
row1_labeled <- cowplot::plot_grid(
  label_top,
  row1,
  ncol = 2,
  rel_widths = c(0.13, 1)
)
row2_labeled <- cowplot::plot_grid(
  label_bottom,
  row2,
  ncol = 2,
  rel_widths = c(0.13, 1)
)
row3_labeled <- cowplot::plot_grid(
  NULL,
  row3,
  ncol = 2,
  rel_widths = c(0.18, 1)
)

p <- cowplot::plot_grid(
  row0_labeled,
  row1_labeled,
  row2_labeled,
  row3_labeled,
  ncol = 1,
  align = "v",
  rel_heights = c(0.2, 1.0, 1.0, 0.2)
)

# Save plot
ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 11.0,
  height = 5.0,
  units = "in",
  dpi = 600,
  bg = "white"
)
