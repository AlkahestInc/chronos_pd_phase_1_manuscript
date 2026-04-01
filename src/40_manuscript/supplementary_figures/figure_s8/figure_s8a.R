## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s8/figure_s8a.R")

input_dataset_files <- list(
  "Chronos-PD" = glue::glue("{get_config()$input_dir}/name=chronos_pd-tornado_plot/version=1/dataset.parquet"),
  "PPMI" = glue::glue("{get_config()$input_dir}/name=PPMI_Curated_Data_Cut_Public/version=20250714/dataset.parquet"),
  "UKB-PD" = glue::glue("{get_config()$input_dir}/name=uk_pd_tornado_data/version=none/dataset.parquet")
)


output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s8/figure_s8a.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

#####
# Load Data
#####
plotting_data <- list()

plotting_data[["Chronos-PD"]] <- input_dataset_files[["Chronos-PD"]] |>
  arrow::open_dataset() |>
  dplyr::select(
    x = DELTA_YEARS,  # Years from disease onset (x-axis)
    y = PERSON_ID,    # Patient identifier (y-axis, already factored with proper order)
  ) |>
  dplyr::mutate(
    group = "Sample"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    y = forcats::fct_reorder(as.factor(y), x, .fun = max)
  )

plotting_data[["UKB-PD"]] <- input_dataset_files[["UKB-PD"]] |>
  arrow::open_dataset() |>
  dplyr::select(
    x = YEARS_SINCE_ONSET,      # Years from disease onset (x-axis)
    y = Participant_ID_hash,    # Patient identifier (y-axis, already factored with proper order)
  ) |>
  dplyr::mutate(
    group = "Sample"
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    y = forcats::fct_reorder(as.factor(y), x, .fun = max)
  )


plotting_data[["PPMI"]] <- input_dataset_files[["PPMI"]] |>
  arrow::open_dataset() |>
  dplyr::select(
    PATNO,
    ageonset,
    agediag,
    visit_date,
    age_at_visit,
    PRIMDIAG,
    EVENT_ID
  ) |>
  dplyr::group_by(PATNO) |>
  dplyr::collect() |>
  dplyr::filter(
    ! EVENT_ID %in% c(
      "PW" # Premature Withdrawal visit
      #"ST" #	Symptomatic Therapy visit
    ),
    # Remove patients that do not have any PRIMDIAG==1
    any(PRIMDIAG == 1)
  ) |>
  dplyr::mutate(
    age_first_primdiag_visit = min(age_at_visit[PRIMDIAG == 1], na.rm = TRUE)
  ) |>
  dplyr::mutate(
    #diff_visit_diag = age_at_visit - agediag,
    diff_visit_diag = age_at_visit - min(age_first_primdiag_visit, agediag, na.rm = TRUE),
    LAST = max(diff_visit_diag),

    HAS_DIAGNOSIS = !is.na(age_first_primdiag_visit)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(
    HAS_DIAGNOSIS
  ) |>
  dplyr::select(
    x = diff_visit_diag,
    y = PATNO,
  ) |>
  dplyr::mutate(
    group = "Sample",
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    y = forcats::fct_reorder(as.factor(y), x, .fun = max)
  )


#####
# Generate tornado plots
#####
layer_config <- tibble::tribble(
  ~group,      ~size,  ~alpha,  ~colour,      ~shape,
  "Sample",    0.25,   0.25,    "#5580B0",    "circle small",
)

plots <- plotting_data |>
  purrr::imap(
    \(data, cohort) {
      cat("Plotting:", cohort, "\n")

      ylab <- switch(
        cohort,
        "Chronos-PD" = "Donors",
        "UKB-PD" = "Participants",
        "PPMI" = "Patients"
      )

      p <- data |>
        rwd_plot(
          layer_config = layer_config,
          title = cohort,
          subtitle = "",
          caption = "",
          xlab = "Years Since Estimated Disease Onset",
          ylab = ylab
        ) +
        ggplot2::coord_cartesian(
          xlim = c(-15L, 15L),
          ylim = c(NA, NA)
        ) +
        ggplot2::theme(
          legend.position = "none"
        )

      return(p)
    }
  )

p <- cowplot::plot_grid(
  plotlist = plots,
  ncol = length(plots)
)


ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  width = 12,
  height = 3.0,
  units = "in",
  dpi = 600,
  bg = "white"
)
