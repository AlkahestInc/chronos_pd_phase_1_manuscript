## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript
# log_dirs <- purrr::map(unlist(strsplit("output/soma_olink_alamar,output/aws_renv,output/aws_conda", ",")), here::here)

here::i_am("scripts/R/run_performance.R")

#####
# Define commandline arguments
#####
option_list <- list(
  optparse::make_option(c("-d", "--dir"), type = "character", help = "A comma separated list of directories containing log files for parsing"),
  optparse::make_option(c("-o", "--out"), type = "character", help = "The filename of the PDF output file"),
  optparse::make_option(c("-w", "--width"), type = "double", default = 8, help = "The width of the output PDF in inches [default: %default]"),
  optparse::make_option(c("-H", "--height"), type = "double", default = 6, help = "The height of the output PDF in inches [default: %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
log_dirs <- purrr::map(unlist(strsplit(opt[["dir"]], ",")), here::here)
#####

# function to extract times from a single log file
extract_data <- function(file) {
  # Setup return object with zero rows
  d <- tibble::tibble(
    file = character(),
    job_id = integer(),
    cores = integer(),
    pair_count = integer(),
    wall_time_start = lubridate::POSIXct(),
    wall_time_end = lubridate::POSIXct()
  )

  # Function for extracting timestamp from line
  extract_time <- function(x) {
    stringr::str_match(x, "^\\s*(?:\\S+) \\[(.*?)\\]")[, 2] |>
      lubridate::parse_date_time(orders = "ymd HMS", quiet = TRUE, tz = "UTC")
  }

  # Extract lines of interest
  lines <- readr::read_lines(file)
  job_line <- stringr::str_subset(lines, "Slurm JOB_ID:")
  file_line <- stringr::str_subset(lines, "File:")
  cores_line <- stringr::str_subset(lines, "Available cores:")
  pair_count_line <- stringr::str_subset(lines, "Number of pairs:")

  # Construct results object for this log file
  r <- tibble::tibble(
    file = basename(file),
    job_id = as.integer(stringr::str_match(
      job_line,
      "Slurm JOB_ID:\\s*(\\d+)"
    )[, 2]),
    cores = as.integer(stringr::str_match(
      cores_line,
      "Available cores:\\s*(\\d+)"
    )[, 2]),
    pair_count = as.integer(stringr::str_match(
      pair_count_line,
      "Number of pairs:\\s*(\\d+)"
    )[, 2]),
    wall_time_start = extract_time(job_line),
    wall_time_end = extract_time(file_line)
  )

  d |>
    dplyr::bind_rows(
      r
    )
}

results <- purrr::map_dfr(log_dirs, function(dir) {
  fs::dir_ls(dir, regexp = "\\.log$", recurse = TRUE) |>
    purrr::map_dfr(extract_data) |>
    dplyr::mutate(
      dir = stringr::str_remove(!!dir, stringr::fixed(paste0(getwd(), "/"))),
      wall_time = lubridate::as.duration(wall_time_end - wall_time_start)
    )
}) |>
  dplyr::mutate(
    cputime_per_ratio = wall_time * cores / pair_count
  )

# counts <- results |>
#   dplyr::filter(pair_count > 0) |>
#   dplyr::count(dir)

p <- results |>
  dplyr::filter(
    pair_count > 0
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = pair_count,
      y = cputime_per_ratio,
      #group = dir,
      colour = dir
    )
  ) +
  ggplot2::geom_point(
    alpha = 0.2
  ) +
  ggplot2::geom_smooth(
    data = ~ dplyr::filter(.x, pair_count >= 10),
    method = "loess",
    formula = "y ~ x",
    se = FALSE
  ) +
  ggplot2::labs(
    x = "Job Size\n(number of ratios)",
    y = "Allocated CPU Time\n(seconds per ratio)",
    #title = ""
  ) +
  ggplot2::scale_x_continuous(
    transform = "log10",
    labels = scales::label_comma(),
    guide = "axis_logticks",
    limits = c(1, 100000),
    #limits = c(5000, 17000)
  ) +
  ggplot2::scale_y_continuous(
    transform = "log10",
    labels = scales::label_comma(),
    guide = "axis_logticks",
    limits = c(0.1, 1000),
    #limits = c(0.3, 1)
  )

p2 <- results |>
  dplyr::filter(
    pair_count > 0
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = pair_count,
      y = wall_time,
      #group = dir,
      colour = dir
    )
  ) +
  ggplot2::geom_point(
    shape = "bullet",
    alpha = 0.1,
    #color = dir
  ) +
  ggplot2::geom_smooth(
    data = ~ dplyr::filter(.x, pair_count >= 10),
    method = "loess",
    formula = "y ~ x",
    se = FALSE,
    #color = dir
  ) +
  # ggplot2::geom_text(
  #   data = counts,
  #   ggplot2::aes(
  #     x = 1500,
  #     y = max(results$wall_time),
  #     label = paste0("n = ", n),
  #     colour = dir
  #   ),
  #   inherit.aes = FALSE,
  #   hjust = 1
  # ) +
  ggplot2::labs(
    x = "Job Size\n(number of ratio)",
    y = "Wall Time\n(seconds)",
    #title = ""
  ) +
  ggplot2::scale_x_continuous(
    #transform = "log10",
    labels = scales::label_comma(),
    #guide = "axis_logticks",
    #limits = c(1, 100000),
    #limits = c(5000, 17000)
  ) +
  ggplot2::scale_y_continuous(
    #transform = "log10",
    labels = scales::label_comma(),
    #guide = "axis_logticks",
    limits = c(0, 1500),
    #limits = c(0.1, 1000),
    #limits = c(0.3, 1)
  ) +
  ggplot2::theme(
    legend.position = "bottom"
  )

pdf(
  here::here(opt[["out"]]),
  width = opt[["width"]],
  height = opt[["height"]]
)
#p
p2
dev.off()

results |>
  readr::write_tsv(
    here::here(stringr::str_replace(opt[["out"]], "\\.pdf$", ".tsv.gz"))
  )

# ggplot2::ggsave(
#   filename = here::here(opt[["out"]]),
#   plot     = p,
#   width    = opt[["width"]],
#   height   = opt[["height"]]
# )
# ggplot2::ggsave(
#   filename = here::here("walltime.pdf"),
#   plot     = p2,
#   width    = 8,
#   height   = 6
# )

# results |>
#   dplyr::group_by(
#     dir
#   ) |>
#   dplyr::summarise(
#     n = dplyr::n(),
#     total_wall_time_hrs = sum(wall_time) / 60 / 60,
#     min_cpu_time_hrs = min(wall_time * cores) / 60 / 60,
#     mean_cpu_time_hrs = mean(wall_time * cores) / 60 / 60,
#     max_cpu_time_hrs = max(wall_time * cores) / 60 / 60,
#     total_cpu_time_hrs = sum(wall_time * cores) / 60 / 60
#   )

# # Full Run: results/-3_-1_middle
# dir <- "output/aws_conda"
# #dir <- "output/aws_renv"
# #dir <- "results/-3_-1_middle"
# loess_fit <- loess(
#   wall_time ~ pair_count,
#   data = results |> dplyr::filter(dir == !!dir)
# )
# newdata <- tibble::tibble(pair_count = 1:max(results$pair_count))
# preds <- newdata |>
#   dplyr::mutate(
#     cores = 8,
#     pred_wall_time = predict(loess_fit, newdata)
#   )
# preds |>
#   dplyr::summarise(
#     n = dplyr::n(),
#     total_wall_time_hrs = sum(pred_wall_time) / 60 / 60,
#     min_cpu_time_hrs = min(pred_wall_time * cores) / 60 / 60,
#     mean_cpu_time_hrs = mean(pred_wall_time * cores) / 60 / 60,
#     max_cpu_time_hrs = max(pred_wall_time * cores) / 60 / 60,
#     total_cpu_time_hrs = sum(pred_wall_time * cores) / 60 / 60
#   ) |>
#   dplyr::mutate(
#     cost = total_wall_time_hrs * 0.308
#   )
