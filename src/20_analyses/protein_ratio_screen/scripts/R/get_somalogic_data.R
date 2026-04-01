## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript
# Usage: scripts/R/get_somalogic_data.R --output assets/olink.parquet

here::i_am("scripts/R/get_somalogic_data.R")

#####
# Define commandline arguments
#####
option_list <- list(
  optparse::make_option(c("-s", "--output"), type = "character", help = "Somalogic output data file (parquet format)"),
  optparse::make_option(c("-l", "--level"), type = "character", default = "INFO", help = "Log level: OFF, FATAL, ERROR, WARN, SUCCESS, INFO, DEBUG, TRACE  [default: INFO])")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

output_parquet_filename <- here::here(opt$output)
log_level <- toupper(opt$level)
#####

input_dir <- "_input_data"

data <- glue::glue(
  "{input_dir}/name={name}/version={version}/dataset.parquet",
  name = "Chronos_phase1_proteomics_somalogic_table",
  version = "none"
)  |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  dplyr::mutate(
    group_name = factor(group_name, levels = c("control", "pd")),
    pair_id = factor(pair_id),
    patient_id = factor(patient_id),
    gender = factor(gender),
    age_at_baseline = as.integer(age_at_baseline)
  )

# Ensure the parent dir of the output files is present
purrr::map(c(output_parquet_filename), ~ fs::dir_create(dirname(.x)))

data |>
  arrow::write_parquet(output_parquet_filename)
