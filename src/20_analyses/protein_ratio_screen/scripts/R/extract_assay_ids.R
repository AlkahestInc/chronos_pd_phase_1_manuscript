## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript

here::i_am("scripts/R/extract_assay_ids.R")

#####
# Define commandline arguments
#####
option_list <- list(
  optparse::make_option(
    c("-o", "--olink"),
    type = "character",
    help = "Olink input data file (parquet format)"
  ),
  optparse::make_option(
    c("-s", "--somalogic"),
    type = "character",
    help = "Somalogic input data file (parquet format)"
  ),
  optparse::make_option(
    c("-a", "--alamar"),
    type = "character",
    help = "Alamar Bio input data file (parquet format)"
  ),
  optparse::make_option(
    c("-f", "--file"),
    type = "character",
    help = "Filename for list of assay IDs (txt format)"
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
#####

input_olink_filename <- here::here(opt$olink)
input_somalogic_filename <- here::here(opt$somalogic)
input_alamar_filename <- here::here(opt$alamar)

#####
# Load Data
#####
olink_data <- arrow::read_parquet(
  input_olink_filename
)
somalogic_data <- arrow::read_parquet(
  input_somalogic_filename
)
alamar_data <- arrow::read_parquet(
  input_alamar_filename
)

assay_ids <- c(
  somalogic_data |>
    dplyr::select(
      years_since_PDonset:tidyselect::last_col(),
      -years_since_PDonset
    ) |>
    colnames(),

  olink_data |>
    dplyr::select(
      years_since_PDonset:tidyselect::last_col(),
      -years_since_PDonset
    ) |>
    colnames(),

  alamar_data |>
    dplyr::select(
      years_since_PDonset:tidyselect::last_col(),
      -years_since_PDonset
    ) |>
    colnames()
) |>
  sort(decreasing = TRUE) |>
  tibble::tibble() |>
  readr::write_tsv(
    here::here(opt$file),
    col_names = FALSE
  )
