## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript
# Usage: scripts/R/initial_split.R --olink assets/olink.parquet --somalogic assets/somalogic.parquet --alamar assets/alamar.parquet --start -2 --end 0 --training assets/training.parquet --testing assets/testing.parquet

here::i_am("scripts/R/initial_split.R")

#####
# Define commandline arguments
#####
option_list <- list(
  optparse::make_option(c("-o", "--olink"), type = "character", help = "Olink input data file (parquet format)"),
  optparse::make_option(c("-s", "--somalogic"), type = "character", help = "Somalogic input data file (parquet format)"),
  optparse::make_option(c("-a", "--alamar"), type = "character", help = "Alamar Bio input data file (parquet format)"),
  optparse::make_option(c("-p", "--training_proportion"), type = "numeric", default = 0.7, help = "Proportion of data for training"),
  optparse::make_option(c("-r", "--training"), type = "character", help = "Training data output file"),
  optparse::make_option(c("-e", "--testing"), type = "character", help = "Testing data output file"),
  optparse::make_option(c("-1", "--window_start"), type = "numeric", help = "Start of date window in years since disease onset"),
  optparse::make_option(c("-2", "--window_end"), type = "numeric", help = "End of date window in years since disease onset"),
  optparse::make_option(c("-t", "--strategy"), type = "character", default = "last", help = "Strategy for seleccting samples from donors [Default: last]"),
  optparse::make_option(c("-l", "--level"), type = "character", default = "INFO", help = "Log level: OFF, FATAL, ERROR, WARN, SUCCESS, INFO, DEBUG, TRACE  [default: INFO])")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

input_olink_filename <- here::here(opt$olink)
input_somalogic_filename <- here::here(opt$somalogic)
input_alamar_filename <- here::here(opt$alamar)
training_proportion <- opt$training_proportion
window_years <- c(opt$window_start, opt$window_end)
output_training_filename <- here::here(opt$training)
output_testing_filename <- here::here(opt$testing)
strategy <- opt$strategy
log_level <- toupper(opt$level)
#####

#####
# Validate input arguments
#####
if (training_proportion < 0 || training_proportion > 1) {
  stop("Error: --training_proportion must be between 0 and 1 (inclusive).")
}

# Setup logger
logger::layout_glue_generator(
  format = "{sprintf('%5s', level)} [{format(time, \"%Y-%m-%d %H:%M:%OS2\")}]: {msg}"
) |>
  logger::log_layout()
logger::log_threshold(log_level)

logger::log_info("Loading datasets")
logger::log_debug("  Loading Olink: {input_olink_filename}")
logger::log_debug("  Loading Somalogic: {input_somalogic_filename}")
logger::log_debug("  Loading Alamar: {input_alamar_filename}")

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

#####
logger::log_info("Joining datasets on PARENT_BARCODE")
# TODO: All three datasets contain metadata need to handle this before merging...preferably also before those datasets were created.
#####
# 3 - olink QC error
# 22 - olink QC warning
combined_data <- somalogic_data |>
  dplyr::select(
    ID.somalogic = ID,
    PARENT_BARCODE,
    years_since_PDonset:tidyselect::last_col(),
    -years_since_PDonset
  ) |>
  dplyr::full_join(
    olink_data |>
      dplyr::select(
        ID.olink = ID,
        PARENT_BARCODE,
        years_since_PDonset:tidyselect::last_col(),
        -years_since_PDonset
      ),
    by = "PARENT_BARCODE"
  ) |>
  dplyr::full_join(
    alamar_data |>
      dplyr::select(
        ID.alamar = ID,
        PARENT_BARCODE,
        years_since_PDonset:tidyselect::last_col(),
        -years_since_PDonset
      ),
    by = "PARENT_BARCODE"
  )

#####
logger::log_info("Excluding samples not measured across all datasets")
logger::log_info(
  "  Number of initial samples: {n_samples}",
  n_samples = combined_data |> dplyr::tally() |> dplyr::pull()
)
#####
# missing <- combined_data |>
#   dplyr::filter(
#     !vctrs::vec_detect_complete(combined_data)
#   )
complete <- combined_data |>
  dplyr::filter(
    vctrs::vec_detect_complete(combined_data)
  )
# combined_data |>
#   tidyr::drop_na()
#####
logger::log_info(
  "  Number of complete samples: {n_samples}",
  n_samples = complete |> dplyr::tally() |> dplyr::pull()
)
#####

#####
logger::log_info("Extracting sample metadata")
#####
sample_metadata <- somalogic_data |>
  dplyr::select(
    "PARENT_BARCODE",
    "pair_id",
    "patient_id",
    "group_name",
    "gender",
    "age_at_baseline",
    "years_since_PDonset"
  )
# a <- sample_metadata |>
#   dplyr::group_by(
#     pair_id,
#     patient_id,
#     group_name
#   ) |>
#   dplyr::summarise(
#     n_samples = dplyr::n_distinct(PARENT_BARCODE)
#   )
# b <- sample_metadata |>
#   dplyr::semi_join(
#     missing,
#     by = dplyr::join_by(PARENT_BARCODE)
#   ) |>
#   dplyr::group_by(
#     pair_id,
#     patient_id,
#     group_name
#   ) |>
#   dplyr::summarise(
#     n_dropped = dplyr::n_distinct(PARENT_BARCODE)
#   )
#
# a |>
#   dplyr::left_join(
#     b,
#     by = dplyr::join_by(
#       pair_id,
#       patient_id,
#       group_name
#     )
#   ) |>
#   dplyr::mutate(
#     n_dropped = dplyr::coalesce(n_dropped, 0),
#     n_remaining = n_samples - n_dropped,
#   ) |>
#   dplyr::arrange(
#     dplyr::desc(n_dropped)
#   )

logger::log_info("Filtering datasets")
logger::log_info("  Window start: {start}", start = window_years[1])
logger::log_info("  Window end: {end}", end = window_years[2])
logger::log_info("  Strategy: {strategy}")

data <- complete |>
  dplyr::left_join(
    sample_metadata,
    by = dplyr::join_by(
      PARENT_BARCODE
    )
  ) |>
  dplyr::relocate(
    ID.somalogic,
    ID.olink,
    ID.alamar,
    tidyselect::all_of(
      sample_metadata |> colnames()
    )
  ) |>
  #####
  # Consider only samples occurring in the defined window
  #####
  dplyr::filter(
    years_since_PDonset >= window_years[1],
    years_since_PDonset < window_years[2]
  ) |>
  #####
  # Remove those donor pairs that are now orphaned
  #####
  dplyr::group_by(
    pair_id
  ) |>
  dplyr::filter(
    dplyr::n_distinct(group_name) == 2
  ) |>
  dplyr::group_by(
    patient_id
  ) |>
  dplyr::arrange(
    years_since_PDonset
  ) |>
  dplyr::mutate(
    selected_sample = dplyr::case_when(
      strategy == "last" ~ dplyr::n(),
      strategy == "first" ~ 1,
      strategy == "random" ~ sample(dplyr::row_number(), 1),
      strategy == "middle" ~
        {
          mid_point <- (min(years_since_PDonset) + max(years_since_PDonset)) / 2
          which.min(abs(years_since_PDonset - mid_point))
        },
      TRUE ~ NA_integer_
    )
  ) |>
  dplyr::filter(selected_sample == dplyr::row_number()) |>
  dplyr::select(-selected_sample) |>
  dplyr::ungroup()

#####
# Deal with situation when no split is required
#####
logger::log_info("Creating initial split")
logger::log_info(
  "  Proportion (training/testing): {training_prop}/{testing_prop} ",
  training_prop = training_proportion,
  testing_prop = 1 - training_proportion
)
logger::log_info(
  "  Training file: {file}",
  file = output_training_filename
)
logger::log_info(
  "  Testing file: {file}",
  file = output_testing_filename
)
if (training_proportion == 1.0) {
  # All data to go into training file
  data |>
    arrow::write_parquet(output_training_filename)
  data |>
    dplyr::slice_head(n = 0) |>
    arrow::write_parquet(output_testing_filename)
  quit(save = "no", status = 0)
} else if (training_proportion == 0.0) {
  # All data to go into testing file
  data |>
    dplyr::slice_head(n = 0) |>
    arrow::write_parquet(output_training_filename)
  data |>
    arrow::write_parquet(output_testing_filename)
  quit(save = "no", status = 0)
}

#####
# Create an initial split of the data for training and testing
#####
data_split <- data |>
  rsample::group_initial_split(
    prop = training_proportion,
    group = pair_id
  )

#####
# Save the subsets to file
#####
logger::log_info("Saving initial split to file(s)")
rsample::training(data_split) |>
  arrow::write_parquet(output_training_filename)
rsample::testing(data_split) |>
  arrow::write_parquet(output_testing_filename)
