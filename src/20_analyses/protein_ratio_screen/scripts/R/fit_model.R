## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#!/usr/bin/env Rscript
# Usage: scripts/R/fit_model.R --training assets/training.parquet --testing assets/testing.parquet --denominator seq.10000.28 --output results/denominator/seq.10000.28.tsv.gz

#####
logger::log_info("START")
#####

here::i_am("scripts/R/fit_model.R")

#####
# Define commandline arguments
#####
option_list <- list(
  optparse::make_option(c("-r", "--training"), type = "character", help = "Training input data file (parquet format)"),
  optparse::make_option(c("-e", "--testing"), type = "character", help = "Testing input data file (parquet format)"),
  optparse::make_option(c("-o", "--output"), type = "character", help = "Pairwise ratio assessment output file (parquet format)"),
  optparse::make_option(c("-d", "--denominator"), type = "character", help = "Proteomics assay identifier"),
  optparse::make_option(c("-n", "--number"), type = "integer", default = 99999999, help = "The maximum number of numerators to process [default: 99999999])"),
  optparse::make_option(c("-v", "--partitions"), type = "integer", default = 3, help = "Number of partitions in V-fold cross-validation [default: 3])"),
  optparse::make_option(c("-p", "--repeats"), type = "integer", default = 3, help = "Number of repeats in V-fold cross-validation [default: 3])"),
  optparse::make_option(c("-m", "--permute"), type = "logical", default = FALSE, help = "Enable permutation of group_name labels within donor pairs [default: FALSE])"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = parallelly::availableCores(), help = "Number of cores [default: parallelly::availableCores()])"),
  optparse::make_option(c("-l", "--level"), type = "character", default = "INFO", help = "Log level: OFF, FATAL, ERROR, WARN, SUCCESS, INFO, DEBUG, TRACE  [default: INFO])")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
#####

log_level <- toupper(opt$level)
cores <- opt$cores
number <- opt$number
partitions <- opt$partitions
repeats <- opt$repeats
denominator <- opt$denominator
output <- here::here(opt$output)
training <- here::here(opt$training)
permute <- opt$permute
testing <- here::here(opt$testing)


calculate_protein_ratios <- function(
  numerator,
  denominator,
  dataset,
  is_logged = TRUE
) {
  logger::log_trace("      Create pairwise combinations")
  comparisons <- expand.grid(
    top = numerator,
    bottom = denominator,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(
      top < bottom
    )

  logger::log_trace("      Calculating ratios")
  op <- if (is_logged) `-` else `/`

  compute_difference <- function(col_name) {
    ratio <- op(dataset[[col_name[["top"]]]], dataset[[col_name[["bottom"]]]])

    return(ratio)
  }
  col_names <- purrr::pmap(comparisons, c)
  res <- purrr::map(col_names, compute_difference)

  logger::log_trace("      Adding names to results list")
  names(res) <- vapply(
    col_names,
    FUN.VALUE = character(1),
    function(x) make.names(paste(x[1], x[2], sep = "_o_"))
  )

  return(res)
}

fit_and_assess.glm <- function(split, protein_pair_data, model_formula) {
  training_data <- rsample::training(split)
  assessment_data <- rsample::assessment(split)

  outcome <- rlang::f_lhs(model_formula)

  logistic.fit <- tryCatch(
    stats::glm(model_formula, training_data, family = binomial),
    error = function(e) {
      NULL
    }
  )
  if (!inherits(logistic.fit, "glm")) {
    result <- tibble::tibble(
      auc_train = as.double(NA_real_),
      auc_assessment = as.double(NA_real_)
    )

    return(result)
  }

  auc_train <- stats::predict(logistic.fit, type = "response") |>
    ROCR::prediction(training_data[[outcome]]) |>
    ROCR::performance(measure = "auc") |>
    slot("y.values") |>
    unlist()

  auc_assessment <- stats::predict(
    logistic.fit,
    newdata = assessment_data,
    type = "response"
  ) |>
    ROCR::prediction(assessment_data[[outcome]]) |>
    ROCR::performance(measure = "auc") |>
    slot("y.values") |>
    unlist()

  result <- tibble::tibble(
    auc_train = as.double(auc_train),
    auc_assessment = as.double(auc_assessment)
  )

  return(result)
}

fit_and_assess.fastglm <- function(split, protein_pair_data, model_formula) {
  training_data <- rsample::training(split)
  assessment_data <- rsample::assessment(split)

  outcome <- rlang::f_lhs(model_formula)

  #####
  # Extract x/y in the formats required for fastglm
  #####
  x_training <- model.matrix(model_formula, data = training_data)
  x_assessment <- model.matrix(model_formula, data = training_data)
  # convert 2-level outcome factor to 0/1
  y_training <- training_data[[outcome]]
  y_training <- as.integer(y_training == levels(y_training)[2]) # use second level as "1"

  logistic.fit <- tryCatch(
    fastglm::fastglm(
      x = x_training,
      y = y_training,
      family = binomial()
    ),
    error = function(e) {
      NULL
    }
  )
  if (!inherits(logistic.fit, "glm")) {
    result <- tibble::tibble(
      auc_train = as.double(NA_real_),
      auc_assessment = as.double(NA_real_)
    )

    return(result)
  }

  # Use model to predict outcome in training data and assessment data
  pred_training <- stats::predict(logistic.fit, x_training, type = "response")
  pred_assessment <- stats::predict(
    logistic.fit,
    x_assessment,
    type = "response"
  )

  # Return the AUC's for the training and assessment data
  result <- tibble::tibble(
    auc_train = as.double(pROC::auc(training_data[[outcome]], pred_training)),
    auc_assessment = as.double(pROC::auc(
      assessment_data[[outcome]],
      pred_assessment
    ))
  )

  return(result)
}

process_protein_pair <- function(protein_pair_id) {
  model_formula <- glue::glue(
    "{outcome} ~ {covariates} + {protein_pair_id}",
    covariates = paste(covariates, collapse = " + ")
  ) |>
    stats::as.formula()

  #####
  logger::log_trace(glue::glue(
    "      {id}: Define minimal training data input",
    id = protein_pair_id
  ))
  #####
  #####
  # TODO: try to avoid creating this object for each core as this should reduce memory requirements
  #####
  training_protein_pair_data <- tibble::tibble(
    #####
    # Get the sample metadata
    #####
    training_data |>
      dplyr::select(
        # Identifiers
        "ID" = "PARENT_BARCODE",
        "pair_id",
        "patient_id",

        # Model variables
        tidyselect::all_of(c(outcome, covariates))
      ),
    #####
    # Add the ratios column for the current protein_pair_id
    #####
    !!rlang::sym(protein_pair_id) := training_ratio_data[[protein_pair_id]]
  )

  #####
  # Perform k-fold cross-validation on the folds
  logger::log_trace(glue::glue(
    "      {id}: Perform V-fold cross-validation",
    id = protein_pair_id
  ))
  logger::log_trace(glue::glue(
    "        {id}: v = {v}",
    id = protein_pair_id,
    v = partitions
  ))
  logger::log_trace(glue::glue(
    "        {id}: repeats = {repeats}",
    id = protein_pair_id,
    repeats = repeats
  ))

  folds <- training_protein_pair_data |>
    rsample::group_vfold_cv(
      v = partitions,
      repeats = repeats,
      group = pair_id,
      # strata = gender
    )

  fold_matrix <- folds |>
    dplyr::bind_cols(
      purrr::map_dfr(
        folds[["splits"]],

        fit_and_assess.glm,
        training_protein_pair_data,
        model_formula
      )
    )

  #####
  logger::log_trace(glue::glue(
    "      {id}: Calculate CV performance metrics",
    id = protein_pair_id
  ))
  #####
  training_cv_results <- fold_matrix |>
    dplyr::summarise(
      n_folds = dplyr::n(),
      n_folds_ok = sum(!is.na(auc_train)),

      auc_train_mean = mean(auc_train, na.rm = TRUE),
      auc_train_sd = sd(auc_train, na.rm = TRUE),

      auc_assessment_mean = mean(auc_assessment, na.rm = TRUE),
      auc_assessment_sd = sd(auc_assessment, na.rm = TRUE)
    )

  logger::log_trace(glue::glue(
    "      {id}: Fit model to whole training set",
    id = protein_pair_id
  ))
  logistic.fit <- stats::glm(
    model_formula,
    training_protein_pair_data,
    family = binomial
  )

  #####
  logger::log_trace(glue::glue(
    "      {id}: Predict case/control in whole training data",
    id = protein_pair_id
  ))
  # We are actually using the full training dataset
  #####
  training_predictions <- tibble::tibble(
    .pred_pd = stats::predict(
      logistic.fit,
      newdata = training_protein_pair_data,
      type = "response"
    ),
    .pred_control = 1 - .pred_pd
  ) |>
    dplyr::bind_cols(
      training_protein_pair_data
    ) |>
    # NOTE: Here is where you can set up an equivocal zone
    dplyr::mutate(
      .pred_class = probably::make_two_class_pred(
        estimate = .pred_control,
        levels = levels(group_name),
        threshold = 0.5,
        buffer = 0.0
      )
    ) |>
    dplyr::relocate(
      .pred_class
    )

  #####
  logger::log_trace(glue::glue(
    "      {id}: Build return object",
    id = protein_pair_id
  ))
  #####
  results <- dplyr::bind_cols(
    tibble::tibble(
      protein_pair = protein_pair_id,
    ),

    training_cv_results,

    training_predictions |>
      class_metrics(
        truth = group_name,
        estimate = .pred_class,
        .pred_control
      ) |>
      dplyr::select(-".estimator") |>
      tidyr::pivot_wider(
        names_from = .metric,
        values_from = .estimate
      )
  ) |>
    dplyr::relocate(
      roc_auc,
      .before = accuracy
    )

  return(results)
}

# Setup logger
logger <- logger::layout_glue_generator(
  format = "{sprintf('%5s', level)} [{format(time, \"%Y-%m-%d %H:%M:%OS2\")}]: {msg}"
)
logger::log_layout(logger)
logger::log_threshold(log_level)

#######################
##### Code Proper #####
#######################

logger::log_info(glue::glue("Setting up for: {denom}", denom = denominator))
if (!is.null(number)) {
  logger::log_info(glue::glue("  Maximum pairs: {n}", n = number))
}
logger::log_info(glue::glue("  Available cores: {n}", n = cores))

future::plan(future::multisession, workers = cores)

logger::log_info("  Loading data files")
#####
logger::log_debug(glue::glue("    training: {file}", file = training))
#####
training_data <- arrow::read_parquet(
  training
)
#####
# Permute/shuffle the case/control labels within each pair if requested
#####
if (permute) {
  training_data <- training_data |>
    dplyr::group_by(pair_id) |>
    dplyr::mutate(group_name = sample(group_name)) |>
    dplyr::ungroup()
}

#####
# Define global modelling parameters
#####
outcome <- "group_name"
covariates <- c("age_at_baseline", "gender")

all_assay_ids <- training_data |>
  dplyr::select(
    years_since_PDonset:tidyselect::last_col(),
    -years_since_PDonset
  ) |>
  colnames() |>
  sort(decreasing = TRUE)

class_metrics <- yardstick::metric_set(
  yardstick::roc_auc,
  yardstick::accuracy,
  yardstick::sensitivity,
  yardstick::specificity,
  yardstick::precision
)

logger::log_info("  Calculating ratios")
#####
logger::log_debug("    training data")
#####
training_ratio_data <- calculate_protein_ratios(
  numerator = all_assay_ids,
  denominator = denominator,
  dataset = training_data,
  is_logged = TRUE
)

#####
logger::log_info("  Model: fitting and assessing")
logger::log_info(glue::glue("    Cores: {cores}", cores = cores))

#####
protein_pair_ids <- names(training_ratio_data) |>
  sort()
if (number >= 1 && number < length(protein_pair_ids)) {
  protein_pair_ids <- tail(protein_pair_ids, number)
}
logger::log_info(glue::glue(
  "    Number of pairs: {n}",
  n = length(protein_pair_ids)
))

results <- future.apply::future_lapply(
  protein_pair_ids,

  function(protein_pair_id) {
    logger::log_trace(glue::glue("    Pair: {protein_pair_id}"))
    result <- process_protein_pair(protein_pair_id)

    return(result)
  },
  future.seed = TRUE
) |>
  dplyr::bind_rows()

if (nrow(results) == 0) {
  # No results, setup empty tibble with headers
  results <- tibble::tibble(
    protein_pair = character(),

    auc_train_mean = double(),
    auc_train_sd = double(),

    auc_assessment_mean = double(),
    auc_assessment_sd = double(),

    roc_auc = double(),
    accuracy = double(),
    sensitivity = double(),
    specificity = double(),
    precision = double()
  )
}
#####
logger::log_info("  Writing result to disk")
logger::log_debug(glue::glue("    File: {f}", f = output))
#####
results |>
  arrow::write_parquet(output)

#####
logger::log_info("DONE")
#####

#####
logger::log_info("  Validating parquet file")
#####
is_valid_parquet <- function(path) {
  tryCatch(
    {
      ds <- arrow::open_dataset(path, format = "parquet")

      return(TRUE)
    },
    error = function(e) {
      logger::log_error(
        "Failed to open parquet file at '{path}': {msg}",
        path = path,
        msg = conditionMessage(e)
      )

      return(FALSE)
    }
  )
}

is_valid <- is_valid_parquet(output)

if (is_valid) {
  logger::log_info("    is valid: true")
} else {
  logger::log_info("    is valid: false")
}
#####
logger::log_info("DONE")
#####

if (!is_valid) {
  quit(save = "no", status = 1)
}
