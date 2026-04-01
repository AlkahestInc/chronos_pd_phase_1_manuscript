## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Stratified Matching for Cohort Data Using MatchIt
#'
#' Performs stratified matching between treatment and control groups
#' based on a set of predictor variables. The function constructs a
#' matching formula, runs `MatchIt::matchit()` using nearest-neighbor
#' Mahalanobis matching (or other user-specified method), and returns
#' both the `matchit` object and the matched dataset.
#'
#' @param data A data frame or tibble containing clinical or cohort-level data.
#' @param predictors Character vector of predictor (covariate) column names
#'   used for sample matching. Default: `c("age_at_visit", "sex", "ad")`.
#' @param response Character. The binary treatment/exposure variable to match on.
#'   Default: `"pd"`.
#' @param tolerance_in_age Numeric. Age caliper tolerance for matching (in years).
#'   Used to restrict matches within a certain age range. Default: `5`.
#' @param ratio Integer. Number of controls to match to each treatment sample.
#'   Default: `1`.
#' @param seed Integer. Random seed for reproducibility. Default: `42`.
#' @param method Character. Matching method passed to `MatchIt`. Default: `"nearest"`.
#' @param distance Character. Distance metric passed to `MatchIt`. Default: `"mahalanobis"`.
#' @param formula Formula or character. Matching formula. If supplied as a character,
#'   it will be parsed; however, this function will *rebuild* a formula using the
#'   `predictors` and `response` supplied. Default: `"pd ~ age_at_visit + sex + ad"`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `m.out` — the `matchit` object returned by `MatchIt::matchit()`
#'   \item `matched_metadata` — a tibble of the matched dataset (from `MatchIt::match.data()`)
#' }
#'
#' @examples
#' \dontrun{
#' matched <- match_stratified_cohort(
#'   data = clinical_data,
#'   predictors = c("age_at_visit", "sex", "ad"),
#'   response = "pd"
#' )
#' }
match_stratified_cohort <- function(
  data,
  predictors = c("age_at_visit", "sex", "ad"),
  response = "pd",
  tolerance_in_age = 5,
  ratio = 1,
  seed = 42,
  method = "nearest",
  distance = "mahalanobis",
  formula = "pd ~ age_at_visit + sex + ad"
) {
  # Validate that required variables exist
  required_vars <- all(all.vars(as.formula(formula)) %in% names(data))
  if (!required_vars) {
    stop("Some variables in the formula are not found in the data.")
  }

  # Construct formula from user-specified predictors + response
  formula <- reformulate(predictors, response = response)

  # Perform matching
  set.seed(seed)
  m.out <- MatchIt::matchit(
    formula,
    data = data,
    method = method,
    distance = distance,
    caliper = c(age_at_visit = tolerance_in_age),
    ratio = ratio
  )

  matched_metadata <- MatchIt::match.data(m.out)

  return(
    list(
      m.out = m.out,
      matched_metadata = matched_metadata
    )
  )
}
