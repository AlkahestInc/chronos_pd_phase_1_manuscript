## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

#' Per-protein linear models versus one metadata predictor (with optional covariates)
#'
#' For each numeric protein column in `protein_intensity`, fit a linear regression
#' of the form \eqn{y \sim x + C}, where \eqn{y} is the protein abundance,
#' \eqn{x} is a single predictor from `metadata` (e.g., age), and \eqn{C} is an
#' optional set of covariates. Returns tidy coefficients (excluding the intercept)
#' along with model diagnostics and BH-adjusted p-values, stacked across proteins.
#'
#' @section Overview:
#' \enumerate{
#' \item Inner-join `protein_intensity` and `metadata` on `id_col`.
#' \item Optionally center/scale numeric predictors and covariates when `scale_predictor = TRUE`.
#' \item For each protein column, drop incomplete rows and fit \code{lm(y ~ value\_col + covariates)}.
#' \item Collect coefficient rows (excluding the intercept), add \eqn{R^2}, adjusted \eqn{R^2}, and \code{n}.
#' \item Adjust p-values across proteins per \code{term} using Benjamini–Hochberg (FDR).
#' }
#'
#' @param protein_intensity A data frame/tibble with the **first column** as sample ID
#'   (`id_col`) and remaining columns as **numeric** protein intensities (responses).
#' @param metadata A data frame/tibble containing `id_col`, the primary predictor
#'   `value_col`, and (optionally) additional columns named in `covariates`.
#' @param value_col Character. Name of the primary predictor column in `metadata`
#'   (e.g., `"age_at_donation"`).
#' @param id_col Character. Join key present in both tables. Defaults to the first
#'   column name of `protein_intensity`.
#' @param covariates Character vector or `NULL`. Optional additional predictors
#'   in `metadata` (e.g., `c("sex", "batch")`). Can include factors or numeric columns.
#' @param scale_predictor Logical. If `TRUE`, center/scale **numeric** predictors/covariates
#'   in `metadata` (excluding `id_col`). Factors are unchanged. Default `FALSE`.
#'
#' @return A tibble with one or more rows per protein (one per kept term), including:
#' \describe{
#'   \item{protein}{Protein column name (response).}
#'   \item{term}{Predictor/covariate name (no intercept).}
#'   \item{estimate}{Coefficient estimate for \code{term}.}
#'   \item{std.error}{Standard error of the estimate.}
#'   \item{statistic}{t-statistic for \code{term}.}
#'   \item{p.value}{Two-sided p-value for \code{term}.}
#'   \item{r.squared}{Model \eqn{R^2} from \code{lm}.}
#'   \item{adj.r.squared}{Adjusted \eqn{R^2} from \code{lm}.}
#'   \item{n}{Number of observations used in that protein's model.}
#'   \item{p.adj}{Benjamini–Hochberg adjusted p-value, computed within each \code{term}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic: protein ~ age
#' res <- regress_proteins_vs_metadata(
#'   protein_intensity = protein_intensity,
#'   metadata = metadata,
#'   value_col = "age_at_donation",
#'   id_col = "ID"
#' )
#'
#' # With covariates and scaling numeric predictors
#' res2 <- regress_proteins_vs_metadata(
#'   protein_intensity = protein_intensity,
#'   metadata = metadata,
#'   value_col = "age_at_donation",
#'   id_col = "ID",
#'   covariates = c("gender"),
#'   scale_predictor = TRUE
#' )
#' }
#' @export
#' @importFrom dplyr inner_join group_by mutate ungroup arrange relocate filter bind_rows
#' @importFrom purrr map_dfr
#' @importFrom broom tidy glance
#' @importFrom stats reformulate lm p.adjust
#' @importFrom tibble tibble
regress_proteins_vs_metadata_for_prot <- function(
  protein_intensity,
  metadata,
  value_col, # e.g. "age_at_donation"
  id_col = colnames(protein_intensity)[1],
  covariates = NULL, # e.g. c("sex","batch")
  scale_predictor = FALSE # center/scale the predictor (and any covariates if numeric)
) {
  # --- QC ---
  if (!id_col %in% names(protein_intensity)) stop("ID column not in protein_intensity: ", id_col)
  if (!id_col %in% names(metadata)) stop("ID column not in metadata: ", id_col)
  if (!value_col %in% names(metadata) & !value_col %in% names(protein_intensity)) stop("value_col not in metadata or protein intensity: ", value_col)
  if (!is.null(covariates) && !all(covariates %in% names(metadata))) {
    missing_cov <- setdiff(covariates, names(metadata))
    stop("Covariates missing in metadata: ", paste(missing_cov, collapse = ", "))
  }

  # Keep only needed metadata columns (ID + value_col + covariates)
  meta_keep <- c(id_col, covariates)
  md <- metadata[, meta_keep, drop = FALSE]

  # Optional scaling of numeric predictors/covariates
  if (scale_predictor) {
    num_cols <- names(md)[sapply(md, is.numeric)]
    num_cols <- setdiff(num_cols, id_col)
    if (length(num_cols)) {
      md[num_cols] <- lapply(md[num_cols], function(z) as.numeric(scale(z)))
    }
  }

  # Merge
  df <- dplyr::inner_join(protein_intensity, md, by = id_col)
  if (nrow(df) < min(nrow(protein_intensity), nrow(metadata))) {
    warning("Merge shrank rows — check for ID mismatches.")
  }

  prot_cols <- setdiff(names(protein_intensity), id_col)
  if (!length(prot_cols)) stop("No protein columns found after excluding ID column.")

  # Helper to fit one model per protein
  fit_one <- function(pcol) {
    dat <- df[, c(pcol, covariates), drop = FALSE]
    # names(dat)[1] <- "y"

    # Drop rows with NA in y or predictors
    dat <- dat[complete.cases(dat), ]
    if (nrow(dat) < 3) {
      return(tibble::tibble(
        protein = pcol, term = value_col,
        estimate = NA_real_, std.error = NA_real_,
        statistic = NA_real_, p.value = NA_real_,
        r.squared = NA_real_, adj.r.squared = NA_real_, n = nrow(dat)
      ))
    }

    # lm: y ~ value_col + covariates
    # f <- reformulate(termlabels = c(value_col, covariates), response = "y")
    f <- reformulate(covariates, response = pcol)
    fit <- lm(f, data = dat)

    tid <- broom::tidy(fit)
    gla <- broom::glance(fit)

    # Return only rows for predictors (exclude intercept)
    out <- tid |>
      dplyr::filter(term != "(Intercept)") |>
      dplyr::mutate(
        protein = pcol,
        r.squared = gla$r.squared,
        adj.r.squared = gla$adj.r.squared,
        n = gla$nobs
      ) |>
      dplyr::relocate(protein)

    # If the term of interest is missing (perfect collinearity etc.), fill NA row
    if (!value_col %in% out$term) {
      out <- dplyr::bind_rows(
        out,
        tibble::tibble(
          protein = pcol, term = value_col,
          estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_,
          r.squared = gla$r.squared, adj.r.squared = gla$adj.r.squared, n = gla$nobs
        )
      )
    }
    out
  }

  res <- purrr::map_dfr(prot_cols, fit_one)

  # For convenience, put the requested predictor row first if multiple terms exist
  res <- res |> dplyr::arrange(match(term, c(value_col, setdiff(unique(term), value_col))))

  res
}
