
# This function reproduces the estimates and standard errors returned by
# aggregate() when applied to a fixest object, i.e. the aggregate.fixest()
# method used for Sun-Abraham aggregations, while also returning the full
# transformed covariance matrix A V A'. The full covariance matrix is not
# directly exposed by aggregate(), so validation is based on:
#   1. matching aggregate() estimates,
#   2. matching aggregate() standard errors,
#   3. direct verification of pairwise covariance entries.
#
# Function base/motivation:
# This function builds on the approach discussed in fixest GitHub Issue #295:
# https://github.com/lrberge/fixest/issues/295

#' Extract aggregated Sun-Abraham coefficients and their covariance matrix
#'
#' @description
#' `sunab_aggregate_vcov()` aggregates non-aggregated Sun-Abraham coefficients
#' from a `fixest::feols()` model and returns both the aggregated coefficient
#' estimates and the full variance-covariance matrix of the aggregated estimates.
#'
#' The function is intended for models estimated with `fixest::sunab(...,
#' no_agg = TRUE)`, especially when the Sun-Abraham terms are interacted with
#' subgroup indicators and the user needs post-estimation inference on aggregated
#' effects. It reproduces the estimates and standard errors returned by
#' `aggregate()` for regex-defined aggregations, while also returning the full
#' transformed covariance matrix:
#'
#' `Sigma_agg = A %*% Sigma %*% t(A)`
#'
#' where `A` is the aggregation matrix and `Sigma` is the covariance matrix of
#' the non-aggregated model coefficients.
#'
#' This is useful for Wald tests, custom contrasts, subgroup comparisons,
#' event-window averages, cohort-bin comparisons, and other post-estimation tests
#' that require the covariance between aggregated coefficients, not just their
#' individual standard errors.
#'
#' @details
#' The aggregation is defined by a regular expression supplied to `agg`.
#' Coefficients whose names match `agg` are selected for aggregation. Capture
#' groups in the regex define the aggregation groups.
#'
#' A key rule is:
#'
#' - To preserve a dimension in the output, capture it with parentheses.
#' - To aggregate over a dimension, match it but do not capture it.
#'
#' For example, with coefficient names like:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' the regex:
#'
#' `"(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"`
#'
#' captures:
#'
#' - `group_1`: the event-time term, e.g. `"year::2"`
#' - `group_2`: the subgroup term, e.g. `"cd_bf"`
#'
#' and aggregates over cohorts within each event-time-by-subgroup cell.
#'
#' A regex such as:
#'
#' `"year::[2-6]:cohort::[0-9]+:(cd_.*)"`
#'
#' captures only the subgroup term and therefore aggregates over both cohorts
#' and event times 2 through 6, producing one estimate per subgroup.
#'
#' For wider event-time windows, non-capturing groups can be useful. For example:
#'
#' `"year::(?:[2-9]|1[0-5]):cohort::[0-9]+:(cd_.*)"`
#'
#' matches event times 2 through 15 but captures only the subgroup, so the output
#' is aggregated over event time and cohort, with one row per subgroup.
#'
#' @section Aggregation weights:
#' By default, `weight_method = "model_matrix"` reconstructs aggregation weights
#' using the same model-matrix logic used by `aggregate()` for `fixest` objects:
#'
#' `colSums(sign(model.matrix(...)))`
#'
#' or, for weighted models:
#'
#' `colSums(weights * sign(model.matrix(...)))`.
#'
#' The resulting coefficient-level weights are normalized within each aggregation
#' group. This means that the function produces model-matrix-weighted averages,
#' not simple equal-weighted averages across coefficient names.
#'
#' The `"model_matrix"` path is the safest reference path, but it can be
#' memory-intensive for large models because it may materialize an `N x K` model
#' matrix, where `N` is the number of observations and `K` is the number of
#' selected non-aggregated coefficients.
#'
#' If `weight_method = "data_count"`, the function avoids calling
#' `model.matrix()` and instead computes aggregation weights from `df_est` by
#' counting, or weighted-counting, observations in each coefficient cell. This is
#' much more memory-efficient, but it assumes coefficient names follow the usual
#' `sunab(..., no_agg = TRUE):dummy` pattern, e.g.:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' In this case, the function counts observations where:
#'
#' - `period_var - cohort_var == 2`;
#' - `cohort_var == 2005`;
#' - `cd_bf != 0`.
#'
#' If the interaction column is signed rather than strictly 0/1, the data-count
#' path uses `sign(interaction_value)`, matching the model-matrix weighting logic.
#'
#' @section Data-count weighting assumptions:
#' The `data_count` path is a fast path for dummy- or signed-dummy-interacted
#' Sun-Abraham designs where each selected coefficient corresponds to a cell
#' defined by:
#'
#' - event time, parsed from the coefficient name;
#' - cohort, parsed from the coefficient name;
#' - an interaction column in `df_est`, parsed from the coefficient name.
#'
#' It assumes that `df_est` is the exact estimation sample used by `feols()`.
#' If `feols()` dropped rows due to missingness, singleton fixed effects,
#' collinearity handling, weights, or other preprocessing, `df_est` must already
#' reflect those dropped rows.
#'
#' If the original `feols()` model was estimated with weights and
#' `weight_method = "data_count"`, `weight_var` should be supplied and should
#' point to the same weight variable used in `feols()`. Otherwise the aggregation
#' weights will be unweighted and may not match `aggregate()`.
#'
#' The `data_count` path may not reproduce the model-matrix weights for
#' continuous interactions, transformed variables, or specifications where the
#' selected coefficient columns are not simple cohort-by-event-time-by-interaction
#' cells.
#'
#' Always inspect `names(sunab_fixest$coefficients)` before using
#' `weight_method = "data_count"`. If the coefficient names differ from the
#' defaults, adjust `event_time_regex`, `cohort_regex`, and `interaction_regex`.
#'
#' The `data_count` path should be validated against `weight_method =
#' "model_matrix"` on a smaller dataset before production use.
#'
#' @section Post-hoc covariance matrices:
#' If `vcov_mat` is supplied, the function uses that covariance matrix instead
#' of `sunab_fixest$cov.scaled`. This allows the same aggregation matrix to be
#' applied to post-hoc covariance matrices, such as alternative clustered,
#' Conley, heteroskedastic, or user-supplied covariance matrices. The covariance
#' matrix must correspond to the same non-aggregated coefficient vector and must
#' have row and column names matching the model coefficient names.
#'
#' @section `group_fun`:
#' `group_fun` can be used to modify the captured grouping variables before
#' aggregation. This is useful for custom aggregations that cannot be expressed
#' with regex capture groups alone, such as recoding cohort years into bins
#' before aggregation.
#'
#' The function passed to `group_fun` receives a data frame containing the
#' captured groups plus a `term` column. It must return a data frame that
#' includes `term` and one or more grouping columns. Rows may be filtered to drop
#' terms from the aggregation, but `group_fun` may not add terms that were not
#' selected by the original `agg` regex.
#'
#' The grouping columns returned by `group_fun`, excluding `term`, define the
#' rows of the aggregated output. For example, returning columns `term`,
#' `event_time`, `cohort_bin`, and `cd` produces one output row per
#' event-time-by-cohort-bin-by-subgroup cell.
#'
#' @section Important:
#' This function should be used with models estimated using
#' `sunab(..., no_agg = TRUE)`. It is not intended for models where `coef()`
#' returns already-aggregated `sunab()` coefficients. In particular,
#' `sunab(..., no_agg = FALSE)` can return an aggregated coefficient view while
#' still retaining an underlying non-aggregated covariance structure, which is
#' not the target use case for this function.
#'
#' @section Dependencies:
#' This function uses `dplyr` internally. The `data_count` path additionally
#' uses `stringr`, `tibble`, and `tidyr`. These packages should be installed and
#' available.
#'
#' The examples use the native R pipe `|>` and additional `dplyr` verbs.
#'
#' @section Limitations:
#' This function reconstructs the model-matrix-weighted aggregation used by
#' `aggregate()` for `fixest` Sun-Abraham coefficients. It is not a general-purpose
#' coefficient-averaging function unless model-matrix weights are the desired
#' estimand.
#'
#' For ordinary TWFE event-study models, equal-weighted aggregation across event
#' periods may sometimes be more appropriate than model-matrix-weighted
#' aggregation.
#'
#' The `model_matrix` path assumes that `model.matrix(sunab_fixest)` can be
#' reconstructed. It may fail for lean model objects or objects where the model
#' matrix/data needed by `model.matrix()` have been removed.
#'
#' The supplied `vcov_mat`, if used, must correspond to the same non-aggregated
#' coefficient vector and must use coefficient names matching
#' `names(sunab_fixest$coefficients)`.
#'
#' `aggregate()` does not directly expose the full aggregated covariance matrix,
#' so the full `sigma` matrix returned here cannot be compared to a built-in
#' full aggregated covariance matrix. Validation should instead check that:
#'
#' - aggregated estimates match `aggregate()`;
#' - aggregated standard errors match `aggregate()`;
#' - pairwise covariance entries satisfy
#'   `sigma[i, j] = A[i, ] %*% V %*% A[j, ]`.
#'
#' @section Attribution:
#' This function builds on the approach discussed in `fixest` GitHub Issue #295,
#' "Post-estimation linear hypothesis testing using sunab()", which raised the
#' need to obtain the variance-covariance matrix for aggregated `sunab()`
#' coefficients in order to conduct post-estimation linear hypothesis tests:
#' <https://github.com/lrberge/fixest/issues/295>.
#'
#' This implementation generalizes that idea by allowing arbitrary regex-defined
#' aggregation groups, optional post-hoc covariance matrices, optional
#' user-defined recoding of aggregation groups through `group_fun`, and an
#' optional data-count weighting path to avoid materializing the model matrix.
#'
#' @param sunab_fixest A `fixest` model object, typically returned by
#'   `fixest::feols()`, estimated with one or more `fixest::sunab(...,
#'   no_agg = TRUE)` terms.
#'
#' @param agg Character string. A Perl-compatible regular expression used to
#'   select and group non-aggregated Sun-Abraham coefficient names. Coefficients
#'   matching `agg` are selected. Capture groups in `agg` define the aggregation
#'   groups. At least one capture group is required.
#'
#' @param vcov_mat Optional numeric covariance matrix. If `NULL`, the function
#'   uses `sunab_fixest$cov.scaled`. If supplied, `vcov_mat` should be a
#'   covariance matrix for the non-aggregated model coefficients, with row and
#'   column names matching `names(sunab_fixest$coefficients)`.
#'
#' @param group_fun Optional function used to modify or recode the captured
#'   grouping variables before aggregation. The function receives a data frame
#'   with columns `group_1`, `group_2`, ..., and `term`. It must return a data
#'   frame containing `term` and at least one grouping column. Rows may be
#'   filtered to drop terms from the aggregation, but `group_fun` may not add
#'   terms that were not selected by the original `agg` regex.
#'
#' @param weight_method Character. Either `"model_matrix"` or `"data_count"`.
#'   `"model_matrix"` is the default and reproduces `aggregate()` weighting
#'   using `model.matrix()`. `"data_count"` computes weights directly from
#'   `df_est` and avoids materializing the model matrix.
#'
#' @param df_est Data frame used when `weight_method = "data_count"`. This
#'   should be the exact estimation sample used by the `fixest` model.
#'
#' @param cohort_var Character. Name of the cohort/treatment-timing variable in
#'   `df_est`, e.g. `"FirstTreat"`. Required for `weight_method = "data_count"`.
#'
#' @param period_var Character. Name of the time-period variable in `df_est`,
#'   e.g. `"year"`. Required for `weight_method = "data_count"`.
#'
#' @param weight_var Optional character. Name of a column in `df_est` containing
#'   estimation weights. If `NULL`, unweighted counts are used in the
#'   `data_count` path. Use this only if the model was estimated with the same
#'   weights.
#'
#' @param event_time_regex Character regex used by the `data_count` path to
#'   parse event time from coefficient names. The first capture group must be
#'   the event-time value. The default assumes coefficient names contain strings
#'   like `"year::2"`.
#'
#' @param cohort_regex Character regex used by the `data_count` path to parse
#'   cohort from coefficient names. The first capture group must be the cohort
#'   value. The default assumes coefficient names contain strings like
#'   `"cohort::2005"`.
#'
#' @param interaction_regex Character regex used by the `data_count` path to
#'   parse the interacted variable from coefficient names. The first capture
#'   group must be the interaction variable name, e.g. `"cd_bf"`. The default
#'   assumes the interaction variable is the final colon-delimited component of
#'   the coefficient name.
#'
#' @returns
#' A list with the following elements:
#'
#' \describe{
#'   \item{beta}{A matrix of aggregated coefficient estimates. Rows correspond
#'   to aggregation groups and the single column is named `"estimate"`.}
#'
#'   \item{sigma}{The full variance-covariance matrix of the aggregated
#'   coefficient estimates, computed as `A %*% V %*% t(A)`.}
#'
#'   \item{transform}{The aggregation matrix `A`. Multiplying `A` by the
#'   non-aggregated coefficient vector produces the aggregated coefficients.}
#'
#'   \item{groups}{A data frame describing the aggregation groups, including
#'   the grouping variables, a `key` column, the aggregated estimate, and its
#'   standard error.}
#'
#'   \item{coef_names}{The original non-aggregated coefficient names used in
#'   the aggregation.}
#'
#'   \item{parsed_terms}{A data frame mapping each selected coefficient term to
#'   its parsed and, if applicable, recoded aggregation group.}
#'
#'   \item{coef_weights}{The unnormalized coefficient-level weights used to
#'   build the rows of `A`. These are useful for debugging and for validating
#'   whether `model_matrix` and `data_count` produce the same aggregation
#'   weights.}
#' }
#'
#' @examples
#' \dontrun{
#' # Disaggregated Sun-Abraham model with subgroup-specific effects
#' est_sunab_dummy <- fixest::feols(
#'   rap_tree ~
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
#'     pt_id + year,
#'   data = test_dats_sn2_small,
#'   cluster = ~ pt_id
#' )
#'
#' # Aggregate over cohorts within event-time-by-subgroup cells.
#' agg_cd_event <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Aggregate over cohorts and event times 2 through 6, separately by subgroup.
#' # The event-time window is matched but not captured.
#' att_2_6_cd <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Use a post-hoc covariance matrix.
#' V_alt <- vcov(
#'   est_sunab_dummy,
#'   vcov = ~ pt_id + year
#' )
#'
#' att_2_6_cd_alt <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)",
#'   vcov_mat = V_alt
#' )
#'
#' # Use the faster data-count path. `df_est` should be the exact estimation
#' # sample used by feols().
#' agg_cd_event_data <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)",
#'   weight_method = "data_count",
#'   df_est = test_dats_sn2_small,
#'   cohort_var = "FirstTreat",
#'   period_var = "year"
#' )
#'
#' # Validate the data-count path against the model-matrix path on a smaller
#' # dataset before using it in production.
#' all.equal(
#'   agg_cd_event$transform,
#'   agg_cd_event_data$transform,
#'   tolerance = 1e-12
#' )
#'
#' all.equal(
#'   agg_cd_event$beta,
#'   agg_cd_event_data$beta,
#'   tolerance = 1e-10
#' )
#'
#' all.equal(
#'   agg_cd_event$sigma,
#'   agg_cd_event_data$sigma,
#'   tolerance = 1e-10
#' )
#'
#' # Example: compare cd_bf versus cd_f for the event-time 2--6 average.
#' g <- att_2_6_cd$groups
#'
#' i_bf <- which(g$key == "cd_bf")
#' i_f  <- which(g$key == "cd_f")
#'
#' L <- matrix(0, nrow = 1, ncol = nrow(g))
#' L[1, i_bf] <- 1
#' L[1, i_f]  <- -1
#'
#' est <- as.numeric(L %*% att_2_6_cd$beta)
#' se <- sqrt(as.numeric(L %*% att_2_6_cd$sigma %*% t(L)))
#' z <- est / se
#' p <- 2 * pnorm(abs(z), lower.tail = FALSE)
#'
#' tibble::tibble(
#'   contrast = "cd_bf - cd_f, event times 2--6",
#'   estimate = est,
#'   se = se,
#'   z = z,
#'   p = p
#' )
#'
#' # Example: aggregate by event time, cohort bin, and subgroup.
#' es_by_cohort_bin <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::([0-9]+):(cd_.*)",
#'   group_fun = function(x) {
#'     x |>
#'       dplyr::mutate(
#'         event_time = group_1,
#'         cohort = as.integer(group_2),
#'         cd = group_3,
#'         cohort_bin = dplyr::case_when(
#'           cohort >= 2000 & cohort <= 2010 ~ "cohort_2000_2010",
#'           cohort >= 2011 & cohort <= 2020 ~ "cohort_2011_2020",
#'           TRUE ~ NA_character_
#'         )
#'       ) |>
#'       dplyr::filter(!is.na(cohort_bin)) |>
#'       dplyr::select(term, event_time, cohort_bin, cd)
#'   }
#' )
#'
#' # Example: aggregate event times 2 through 6 by cohort bin and subgroup.
#' # The event-time window is matched but not captured, so the output rows are
#' # cohort_bin x subgroup rather than event_time x cohort_bin x subgroup.
#' att_2_6_by_cohort_bin <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::([0-9]+):(cd_.*)",
#'   group_fun = function(x) {
#'     x |>
#'       dplyr::mutate(
#'         cohort = as.integer(group_1),
#'         cd = group_2,
#'         cohort_bin = dplyr::case_when(
#'           cohort >= 2000 & cohort <= 2010 ~ "cohort_2000_2010",
#'           cohort >= 2011 & cohort <= 2020 ~ "cohort_2011_2020",
#'           TRUE ~ NA_character_
#'         )
#'       ) |>
#'       dplyr::filter(!is.na(cohort_bin)) |>
#'       dplyr::select(term, cohort_bin, cd)
#'   }
#' )
#'
#' # Validate estimates and SEs against aggregate().
#' agg_check <- aggregate(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"
#' )
#'
#' all.equal(
#'   as.numeric(agg_check[, "Estimate"]),
#'   as.numeric(agg_cd_event$beta),
#'   tolerance = 1e-8
#' )
#'
#' all.equal(
#'   as.numeric(agg_check[, "Std. Error"]),
#'   sqrt(diag(agg_cd_event$sigma)),
#'   tolerance = 1e-8
#' )
#'
#' # Validate selected pairwise covariance entries.
#' A <- agg_cd_event$transform
#' V <- est_sunab_dummy$cov.scaled[colnames(A), colnames(A), drop = FALSE]
#'
#' i <- 1
#' j <- 2
#'
#' direct_cov_ij <- as.numeric(
#'   A[i, , drop = FALSE] %*% V %*% t(A[j, , drop = FALSE])
#' )
#'
#' stored_cov_ij <- agg_cd_event$sigma[i, j]
#'
#' all.equal(direct_cov_ij, stored_cov_ij, tolerance = 1e-8)
#' }
#'
#' @seealso
#' [fixest::sunab()], [aggregate()], [fixest::vcov.fixest()]
#'
#' @export
sunab_aggregate_vcov <- function(
    sunab_fixest,
    agg,
    vcov_mat = NULL,
    group_fun = NULL,
    weight_method = c("model_matrix", "data_count"),
    df_est = NULL,
    cohort_var = NULL,
    period_var = NULL,
    weight_var = NULL,
    event_time_regex = "year::(-?[0-9]+)",
    cohort_regex = "cohort::([0-9]+)",
    interaction_regex = ":([^:]+)$"
) {
  
  weight_method <- match.arg(weight_method)
  
  coef_vec <- sunab_fixest$coefficients
  coef_names <- names(coef_vec)
  
  if (is.null(coef_names)) {
    stop("The model coefficient vector must have names.")
  }
  
  is_match <- grepl(agg, coef_names, perl = TRUE)
  agg_coef_names <- coef_names[is_match]
  
  if (length(agg_coef_names) == 0) {
    stop("No coefficients matched the supplied `agg` regex.")
  }
  
  matches <- regexec(agg, agg_coef_names, perl = TRUE)
  captures <- regmatches(agg_coef_names, matches)
  
  n_captures <- length(captures[[1]]) - 1
  
  if (n_captures < 1) {
    stop("`agg` must contain at least one capture group.")
  }
  
  capture_mat <- do.call(
    rbind,
    lapply(captures, function(x) x[-1])
  )
  
  groups_raw <- as.data.frame(capture_mat, stringsAsFactors = FALSE)
  names(groups_raw) <- paste0("group_", seq_len(n_captures))
  groups_raw$term <- agg_coef_names
  
  if (!is.null(group_fun)) {
    groups <- group_fun(groups_raw)
  } else {
    groups <- groups_raw
  }
  
  if (!"term" %in% names(groups)) {
    stop("`group_fun` must return a data frame that includes the `term` column.")
  }
  
  groups <- groups |>
    dplyr::filter(!is.na(term))
  
  if (nrow(groups) == 0) {
    stop("No terms remain after applying `group_fun`.")
  }
  
  if (anyDuplicated(groups$term)) {
    stop("`group_fun` returned duplicated coefficient terms.")
  }
  
  invalid_terms <- setdiff(groups$term, groups_raw$term)
  
  if (length(invalid_terms) > 0) {
    stop(
      "`group_fun` returned terms that were not selected by the original `agg` regex: ",
      paste(utils::head(invalid_terms, 5), collapse = ", "),
      if (length(invalid_terms) > 5) " ..."
    )
  }
  
  missing_terms <- setdiff(groups$term, coef_names)
  
  if (length(missing_terms) > 0) {
    stop(
      "`group_fun` returned terms that are not model coefficients: ",
      paste(utils::head(missing_terms, 5), collapse = ", "),
      if (length(missing_terms) > 5) " ..."
    )
  }
  
  agg_coef_names <- groups$term
  
  group_cols <- setdiff(names(groups), "term")
  
  if (length(group_cols) == 0) {
    stop("No grouping columns remain after applying `group_fun`.")
  }
  
  group_key <- apply(
    groups[, group_cols, drop = FALSE],
    1,
    paste,
    collapse = "::"
  )
  
  unique_groups <- unique(groups[, group_cols, drop = FALSE])
  unique_groups$key <- apply(
    unique_groups,
    1,
    paste,
    collapse = "::"
  )
  
  if (weight_method == "model_matrix") {
    
    mm <- model.matrix(sunab_fixest)[, agg_coef_names, drop = FALSE]
    
    if (!is.null(sunab_fixest$weights)) {
      coef_wgt <- colSums(sunab_fixest$weights * sign(mm))
    } else {
      coef_wgt <- colSums(sign(mm))
    }
    
  } else if (weight_method == "data_count") {
    
    coef_wgt <- compute_sunab_coef_weights_from_data(
      coef_names = agg_coef_names,
      df_est = df_est,
      cohort_var = cohort_var,
      period_var = period_var,
      weight_var = weight_var,
      event_time_regex = event_time_regex,
      cohort_regex = cohort_regex,
      interaction_regex = interaction_regex
    )
  }
  
  A <- matrix(
    0,
    nrow = nrow(unique_groups),
    ncol = length(coef_vec),
    dimnames = list(unique_groups$key, coef_names)
  )
  
  for (i in seq_len(nrow(unique_groups))) {
    idx_local <- which(group_key == unique_groups$key[i])
    idx_global <- match(agg_coef_names[idx_local], coef_names)
    
    w <- coef_wgt[idx_local]
    
    if (sum(w) == 0) {
      stop("Aggregation weights sum to zero for group: ", unique_groups$key[i])
    }
    
    A[i, idx_global] <- w / sum(w)
  }
  
  if (is.null(vcov_mat)) {
    vcov_mat <- sunab_fixest$cov.scaled
  }
  
  if (is.null(vcov_mat)) {
    stop("No covariance matrix found. Provide `vcov_mat` or use a model with `cov.scaled`.")
  }
  
  if (is.null(rownames(vcov_mat)) || is.null(colnames(vcov_mat))) {
    stop("`vcov_mat` must have rownames and colnames matching coefficient names.")
  }
  
  missing_vcov_names <- setdiff(coef_names, rownames(vcov_mat))
  
  if (length(missing_vcov_names) > 0) {
    stop(
      "`vcov_mat` is missing model coefficients: ",
      paste(utils::head(missing_vcov_names, 5), collapse = ", "),
      if (length(missing_vcov_names) > 5) " ..."
    )
  }
  
  V <- vcov_mat[coef_names, coef_names, drop = FALSE]
  
  beta <- A %*% cbind(coef_vec)
  sigma <- A %*% V %*% t(A)
  
  colnames(beta) <- "estimate"
  
  unique_groups$estimate <- as.numeric(beta)
  unique_groups$se <- sqrt(diag(sigma))
  
  
  se <- sqrt(diag(sigma))
  t_value <- as.numeric(beta) / se
  
  # Try to mimic the p-value behavior of fixest aggregate().
  # For most clustered / robust fixest outputs, this is a t reference distribution.
  # If df extraction fails, fall back to normal approximation.
  t_df <- tryCatch(
    {
      if ("fixest" %in% loadedNamespaces() &&
          exists("degrees_freedom", where = asNamespace("fixest"), inherits = FALSE)) {
        fixest::degrees_freedom(sunab_fixest, type = "t")
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(t_df) || length(t_df) != 1 || is.na(t_df) || is.infinite(t_df)) {
    p_value <- 2 * stats::pnorm(abs(t_value), lower.tail = FALSE)
  } else {
    p_value <- 2 * stats::pt(abs(t_value), df = t_df, lower.tail = FALSE)
  }
  
  feols_structure <- cbind(
    Estimate = as.numeric(beta),
    `Std. Error` = se,
    `t value` = t_value,
    `Pr(>|t|)` = p_value
  )
  
  # Match aggregate.fixest()-style row names:
  # rownames are the aggregated group keys, with names corresponding to
  # the position of the first contributing coefficient in the original coef vector.
  feols_row_names <- unique_groups$key
  names(feols_row_names) <- vapply(
    unique_groups$key,
    function(k) {
      first_term <- agg_coef_names[which(group_key == k)[1]]
      as.character(match(first_term, coef_names))
    },
    character(1)
  )
  
  rownames(feols_structure) <- feols_row_names
  
  unique_groups$estimate <- as.numeric(beta)
  unique_groups$se <- se
  
  list(
    beta = beta,
    sigma = sigma,
    transform = A,
    groups = unique_groups,
    coef_names = agg_coef_names,
    parsed_terms = groups,
    coef_weights = coef_wgt,
    feols_structure = feols_structure
  )
  
  # list(
  #   beta = beta,
  #   sigma = sigma,
  #   transform = A,
  #   groups = unique_groups,
  #   coef_names = agg_coef_names,
  #   parsed_terms = groups,
  #   coef_weights = coef_wgt
  # )
}


#' Compute Sun-Abraham coefficient aggregation weights from estimation data
#'
#' @description
#' Internal helper used by `sunab_aggregate_vcov()` when
#' `weight_method = "data_count"`. It parses event time, cohort, and interaction
#' variable names from non-aggregated Sun-Abraham coefficient names, then
#' computes coefficient-level weights directly from the estimation data.
#'
#' @details
#' This helper is intended for coefficient names like:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' and data columns like:
#'
#' - `FirstTreat`, the cohort/treatment-timing variable;
#' - `year`, the period variable;
#' - `cd_bf`, the interacted dummy/signed-dummy variable.
#'
#' The returned weights are designed to match:
#'
#' `colSums(sign(model.matrix(model)[, coef_names]))`
#'
#' or, when `weight_var` is supplied:
#'
#' `colSums(weights * sign(model.matrix(model)[, coef_names]))`.
#'
#' This helper should generally not be called directly by users.
#'
#' @param coef_names Character vector of coefficient names.
#' @param df_est Data frame containing the exact estimation sample.
#' @param cohort_var Character name of the cohort/treatment-timing variable.
#' @param period_var Character name of the time-period variable.
#' @param weight_var Optional character name of the model weight variable.
#' @param event_time_regex Regex with first capture group identifying event time.
#' @param cohort_regex Regex with first capture group identifying cohort.
#' @param interaction_regex Regex with first capture group identifying the
#'   interaction variable.
#'
#' @returns
#' A named numeric vector of unnormalized coefficient-level aggregation weights.
compute_sunab_coef_weights_from_data <- function(
    coef_names,
    df_est,
    cohort_var,
    period_var,
    weight_var = NULL,
    event_time_regex = "year::(-?[0-9]+)",
    cohort_regex = "cohort::([0-9]+)",
    interaction_regex = ":([^:]+)$"
) {
  
  if (is.null(df_est)) {
    stop("`df_est` must be supplied when `weight_method = 'data_count'`.")
  }
  
  if (is.null(cohort_var) || is.null(period_var)) {
    stop(
      "`cohort_var` and `period_var` must be supplied when ",
      "`weight_method = 'data_count'`."
    )
  }
  
  required_vars <- c(cohort_var, period_var)
  
  if (!is.null(weight_var)) {
    required_vars <- c(required_vars, weight_var)
  }
  
  missing_required <- setdiff(required_vars, names(df_est))
  
  if (length(missing_required) > 0) {
    stop(
      "`df_est` is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  event_time <- stringr::str_match(coef_names, event_time_regex)[, 2]
  cohort <- stringr::str_match(coef_names, cohort_regex)[, 2]
  interaction_var <- stringr::str_match(coef_names, interaction_regex)[, 2]
  
  parsed <- tibble::tibble(
    term = coef_names,
    event_time = as.integer(event_time),
    cohort = as.integer(cohort),
    interaction_var = interaction_var
  )
  
  if (anyNA(parsed$event_time)) {
    bad <- parsed$term[is.na(parsed$event_time)]
    stop(
      "Could not parse event time from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  if (anyNA(parsed$cohort)) {
    bad <- parsed$term[is.na(parsed$cohort)]
    stop(
      "Could not parse cohort from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  if (anyNA(parsed$interaction_var)) {
    bad <- parsed$term[is.na(parsed$interaction_var)]
    stop(
      "Could not parse interaction variable from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  interaction_vars <- unique(parsed$interaction_var)
  missing_interactions <- setdiff(interaction_vars, names(df_est))
  
  if (length(missing_interactions) > 0) {
    stop(
      "`df_est` is missing interaction/dummy columns parsed from coefficient names: ",
      paste(utils::head(missing_interactions, 10), collapse = ", "),
      if (length(missing_interactions) > 10) " ..."
    )
  }
  
  df_counts <- df_est |>
    dplyr::mutate(
      .event_time = .data[[period_var]] - .data[[cohort_var]],
      .cohort = .data[[cohort_var]],
      .row_weight = if (is.null(weight_var)) 1 else .data[[weight_var]]
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(interaction_vars),
      names_to = "interaction_var",
      values_to = ".interaction_value"
    ) |>
    dplyr::filter(.interaction_value != 0, !is.na(.interaction_value)) |>
    dplyr::group_by(.event_time, .cohort, interaction_var) |>
    dplyr::summarise(
      coef_wgt = sum(.row_weight * sign(.interaction_value), na.rm = TRUE),
      .groups = "drop"
    )
  
  weights <- parsed |>
    dplyr::left_join(
      df_counts,
      by = c(
        "event_time" = ".event_time",
        "cohort" = ".cohort",
        "interaction_var" = "interaction_var"
      )
    ) |>
    dplyr::mutate(
      coef_wgt = tidyr::replace_na(coef_wgt, 0)
    ) |>
    dplyr::pull(coef_wgt)
  
  stats::setNames(weights, coef_names)
}