# weight_dydid_pipeline_v7.R
# Portable Sun-Abraham / TWFE DiD estimation pipeline — v7
# -------------------------------------------------------
# Key changes from v6:
#   - Unified dummy model: all treatment groups in one feols() call via
#     sunab():cd_group interactions
#   - no_agg = TRUE: cohort-specific coefficients as primary model output
#   - sunab_aggregate_vcov() replaces manual A-matrix aggregation; uses the
#     data_count weight path to avoid materialising the full model matrix
#   - agg_specs: list of aggregation specs (regex + group_fun), each run per
#     (run x vcov_spec); produces agg_{id}.parquet + agg_{id}__{vcov_id}.rds
#   - model is discarded (not saved) after vcov extraction and aggregation
#   - Registry: one row per dummy_group
#
# External dependency: sunab_aggregate_vcov.R must be sourced before this file.
#
# Output files per run (tables/by_run/{run_stub}/):
#   coef.parquet                raw disaggregated coeftable (all vcov_ids)
#   agg_{agg_id}.parquet        aggregated estimates (all vcov_ids; Script 4)
#   agg_{agg_id}__{vcov_id}.rds list(coef,vcov,groups) (Script 3 inference)
#   dummy_group.parquet         N stats per dummy_group
#   support.parquet             event-time support per dummy_group
#   registry.parquet            one row per dummy_group
#   run_spec.rds                full run spec (completion sentinel)
# -------------------------------------------------------


# ── Required packages ────────────────────────────────────────────────────────

.required_pkgs <- c(
  "arrow", "dplyr", "fixest", "glue", "purrr",
  "readr", "rlang", "stringr", "tibble", "tidyr", "jsonlite",
  "WeightIt"
)
.missing_pkgs <- .required_pkgs[
  !vapply(.required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(.missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(.missing_pkgs, collapse = ", "), call. = FALSE)
}
rm(.required_pkgs, .missing_pkgs)

# SECTION 1 — Dataset spec
# ══════════════════════════════════════════════════════════════════════════════

#' Build a dataset spec
#'
#' All column-name mappings that are properties of the data, not the model.
#' One spec per experiment call.
make_dataset_spec <- function(unit_id,
                              time_var,
                              trt_col,
                              cohort_var = NA_character_,
                              event_id   = NA_character_) {
  stopifnot(
    is.character(unit_id)  && length(unit_id)  == 1,
    is.character(time_var) && length(time_var) == 1,
    is.character(trt_col)  && length(trt_col)  == 1
  )
  list(
    unit_id    = unit_id,
    time_var   = time_var,
    trt_col    = trt_col,
    cohort_var = normalize_optional_colname(cohort_var),
    event_id   = normalize_optional_colname(event_id)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Analysis subset specs
# ══════════════════════════════════════════════════════════════════════════════

make_analysis_subset_spec <- function(subset_id,
                                      long_data_source,
                                      data_filter       = NULL,
                                      short_data_source = NULL) {
  if (!is.character(subset_id) || length(subset_id) != 1 || is.na(subset_id)) {
    stop("subset_id must be a length-1 non-NA character string.")
  }
  validate_data_source(long_data_source,
                       label = glue::glue("subset '{subset_id}' long_data_source"))
  if (!is.null(short_data_source)) {
    validate_data_source(short_data_source,
                         label = glue::glue("subset '{subset_id}' short_data_source"))
  }
  if (!is.null(data_filter) && !inherits(data_filter, "formula")) {
    stop("data_filter must be NULL or a one-sided formula, e.g. ~ year >= 1997")
  }
  tibble::tibble(
    subset_id         = subset_id,
    long_data_source  = list(long_data_source),
    short_data_source = list(short_data_source),
    data_filter       = list(data_filter)
  )
}


expand_analysis_subset_specs_by_col <- function(long_data_source,
                                                split_col,
                                                id_prefix         = split_col,
                                                base_filter       = NULL,
                                                values            = NULL,
                                                short_data_source = NULL,
                                                check_all_files   = FALSE) {

  if (isTRUE(check_all_files)) {
    parquet_files <- if (length(long_data_source) == 1 && dir.exists(long_data_source)) {
      list.files(long_data_source, pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)
    } else {
      long_data_source
    }
    if (length(parquet_files) == 0) stop("No parquet files found in long_data_source.")
    unique_vals <- purrr::map(parquet_files, \(f) {
      arrow::read_parquet(f, col_select = dplyr::all_of(split_col)) |>
        dplyr::pull(dplyr::all_of(split_col)) |>
        as.character() |>
        unique()
    }) |>
      unlist(use.names = FALSE) |>
      stats::na.omit() |>
      unique() |>
      sort()
  } else {
    unique_vals <- open_arrow_source(long_data_source) |>
      dplyr::select(dplyr::all_of(split_col)) |>
      dplyr::distinct() |>
      dplyr::collect() |>
      dplyr::pull(1) |>
      stats::na.omit() |>
      sort()
  }

  if (!is.null(values)) {
    unique_vals <- intersect(unique_vals, values)
    if (length(unique_vals) == 0) {
      stop("No overlap between provided values and unique values of '", split_col, "'.")
    }
  }

  purrr::map_dfr(unique_vals, function(val) {
    make_analysis_subset_spec(
      subset_id         = as.character(glue::glue("{id_prefix}_{safe_path_component(val)}")),
      long_data_source  = long_data_source,
      data_filter       = combine_filters(base_filter, build_equality_filter(split_col, val)),
      short_data_source = short_data_source
    )
  })
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Treatment group specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build one treatment group spec row
#'
#' group_fun contract (v7):
#'   Receives: df, group_col (injected), plus everything in group_args
#'             (including include_control = TRUE when called from weighting).
#'   Must produce:
#'     - group_col: multi-valued character column (f/bf/df/bdf/control) for WeightIt
#'     - all columns named in dummy_cols: integer 0/1 for feols interaction terms.
#'       NA values in dummies must be replaced with 0L before return.
#'
#' @param group_id    Short identifier.
#' @param group_col   Name of the multi-valued group column (for WeightIt).
#' @param dummy_cols  Character vector of binary dummy column names (for feols).
#' @param group_fun   Function with signature f(df, group_col, ...).
#' @param group_args  Named list forwarded to group_fun.
make_treatment_group_spec <- function(group_id,
                                      group_col,
                                      dummy_cols,
                                      group_fun,
                                      group_args = list()) {
  if (!is.character(group_id)   || length(group_id)   != 1) stop("group_id must be length-1 character.")
  if (!is.character(group_col)  || length(group_col)  != 1) stop("group_col must be length-1 character.")
  if (!is.character(dummy_cols) || length(dummy_cols) == 0)  stop("dummy_cols must be a non-empty character vector.")
  if (!is.function(group_fun))                               stop("group_fun must be a function.")
  if ("group_col" %in% names(group_args)) {
    stop("group_args must not contain 'group_col'.")
  }

  tibble::tibble(
    group_id   = group_id,
    group_col  = group_col,
    dummy_cols = list(dummy_cols),
    group_fun  = list(group_fun),
    group_args = list(group_args)
  )
}


# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3.5 — Aggregation specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build one aggregation spec
#'
#' Each spec runs sunab_aggregate_vcov() once per (run x vcov_spec), producing:
#'   agg_{id}.parquet          estimates + CIs for all vcov_ids (Script 4 plots)
#'   agg_{id}__{vcov_id}.rds   list(coef, vcov, groups) (Script 3 inference)
#'
#' group_fun contract:
#'   Input  — data frame with columns group_1, group_2, ..., term
#'            (capture groups from the agg regex, plus the matched term name)
#'   Output — data frame with column `term` plus AT MINIMUM:
#'              event_time  (integer) — required by Script 3 inference functions
#'              dummy_group (character, e.g. "cd_f") — required by Script 3
#'            Additional columns (e.g. cohort_bin) define finer aggregations.
#'            Rows with NA in any grouping column are dropped before aggregation.
#'
#' @param id       Short identifier used in output filenames; no spaces.
#' @param agg      Perl-compatible regex for sunab_aggregate_vcov(). At least
#'                 one capture group required.
#' @param group_fun Function transforming captured groups; see contract above.
#'                  NULL keeps group_1, group_2, ... as column names, which will
#'                  not satisfy Script 3 expectations — always supply group_fun.
#' @param label    Human-readable description stored in output parquets.
make_agg_spec <- function(id, agg, group_fun = NULL, label = NA_character_) {
  if (!is.character(id)  || length(id)  != 1) stop("id must be length-1 character.")
  if (!is.character(agg) || length(agg) != 1) stop("agg must be length-1 character.")
  if (!is.null(group_fun) && !is.function(group_fun)) stop("group_fun must be NULL or a function.")
  list(id = id, agg = agg, group_fun = group_fun, label = label %||% NA_character_)
}


validate_agg_specs <- function(agg_specs) {
  if (!is.list(agg_specs) || length(agg_specs) == 0) {
    stop("agg_specs must be a non-empty list of make_agg_spec() outputs.")
  }
  for (i in seq_along(agg_specs)) {
    s <- agg_specs[[i]]
    if (!all(c("id", "agg") %in% names(s))) {
      stop("agg_spec [[", i, "]] is missing required fields: id and/or agg.")
    }
    if (!is.character(s$id)  || length(s$id)  != 1) stop("agg_spec[[", i, "]]$id must be length-1 character.")
    if (!is.character(s$agg) || length(s$agg) != 1) stop("agg_spec[[", i, "]]$agg must be length-1 character.")
  }
  ids <- purrr::map_chr(agg_specs, "id")
  if (anyDuplicated(ids)) {
    stop("Duplicate agg_spec ids: ", paste(ids[duplicated(ids)], collapse = ", "))
  }
  invisible(TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Model specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build one model spec row
#'
#' @param model_id         Short identifier.
#' @param formula_template Glue string with {outcome} placeholder. Must contain
#'                         the full sunab():dummy_col interaction structure.
#' @param estimator_type   "sunab" or "twfe".
#' @param term_pattern     Regex matching event-time coefficient names.
#' @param weights_col      Name of the weight column. NA_character_ = unweighted.
#' @param feols_args       Named list of extra args passed to feols().
#'                         mem.clean = TRUE is recommended for large models.
#'                         lean = TRUE is removed with a warning — the pipeline
#'                         needs scores for vcov extraction before discarding
#'                         the model.
make_model_spec <- function(model_id,
                            formula_template,
                            estimator_type  = c("sunab", "twfe"),
                            term_pattern    = ".*",
                            weights_col     = NA_character_,
                            feols_args      = list()) {
  estimator_type <- match.arg(estimator_type)
  if (!is.character(model_id)         || length(model_id)         != 1) stop("model_id must be length-1 character.")
  if (!is.character(formula_template) || length(formula_template) != 1) stop("formula_template must be length-1 character.")
  if ("lean" %in% names(feols_args) && isTRUE(feols_args$lean)) {
    warning(
      "lean = TRUE in feols_args prevents vcov extraction. ",
      "The pipeline extracts all vcovs before discarding the model. Removing lean = TRUE."
    )
    feols_args$lean <- NULL
  }

  tibble::tibble(
    model_id         = model_id,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    term_pattern     = term_pattern,
    weights_col      = normalize_optional_colname(weights_col),
    feols_args       = list(feols_args)
  )
}

# SECTION 5 — Weighting specs
# ══════════════════════════════════════════════════════════════════════════════

make_weighting_spec <- function(weighting_id,
                                weight_formula,
                                method         = "glm",
                                estimand       = "ATO",
                                weighting_name = NA_character_) {
  if (!is.character(weighting_id) || length(weighting_id) != 1) {
    stop("weighting_id must be a length-1 character string.")
  }
  if (!inherits(weight_formula, "formula")) {
    stop("weight_formula must be a formula, e.g. ~ chili + def + aet")
  }
  if (length(attr(stats::terms(weight_formula), "term.labels")) == 0) {
    stop("weight_formula must contain at least one covariate on the RHS.")
  }
  tibble::tibble(
    weighting_id   = weighting_id,
    weight_formula = list(weight_formula),
    method         = method,
    estimand       = estimand,
    weighting_name = normalize_optional_colname(weighting_name)
  )
}

make_weight_col_name <- function(weighting_id, weighting_name) {
  if (!is.na(weighting_name) && nchar(weighting_name) > 0) {
    paste(weighting_id, weighting_name, "weights", sep = "_")
  } else {
    paste(weighting_id, "weights", sep = "_")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Weighting experiment runner
# ══════════════════════════════════════════════════════════════════════════════

run_weighting_experiment <- function(dataset_spec,
                                     analysis_subset_specs,
                                     treatment_group_specs,
                                     weighting_specs,
                                     dir_out,
                                     skip_existing  = TRUE,
                                     verbose_timing = FALSE,
                                     .progress      = FALSE) {

  for (pkg in c("WeightIt", "cobalt", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for weighting. Install with install.packages('", pkg, "').")
    }
  }

  validate_dataset_spec(dataset_spec)
  validate_spec_table(analysis_subset_specs, c("subset_id", "long_data_source", "short_data_source", "data_filter"), "analysis_subset_specs")
  validate_spec_table(treatment_group_specs,  c("group_id", "group_col", "dummy_cols", "group_fun", "group_args"),    "treatment_group_specs")
  validate_spec_table(weighting_specs,         c("weighting_id", "weight_formula", "method", "estimand"),              "weighting_specs")

  missing_short <- analysis_subset_specs |>
    dplyr::filter(purrr::map_lgl(short_data_source, is.null)) |>
    dplyr::pull(subset_id)
  if (length(missing_short) > 0) {
    stop("These subsets have no short_data_source:\n", paste(missing_short, collapse = "\n"))
  }

  dir_ensure_local(c(
    dir_out,
    file.path(dir_out, "weights", "by_run"),
    file.path(dir_out, "metadata")
  ))

  run_grid <- tidyr::crossing(analysis_subset_specs, treatment_group_specs, weighting_specs) |>
    dplyr::mutate(
      weight_run_id = as.character(glue::glue("{subset_id}__{group_id}__{weighting_id}")),
      weight_col    = purrr::map2_chr(weighting_id, weighting_name, make_weight_col_name)
    )

  ts <- timestamp_now()
  saveRDS(
    list(dataset_spec = dataset_spec, analysis_subset_specs = analysis_subset_specs,
         treatment_group_specs = treatment_group_specs, weighting_specs = weighting_specs,
         skip_existing = skip_existing),
    file.path(dir_out, "metadata", glue::glue("weighting_specs_snapshot__{ts}.rds"))
  )

  if (verbose_timing) message(glue::glue("[{timestamp_now()}] Weighting: {nrow(run_grid)} combinations."))

  results <- purrr::pmap(
    list(
      short_data_source = run_grid$short_data_source,
      data_filter       = run_grid$data_filter,
      subset_id         = run_grid$subset_id,
      group_id          = run_grid$group_id,
      group_col         = run_grid$group_col,
      group_fun         = run_grid$group_fun,
      group_args        = run_grid$group_args,
      weighting_id      = run_grid$weighting_id,
      weight_formula    = run_grid$weight_formula,
      method            = run_grid$method,
      estimand          = run_grid$estimand,
      weight_col        = run_grid$weight_col,
      weight_run_id     = run_grid$weight_run_id
    ),
    function(short_data_source, data_filter, subset_id, group_id, group_col,
             group_fun, group_args, weighting_id, weight_formula, method,
             estimand, weight_col, weight_run_id) {

      result <- list(weight_run_id = weight_run_id, error = NULL, skipped = FALSE)
      t_start <- proc.time()

      tryCatch({
        result <- run_one_weighting(
          short_data_source = short_data_source, data_filter = data_filter,
          dataset_spec = dataset_spec, subset_id = subset_id, group_id = group_id,
          group_col = group_col, group_fun = group_fun, group_args = group_args,
          weighting_id = weighting_id, weight_formula = weight_formula,
          method = method, estimand = estimand, weight_col = weight_col,
          weight_run_id = weight_run_id, dir_out = dir_out, skip_existing = skip_existing
        )
      }, error = function(e) {
        result$error <<- conditionMessage(e)
        message(glue::glue("[ERROR] weight_run_id={weight_run_id}: {conditionMessage(e)}"))
      })

      if (verbose_timing) {
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        status  <- if (!is.null(result$error)) "FAILED" else if (isTRUE(result$skipped)) "SKIPPED" else "OK"
        message(glue::glue("[{timestamp_now()}] {weight_run_id} | {status} | {round(elapsed, 1)}s"))
      }
      result
    },
    .progress = .progress
  )

  failed <- purrr::keep(results, \(r) !is.null(r$error))
  if (length(failed) > 0) {
    message(glue::glue("\n{length(failed)} weighting run(s) failed:\n",
                       paste(purrr::map_chr(failed, "weight_run_id"), collapse = "\n")))
  }

  invisible(list(run_grid = run_grid, run_results = results))
}


run_one_weighting <- function(short_data_source, data_filter, dataset_spec,
                               subset_id, group_id, group_col, group_fun, group_args,
                               weighting_id, weight_formula, method, estimand,
                               weight_col, weight_run_id, dir_out, skip_existing = TRUE) {

  run_stub <- glue::glue("{safe_path_component(subset_id)}__{safe_path_component(group_id)}__{safe_path_component(weighting_id)}")
  run_dir  <- file.path(dir_out, "weights", "by_run", run_stub)
  dir_ensure_local(run_dir)

  weights_file              <- file.path(run_dir, "weights.parquet")
  registry_file             <- file.path(run_dir, "registry.parquet")
  bal_rds_file              <- file.path(run_dir, "balance_objects.rds")
  obs_file                  <- file.path(run_dir, "obs.csv")
  love_plot_file            <- file.path(run_dir, "love_plot.png")
  love_plot_fullpairwise_file <- file.path(run_dir, "love_plot_fullpairwise.png")
  run_spec_file             <- file.path(run_dir, "weight_run_spec.rds")
  weightit_file <- file.path(run_dir, "weightit_object.rds")
  

  if (skip_existing && all(file.exists(c(weights_file, registry_file, bal_rds_file, love_plot_file, run_spec_file, weightit_file)))) {
    return(list(weight_run_id = weight_run_id, weights_file = weights_file,
                registry_file = registry_file, skipped = TRUE, error = NULL))
  }

  run_started <- Sys.time()
  df_short    <- load_arrow_data(short_data_source, data_filter)

  df_grouped <- apply_treatment_grouping(
    df         = df_short,
    group_fun  = group_fun,
    group_col  = group_col,
    group_args = c(group_args, list(include_control = TRUE))
  ) |>
    dplyr::filter(!is.na(.data[[group_col]]))

  if (nrow(df_grouped) == 0) stop("No rows remain after grouping for ", weight_run_id)

  rhs_chr   <- sub("^~\\s*", "", paste(deparse(weight_formula), collapse = " "))
  w_formula <- stats::as.formula(paste(group_col, "~", rhs_chr))

  w_out <- WeightIt::weightit(formula = w_formula, data = df_grouped,
                               method = method, estimand = estimand)

  unit_id     <- dataset_spec$unit_id
  weights_tbl <- tibble::tibble(
    !!unit_id     := df_grouped[[unit_id]],
    !!weight_col  := w_out$weights,
    subset_id     = subset_id,
    group_id      = group_id,
    weighting_id  = weighting_id,
    weight_run_id = weight_run_id
  )

  bal_obj <- cobalt::bal.tab(w_out, un = TRUE, pairwise = TRUE, abs = TRUE,
                              stats = c("mean.diffs", "variance.ratios"),
                              thresholds = c(m = 0.1, v = 2))

  readr::write_csv(bal_obj$Observations, obs_file)

  love_p <- cobalt::love.plot(w_out, stats = "mean.diffs", abs = TRUE, pairwise = TRUE,
                               thresholds = c(m = 0.1), var.order = "unadjusted", stars = "raw") +
    ggplot2::labs(title = glue::glue("{subset_id} | {group_id} | {weighting_id} ({method}, {estimand})"))
  ggplot2::ggsave(love_plot_file, plot = love_p, width = 2500, height = 2000, units = "px", bg = "white")

  love_p_fp <- cobalt::love.plot(w_out, stats = "mean.diffs", abs = TRUE, pairwise = TRUE,
                                  thresholds = c(m = 0.1), var.order = "unadjusted", stars = "raw",
                                  which.treat = .all) +
    ggplot2::labs(title = glue::glue("{subset_id} | {group_id} | {weighting_id} ({method}, {estimand})"))
  ggplot2::ggsave(love_plot_fullpairwise_file, plot = love_p_fp, width = 5000, height = 2000, units = "px", bg = "white")

  run_finished  <- Sys.time()
  registry_tbl  <- tibble::tibble(
    weight_run_id     = weight_run_id, subset_id = subset_id, group_id = group_id,
    weighting_id      = weighting_id, unit_id = unit_id, group_col = group_col,
    weight_col        = weight_col, method = method, estimand = estimand,
    weight_formula    = paste(deparse(w_formula), collapse = " "),
    short_data_source = paste(unlist(short_data_source), collapse = " | "),
    data_filter       = data_filter_to_chr(data_filter),
    n_units           = nrow(df_grouped),
    n_groups          = dplyr::n_distinct(df_grouped[[group_col]]),
    group_levels      = paste(sort(unique(df_grouped[[group_col]])), collapse = ","),
    run_started       = run_started, run_finished = run_finished,
    weights_file      = weights_file, registry_file = registry_file,
    bal_rds_file      = bal_rds_file, love_plot_file = love_plot_file, run_spec_file = run_spec_file,
    weightit_file = weightit_file
  )

  run_spec <- list(weight_run_id = weight_run_id, run_stub = run_stub,
                   subset_id = subset_id, group_id = group_id, weighting_id = weighting_id,
                   group_col = group_col, weight_col = weight_col, dataset_spec = dataset_spec,
                   weight_formula = w_formula, method = method, estimand = estimand,
                   short_data_source = short_data_source, data_filter = data_filter,
                   group_fun = group_fun, group_args = group_args)

  arrow::write_parquet(weights_tbl, weights_file)
  arrow::write_parquet(registry_tbl, registry_file)
  saveRDS(bal_obj, bal_rds_file)
  saveRDS(run_spec, run_spec_file)
  saveRDS(w_out, weightit_file)
  
  rm(df_short, df_grouped, weights_tbl, registry_tbl, run_spec, w_out, bal_obj); gc()

  list(weight_run_id = weight_run_id, weights_file = weights_file,
       registry_file = registry_file, skipped = FALSE, error = NULL)
}


rebuild_weighting_tables <- function(dir_out, write_csv = TRUE) {
  dir_by_run <- file.path(dir_out, "weights", "by_run")
  dir_all    <- file.path(dir_out, "weights", "all")
  if (!dir.exists(dir_by_run)) stop("Weights by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)

  registry_files <- list.files(dir_by_run, pattern = "^registry\\.parquet$",
                                recursive = TRUE, full.names = TRUE)
  if (length(registry_files) == 0) stop("No weight registry.parquet files found under: ", dir_by_run)

  registry_tbl <- purrr::map_dfr(registry_files,
    \(f) arrow::read_parquet(f) |> dplyr::mutate(.mtime = file.mtime(f))
  ) |>
    dplyr::group_by(weight_run_id) |>
    dplyr::filter(.mtime == max(.mtime)) |>
    dplyr::ungroup() |>
    dplyr::select(-.mtime)

  pq <- file.path(dir_all, "weight_registry.parquet")
  arrow::write_parquet(registry_tbl, pq)
  if (write_csv) readr::write_csv(registry_tbl, file.path(dir_all, "weight_registry.csv"))
  invisible(list(weight_registry = registry_tbl, files = list(parquet = pq)))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Weights path resolver
# ══════════════════════════════════════════════════════════════════════════════

resolve_weights_parquet <- function(weights_col, subset_id, group_id, dir_out) {
  if (is.na(weights_col)) return(NULL)

  dir_weights_run <- file.path(dir_out, "weights", "by_run")
  if (!dir.exists(dir_weights_run)) {
    stop(glue::glue("Model requires weights_col='{weights_col}' but no weights/by_run exists. ",
                    "Run run_weighting_experiment() first."))
  }

  registry_files <- list.files(dir_weights_run, pattern = "^registry\\.parquet$",
                                recursive = TRUE, full.names = TRUE)
  if (length(registry_files) == 0) stop("No weight registry files found.")

  registry_tbl <- purrr::map_dfr(registry_files, arrow::read_parquet)

  match_row <- registry_tbl |>
    dplyr::filter(.data$subset_id == !!subset_id,
                  .data$group_id  == !!group_id,
                  .data$weight_col == !!weights_col)

  if (nrow(match_row) == 0) {
    stop(glue::glue("No weights for subset_id='{subset_id}', group_id='{group_id}', ",
                    "weight_col='{weights_col}'."))
  }

  match_row <- match_row |>
    dplyr::mutate(.mtime = purrr::map_dbl(weights_file, \(f) as.numeric(file.mtime(f)))) |>
    dplyr::slice_max(.mtime, n = 1, with_ties = FALSE)

  match_row$weights_file[[1]]
}


# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Top-level estimation experiment runner
# ══════════════════════════════════════════════════════════════════════════════

#' Run a full DiD experiment (v7: unified dummy model)
#'
#' Crosses analysis_subset_specs x outcome_specs x treatment_group_specs x
#' model_specs. Each combination fits ONE feols() model covering all dummy
#' groups, then runs all agg_specs x vcov_specs aggregations.
run_experiment <- function(dataset_spec,
                           analysis_subset_specs,
                           outcome_specs,
                           treatment_group_specs,
                           model_specs,
                           vcov_specs,
                           agg_specs,
                           dir_out,
                           group_palette    = NULL,
                           ci_level         = 0.95,
                           run_estimation   = TRUE,
                           run_descriptive  = FALSE,
                           descriptive_args = list(),
                           skip_existing    = TRUE,
                           verbose_timing   = FALSE,
                           .progress        = FALSE) {

  validate_dataset_spec(dataset_spec)
  validate_spec_table(analysis_subset_specs, c("subset_id", "long_data_source", "data_filter"),                  "analysis_subset_specs")
  validate_spec_table(outcome_specs,          c("outcome"),                                                        "outcome_specs")
  validate_spec_table(treatment_group_specs,  c("group_id", "group_col", "dummy_cols", "group_fun", "group_args"), "treatment_group_specs")
  validate_spec_table(model_specs,            c("model_id", "formula_template"),                                   "model_specs")
  validate_spec_table(vcov_specs,             c("vcov_id", "vcov", "vcov_label"),                                  "vcov_specs")
  validate_agg_specs(agg_specs)

  if (!"vcov_vars"   %in% names(vcov_specs))  vcov_specs$vcov_vars  <- replicate(nrow(vcov_specs), character(0), simplify = FALSE)
  if (!"weights_col" %in% names(model_specs)) model_specs$weights_col <- NA_character_
  if (!"feols_args"  %in% names(model_specs)) model_specs$feols_args  <- replicate(nrow(model_specs), list(), simplify = FALSE)

  if (!run_estimation && !run_descriptive) stop("At least one of run_estimation or run_descriptive must be TRUE.")

  desc_args <- list(treated_year_var = dataset_spec$time_var, control_year_var = dataset_spec$time_var)
  desc_args[names(descriptive_args)] <- descriptive_args

  dir_ensure_local(c(dir_out, file.path(dir_out, "tables", "by_run"),
                     file.path(dir_out, "descriptive", "by_run"), file.path(dir_out, "metadata")))

  run_grid <- tidyr::crossing(analysis_subset_specs, outcome_specs, treatment_group_specs, model_specs) |>
    dplyr::mutate(
      run_id   = as.character(glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}")),
      run_stub = purrr::pmap_chr(list(subset_id, outcome, group_id, model_id),
                                  \(s, o, g, m) make_run_stub(s, o, g, m))
    )

  ts <- timestamp_now()
  saveRDS(run_grid, file.path(dir_out, "metadata", glue::glue("run_grid__{ts}.rds")))
  saveRDS(list(dataset_spec = dataset_spec, analysis_subset_specs = analysis_subset_specs,
               outcome_specs = outcome_specs, treatment_group_specs = treatment_group_specs,
               model_specs = model_specs, vcov_specs = vcov_specs, agg_specs = agg_specs,
               group_palette = group_palette, ci_level = ci_level, skip_existing = skip_existing),
          file.path(dir_out, "metadata", glue::glue("specs_snapshot__{ts}.rds")))

  if (verbose_timing) message(glue::glue("[{timestamp_now()}] Starting: {nrow(run_grid)} runs."))

  run_results <- purrr::pmap(
    list(
      long_data_source = run_grid$long_data_source,
      data_filter      = run_grid$data_filter,
      subset_id        = run_grid$subset_id,
      outcome          = run_grid$outcome,
      group_id         = run_grid$group_id,
      group_col        = run_grid$group_col,
      dummy_cols       = run_grid$dummy_cols,
      group_fun        = run_grid$group_fun,
      group_args       = run_grid$group_args,
      model_id         = run_grid$model_id,
      formula_template = run_grid$formula_template,
      estimator_type   = run_grid$estimator_type,
      term_pattern     = run_grid$term_pattern,
      weights_col      = run_grid$weights_col,
      run_id           = run_grid$run_id,
      run_stub         = run_grid$run_stub,
      feols_args       = run_grid$feols_args
    ),
    function(long_data_source, data_filter, subset_id, outcome, group_id, group_col,
             dummy_cols, group_fun, group_args, model_id, formula_template, estimator_type,
             term_pattern, weights_col, run_id, run_stub, feols_args) {

      result  <- list(run_id = run_id, estimation = NULL, descriptive = NULL,
                      error = NULL, skipped = FALSE)
      t_start <- proc.time()

      tryCatch({
        weights_parquet_path <- resolve_weights_parquet(
          weights_col = weights_col, subset_id = subset_id,
          group_id = group_id, dir_out = dir_out
        )

        if (run_estimation) {
          result$estimation <- run_one_estimation(
            data_source          = long_data_source,
            data_filter          = data_filter,
            dataset_spec         = dataset_spec,
            subset_id            = subset_id,
            outcome              = outcome,
            group_id             = group_id,
            group_col            = group_col,
            dummy_cols           = dummy_cols,
            group_fun            = group_fun,
            group_args           = group_args,
            model_id             = model_id,
            formula_template     = formula_template,
            estimator_type       = estimator_type,
            term_pattern         = term_pattern,
            weights_col          = weights_col,
            weights_parquet_path = weights_parquet_path,
            vcov_specs           = vcov_specs,
            agg_specs            = agg_specs,
            dir_out              = dir_out,
            run_id               = run_id,
            run_stub             = run_stub,
            group_palette        = group_palette,
            ci_level             = ci_level,
            skip_existing        = skip_existing,
            feols_args           = feols_args
          )
          result$skipped <- isTRUE(result$estimation$skipped_existing)
        }

        if (run_descriptive) {
          result$descriptive <- run_one_descriptive(
            data_source = long_data_source, data_filter = data_filter,
            dataset_spec = dataset_spec, subset_id = subset_id, outcome = outcome,
            group_id = group_id, group_col = group_col, dummy_cols = dummy_cols,
            group_fun = group_fun, group_args = group_args, model_id = model_id,
            formula_template = formula_template, dir_out = dir_out,
            run_id = run_id, run_stub = run_stub, group_palette = group_palette,
            treated_year_var = desc_args$treated_year_var,
            control_year_var = desc_args$control_year_var, skip_existing = skip_existing
          )
        }

      }, error = function(e) {
        result$error <<- conditionMessage(e)
        message(glue::glue("[ERROR] run_id={run_id}: {conditionMessage(e)}"))
      })

      if (verbose_timing) {
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        status  <- if (!is.null(result$error)) "FAILED" else if (result$skipped) "SKIPPED" else "OK"
        message(glue::glue("[{timestamp_now()}] {run_id} | {status} | {round(elapsed, 1)}s"))
      }
      result
    },
    .progress = .progress
  )

  failed <- purrr::keep(run_results, \(r) !is.null(r$error))
  if (length(failed) > 0) {
    message(glue::glue("\n{length(failed)} run(s) failed:\n",
                       paste(purrr::map_chr(failed, "run_id"), collapse = "\n")))
  }

  info_file <- file.path(dir_out, "metadata", glue::glue("experiment_info__{ts}.rds"))
  saveRDS(list(dir_out = dir_out, timestamp = ts, n_runs_planned = nrow(run_grid),
               n_runs_ok      = sum(purrr::map_lgl(run_results, \(r) is.null(r$error))),
               n_runs_failed  = length(failed),
               n_runs_skipped = sum(purrr::map_lgl(run_results, \(r) isTRUE(r$skipped))),
               failed_run_ids = purrr::map_chr(failed, "run_id")), info_file)

  invisible(list(run_grid = run_grid, run_results = run_results))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Single-run estimation worker
# ══════════════════════════════════════════════════════════════════════════════

run_one_estimation <- function(data_source,
                               data_filter,
                               dataset_spec,
                               subset_id,
                               outcome,
                               group_id,
                               group_col,
                               dummy_cols,
                               group_fun,
                               group_args,
                               model_id,
                               formula_template,
                               estimator_type,
                               term_pattern,
                               weights_col,
                               weights_parquet_path,
                               vcov_specs,
                               agg_specs,
                               dir_out,
                               run_id,
                               run_stub,
                               group_palette  = NULL,
                               ci_level       = 0.95,
                               skip_existing  = TRUE,
                               feols_args     = list()) {
  
  on.exit({
    suppressWarnings(rm(df, df_grouped, df_est, model, V_list,
                        coef_rows, pq_rows, agg_pq_tbl,
                        coef_tbl, dummy_group_tbl, support_tbl))
    gc()
  }, add = TRUE)
  
  # ── Logging ──────────────────────────────────────────────────────────────────
  status_file <- file.path(here::here(), "status_file.txt")
  log_status <- function(msg) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", run_id, " | ", msg)
    cat(line, "\n")
    write(line, file = status_file, append = TRUE)
  }
  log_mem <- function(msg) {
    # full = TRUE forces a complete collection before reading heap size
    g       <- gc(reset = FALSE, verbose = FALSE, full = TRUE)
    used_mb <- round(sum(g[, 2]) / 1024, 1)
    max_mb  <- round(sum(g[, 4]) / 1024, 1)
    log_status(glue::glue("{msg} | heap_used={used_mb}MB heap_max={max_mb}MB"))
  }
  
  # ── Output paths ─────────────────────────────────────────────────────────────
  run_dir <- file.path(dir_out, "tables", "by_run", run_stub)
  dir_ensure_local(run_dir)
  log_status(glue::glue("START run_stub={run_stub} fixest_version={packageVersion('fixest')}"))
  
  coef_file        <- file.path(run_dir, "coef.parquet")
  dummy_group_file <- file.path(run_dir, "dummy_group.parquet")
  support_file     <- file.path(run_dir, "support.parquet")
  registry_file    <- file.path(run_dir, "registry.parquet")
  run_spec_file    <- file.path(run_dir, "run_spec.rds")
  
  if (skip_existing && file.exists(run_spec_file)) {
    log_status("SKIP run_spec.rds exists")
    return(list(coef_file = coef_file, dummy_group_file = dummy_group_file,
                support_file = support_file, registry_file = registry_file,
                run_spec_file = run_spec_file, skipped_existing = TRUE))
  }
  
  fit_started <- Sys.time()
  
  # ── Load data ────────────────────────────────────────────────────────────────
  log_mem("before load_arrow_data")
  df <- load_arrow_data(data_source, data_filter)
  n_rows_read <- nrow(df)
  log_mem(glue::glue("after load n_rows={n_rows_read}"))
  
  # ── Join weights ─────────────────────────────────────────────────────────────
  if (!is.na(weights_col) && !is.null(weights_parquet_path)) {
    log_status(glue::glue("joining weights col={weights_col}"))
    weights_tbl <- arrow::read_parquet(weights_parquet_path) |>
      dplyr::select(dplyr::all_of(c(dataset_spec$unit_id, weights_col)))
    df <- dplyr::left_join(df, weights_tbl, by = dataset_spec$unit_id)
    n_na_w <- sum(is.na(df[[weights_col]]))
    if (n_na_w > 0) warning(glue::glue("[{run_id}] {n_na_w} rows have NA weights."))
    log_mem("after weights join")
  }
  
  # ── Apply grouping ────────────────────────────────────────────────────────────
  log_mem("before apply_treatment_grouping")
  df_grouped <- apply_treatment_grouping(df = df, group_fun = group_fun,
                                         group_col = group_col, group_args = group_args)
  rm(df)
  log_mem("after grouping rm(df)")
  
  missing_dummy <- setdiff(dummy_cols, names(df_grouped))
  if (length(missing_dummy) > 0) {
    stop(glue::glue("group_fun did not produce expected dummy columns: ",
                    paste(missing_dummy, collapse = ", ")))
  }
  
  # ── Drop treated units assigned to no group ──────────────────────────────────
  trt_col       <- dataset_spec$trt_col
  all_zero_mask <- df_grouped[[trt_col]] == 1 &
    rowSums(as.matrix(dplyr::select(df_grouped, dplyr::all_of(dummy_cols)))) == 0
  n_orphan <- sum(all_zero_mask, na.rm = TRUE)
  if (n_orphan > 0) {
    warning(glue::glue("[{run_id}] Dropping {n_orphan} treated units with all dummy_cols == 0."))
    df_grouped <- df_grouped[!all_zero_mask, ]
  }
  log_status(glue::glue("orphan treated units dropped: {n_orphan}"))
  
  # ── Build formula and prune columns ──────────────────────────────────────────
  model_formula <- build_model_formula(formula_template, outcome)
  needed_cols   <- get_needed_columns(
    formula     = model_formula, trt_col = trt_col, group_col = group_col,
    dummy_cols  = dummy_cols,    unit_id = dataset_spec$unit_id,
    event_id    = dataset_spec$event_id, vcov_vars = vcov_specs$vcov_vars,
    weights_col = weights_col
  )
  check_required_columns(df_grouped, needed_cols, context = run_id)
  df_grouped <- df_grouped |> dplyr::select(dplyr::all_of(needed_cols))
  log_status(glue::glue("columns pruned [{ncol(df_grouped)}]"))
  
  feols_weights <- if (!is.na(weights_col) && weights_col %in% names(df_grouped)) {
    stats::as.formula(paste("~", weights_col))
  } else NULL
  
  # ── Fit model ─────────────────────────────────────────────────────────────────
  log_mem("before feols")
  model <- do.call(
    fixest::feols,
    c(list(fml = model_formula, data = df_grouped, weights = feols_weights), feols_args)
  )
  log_mem("after feols")
  
  # ── Recover exact estimation sample ──────────────────────────────────────────
  df_est <- df_grouped[fixest::obs(model), , drop = FALSE]
  rm(df_grouped)
  log_mem("after rm(df_grouped)")
  
  # ── Precompute all vcov matrices (requires scores; must precede rm(model)) ────
  log_status("precomputing vcov matrices")
  V_list <- vector("list", nrow(vcov_specs))
  names(V_list) <- vcov_specs$vcov_id
  
  for (vi in seq_len(nrow(vcov_specs))) {
    v_id   <- vcov_specs$vcov_id[[vi]]
    v_spec <- vcov_specs$vcov[[vi]]
    V_list[[v_id]] <- tryCatch(
      if (is.null(v_spec)) stats::vcov(model) else stats::vcov(model, vcov = v_spec),
      error = function(e) {
        warning(glue::glue("[{run_id}] vcov_id={v_id}: {e$message}. Skipping."))
        NULL
      }
    )
    log_status(glue::glue("vcov {v_id}: {if (!is.null(V_list[[v_id]])) 'OK' else 'FAILED'}"))
  }
  
  # ── Raw coef tables — write immediately and clear before aggregation loop ─────
  # coef_rows is not needed by the aggregation loop; clearing it here frees
  # several hundred MB before the most memory-intensive phase begins
  log_status("extracting coef tables")
  coef_rows <- vector("list", nrow(vcov_specs))
  
  for (vi in seq_len(nrow(vcov_specs))) {
    v_id    <- vcov_specs$vcov_id[[vi]]
    v_label <- vcov_specs$vcov_label[[vi]]
    V       <- V_list[[v_id]]
    if (is.null(V)) next
    
    ct <- tryCatch(
      fixest::coeftable(model, vcov = V),
      error = function(e) {
        warning(glue::glue("[{run_id}] coeftable vcov_id={v_id}: {e$message}"))
        NULL
      }
    )
    if (is.null(ct)) next
    
    coef_rows[[vi]] <- coeftable_to_tbl(
      ct           = ct,
      vcov_id      = v_id,
      vcov_label   = v_label,
      term_pattern = term_pattern,
      ci_level     = ci_level,
      meta         = list(subset_id = subset_id, outcome = outcome, group_id = group_id,
                          group_col = group_col, model_id = model_id,
                          formula_template = formula_template, estimator_type = estimator_type,
                          weights_col = weights_col, group_palette = group_palette),
      run_id       = run_id
    )
  }
  
  coef_tbl <- dplyr::bind_rows(purrr::compact(coef_rows))
  arrow::write_parquet(coef_tbl, coef_file)
  rm(coef_rows, coef_tbl)
  gc()
  log_mem("after coef write rm(coef_rows coef_tbl)")
  
  # ── Aggregation loop: agg_spec x vcov_spec ────────────────────────────────────
  log_status(glue::glue(
    "running aggregations: {length(agg_specs)} agg_specs x {nrow(vcov_specs)} vcov_specs"
  ))
  
  for (ai in seq_along(agg_specs)) {
    agg_spec    <- agg_specs[[ai]]
    agg_id_safe <- safe_path_component(agg_spec$id)
    pq_rows     <- vector("list", nrow(vcov_specs))
    
    for (vi in seq_len(nrow(vcov_specs))) {
      v_id    <- vcov_specs$vcov_id[[vi]]
      v_label <- vcov_specs$vcov_label[[vi]]
      V       <- V_list[[v_id]]
      if (is.null(V)) next
      
      result <- tryCatch(
        sunab_aggregate_vcov(
          sunab_fixest  = model,
          agg           = agg_spec$agg,
          vcov_mat      = V,
          group_fun     = agg_spec$group_fun,
          weight_method = "data_count",
          df_est        = df_est,
          cohort_var    = dataset_spec$cohort_var,
          period_var    = dataset_spec$time_var,
          weight_var    = if (!is.na(weights_col) && weights_col %in% names(df_est)) weights_col else NULL
        ),
        error = function(e) {
          warning(glue::glue("[{run_id}] agg_id={agg_spec$id} vcov_id={v_id}: {e$message}"))
          NULL
        }
      )
      if (is.null(result)) next
      
      agg_obj <- list(
        coef   = as.numeric(result$beta),
        vcov   = result$sigma,
        groups = result$groups
      )
      agg_rds_file <- file.path(run_dir,
                                glue::glue("agg_{agg_id_safe}__{safe_path_component(v_id)}.rds"))
      saveRDS(agg_obj, agg_rds_file)
      log_status(glue::glue("saved {basename(agg_rds_file)}"))
      
      z_val <- stats::qnorm(1 - (1 - ci_level) / 2)
      pq_rows[[vi]] <- result$groups |>
        dplyr::mutate(
          ci_lower   = estimate - z_val * se,
          ci_upper   = estimate + z_val * se,
          vcov_id    = v_id,
          vcov_label = v_label,
          agg_id     = agg_spec$id,
          agg_label  = agg_spec$label %||% NA_character_,
          subset_id  = subset_id,
          outcome    = outcome,
          group_id   = group_id,
          model_id   = model_id,
          run_id     = run_id
        )
    }
    
    agg_pq_tbl <- dplyr::bind_rows(purrr::compact(pq_rows))
    if (nrow(agg_pq_tbl) > 0) {
      agg_pq_file <- file.path(run_dir, glue::glue("agg_{agg_id_safe}.parquet"))
      arrow::write_parquet(agg_pq_tbl, agg_pq_file)
      log_status(glue::glue("saved {basename(agg_pq_file)} [{nrow(agg_pq_tbl)} rows]"))
    }
    
    # clear per-agg-spec objects before next iteration
    rm(pq_rows, agg_pq_tbl)
    gc()
    log_mem(glue::glue("after agg_spec={agg_spec$id}"))
  }
  
  # ── Discard model — scores no longer needed ───────────────────────────────────
  rm(model)
  V_list <- NULL
  gc()
  log_mem("after rm(model) gc()")
  
  fit_finished <- Sys.time()
  
  # ── N stats and event-time support ───────────────────────────────────────────
  dummy_group_tbl <- make_dummy_group_summary(
    df_est = df_est, dummy_cols = dummy_cols, trt_col = trt_col,
    unit_id = dataset_spec$unit_id, event_id = dataset_spec$event_id,
    meta = list(subset_id = subset_id, outcome = outcome, group_id = group_id,
                group_col = group_col, model_id = model_id,
                formula_template = formula_template, estimator_type = estimator_type,
                weights_col = weights_col, group_palette = group_palette),
    run_id = run_id
  )
  
  support_tbl <- make_event_time_support_unified(
    df_est     = df_est,     dummy_cols = dummy_cols, trt_col  = trt_col,
    unit_id    = dataset_spec$unit_id, event_id = dataset_spec$event_id,
    cohort_var = dataset_spec$cohort_var,    time_var = dataset_spec$time_var
  ) |>
    dplyr::mutate(subset_id = subset_id, outcome = outcome,
                  group_id  = group_id,  model_id = model_id, run_id = run_id)
  
  rm(df_est); gc()
  log_mem("after rm(df_est) gc()")
  
  # ── Registry and run spec ─────────────────────────────────────────────────────
  registry_tbl <- make_estimation_registry(
    run_id               = run_id,
    subset_id            = subset_id,
    outcome              = outcome,
    group_id             = group_id,
    model_id             = model_id,
    group_col            = group_col,
    dummy_cols           = dummy_cols,
    dataset_spec         = dataset_spec,
    formula_template     = formula_template,
    estimator_type       = estimator_type,
    term_pattern         = term_pattern,
    weights_col          = weights_col,
    weights_parquet_path = weights_parquet_path,
    data_filter          = data_filter,
    data_source          = data_source,
    vcov_specs           = vcov_specs,
    agg_specs            = agg_specs,
    group_args           = group_args,
    n_rows_read          = n_rows_read,
    dummy_group_tbl      = dummy_group_tbl,
    fit_started          = fit_started,
    fit_finished         = fit_finished,
    coef_file            = coef_file,
    dummy_group_file     = dummy_group_file,
    support_file         = support_file,
    registry_file        = registry_file,
    run_spec_file        = run_spec_file
  )
  
  run_spec <- list(
    run_id = run_id, run_stub = run_stub, subset_id = subset_id, outcome = outcome,
    group_id = group_id, model_id = model_id, group_col = group_col, dummy_cols = dummy_cols,
    dataset_spec = dataset_spec, formula_template = formula_template,
    estimator_type = estimator_type, term_pattern = term_pattern,
    weights_col = weights_col, weights_parquet_path = weights_parquet_path,
    data_filter = data_filter, data_source = data_source,
    group_fun = group_fun, group_args = group_args,
    vcov_specs = vcov_specs, agg_specs = agg_specs,
    group_palette = group_palette, ci_level = ci_level
  )
  
  log_status("writing final parquet outputs")
  arrow::write_parquet(dummy_group_tbl, dummy_group_file)
  arrow::write_parquet(support_tbl,     support_file)
  arrow::write_parquet(registry_tbl,    registry_file)
  saveRDS(run_spec, run_spec_file)
  
  rm(dummy_group_tbl, support_tbl, registry_tbl, run_spec)
  gc()
  log_mem("END outputs written")
  
  list(coef_file = coef_file, dummy_group_file = dummy_group_file,
       support_file = support_file, registry_file = registry_file,
       run_spec_file = run_spec_file, skipped_existing = FALSE)
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — Aggregation note
# ══════════════════════════════════════════════════════════════════════════════
#
# sunab_aggregate_vcov() and compute_sunab_coef_weights_from_data() are defined
# in sunab_aggregate_vcov.R. Source that file before sourcing this pipeline.
# They replace the old make_lean(), build_agg_eventstudy_from_bV(), and
# build_agg_eventstudy_cohortgroup_from_bV(), which are no longer used.

# SECTION 11 — Inference helpers (called from Script 3)
# ══════════════════════════════════════════════════════════════════════════════

#' Joint Wald test that all pre-trend estimates are zero (flatness test)
#'
#' @param agg_obj     Output of build_agg_eventstudy_from_bV.
#' @param dummy_group Name of the dummy group to test (e.g. "cd_f").
#' @param pre_years   Integer vector of event-time periods to include.
#' @return One-row tibble with test statistics and p-value.
wald_pretrend_flat <- function(agg_obj, dummy_group, pre_years) {
  idx <- which(agg_obj$groups$dummy_group == dummy_group &
                 agg_obj$groups$event_time  %in% pre_years)
  idx <- idx[order(agg_obj$groups$event_time[idx])]
  
  if (length(idx) == 0) {
    warning("No pre-trend periods found for dummy_group=", dummy_group)
    return(tibble::tibble(
      dummy_group = dummy_group, stat = NA_real_, df = NA_integer_, p = NA_real_,
      n_periods_used = 0L, n_periods_theoretical = length(pre_years),
      frac_periods_used = 0
    ))
  }
  
  b    <- agg_obj$coef[idx]
  V    <- agg_obj$vcov[idx, idx, drop = FALSE]
  stat <- tryCatch(
    as.numeric(t(b) %*% solve(V) %*% b),
    error = function(e) { warning("Singular vcov in wald_pretrend_flat: ", e$message); NA_real_ }
  )
  
  tibble::tibble(
    dummy_group           = dummy_group,
    stat                  = stat,
    df                    = length(idx),
    p                     = if (!is.na(stat)) pchisq(stat, df = length(idx), lower.tail = FALSE) else NA_real_,
    n_periods_used        = length(idx),
    n_periods_theoretical = length(pre_years),
    frac_periods_used     = length(idx) / length(pre_years)
  )
}


#' GLS test for a linear pre-trend slope
#'
#' Fits intercept + slope via GLS using the aggregated vcov as the error
#' covariance matrix and tests whether slope == 0.
wald_pretrend_slope <- function(agg_obj, dummy_group, pre_years) {
  idx <- which(agg_obj$groups$dummy_group == dummy_group &
                 agg_obj$groups$event_time  %in% pre_years)
  idx <- idx[order(agg_obj$groups$event_time[idx])]
  
  na_row <- tibble::tibble(
    dummy_group = dummy_group, slope = NA_real_, slope_se = NA_real_,
    stat = NA_real_, df = 1L, p = NA_real_,
    n_periods_used = length(idx), n_periods_theoretical = length(pre_years),
    frac_periods_used = length(idx) / length(pre_years)
  )
  
  if (length(idx) < 2) return(na_row)
  
  y <- agg_obj$coef[idx]; V <- agg_obj$vcov[idx, idx, drop = FALSE]
  t <- agg_obj$groups$event_time[idx]; X <- cbind(intercept = 1, slope = t)
  
  Vinv <- tryCatch(solve(V), error = function(e) {
    warning("Singular V in wald_pretrend_slope (", dummy_group, "): ", e$message); NULL
  })
  if (is.null(Vinv)) return(na_row)
  
  XVX      <- t(X) %*% Vinv %*% X
  beta_hat <- tryCatch(solve(XVX) %*% t(X) %*% Vinv %*% y, error = function(e) {
    warning("Singular XVX in wald_pretrend_slope (", dummy_group, "): ", e$message); NULL
  })
  if (is.null(beta_hat)) return(na_row)
  
  vcov_beta <- solve(XVX)
  slope_est <- beta_hat["slope", 1]
  slope_se  <- sqrt(vcov_beta["slope", "slope"])
  stat      <- (slope_est / slope_se)^2
  
  tibble::tibble(
    dummy_group           = dummy_group,
    slope                 = slope_est,
    slope_se              = slope_se,
    stat                  = stat,
    df                    = 1L,
    p                     = pchisq(stat, df = 1, lower.tail = FALSE),
    n_periods_used        = length(idx),
    n_periods_theoretical = length(pre_years),
    frac_periods_used     = length(idx) / length(pre_years)
  )
}


#' Average treatment effect (ATT) over a post-treatment window for one group
att_window <- function(agg_obj, dummy_group, years) {
  idx <- which(agg_obj$groups$dummy_group == dummy_group &
                 agg_obj$groups$event_time  %in% years)
  
  if (length(idx) == 0) {
    return(tibble::tibble(
      dummy_group = dummy_group, estimate = NA_real_, se = NA_real_,
      stat = NA_real_, df = 1L, p = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      n_periods_used = 0L, n_periods_theoretical = length(years),
      frac_periods_used = 0
    ))
  }
  
  L      <- numeric(length(agg_obj$coef))
  L[idx] <- 1 / length(idx)
  est    <- sum(L * agg_obj$coef)
  se     <- sqrt(as.numeric(t(L) %*% agg_obj$vcov %*% L))
  stat   <- (est / se)^2
  
  tibble::tibble(
    dummy_group           = dummy_group,
    estimate              = est,
    se                    = se,
    stat                  = stat,
    df                    = 1L,
    p                     = pchisq(stat, df = 1, lower.tail = FALSE),
    ci_lower              = est - 1.96 * se,
    ci_upper              = est + 1.96 * se,
    n_periods_used        = length(idx),
    n_periods_theoretical = length(years),
    frac_periods_used     = length(idx) / length(years)
  )
}


#' GLS slope estimate over a post-treatment window for one group
gls_slope_one_group <- function(agg_obj, dummy_group, years) {
  idx <- which(agg_obj$groups$dummy_group == dummy_group &
                 agg_obj$groups$event_time  %in% years)
  idx <- idx[order(agg_obj$groups$event_time[idx])]
  
  na_row <- tibble::tibble(
    dummy_group = dummy_group, slope = NA_real_, slope_se = NA_real_,
    stat = NA_real_, p = NA_real_,
    n_periods_used = length(idx), n_periods_theoretical = length(years),
    frac_periods_used = length(idx) / length(years)
  )
  
  if (length(idx) < 2) return(na_row)
  
  y <- agg_obj$coef[idx]; V <- agg_obj$vcov[idx, idx, drop = FALSE]
  t <- agg_obj$groups$event_time[idx]; X <- cbind(intercept = 1, slope = t)
  
  Vinv <- tryCatch(solve(V), error = function(e) {
    warning("Singular V in gls_slope_one_group (", dummy_group, "): ", e$message); NULL
  })
  if (is.null(Vinv)) return(na_row)
  
  XVX      <- t(X) %*% Vinv %*% X
  beta_hat <- tryCatch(solve(XVX) %*% t(X) %*% Vinv %*% y, error = function(e) {
    warning("Singular XVX in gls_slope_one_group (", dummy_group, "): ", e$message); NULL
  })
  if (is.null(beta_hat)) return(na_row)
  
  vcov_beta <- solve(XVX)
  slope_est <- beta_hat["slope", 1]
  slope_se  <- sqrt(vcov_beta["slope", "slope"])
  stat      <- (slope_est / slope_se)^2
  
  tibble::tibble(
    dummy_group           = dummy_group,
    slope                 = slope_est,
    slope_se              = slope_se,
    stat                  = stat,
    p                     = pchisq(stat, df = 1, lower.tail = FALSE),
    n_periods_used        = length(idx),
    n_periods_theoretical = length(years),
    frac_periods_used     = length(idx) / length(years)
  )
}


#' Pairwise contrast of average ATTs over a window
wald_compare_att <- function(agg_obj, group_a, group_b, years) {
  idx_a <- which(agg_obj$groups$dummy_group == group_a & agg_obj$groups$event_time %in% years)
  idx_b <- which(agg_obj$groups$dummy_group == group_b & agg_obj$groups$event_time %in% years)
  
  # n_periods_used is the overlap — both groups must contribute to the contrast
  n_used <- min(length(idx_a), length(idx_b))
  
  if (length(idx_a) == 0 || length(idx_b) == 0) {
    return(tibble::tibble(
      contrast = paste(group_a, "vs", group_b), group_a = group_a, group_b = group_b,
      estimate = NA_real_, se = NA_real_, stat = NA_real_, df = 1L, p = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      n_periods_used = n_used, n_periods_theoretical = length(years),
      frac_periods_used = n_used / length(years)
    ))
  }
  
  L        <- numeric(length(agg_obj$coef))
  L[idx_a] <-  1 / length(idx_a)
  L[idx_b] <- -1 / length(idx_b)
  est  <- sum(L * agg_obj$coef)
  se   <- sqrt(as.numeric(t(L) %*% agg_obj$vcov %*% L))
  stat <- (est / se)^2
  
  tibble::tibble(
    contrast              = paste(group_a, "vs", group_b),
    group_a               = group_a,
    group_b               = group_b,
    estimate              = est,
    se                    = se,
    stat                  = stat,
    df                    = 1L,
    p                     = pchisq(stat, df = 1, lower.tail = FALSE),
    ci_lower              = est - 1.96 * se,
    ci_upper              = est + 1.96 * se,
    n_periods_used        = n_used,
    n_periods_theoretical = length(years),
    frac_periods_used     = n_used / length(years)
  )
}

# a version of the att-comparison, but normalized by a reference year (e.g. year prior to burn)
wald_compare_att_normalized <- function(agg_obj, group_a, group_b, years,
                                        ref_year = -1L, ci_level = 0.95) {
  groups <- agg_obj$groups
  coef   <- agg_obj$coef
  sigma  <- agg_obj$vcov
  
  idx_a_win <- which(groups$dummy_group == group_a & groups$event_time %in% years)
  idx_a_ref <- which(groups$dummy_group == group_a & groups$event_time == ref_year)
  idx_b_win <- which(groups$dummy_group == group_b & groups$event_time %in% years)
  
  if (any(c(length(idx_a_win), length(idx_b_win)) == 0) || length(idx_a_ref) != 1) {
    warning(glue::glue("Missing support: {group_a} window({length(idx_a_win)}) ",
                       "{group_a} ref({length(idx_a_ref)}) {group_b} window({length(idx_b_win)})"))
    return(tibble::tibble())
  }
  
  # contrast: (mean(group_a[window]) - group_a[ref_year]) - mean(group_b[window])
  c_vec            <- numeric(length(coef))
  c_vec[idx_a_win] <-  1 / length(idx_a_win)
  c_vec[idx_a_ref] <- -1                        # single year, not averaged
  c_vec[idx_b_win] <- -1 / length(idx_b_win)
  
  est  <- as.numeric(c_vec %*% coef)
  se   <- sqrt(as.numeric(c_vec %*% sigma %*% c_vec))
  z    <- est / se
  p    <- 2 * pnorm(-abs(z))
  crit <- qnorm(1 - (1 - ci_level) / 2)
  
  tibble::tibble(
    contrast = glue::glue("{group_a}[win-ref{ref_year}]_vs_{group_b}[win]"),
    group_a  = group_a,
    group_b  = group_b,
    ref_year = ref_year,
    estimate = est,
    se       = se,
    z        = z,
    p        = p,
    ci_lower    = est - crit * se,
    ci_upper    = est + crit * se
  )
}


wald_compare_att_multi <- function(agg_obj, group_a, group_b, group_c, years,
                                   ci_level = 0.95) {
  groups <- agg_obj$groups
  coef   <- agg_obj$coef
  sigma  <- agg_obj$vcov
  
  idx_a <- which(groups$dummy_group == group_a & groups$event_time %in% years)
  idx_b <- which(groups$dummy_group == group_b & groups$event_time %in% years)
  idx_c <- which(groups$dummy_group == group_c & groups$event_time %in% years)
  
  if (any(c(length(idx_a), length(idx_b), length(idx_c)) == 0)) {
    warning(glue::glue("Missing window support for one or more groups: ",
                       "{group_a}({length(idx_a)}) {group_b}({length(idx_b)}) ",
                       "{group_c}({length(idx_c)})"))
    return(tibble::tibble())
  }
  
  # contrast vector: ATT_a + ATT_b - ATT_c
  # each group is averaged separately over its own window indices
  c_vec          <- numeric(length(coef))
  c_vec[idx_a]   <-  1 / length(idx_a)
  c_vec[idx_b]   <-  1 / length(idx_b)
  c_vec[idx_c]   <- -1 / length(idx_c)
  
  est    <- as.numeric(c_vec %*% coef)
  var_e  <- as.numeric(c_vec %*% sigma %*% c_vec)
  se     <- sqrt(var_e)
  z      <- est / se
  p      <- 2 * pnorm(-abs(z))
  crit   <- qnorm(1 - (1 - ci_level) / 2)
  
  tibble::tibble(
    contrast = glue::glue("{group_a}+{group_b}_vs_{group_c}"),
    group_a  = group_a,
    group_b  = group_b,
    group_c  = group_c,
    estimate = est,    # positive = (a+b) > c
    se       = se,
    z        = z,
    p        = p,
    ci_lower    = est - crit * se,
    ci_upper    = est + crit * se
  )
}




#' Pairwise contrast of post-treatment GLS slopes
compare_gls_slopes <- function(agg_obj, group_a, group_b, years) {
  idx_a <- which(agg_obj$groups$dummy_group == group_a & agg_obj$groups$event_time %in% years)
  idx_b <- which(agg_obj$groups$dummy_group == group_b & agg_obj$groups$event_time %in% years)
  
  n_used <- min(length(idx_a), length(idx_b))
  
  na_row <- tibble::tibble(
    contrast = paste(group_a, "vs", group_b), group_a = group_a, group_b = group_b,
    slope_diff = NA_real_, slope_se = NA_real_, stat = NA_real_, p = NA_real_,
    n_periods_used = n_used, n_periods_theoretical = length(years),
    frac_periods_used = n_used / length(years)
  )
  
  if (length(idx_a) < 2 || length(idx_b) < 2) return(na_row)
  
  idx    <- c(idx_a, idx_b)
  df_sub <- agg_obj$groups[idx, ] |>
    dplyr::mutate(estimate = agg_obj$coef[idx],
                  group    = factor(dummy_group, levels = c(group_a, group_b)))
  V <- agg_obj$vcov[idx, idx, drop = FALSE]
  X <- model.matrix(~ event_time * group, data = df_sub)
  
  Vinv <- tryCatch(solve(V), error = function(e) {
    warning("Singular V in compare_gls_slopes (", group_a, " vs ", group_b, "): ", e$message); NULL
  })
  if (is.null(Vinv)) return(na_row)
  
  XVX      <- t(X) %*% Vinv %*% X
  beta_hat <- tryCatch(solve(XVX) %*% t(X) %*% Vinv %*% df_sub$estimate, error = function(e) {
    warning("Singular XVX in compare_gls_slopes (", group_a, " vs ", group_b, "): ", e$message); NULL
  })
  if (is.null(beta_hat)) return(na_row)
  
  vcov_beta      <- solve(XVX)
  interaction_nm <- grep("^event_time:group", rownames(vcov_beta), value = TRUE)
  
  if (length(interaction_nm) == 0) {
    warning("No interaction term in compare_gls_slopes (", group_a, " vs ", group_b, ")")
    return(na_row)
  }
  
  slope_diff <- beta_hat[interaction_nm, 1]
  slope_se   <- sqrt(vcov_beta[interaction_nm, interaction_nm])
  stat       <- (slope_diff / slope_se)^2
  
  tibble::tibble(
    contrast              = paste(group_a, "vs", group_b),
    group_a               = group_a,
    group_b               = group_b,
    slope_diff            = slope_diff,
    slope_se              = slope_se,
    stat                  = stat,
    p                     = pchisq(stat, df = 1, lower.tail = FALSE),
    n_periods_used        = n_used,
    n_periods_theoretical = length(years),
    frac_periods_used     = n_used / length(years)
  )
}


# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — Rebuild (merge across runs)
# ══════════════════════════════════════════════════════════════════════════════

rebuild_estimation_tables <- function(dir_out, write_csv = TRUE, recursive = TRUE) {
  dir_by_run <- file.path(dir_out, "tables", "by_run")
  dir_all    <- file.path(dir_out, "tables", "all")
  if (!dir.exists(dir_by_run)) stop("by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)

  find_pq <- function(pattern) {
    list.files(dir_by_run, pattern = pattern, recursive = recursive, full.names = TRUE)
  }

  coef_files        <- find_pq("^coef\\.parquet$")
  dummy_group_files <- find_pq("^dummy_group\\.parquet$")
  support_files     <- find_pq("^support\\.parquet$")
  registry_files    <- find_pq("^registry\\.parquet$")
  run_spec_files    <- list.files(dir_by_run, "^run_spec\\.rds$",
                                   recursive = recursive, full.names = TRUE)

  # agg parquets: one per agg_id per run, named agg_{agg_id}.parquet
  # (RDS files match agg_*__*.rds — the double underscore separates agg_id
  # from vcov_id; they are NOT picked up here)
  agg_pq_files <- find_pq("^agg_[^.]+\\.parquet$")

  for (nm in c("coef_files", "dummy_group_files", "support_files", "registry_files")) {
    if (length(get(nm)) == 0) {
      message("No ", sub("_files", ".parquet", nm), " found under: ", dir_by_run)
    }
  }

  read_and_dedup <- function(files, id_col) {
    if (length(files) == 0) return(NULL)
    purrr::map(files, \(f) {
      x <- arrow::read_parquet(f)
      x[] <- lapply(x, \(col) if (is.list(col)) as.character(col) else col)
      x$.mtime <- file.mtime(f)
      x
    }) |>
      dplyr::bind_rows() |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }

  write_pair <- function(tbl, stem) {
    if (is.null(tbl) || nrow(tbl) == 0) return(list(parquet = NULL, csv = NULL))
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }

  coef_tbl        <- read_and_dedup(coef_files,        c("run_id", "term", "vcov_id"))
  dummy_group_tbl <- read_and_dedup(dummy_group_files, c("run_id", "dummy_group"))
  support_tbl     <- read_and_dedup(support_files,     c("run_id", "dummy_group", "event_time"))
  registry_tbl    <- read_and_dedup(registry_files,    c("run_id", "dummy_group"))

  # merge agg parquets per agg_id; different agg_ids may have different group columns
  agg_stems <- unique(basename(agg_pq_files))
  agg_tbls  <- purrr::set_names(
    purrr::map(agg_stems, function(stem) {
      stem_files <- agg_pq_files[basename(agg_pq_files) == stem]
      read_and_dedup(stem_files, c("run_id", "vcov_id", "key"))
    }),
    tools::file_path_sans_ext(agg_stems)
  )

  spec_index_tbl <- if (length(run_spec_files) > 0) {
    purrr::map_dfr(run_spec_files, \(f) {
      spec <- readRDS(f)
      tibble::tibble(
        run_id           = spec$run_id           %||% NA_character_,
        subset_id        = spec$subset_id         %||% NA_character_,
        outcome          = spec$outcome           %||% NA_character_,
        group_id         = spec$group_id          %||% NA_character_,
        model_id         = spec$model_id          %||% NA_character_,
        group_col        = spec$group_col         %||% NA_character_,
        dummy_cols       = paste(unlist(spec$dummy_cols), collapse = ","),
        estimator_type   = spec$estimator_type    %||% NA_character_,
        formula_template = spec$formula_template  %||% NA_character_,
        weights_col      = spec$weights_col       %||% NA_character_,
        data_filter      = data_filter_to_chr(spec$data_filter),
        data_source      = paste(unlist(spec$data_source), collapse = " | "),
        vcov_ids         = paste(spec$vcov_specs$vcov_id, collapse = " | "),
        agg_ids          = paste(purrr::map_chr(spec$agg_specs, "id"), collapse = " | "),
        run_spec_file    = f,
        .mtime           = file.mtime(f)
      )
    }) |>
      dplyr::group_by(run_id) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  } else NULL

  files <- list(
    coef               = write_pair(coef_tbl,        "coef"),
    dummy_group        = write_pair(dummy_group_tbl, "dummy_group_summary"),
    event_time_support = write_pair(support_tbl,     "event_time_support"),
    run_registry       = write_pair(registry_tbl,    "run_registry"),
    run_spec_index     = write_pair(spec_index_tbl,  "run_spec_index")
  )

  # write one merged parquet per agg_id into tables/all/
  agg_files <- purrr::imap(agg_tbls, function(tbl, stem) write_pair(tbl, stem))

  invisible(list(
    coef_tbl           = coef_tbl,
    dummy_group_tbl    = dummy_group_tbl,
    event_time_support = support_tbl,
    run_registry       = registry_tbl,
    run_spec_index     = spec_index_tbl,
    agg_tbls           = agg_tbls,
    files              = c(files, agg_files)
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13 — Extraction internals
# ══════════════════════════════════════════════════════════════════════════════

#' Convert a fixest coeftable to a tidy tibble with metadata
coeftable_to_tbl <- function(ct, vcov_id, vcov_label, term_pattern,
                              ci_level, meta, run_id) {
  ct_tbl <- tibble::as_tibble(as.data.frame(ct), rownames = "term")
  se_col <- if ("Std. Error" %in% names(ct_tbl)) "Std. Error" else names(ct_tbl)[3]
  t_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^t value$|^t-value$")]
  p_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^Pr\\(>\\|t\\|\\)$|^p-value$")]
  z_val  <- stats::qnorm(1 - (1 - ci_level) / 2)

  out <- ct_tbl |>
    dplyr::mutate(
      subset_id            = meta$subset_id,
      outcome              = meta$outcome,
      group_id             = meta$group_id,
      group_col            = meta$group_col,
      model_id             = meta$model_id,
      estimator_type       = meta$estimator_type,
      weights_col          = meta$weights_col %||% NA_character_,
      vcov_id              = vcov_id,
      vcov_label           = vcov_label,
      formula_template     = meta$formula_template,
      estimate             = .data[["Estimate"]],
      std_error            = .data[[se_col]],
      ci_lower             = estimate - z_val * std_error,
      ci_upper             = estimate + z_val * std_error,
      term_matches_pattern = stringr::str_detect(term, term_pattern),
      term_value           = extract_first_number(term),
      run_id               = run_id
    )

  if (length(t_col) == 1) out$t_value <- ct_tbl[[t_col]] else out$t_value <- NA_real_
  if (length(p_col) == 1) out$p_value <- ct_tbl[[p_col]] else out$p_value <- NA_real_

  out |>
    dplyr::select(
      subset_id, outcome, group_id, group_col, model_id, estimator_type, weights_col,
      vcov_id, vcov_label, formula_template,
      term, estimate, std_error, t_value, p_value, ci_lower, ci_upper,
      term_matches_pattern, term_value, run_id
    )
}


#' N statistics per dummy_group from the estimation sample
make_dummy_group_summary <- function(df_est, dummy_cols, trt_col, unit_id,
                                      event_id, meta, run_id) {
  n_control <- df_est |>
    dplyr::filter(.data[[trt_col]] == 0) |>
    dplyr::summarise(n = dplyr::n_distinct(.data[[unit_id]])) |>
    dplyr::pull(n)

  purrr::map_dfr(dummy_cols, function(dc) {
    d_group <- df_est |> dplyr::filter(.data[[dc]] == 1)
    tibble::tibble(
      subset_id         = meta$subset_id,
      outcome           = meta$outcome,
      group_id          = meta$group_id,
      group_col         = meta$group_col,
      dummy_group       = dc,
      model_id          = meta$model_id,
      formula_template  = meta$formula_template,
      estimator_type    = meta$estimator_type,
      weights_col       = meta$weights_col %||% NA_character_,
      n_treated_units   = dplyr::n_distinct(d_group[[unit_id]]),
      n_treated_events  = compute_n_distinct_optional(d_group, event_id),
      n_control_units   = n_control,
      n_total_units     = dplyr::n_distinct(df_est[[unit_id]]),
      n_rows_model_data = nrow(df_est),
      dummy_group_color = get_group_color(dc, meta$group_palette),
      dummy_group_run_id = make_dummy_group_run_id(meta, dc),
      run_id            = run_id
    )
  })
}


#' Event-time support (N stats per event_time x dummy_group)
make_event_time_support_unified <- function(df_est, dummy_cols, trt_col,
                                             unit_id, event_id, cohort_var, time_var) {
  event_id   <- normalize_optional_colname(event_id)
  cohort_var <- normalize_optional_colname(cohort_var)

  if (is.na(cohort_var) || is.na(time_var)) {
    warning("cohort_var or time_var is NA; event_time_support will be empty.")
    return(tibble::tibble(dummy_group = character(), event_time = numeric(),
                          n_ptids = integer(), n_fireids = integer(), n_rows_treated = integer()))
  }

  purrr::map_dfr(dummy_cols, function(dc) {
    d_treated <- df_est |>
      dplyr::filter(.data[[dc]] == 1) |>
      dplyr::mutate(event_time = suppressWarnings(
        as.numeric(.data[[time_var]]) - as.numeric(.data[[cohort_var]])
      )) |>
      dplyr::filter(!is.na(event_time), is.finite(event_time))

    if (nrow(d_treated) == 0) return(tibble::tibble())

    if (is.na(event_id)) {
      d_treated |>
        dplyr::group_by(event_time) |>
        dplyr::summarise(n_ptids = dplyr::n_distinct(.data[[unit_id]]),
                         n_fireids = NA_integer_,
                         n_rows_treated = dplyr::n(), .groups = "drop") |>
        dplyr::mutate(dummy_group = dc)
    } else {
      d_treated |>
        dplyr::group_by(event_time) |>
        dplyr::summarise(n_ptids = dplyr::n_distinct(.data[[unit_id]]),
                         n_fireids = dplyr::n_distinct(.data[[event_id]]),
                         n_rows_treated = dplyr::n(), .groups = "drop") |>
        dplyr::mutate(dummy_group = dc)
    }
  })
}


make_estimation_registry <- function(run_id, subset_id, outcome, group_id, model_id,
                                      group_col, dummy_cols, dataset_spec,
                                      formula_template, estimator_type, term_pattern,
                                      weights_col, weights_parquet_path,
                                      data_filter, data_source, vcov_specs, agg_specs,
                                      group_args, n_rows_read, dummy_group_tbl,
                                      fit_started, fit_finished,
                                      coef_file, dummy_group_file, support_file,
                                      registry_file, run_spec_file) {
  model_formula <- build_model_formula(formula_template, outcome)

  dummy_group_tbl |>
    dplyr::select(dummy_group, n_treated_units, n_treated_events,
                  n_control_units, n_total_units, n_rows_model_data) |>
    dplyr::mutate(
      run_id               = run_id,
      subset_id            = subset_id,
      outcome              = outcome,
      group_id             = group_id,
      model_id             = model_id,
      group_col            = group_col,
      dummy_cols           = paste(dummy_cols, collapse = ","),
      unit_id              = dataset_spec$unit_id,
      event_id             = dataset_spec$event_id,
      trt_col              = dataset_spec$trt_col,
      cohort_var           = dataset_spec$cohort_var,
      time_var             = dataset_spec$time_var,
      formula_template     = formula_template,
      formula_resolved     = paste(deparse(model_formula), collapse = " "),
      estimator_type       = estimator_type,
      term_pattern         = term_pattern,
      weights_col          = weights_col          %||% NA_character_,
      weights_parquet_path = weights_parquet_path %||% NA_character_,
      data_filter          = data_filter_to_chr(data_filter),
      input_source         = paste(unlist(data_source), collapse = " | "),
      group_args_json      = as.character(serialize_object_json(group_args)),
      vcov_ids             = paste(vcov_specs$vcov_id, collapse = " | "),
      vcov_labels          = paste(vcov_specs$vcov_label, collapse = " | "),
      agg_ids              = paste(purrr::map_chr(agg_specs, "id"), collapse = " | "),
      n_rows_read          = n_rows_read,
      fit_started          = fit_started,
      fit_finished         = fit_finished,
      coef_file            = coef_file,
      dummy_group_file     = dummy_group_file,
      support_file         = support_file,
      registry_file        = registry_file,
      run_spec_file        = run_spec_file
    )
}

# SECTION 14 — Descriptive pipeline worker
# ══════════════════════════════════════════════════════════════════════════════

run_one_descriptive <- function(data_source, data_filter, dataset_spec,
                                subset_id, outcome, group_id, group_col, dummy_cols,
                                group_fun, group_args, model_id, formula_template,
                                dir_out, run_id, run_stub, group_palette = NULL,
                                treated_year_var, control_year_var, skip_existing = TRUE) {

  run_dir         <- file.path(dir_out, "descriptive", "by_run", run_stub)
  dir_ensure_local(run_dir)
  trajectory_file <- file.path(run_dir, "event_time_trajectory.parquet")
  registry_file   <- file.path(run_dir, "registry.parquet")
  run_spec_file   <- file.path(run_dir, "descriptive_spec.rds")

  if (skip_existing && all(file.exists(c(trajectory_file, registry_file, run_spec_file)))) {
    return(list(trajectory_file = trajectory_file, registry_file = registry_file,
                run_spec_file = run_spec_file, skipped_existing = TRUE))
  }

  run_started  <- Sys.time()
  df           <- load_arrow_data(data_source, data_filter)
  df_grouped   <- apply_treatment_grouping(df = df, group_fun = group_fun,
                                            group_col = group_col, group_args = group_args)

  traj_tbl <- purrr::map_dfr(dummy_cols, function(dc) {
    make_descriptive_one_dummy_group(
      df               = df_grouped,
      dummy_col        = dc,
      trt_col          = dataset_spec$trt_col,
      outcome          = outcome,
      unit_id          = dataset_spec$unit_id,
      event_id         = dataset_spec$event_id,
      time_var         = dataset_spec$time_var,
      treated_year_var = treated_year_var,
      control_year_var = control_year_var
    ) |>
      dplyr::mutate(
        subset_id        = subset_id, outcome = outcome, group_id = group_id,
        group_col        = group_col, dummy_group = dc, model_id = model_id,
        formula_template = formula_template, treated_year_var = treated_year_var,
        control_year_var = control_year_var, time_var = dataset_spec$time_var,
        dummy_group_color = get_group_color(dc, group_palette),
        run_id = run_id
      )
  })

  run_finished <- Sys.time()
  registry_tbl <- tibble::tibble(
    run_id = run_id, subset_id = subset_id, outcome = outcome, group_id = group_id,
    model_id = model_id, group_col = group_col, dummy_cols = paste(dummy_cols, collapse = ","),
    unit_id = dataset_spec$unit_id, trt_col = dataset_spec$trt_col,
    formula_template = formula_template, data_filter = data_filter_to_chr(data_filter),
    treated_year_var = treated_year_var, control_year_var = control_year_var,
    n_rows_read = nrow(df), n_dummy_groups = length(dummy_cols),
    run_started = run_started, run_finished = run_finished,
    trajectory_file = trajectory_file, registry_file = registry_file, run_spec_file = run_spec_file
  )

  run_spec <- list(run_id = run_id, run_stub = run_stub, subset_id = subset_id,
                   outcome = outcome, group_id = group_id, model_id = model_id,
                   group_col = group_col, dummy_cols = dummy_cols,
                   dataset_spec = dataset_spec, formula_template = formula_template,
                   data_filter = data_filter, data_source = data_source,
                   group_fun = group_fun, group_args = group_args,
                   treated_year_var = treated_year_var, control_year_var = control_year_var)

  arrow::write_parquet(traj_tbl,    trajectory_file)
  arrow::write_parquet(registry_tbl, registry_file)
  saveRDS(run_spec, run_spec_file)
  rm(df, df_grouped, traj_tbl, registry_tbl, run_spec); gc()

  list(trajectory_file = trajectory_file, registry_file = registry_file,
       run_spec_file = run_spec_file, skipped_existing = FALSE)
}


make_descriptive_one_dummy_group <- function(df, dummy_col, trt_col, outcome,
                                              unit_id, event_id, time_var,
                                              treated_year_var, control_year_var) {
  event_id <- normalize_optional_colname(event_id)

  d_treated <- df |>
    dplyr::filter(.data[[trt_col]] == 1, .data[[dummy_col]] == 1) |>
    dplyr::mutate(event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[treated_year_var]]),
                  series = "treated") |>
    dplyr::filter(!is.na(event_time), is.finite(event_time))

  d_control <- df |>
    dplyr::filter(.data[[trt_col]] == 0) |>
    dplyr::mutate(event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[control_year_var]]),
                  series = "control_mock") |>
    dplyr::filter(!is.na(event_time), is.finite(event_time))

  d_plot   <- dplyr::bind_rows(d_treated, d_control)
  grp_vars <- c("series", "event_time")

  if (nrow(d_plot) == 0) {
    return(tibble::tibble(series = character(), event_time = numeric(),
                          mean_outcome = numeric(), sd_outcome = numeric(),
                          n_rows = integer(), n_ptids = integer(), n_fireids = integer()))
  }

  if (is.na(event_id)) {
    d_plot |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) |>
      dplyr::summarise(mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
                       sd_outcome   = stats::sd(.data[[outcome]], na.rm = TRUE),
                       n_rows       = dplyr::n(),
                       n_ptids      = dplyr::n_distinct(.data[[unit_id]]),
                       n_fireids    = NA_integer_, .groups = "drop")
  } else {
    d_plot |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) |>
      dplyr::summarise(mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
                       sd_outcome   = stats::sd(.data[[outcome]], na.rm = TRUE),
                       n_rows       = dplyr::n(),
                       n_ptids      = dplyr::n_distinct(.data[[unit_id]]),
                       n_fireids    = dplyr::n_distinct(.data[[event_id]]), .groups = "drop")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 15 — Sun-Abraham support helpers
# ══════════════════════════════════════════════════════════════════════════════

parse_sunab_vars <- function(formula) {
  txt <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  m   <- stringr::str_match(txt, "sunab\\s*\\(\\s*([^,]+?)\\s*,\\s*([^,\\)]+?)\\s*(?:,|\\))")
  if (length(m) == 0 || all(is.na(m[1, ]))) {
    warning("Could not parse sunab() vars from formula: ", txt)
    return(list(cohort_var = NA_character_, time_var = NA_character_))
  }
  list(cohort_var = trimws(m[1, 2]), time_var = trimws(m[1, 3]))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 16 — Arrow data loading helpers
# ══════════════════════════════════════════════════════════════════════════════

open_arrow_source <- function(data_source) {
  arrow::open_dataset(unlist(data_source), format = "parquet")
}

load_arrow_data <- function(data_source, data_filter = NULL) {
  ds <- open_arrow_source(data_source)
  if (is.null(data_filter)) return(dplyr::collect(ds))

  filter_expr <- rlang::f_rhs(data_filter)
  tryCatch(
    ds |> dplyr::filter(!!filter_expr) |> dplyr::collect(),
    error = function(e) {
      message("[load_arrow_data] Arrow could not push down filter; collecting full dataset. Filter: ",
              data_filter_to_chr(data_filter))
      dplyr::collect(ds) |> dplyr::filter(!!filter_expr)
    }
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 17 — Treatment grouping helper
# ══════════════════════════════════════════════════════════════════════════════

apply_treatment_grouping <- function(df, group_fun, group_col, group_args) {
  df_out <- do.call(group_fun, c(list(df = df, group_col = group_col), group_args))
  if (!group_col %in% names(df_out)) {
    stop(glue::glue("group_fun did not produce expected column '{group_col}'."))
  }
  df_out
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 18 — Utility functions
# ══════════════════════════════════════════════════════════════════════════════

build_model_formula <- function(formula_template, outcome) {
  if (!is.character(formula_template) || length(formula_template) != 1) {
    stop("formula_template must be a length-1 character string.")
  }
  stats::as.formula(glue::glue(formula_template, outcome = outcome))
}

normalize_optional_colname <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  if (length(x) > 1) { warning("normalize_optional_colname: using first element only."); x <- x[1] }
  if (is.na(x)) return(NA_character_)
  as.character(x)
}

get_needed_columns <- function(formula, trt_col, group_col, dummy_cols,
                                unit_id, event_id = NA_character_,
                                vcov_vars = NULL, weights_col = NA_character_) {
  cols <- c(
    all.vars(formula), trt_col, group_col, dummy_cols, unit_id, event_id,
    unlist(vcov_vars, use.names = FALSE),
    if (!is.na(weights_col)) weights_col else NULL
  )
  cols |> unique() |> stats::na.omit() |> as.character()
}

check_required_columns <- function(df, required, context = "") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(if (nchar(context) > 0) paste0("[", context, "] ") else "",
         "Missing required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

compute_n_distinct_optional <- function(df, col_nm) {
  col_nm <- normalize_optional_colname(col_nm)
  if (is.na(col_nm)) return(NA_integer_)
  if (!col_nm %in% names(df)) stop(glue::glue("Column '{col_nm}' not found in data."))
  dplyr::n_distinct(stats::na.omit(df[[col_nm]]))
}

safe_path_component <- function(x) gsub("[^[:alnum:]_\\-]+", "_", as.character(x))

make_run_stub <- function(subset_id, outcome, group_id, model_id) {
  paste(safe_path_component(subset_id), safe_path_component(outcome),
        safe_path_component(group_id),  safe_path_component(model_id), sep = "__")
}

make_dummy_group_run_id <- function(meta, dummy_group) {
  as.character(glue::glue("{meta$subset_id}__{meta$outcome}__{meta$group_id}__{meta$model_id}__{dummy_group}"))
}

extract_first_number <- function(x) suppressWarnings(as.numeric(stringr::str_extract(x, "-?\\d+")))

dir_ensure_local <- function(paths) {
  purrr::walk(paths, \(p) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p)) stop("Failed to create directory: ", p)
  })
}

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

serialize_object_json <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, pretty = FALSE))
}

data_filter_to_chr <- function(data_filter = NULL) {
  if (is.null(data_filter))             return(NA_character_)
  if (inherits(data_filter, "formula")) return(paste(deparse(data_filter), collapse = " "))
  as.character(data_filter)
}

get_group_color <- function(group_value, group_palette) {
  if (!is.null(group_palette) && group_value %in% names(group_palette)) {
    return(unname(group_palette[[group_value]]))
  }
  NA_character_
}

validate_data_source <- function(src, label = "data_source") {
  if (!is.character(src) || length(src) == 0) stop(label, " must be a non-empty character vector.")
  if (length(src) == 1 && dir.exists(src)) return(invisible(TRUE))
  missing_files <- src[!file.exists(src)]
  if (length(missing_files) > 0) warning(label, ": these files do not exist yet:\n", paste(missing_files, collapse = "\n"))
  invisible(TRUE)
}

validate_dataset_spec <- function(spec) {
  required <- c("unit_id", "time_var", "trt_col", "cohort_var", "event_id")
  missing  <- setdiff(required, names(spec))
  if (length(missing) > 0) stop("dataset_spec is missing fields: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

validate_spec_table <- function(tbl, required_cols, name) {
  if (!is.data.frame(tbl)) stop(name, " must be a data frame / tibble.")
  missing <- setdiff(required_cols, names(tbl))
  if (length(missing) > 0) stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

build_equality_filter <- function(col, value) {
  rlang::new_formula(lhs = NULL, rhs = call("==", as.name(col), value), env = baseenv())
}

combine_filters <- function(filter_a, filter_b) {
  if (is.null(filter_a)) return(filter_b)
  if (is.null(filter_b)) return(filter_a)
  rlang::new_formula(lhs = NULL, rhs = call("&", rlang::f_rhs(filter_a), rlang::f_rhs(filter_b)),
                     env = rlang::f_env(filter_a))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


#### Post-hoc relative path helper

relativize_result_paths <- function(dir_results, base = here::here()) {
  all_files <- c(
    fs::dir_ls(dir_results, recurse = TRUE, glob = "*.parquet"),
    fs::dir_ls(dir_results, recurse = TRUE, glob = "*.csv")
  )
  
  path_cols <- c("coef_file", "dummy_group_file", "support_file", "registry_file",
                 "run_spec_file", "agg_eventstudy_file", "weights_file", "bal_rds_file",
                 "love_plot_file", "weight_run_spec_file", "weights_parquet_path")
  
  purrr::walk(all_files, function(f) {
    is_parquet <- fs::path_ext(f) == "parquet"
    
    tbl <- if (is_parquet) {
      arrow::read_parquet(f) |> tibble::as_tibble()
    } else {
      readr::read_csv(f, show_col_types = FALSE)
    }
    
    cols_present <- intersect(path_cols, names(tbl))
    if (length(cols_present) == 0) return(invisible(NULL))
    
    tbl <- tbl |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(cols_present),
        \(x) as.character(fs::path_rel(x, start = base))
      ))
    
    # force Arrow to release the memory-mapped file handle before writing
    gc()
    
    tmp <- fs::path(fs::path_dir(f), paste0(".tmp_", fs::path_file(f)))
    if (is_parquet) arrow::write_parquet(tbl, tmp) else readr::write_csv(tbl, tmp)
    fs::file_move(tmp, f)
  })
  
  invisible(NULL)
}





rebuild_descriptive_tables <- function(dir_out, write_csv = TRUE, recursive = TRUE) {
  dir_by_run <- file.path(dir_out, "descriptive", "by_run")
  dir_all    <- file.path(dir_out, "descriptive", "all")
  if (!dir.exists(dir_by_run)) stop("descriptive by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  find_pq <- function(pattern) {
    list.files(dir_by_run, pattern = pattern, recursive = recursive, full.names = TRUE)
  }
  
  traj_files     <- find_pq("^event_time_trajectory\\.parquet$")
  registry_files <- find_pq("^registry\\.parquet$")
  
  if (length(traj_files) == 0) {
    message("No event_time_trajectory.parquet files found under: ", dir_by_run)
  }
  
  # mirrors the read_and_dedup pattern from rebuild_estimation_tables();
  # keeps the most recently written file when the same run appears twice
  read_and_dedup <- function(files, id_cols) {
    if (length(files) == 0) return(NULL)
    purrr::map(files, \(f) {
      x <- arrow::read_parquet(f)
      x[] <- lapply(x, \(col) if (is.list(col)) as.character(col) else col)
      x$.mtime <- file.mtime(f)
      x
    }) |>
      dplyr::bind_rows() |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }
  
  write_pair <- function(tbl, stem) {
    if (is.null(tbl) || nrow(tbl) == 0) return(list(parquet = NULL, csv = NULL))
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  # trajectory dedup key: one row per run x dummy_group x series x event_time
  traj_tbl     <- read_and_dedup(traj_files,     c("run_id", "dummy_group", "series", "event_time"))
  registry_tbl <- read_and_dedup(registry_files, c("run_id"))
  
  files <- list(
    event_time_trajectory = write_pair(traj_tbl,     "event_time_trajectory"),
    run_registry          = write_pair(registry_tbl, "descriptive_run_registry")
  )
  
  if (!is.null(traj_tbl)) {
    message(glue::glue(
      "Descriptive tables merged: {nrow(traj_tbl)} trajectory rows, ",
      "{dplyr::n_distinct(traj_tbl$run_id)} runs"
    ))
  }
  
  invisible(list(
    traj_tbl     = traj_tbl,
    registry_tbl = registry_tbl,
    files        = files
  ))
}
