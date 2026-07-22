# perform estimation using dydid sun-ab pipeline from other compound disturbance work

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/03_dydid_estimation_v2.R")

required_pkgs <- c(
  "dplyr", "ggplot2", "tidyr", "readr", "purrr", "tibble", "stringr",
  "forcats", "fixest", "arrow", "glue", "here", "WeightIt", "tictoc"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr); library(ggplot2); library(tidyr);  library(readr)
library(tibble); library(purrr);  library(stringr); library(forcats)
library(fixest); library(arrow);  library(glue);    library(WeightIt)
library(tictoc)

# sunab_aggregate_vcov must be sourced before the pipeline file
source(here::here("code", "sunab_aggregate_vcov.R"))
source(here::here("code", "weight_dydid_pipeline_v7.R"))

seed <- 1234
set.seed(seed)

# Set number of cores to use in FEOLS call
fixest::setFixest_nthreads(48)

# log any unhandled errors to the status file before R exits
options(error = function() {
  err_msg <- geterrmessage()
  line    <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    " | FATAL | ", err_msg)
  write(line, file = here::here("logs/status_file.txt"), append = TRUE)
})

# ── Directory layout ──────────────────────────────────────────────────────────

version <- "v2"

dir_data    <- here::here("data", "derived")

dir_results <- here::here("results", version)

dir_long  <- file.path(dir_data, "parquet_long")

dir_ensure_local(c(dir_data, dir_long, dir_results))

#x <- arrow::open_dataset(dir_long) |> collect()


# ══════════════════════════════════════════════════════════════════════════════
# 2. Dataset spec ----
# ══════════════════════════════════════════════════════════════════════════════

dataset_spec <- make_dataset_spec(
  unit_id    = "xy",
  time_var   = "year",
  trt_col    = "treated",
  cohort_var = "FirstTreat",
  event_id   = NA_character_ #pass NA_character_ to get NAs for event_time_support instead of error if needed
)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Analysis subset specs ----
# ══════════════════════════════════════════════════════════════════════════════

typ_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_long,
  split_col         = "typ",
  id_prefix         = "typc1",
  check_all_files   = TRUE#,
  #base_filter       = ~ (treated == 1 & treated_subset %in% 1:10) | (treated == 0 & control_subset %in% 1:5)
)


all_subset_specs <- make_analysis_subset_spec(
  subset_id = "alldatac1",
  long_data_source = dir_long#,
  #data_filter = ~ (treated == 1 & treated_subset %in% 1:10) | (treated == 0 & control_subset %in% 1:5)
)



# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs ----
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(outcome = c("tec_flux_v2"))
                                


# ══════════════════════════════════════════════════════════════════════════════
# 5. Treatment group specs ----
# ══════════════════════════════════════════════════════════════════════════════

# function unneeded;cd_group set in prior script - provide empty function to avoid errors
set_cd_groups <- function(df,
                          group_col,
                          dummy_cols = c("cd_tc", "cd_d1_tc", "cd_d2p_tc"),
                          include_control = FALSE) {
  df
}


cd_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "v1",
    group_col  = "cd_group",
    dummy_cols = c("cd_tc", "cd_d1_tc", "cd_d2p_tc"),
    group_fun  = set_cd_groups
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# 6. Weighting specs ----
# ══════════════════════════════════════════════════════════════════════════════

# NONE
                                
# ══════════════════════════════════════════════════════════════════════════════
# 7. Model specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# formula_template contains the full unified dummy interaction structure.
# no_agg = TRUE gives cohort-specific coefficients required by agg_specs.
# mem.clean = TRUE is recommended for models of this size.

sunab_formula_nocovar <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_tc + ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d1_tc + ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d2p_tc",
  " | xy + year"
)



# sunab_formula_somecovar <- paste0(
#   "{outcome} ~ ",
#   "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_tc + ",
#   "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d1_tc + ",
#   "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d2p_tc +",
#   "min_t2m",
#   " | xy + year"
# )


sunab_formula_allcovar <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_tc + ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d1_tc + ",
  "sunab(FirstTreat, year, ref.p = -4, no_agg = TRUE):cd_d2p_tc + ",
  "min_t2m + mean_t2m + max_t2m + min_tp + mean_tp + max_tp",
  " | xy + year"
)

model_specs <- dplyr::bind_rows(
  
  make_model_spec(
    model_id         = "sunab_nocovar",
    formula_template = sunab_formula_nocovar,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  ),
  
  make_model_spec(
    model_id         = "sunab_allcovar",
    formula_template = sunab_formula_allcovar,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  )
  
)



# ══════════════════════════════════════════════════════════════════════════════
# 8. Vcov specs ----
# ══════════════════════════════════════════════════════════════════════════════

vcov_specs <- tibble::tibble(
  vcov_id = c(
    "conley_75km_5km"
  ),
  vcov_label = c(
    "Conley SEs: 75 km cutoff, 5 km pixel"
  ),
  vcov = list(
    fixest::vcov_conley(lat = "lat", lon = "long", cutoff = 75, pixel = 5, distance = "triangular")
  ),
  vcov_vars = list(c("lat", "long"))
)


# ══════════════════════════════════════════════════════════════════════════════
# 9. Aggregation specs ----
# ══════════════════════════════════════════════════════════════════════════════

agg_specs <- list(
  
  # Standard event study: cohort-averaged ATT per event_time x dummy_group.
  # Primary aggregation for pre-trend tests, ATT windows, and pairwise comparisons.
  make_agg_spec(
    id    = "event_study",
    agg   = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)",
    group_fun = function(x) {
      x |>
        dplyr::mutate(
          event_time  = as.integer(stringr::str_extract(group_1, "-?[0-9]+")),
          dummy_group = group_2
        ) |>
        dplyr::select(term, event_time, dummy_group)
    },
    label = "Cohort-averaged event study by dummy group"
  )
)                          





preview_run_grid <- function(subset_specs, outcome_specs, treatment_group_specs,
                             model_specs, vcov_specs, agg_specs) {
  grid <- tidyr::crossing(
    subset_specs |> dplyr::select(subset_id),
    outcome_specs,
    treatment_group_specs |> dplyr::select(group_id, dummy_cols),
    model_specs |> dplyr::select(model_id, estimator_type, weights_col)
  ) |>
    dplyr::mutate(run_id     = glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}"),
                  dummy_cols = purrr::map_chr(dummy_cols, \(dc) paste(dc, collapse = ",")))
  
  message(glue::glue("Total runs planned: {nrow(grid)}"))
  message(glue::glue("Vcov specs per run: {paste(vcov_specs$vcov_id, collapse = ', ')}"))
  message(glue::glue("Agg specs per run:  {paste(purrr::map_chr(agg_specs, 'id'), collapse = ', ')}"))
  print(grid |> dplyr::select(subset_id, outcome, group_id, model_id,
                              weights_col, dummy_cols))
}



preview <- preview_run_grid(typ_subset_specs,
                            outcome_specs,
                            cd_specs,
                            model_specs,
                            vcov_specs,
                            agg_specs)



tic('sunab estimation')
results_sunab_typ <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = typ_subset_specs,
  outcome_specs         = outcome_specs,
  treatment_group_specs = cd_specs,
  model_specs           = model_specs,
  vcov_specs            = vcov_specs,
  agg_specs             = agg_specs,
  dir_out               = dir_results,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = FALSE,
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_typ <- purrr::keep(results_sunab_typ$run_results, \(r) !is.null(r$error))
if (length(failed_typ) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_typ, \(r) message("  ", r$run_id, ": ", r$error))
}

all_estimation_tables <- rebuild_estimation_tables(dir_out = dir_results, write_csv = TRUE)

message(glue::glue(
  "Coef rows:       {nrow(all_estimation_tables$coef_tbl)}\n",
  "Agg specs found: {paste(names(all_estimation_tables$agg_tbls), collapse = ', ')}\n",
  "Unique run_ids:  {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}"
))

#relativize_result_paths(dir_results = dir_results, base = here::here())


