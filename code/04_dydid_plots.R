# 08_dydid_plots_v7.R
# Plotting — skeleton
# -------------------------------------------------------
# Reads from:
#   tables/all/agg_eventstudy.parquet          — aggregated event studies
#   tables/inference/all/*.parquet             — ATT windows, comparisons, etc.
#   tables/all/dummy_group_summary.parquet     — N by dummy_group
#   tables/all/event_time_support.parquet      — N by event_time x dummy_group
# -------------------------------------------------------

rm(list = ls())

library(here)
here::i_am("code/04_dydid_plots.R")

library(dplyr); library(ggplot2); library(readr); library(purrr)
library(tibble); library(stringr); library(arrow); library(glue)
library(patchwork)

source(here::here("code", "weight_dydid_pipeline_v7.R"))
source(here("code", "functions.R"))

version     <- "v2"
dir_results <- here::here("results", version)
dir_figs    <- here::here("figs",    version)

dir_ensure(dir_figs)



# # Make sure all tables have been built if downloaded incomplete set from exo
# all_estimation_tables <- rebuild_estimation_tables(dir_out = dir_results, write_csv = TRUE)
# all_descriptive_tables <- rebuild_descriptive_tables(dir_out = dir_results, write_csv = TRUE)
#
# message(glue::glue(
#   "Coef rows:       {nrow(all_estimation_tables$coef_tbl)}\n",
#   "Agg specs found: {paste(names(all_estimation_tables$agg_tbls), collapse = ', ')}\n",
#   "Unique run_ids:  {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}"
# ))
#
# relativize_result_paths(dir_results = dir_results, base = here::here())

# ══════════════════════════════════════════════════════════════════════════════
# Load merged tables ----
# ══════════════════════════════════════════════════════════════════════════════

dir_all      <- file.path(dir_results, "tables", "all")

agg_es_tbl  <- arrow::read_parquet(file.path(dir_all, "agg_event_study.parquet"))
support <- arrow::read_parquet(file.path(dir_all, "event_time_support.parquet"))


support_summary <- support


# ══════════════════════════════════════════════════════════════════════════════
# Plot helpers ----
# ══════════════════════════════════════════════════════════════════════════════

# Dummy group labels for axes
dummy_group_labels <- c(
  "cd_tc"   = "Tropical cyclone",
  "cd_d1_tc"  = "Tropical cyclone + 1 year of drought",
  "cd_d2p_tc"  = "Tropical cyclone + 2 years of drought"
)

# Default palette (override via group_palette from Script 6)
dummy_group_palette <- c(
  "cd_tc"   = "#0072B2",
  "cd_d1_tc"  = "goldenrod2",
  "cd_d2p_tc"  = "#7A0177"
)



subset_read_codes <- tibble(
  subset_id = c("nfg_aspen_birch",
                "nfg_california_mixed_conifer",
                "nfg_douglas_fir",
                "nfg_fir_spruce_mountain_hemlock",
                "nfg_lodgepole_pine",
                "nfg_pinyon_juniper",
                "nfg_ponderosa_pine",
                "nfg_western_oak"),
  subset_label = c("Aspen/birch",
                   "California mixed conifer",
                   "Douglas-fir",
                   "Fir/spruce/mountain hemlock",
                   "Lodgepole pine",
                   "Pinyon juniper",
                   "Ponderosa pine",
                   "Western oak")
)


# ══════════════════════════════════════════════════════════════════════════════
# Event study plot ----
# ══════════════════════════════════════════════════════════════════════════════

plot_event_study <- function(agg_es,
                             subset_id_filter    = NULL,
                             outcome_filter,
                             group_id_filter,
                             model_id_filter,
                             vcov_id_filter      = NULL,
                             event_time_range    = c(-15, 20),
                             ref_period          = -6,
                             palette             = dummy_group_palette,
                             group_labels        = dummy_group_labels,
                             title               = NULL,
                             facet_by_dummy      = TRUE,
                             facet               = NULL,
                             free_y              = FALSE,
                             support             = NULL,
                             min_n_events        = NULL,
                             min_n_points        = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- agg_es |>
    dplyr::filter(
      if (is.null(subset_id_filter)) TRUE else subset_id %in% subset_id_filter,
      outcome    %in% outcome_filter,
      group_id   %in% group_id_filter,
      model_id   %in% model_id_filter,
      event_time >= event_time_range[1],
      event_time <= event_time_range[2]
    )
  
  if (!is.null(vcov_id_filter)) {
    d <- dplyr::filter(d, vcov_id %in% vcov_id_filter)
  }
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check subset/outcome/group/model filters.")
    return(NULL)
  }
  
  d <- d |>
    dplyr::mutate(
      dummy_group_label = factor(
        dplyr::recode(dummy_group, !!!group_labels),
        levels = unname(group_labels)
      )
    )
  
  # ── Support threshold filtering ───────────────────────────────────────────
  d_plot <- d
  
  if (!is.null(support) && (!is.null(min_n_events) || !is.null(min_n_points))) {
    support_sub <- support |>
      dplyr::filter(
        if (is.null(subset_id_filter)) TRUE else subset_id %in% subset_id_filter,
        outcome    %in% outcome_filter,
        group_id   %in% group_id_filter,
        model_id   %in% model_id_filter,
        event_time >= event_time_range[1],
        event_time <= event_time_range[2]
      ) |>
      dplyr::select(
        subset_id, outcome, group_id, model_id,
        event_time, dummy_group, n_fireids, n_ptids
      )
    
    d_plot <- d |>
      dplyr::left_join(
        support_sub,
        by = c(
          "subset_id", "outcome", "group_id", "model_id",
          "event_time", "dummy_group"
        )
      )
    
    if (!is.null(min_n_events)) {
      d_plot <- dplyr::filter(d_plot, !is.na(n_fireids), n_fireids >= min_n_events)
    }
    
    if (!is.null(min_n_points)) {
      d_plot <- dplyr::filter(d_plot, !is.na(n_ptids), n_ptids >= min_n_points)
    }
  }
  
  if (nrow(d_plot) == 0) {
    warning("No rows remain after support filtering. Try relaxing min_n_events or min_n_points.")
    return(NULL)
  }
  
  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    d_plot,
    ggplot2::aes(
      x = event_time,
      y = estimate,
      color = dummy_group,
      fill = dummy_group,
      group = dummy_group
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = ref_period, linetype = "dotted", color = "grey60") +
    ggplot2::annotate(
      "rect",
      xmin = min(d_plot$event_time),
      xmax = -0.5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.03,
      fill = "grey50"
    ) +
    ggplot2::scale_color_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::scale_fill_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::labs(
      x     = "Event time (years since TC)",
      y     = "Effect on tec_flux_v2",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # ── Faceting ──────────────────────────────────────────────────────────────
  facet_vars <- c(
    if (facet_by_dummy) "dummy_group_label" else NULL,
    facet
  )
  
  if (length(facet_vars) > 0) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(facet_vars)),
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  
  p
}


# USE FUNCTIONS ----

subsets <- agg_es_tbl |>
  pull(subset_id) |>
  unique()


plot_combinations <- agg_es_tbl |>
  dplyr::filter(
    subset_id %in% subsets,
    !is.na(model_id),
    !is.na(group_id),
    !is.na(outcome),
    !is.na(vcov_id),
    !is.na(subset_id)
  ) |>
  dplyr::distinct(
    model_id,
    group_id,
    outcome,
    vcov_id,
    subset_id
  )

# summary_plot_combinations <- agg_es_tbl |>
#   select(-subset_id) |>
#   dplyr::distinct(
#     model_id,
#     group_id,
#     outcome,
#     vcov_id
#   )

# PER-GROUP PLOTS ----

#for(i in 1:20) {
for(i in seq_len(nrow(plot_combinations))) {

  pars <- plot_combinations[i, ]

  dir_figs_subset <- here(
    dir_figs,
    pars$model_id,
    pars$group_id,
    pars$outcome,
    pars$vcov_id
  )

  dir_ensure(dir_figs_subset)

  # Normal event study
  p <- plot_event_study(
    agg_es           = agg_es_tbl,
    subset_id_filter = pars$subset_id,
    outcome_filter   = pars$outcome,
    group_id_filter  = pars$group_id,
    model_id_filter  = pars$model_id,
    vcov_id_filter   = pars$vcov_id,
    support          = support,
    facet_by_dummy   = FALSE,,
    ref_period = -3,
    #min_n_events     = 3,
    min_n_points     = 50,
    title = glue::glue(
      "{pars$model_id}: {pars$group_id}\n",
      "{pars$subset_id}\n",
      "{pars$outcome}-{pars$vcov_id}"
    )
  ) +
    coord_cartesian(xlim = c(-15, 15))

  ggsave(
    plot = p,
    filename = here(
      dir_figs_subset,
      glue::glue("event_study_{pars$subset_id}_.png")
    )
  )

}



# 
# 
# 
# ## SUMMARY PLOTS ----
# 
# #for(i in 1:5) {
# for(i in seq_len(nrow(summary_plot_combinations))) {
#   
#   pars <- summary_plot_combinations[i, ]
#   
#   dir_fig_summary <- here(
#     dir_figs,
#     "summary_figs"
#     # pars$model_id,
#     # pars$group_id,
#     # pars$outcome,
#     # pars$vcov_id
#   )
#   
#   dir_ensure(dir_fig_summary)
#   
#   # Raw estimates
#   
#   dir_fig_summary_one_window <- here(dir_fig_summary, "one_window")
#   dir_ensure(dir_fig_summary_one_window)
#   
#   p_one_window <- plot_att_windows(
#     att_windows = att_windows |>
#       filter(grepl("nfg", subset_id)) |>
#       filter(dummy_group %in% c("cd_f", "cd_bf", "cd_df", "cd_bdf")),
#     outcome_filter  = pars$outcome,
#     group_id_filter = pars$group_id,
#     model_id_filter = pars$model_id,
#     vcov_id_filter  = pars$vcov_id,
#     window_id_filter = "years_2_16",
#     dodge_width     = 0.6,
#     sig_threshold = 0.95,
#     free_x = FALSE,
#     title = glue::glue(
#       "2-16 year post-fire aggregated estimates"
#     )#,
#     #xlim = c(-15, 20)
#   )
#   ggsave(plot = p_one_window,
#          here(dir_fig_summary_one_window, glue("one_window_summary_nfgs_{pars$model_id}_{pars$group_id}_{pars$outcome}_{pars$vcov_id}.png")),
#          units = "px",
#          width = 2000,
#          height = 1500)
#   
#   
#   dir_fig_summary_multi_window <- here(dir_fig_summary, "multi_window")
#   dir_ensure(dir_fig_summary_multi_window)
#   
#   p_multi_window <- plot_att_windows(
#     att_windows = att_windows |>
#       filter(grepl("nfg", subset_id)) |>
#       filter(dummy_group %in% c("cd_f", "cd_bf", "cd_df", "cd_bdf")),
#     outcome_filter  = pars$outcome,
#     group_id_filter = pars$group_id,
#     model_id_filter = pars$model_id,
#     vcov_id_filter  = pars$vcov_id,
#     window_id_filter = c("years_2_6", "years_7_11", "years_12_16"),
#     dodge_width     = 0.6,
#     sig_threshold = 0.95,
#     free_x = FALSE,
#     title = glue::glue(
#       "Post-fire aggregated estimates"
#     )#,
#     #xlim = c(-15, 20)
#   )
#   ggsave(plot = p_multi_window,
#          here(dir_fig_summary_multi_window, glue("multi_window_summary_nfgs_{pars$model_id}_{pars$group_id}_{pars$outcome}_{pars$vcov_id}.png")),
#          units = "px",
#          width = 3000,
#          height = 1500)
#   
#   # Comparison to fire
#   dir_fig_summary_fire_compare <- here(dir_fig_summary, "fire_compare")
#   dir_ensure(dir_fig_summary_fire_compare)
#   
#   p_vs_f <- plot_att_comparisons_vs_f_windows(
#     att_comps       = att_comps |>
#       filter(grepl("nfg", subset_id)),
#     outcome_filter  = pars$outcome,
#     group_id_filter = pars$group_id,
#     model_id_filter = pars$model_id,
#     vcov_id_filter  = pars$vcov_id,
#     dodge_width     = 0.6,
#     title = glue::glue(
#       "Formal comparison to fire-only group"
#     ),
#     xlim = c(-20, 15)
#   )
#   ggsave(plot = p_vs_f,
#          here(dir_fig_summary_fire_compare, glue("fire_comparison_summary_nfgs_{pars$model_id}_{pars$group_id}_{pars$outcome}_{pars$vcov_id}.png")),
#          units = "px",
#          width = 2500,
#          height = 1750)
#   
#   
#   # Normalized comparison to fire
#   
#   dir_fig_summary_fire_compare_normalized <- here(dir_fig_summary, "fire_compare_normalized")
#   dir_ensure(dir_fig_summary_fire_compare_normalized)
#   
#   p_vs_f <- plot_att_comparisons_vs_f_windows(
#     att_comps       = att_comps_normalized |>
#       filter(grepl("nfg", subset_id)),
#     outcome_filter  = pars$outcome,
#     group_id_filter = pars$group_id,
#     model_id_filter = pars$model_id,
#     vcov_id_filter  = pars$vcov_id,
#     dodge_width     = 0.6,
#     title = glue::glue(
#       "Normalized comparison to fire-only group"
#     ),
#     xlim = c(-20, 15)
#   )
#   ggsave(plot = p_vs_f,
#          here(dir_fig_summary_fire_compare_normalized, glue("fire_comparison_normalized_summary_nfgs_{pars$model_id}_{pars$group_id}_{pars$outcome}_{pars$vcov_id}.png")),
#          units = "px",
#          width = 2500,
#          height = 1750)
#   
#   
#   
#   # Multi-att comparison
#   
#   dir_fig_summary_fire_compare_multi <- here(dir_fig_summary, "fire_compare_multi")
#   dir_ensure(dir_fig_summary_fire_compare_multi)
#   
#   if(grepl("ext", pars$model_id)) {
#     
#     p_vs_f <- plot_att_comparisons_vs_f_windows(
#       att_comps       = att_comps_multiple |>
#         filter(grepl("nfg", subset_id)),
#       outcome_filter  = pars$outcome,
#       group_id_filter = pars$group_id,
#       model_id_filter = pars$model_id,
#       vcov_id_filter  = pars$vcov_id,
#       dodge_width     = 0.6,
#       title = glue::glue(
#         "Multi-comparison to fire-only group"
#       ),
#       xlim = c(-20, 15)
#     )
#     ggsave(plot = p_vs_f,
#            here(dir_fig_summary_fire_compare_multi, glue("fire_comparison_multi_summary_nfgs_{pars$model_id}_{pars$group_id}_{pars$outcome}_{pars$vcov_id}.png")),
#            units = "px",
#            width = 2500,
#            height = 1750)
#     
#   }
#   
#   
#   
# }



