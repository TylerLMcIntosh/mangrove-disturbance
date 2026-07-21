
rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("code", "functions.R"))
source(here::here("code", "sunab_aggregate_vcov.R"))

install_and_load_packages(c("glue",
                            "tidyverse",
                            "sf",
                            "ggspatial",
                            "prettymapr",
                            "scico",
                            "fixest",
                            "fect",
                            "data.table",
                            "arrow"))

dir_raw <- here::here("data", "raw")
dir_derived <- here::here("data", "derived")
dir_figs <- here::here("figs")
dir_ensure(c(dir_raw,
             dir_derived,
             dir_figs))

cd_group_colors <- c(
  "tc"      = "#7F77DD",
  "d1_tc"    = "#EF9F27",
  "d2p_tc" = "red2",
  "control" = "#639922"
)


set.seed(seed = 1234)


# 19.8 million time-space points, 1,028,686; 2002-2022 (~20 years)
dats_filtered_did <- arrow::read_parquet(here::here(dir_derived, "dir_long", "did_ready_every_third_subsample.parquet"))
length(unique(dats_filtered_did$xy))

glimpse(dats_filtered_did)

ecor_counts <- dats_filtered_did |>
  count(ECOREGION)

typ_counts <- dats_filtered_did |>
  count(typ)

trt_counts <- dats_filtered_did |>
  count(treated, FirstTreat)


sample_xys <- sample(unique(dats_filtered_did$xy), 1/100 * nrow(dats_filtered_did))

dats_filtered_did <- dats_filtered_did |> filter(xy %in% sample_xys)

# Check data pre-estimation

# distribution of dummies
dats_filtered_did[, .(n = .N), by = .(cd_tc, cd_d1_tc, cd_d2p_tc)]

# FirstTreat for treated units
checks <- dats_filtered_did[cd_tc == 1 | cd_d1_tc == 1 | cd_d2p_tc == 1, 
                            .(n = .N, 
                              n_na_FirstTreat = sum(is.na(FirstTreat)),
                              n_1000 = sum(FirstTreat == 1000),
                              unique_cohorts = uniqueN(FirstTreat)),
                            by = .(cd_tc, cd_d1_tc, cd_d2p_tc)]

sum(checks$n)

hist(dats_filtered_did[FirstTreat > 1000]$FirstTreat)


# Estimate

est_sunab_dummy <- feols(
  tec_flux_v2 ~
    sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_tc +
    sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_d1_tc +
    sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_d2p_tc |
    xy + year,
  data = dats_filtered_did,
  cluster = ~ xy
)





# see dummy coefficient names
cn_dummy <- names(coef(est_sunab_dummy, agg = FALSE))
cat(cn_dummy, sep = "\n")


# Aggregate using new function
aggregation_code <- "(year::-?[0-9]+):cohort::[0-9]+:(cd.*)"

aggnew_cd_event_dummy <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = aggregation_code,
  weight_method = "model_matrix"
)



# Plot


create_plot <- function(model, version) {
  
  es_df_dummy <- as.data.frame(model) |>
    rownames_to_column("term") |>
    rename(
      estimate = Estimate,
      se = `Std. Error`
    ) |>
    mutate(
      event_time = as.integer(str_extract(term, "(?<=year::)-?[0-9]+")),
      cd = str_extract(term, "cd_.*$"),
      conf_low = estimate - 1.96 * se,
      conf_high = estimate + 1.96 * se
    )
  
  p1 <- ggplot(es_df_dummy, aes(x = event_time, y = estimate, color = cd)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(
      aes(ymin = conf_low, ymax = conf_high, fill = cd),
      alpha = 0.15,
      color = NA
    ) +
    geom_line() +
    geom_point() +
    labs(
      x = "Years since treatment",
      y = "Estimated effect on tec_flux_v2",
      color = "Group",
      fill = "Group",
      title = version,
      caption = "cd_d1_tc = at least one year of drought in 3 years before cyclone\ncd_d2p_tc = two or more years of drought in 3 years before cyclone\ncd_tc = zero years of drought pre-cyclone"
    ) +
    #scale_color_manual(values = cd_group_colors) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlim(-10, 15)
  
  return(p1)  
  
}

p1 <- create_plot(aggnew_cd_event_dummy$feols_structure, version = "Tropical cyclone impact on TEC_FLUX_V2 (1/100th of sub-sample)")


ggsave(filename = here(dir_figs, "test_plot_1perc.png"),
       plot = p1)
