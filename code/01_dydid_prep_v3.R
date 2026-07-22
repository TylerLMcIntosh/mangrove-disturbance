
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
dir_long <- here::here(dir_derived, "parquet_long_full")
dir_ensure(c(dir_raw,
             dir_derived,
             dir_figs))

cd_group_colors <- c(
  "tc"      = "#7F77DD",
  "d1_tc"    = "#EF9F27",
  "d2p_tc" = "red2",
  "control" = "#639922"
)

figs = FALSE


set.seed(seed = 1234)


parquet_files <- list.files(
  here(dir_raw, "drought_and_TC"),
  pattern    = "\\.parquet$",
  full.names = TRUE
)


dats <- arrow::open_dataset(parquet_files) |>
  collect() |>
  data.table::as.data.table()

#add lat/long to data for conley vcov
dats[
  ,
  c("long", "lat") := tstrsplit(
    xy,
    "_",
    fixed = TRUE,
    type.convert = TRUE
  )
]
  



# -- New 4-year rolling window -------------------------------------------------

setorder(dats, xy, year) # set order of rows for frollsum

dats[, drought_rc4 := frollsum(drought_c, n = 4, align = "right", fill = NA), by = xy]

# Do NOT shift by 1, maintain inclustion of T itself
## frollsum at position t includes t itself, so shift forward by 1
#dats[, drought_rc4 := shift(drought_rc4, n = 1, type = "lag"), by = xy]


dats <- dats[!is.na(drought_rc4)] # filter out data without rc3


# ── TC and drought counts ─────────────────────────────────────────────────────
tc_counts <- dats[, .(n_tc = sum(TC_c)), by = xy]

tc_counts_counts <- tc_counts |>
  count(n_tc) |>
  mutate(perc_dats = n / sum(n) * 100)

if(figs) {
  ggplot(tc_counts) +
    geom_histogram(aes(x = n_tc)) +
    labs(title = "Number of TCs at unit (>=2 removed)")
  ggsave(filename = here(dir_figs, "n_tc.png"))
}

d_counts <- dats[, .(n_d = sum(drought_c)), by = xy]
hist(d_counts$n_d)

# ── Filter to units with at most one TC event ─────────────────────────────────
no_multi_tc <- tc_counts[n_tc <= 1, .(xy)]

dats_filtered <- dats[xy %in% no_multi_tc$xy]


# ── Build short-form data ─────────────────────────────────────────────────────
xy_single_tc <- tc_counts[n_tc == 1, unique(xy)]
xy_no_tc <- tc_counts[n_tc == 0, unique(xy)]
xy_multi_tc <- tc_counts[n_tc > 1, unique(xy)]

length(xy_multi_tc) / length(tc_counts$n_tc)
length(xy_no_tc) / length(tc_counts$n_tc)
length(xy_single_tc) / length(tc_counts$n_tc)


treated_sample <- dats_filtered[TC_c == 1]


dats_filtered_short <- data.table::rbindlist(list(
  dats_filtered[TC_c == 1],
  unique(dats_filtered[TC_c == 0 & !xy %in% xy_single_tc], by = "xy")
))[, drought_years_pre_tc := drought_rc4]

unique(dats_filtered_short$drought_years_pre_tc)


# ── Group assignment ──────────────────────────────────────────────────────────
dats_filtered_short[, cd_group := data.table::fcase(
  TC_c == 1 & drought_years_pre_tc >= 2, "d2p_tc",
  TC_c == 1 & drought_years_pre_tc == 1, "d1_tc",
  TC_c == 1 & drought_years_pre_tc == 0, "tc",
  TC_c == 0,                              "control",
  default = NA_character_
)]

dats_filtered_short[, treated := data.table::fcase(
  TC_c == 1, 1L,
  TC_c == 0, 0L,
  default   = NA_integer_
)]

dats_filtered_short[, cd_group := relevel(factor(cd_group), ref = "tc")]

grouping_summary <- dats_filtered_short |>
  group_by(cd_group) |>
  summarize(n = n())


# ── Binary dummy columns ──────────────────────────────────────────────────────
level_names <- c("tc", "d1_tc", "d2p_tc")
dummy_cols  <- c("cd_tc", "cd_d1_tc", "cd_d2p_tc")

for (i in seq_along(dummy_cols)) {
  lv <- level_names[i]
  dats_filtered_short[, (dummy_cols[i]) := as.integer(
    !is.na(cd_group) & as.character(cd_group) == lv
  )]
  dats_filtered_short[is.na(get(dummy_cols[i])), (dummy_cols[i]) := 0L]
}

#──────────────────────────────────────────────────────
## Summary stats and descriptive figs ----

if(figs) {
  tc_year_p <- ggplot(treated_sample) +
    geom_bar(aes(x = year)) +
    theme_minimal() +
    labs(title ="Year of TC")
  tc_year_p
  ggsave(tc_year_p, filename = here(dir_figs, "year_of_tc.png"))
  
  
  # Show drought years
  d_tc_p <- ggplot(treated_sample) +
    geom_bar(aes(x = drought_rc4)) +
    theme_minimal() +
    labs(title ="Drought years prior to TC (single-TC only)")
  d_tc_p
  ggsave(d_tc_p, filename = here(dir_figs, "n_d_prior_to_tc.png"))
  
  grouping_summary |> filter(cd_group == "tc") |> pull(n) / sum(grouping_summary |> filter(cd_group == "tc" | cd_group == "d_tc") |> pull(n))
  
  prop <- ggplot(grouping_summary) +
    geom_bar(aes(x = "", y = n, fill = cd_group), stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = cd_group_colors) +
    labs(title = "Proportions of final dataset")
  
  ggsave(prop, filename = here(dir_figs, "proportions_fig.png"))
  
  
  group_by_year <- ggplot(dats_filtered_short |> filter(cd_group != "control")) +
    geom_bar(aes(x = year, fill = cd_group)) +
    scale_fill_manual(values =cd_group_colors) +
    theme_minimal() +
    labs(title ="Year of TC, split by treatment type")
  group_by_year
  ggsave(group_by_year, filename = here(dir_figs, "year_of_tc_by_trt.png"))
  
}






# MAP DATASET
dats_filtered_short[, c("lon", "lat") := tstrsplit(xy, "_", fixed = TRUE, type.convert = TRUE)]

bb <- c(
  xmin = min(dats_filtered_short$lon),
  xmax = max(dats_filtered_short$lon),
  ymin = min(dats_filtered_short$lat),
  ymax = max(dats_filtered_short$lat)
)

# convert to sf only after parsing — avoids doing it on the full dataset
pts_sf <- sf::st_as_sf(dats_filtered_short, coords = c("lon", "lat"), crs = 4326)

if(figs) {
  d_map <- ggplot() +
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 5) +
    ggspatial::layer_spatial(
      data = pts_sf,
      mapping = aes(color = cd_group),
      size = 0.2,
      alpha = 0.6
    ) +
    scale_color_manual(values = cd_group_colors) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      crs = sf::st_crs(4326),
      expand = FALSE
    ) +
    theme_minimal() +
    labs(title = "Drought group distribution", color = "cd_group")
  
  ggsave(
    filename = here::here(dir_figs, "drought_cd_group.png"),
    plot     = d_map,
    units    = "px",
    width    = 6000,
    height   = 4000
  )
}




# DID PREP
# Test new get.cohort() - re-implementation of fect() version for data.table
get_cohort_dt <- function(dt, D, index) {
  unit_var <- index[1]
  time_var <- index[2]
  
  cohort_dt <- dt[
    get(D) == 1 & !is.na(get(time_var)),
    .(FirstTreat = min(get(time_var))),
    by = unit_var
  ]
  
  dt <- merge(
    dt,
    cohort_dt,
    by = unit_var,
    all.x = TRUE,
    sort = FALSE
  )
  
  dt[, rel_time := get(time_var) - FirstTreat]
  
  dt[]
}

arrow::write_parquet(
  dats_filtered,
  sink = here::here(
    dir_derived,
    "test5.parquet"
  ),
  compression = "snappy"
)


dats_filtered_did <- data.table::as.data.table(dats_filtered) |>
  merge(
    data.table::as.data.table(dats_filtered_short)[, .(xy, drought_years_pre_tc, cd_group, cd_tc, cd_d1_tc, cd_d2p_tc, treated)],
    by = "xy",
    all.x = TRUE
  )


dats_filtered_did <- get_cohort_dt(dats_filtered_did, D = "TC_c", index = c("xy", "year"))

dats_filtered_did[is.na(FirstTreat), FirstTreat := 1000]

# set random control group set
set.seed(1234)

control_subsets <- unique(
  dats_filtered_did[treated == 0, .(xy)]
)[
  , control_subset := sample(rep(1:10, length.out = .N))
]

treated_subsets <- unique(
  dats_filtered_did[treated == 1, .(xy)]
)[
  , treated_subset := sample(rep(1:10, length.out = .N))
]

dats_filtered_did[
  control_subsets,
  on = "xy",
  control_subset := i.control_subset
]

dats_filtered_did[
  treated_subsets,
  on = "xy",
  treated_subset := i.treated_subset
]


# get counts
typ_counts <- dats_filtered_did |>
  count(typ)

typ_counts_cd <- dats_filtered_did |>
  count(typ, cd_group) |>
  mutate(n_distinct_ids = n / 18)

write_csv(typ_counts, here::here(dir_long, "all_samples_typ_counts.csv"))
write_csv(typ_counts_cd, here::here(dir_long, "all_samples_typ_counts_cd.csv"))



arrow::write_parquet(as.data.frame(dats_filtered_did),
                     sink = here::here(dir_long, "did_ready_all_samples.parquet"))

message(
  "Finished dats_filtered_did write: ",
  round(file.size(here::here(dir_long, "did_ready_all_samples.parquet")) / 1024^2, 1),
  " MB"
)



