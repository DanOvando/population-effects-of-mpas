#' study mpas
#'
#' study mpas calculates density ratios, BACI experiments, etc. from simulated MPAs
#'
#' @param results results from MPA experiments generated as part of zissou
#'
#' @return a list object with MPA studies appended
#' @export
#'
study_mpa <- function(results){


calc_mpa_effect <- function(outcomes) {
  mpa_effect <- outcomes %>%
    group_by(year) %>%
    mutate(mpa_size = max(percent_mpa)) %>%
    ungroup() %>%
    select(year, experiment, biomass, mpa_size) %>%
    spread(experiment, biomass) %>%
    mutate(mpa_effect = `with-mpa` / `no-mpa` - 1) # %>%
  # select(year, mpa_size, mpa_effect)
}


calc_abs_mpa_effect <- function(outcomes) {
  mpa_effect <- outcomes %>%
    group_by(year) %>%
    mutate(mpa_size = max(percent_mpa)) %>%
    ungroup() %>%
    select(year, experiment, biomass, mpa_size) %>%
    spread(experiment, biomass) %>%
    mutate(mpa_effect = `with-mpa` - `no-mpa`) # %>%
  # select(year, mpa_size, mpa_effect)
}


calc_mpa_fishery_effect <- function(outcomes) {
  mpa_effect <- outcomes %>%
    group_by(year) %>%
    mutate(mpa_size = max(percent_mpa)) %>%
    ungroup() %>%
    select(year, experiment, catch, mpa_size) %>%
    spread(experiment, catch) %>%
    mutate(mpa_effect = `with-mpa` / `no-mpa` - 1)
}


patch_weight <-  results$mpa_experiment[[1]]$raw_outcomes %>%
  filter(year == max(year), experiment == "with-mpa") %>%
  group_by(patch) %>%
  summarise(mpa = unique(mpa))

# weight fished mpa scenario densities by distance from nearest MPA
distance_grid <-
  expand.grid(from = 1:num_patches, to = 1:num_patches) %>%
  as.data.frame() %>%
  mutate(distance = purrr::map2_dbl(from, to, ~ min(
    c(abs(.x - .y),
      .x + num_patches - .y,
      num_patches - .x + .y)
  ))) %>%
  left_join(patch_weight, by = c("to" = "patch")) %>%
  group_by(from) %>%
  summarise(min_dist = min(distance[mpa == TRUE]),
            min_fished_dist = min(distance[mpa == FALSE]))

# calculate density ratios
density_ratio <- results$mpa_experiment[[1]]$raw_outcomes %>%
  left_join(distance_grid, by = c("patch" = "from")) %>%
  group_by(year) %>%
  summarise(
    fished_density = weighted.mean(biomass[eventual_mpa == FALSE &
                                             experiment == "with-mpa"],
                                   w = min_dist[eventual_mpa == FALSE &
                                                  experiment == "with-mpa"]),
    mpa_density = weighted.mean(biomass[eventual_mpa == TRUE &
                                          experiment == "with-mpa"],
                                w = min_fished_dist[eventual_mpa == TRUE &
                                                      experiment == "with-mpa"]),
    nompa_density = mean(biomass[experiment == "no-mpa"]),
    mpa_period = any(mpa[experiment == "with-mpa"] == TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    biased_density_ratio = mpa_density / fished_density,
    true_density_ratio = mpa_density / nompa_density
  )

# calculate BACI
results <- results %>%
  mutate(
    mpa_effect = map(map(mpa_experiment, "outcomes"), calc_mpa_effect),
    absolute_mpa_effect = map(map(mpa_experiment, "outcomes"), calc_abs_mpa_effect),
    fishery_effect = map(map(mpa_experiment, "outcomes"), calc_mpa_fishery_effect),
    density_ratio = list(density_ratio),
    baci = map(map(mpa_experiment, "raw_outcomes"), calculate_baci, distance = distance_grid, type = "biased")
  )
}