calculate_mpa_outcomes <- function(fish,
                     fleet,
                     year_mpa,
                     mpa_size,
                     sim_years,
                     num_patches,
                     burn_years,
                     run_id,
                     num_runs) {
  write(
    paste('run', run_id, 'out of', num_runs),
    file = 'simprog.txt',
    append = T,
    ncolumns = 1
  )

  no_mpa <-
    sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(year_mpa = year_mpa, mpa_size = 0),
      sim_years = sim_years,
      num_patches = num_patches,
      burn_years = burn_years
    ) %>%
    mutate(experiment = 'no-mpa')

  wi_mpa <-
    sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(year_mpa = year_mpa, mpa_size = mpa_size),
      sim_years = sim_years,
      num_patches = num_patches,
      burn_years = burn_years
    ) %>%
    mutate(experiment = 'with-mpa')

  outcomes <- no_mpa %>%
    bind_rows(wi_mpa) %>%
    group_by(year, experiment) %>%
    summarise(
      ssb = sum(ssb),
      percent_mpa = mean(mpa),
      catch = sum(biomass_caught),
      profits = sum(profits),
      effort = sum(effort)
    )



}
