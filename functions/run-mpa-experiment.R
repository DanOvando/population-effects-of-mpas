run_mpa_experiment <- function(fish,
                     fleet,
                     year_mpa,
                     mpa_size,
                     sim_years,
                     num_patches,
                     burn_years,
                     sprinkler = FALSE,
                     mpa_habfactor = 1,
                     enviro = NA,
                     enviro_strength = NA,
                     rec_driver = 'stochastic',
                     random_mpa = FALSE,
                     min_size = 1,
                     simseed = 42) {


  set.seed(simseed)
  wi_mpa <-
    sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(year_mpa = year_mpa, mpa_size = mpa_size),
      sim_years = sim_years,
      num_patches = num_patches,
      burn_years = burn_years,
      enviro = enviro,
      enviro_strength = enviro_strength,
      rec_driver = rec_driver,
      sprinkler = sprinkler,
      mpa_habfactor = mpa_habfactor,
      tune_costs = FALSE,
      est_msy = FALSE,
      random_mpa = random_mpa,
      min_size = min_size
    ) %>%
    mutate(experiment = 'with-mpa')

  mpa_locations <- unique(wi_mpa$patch[wi_mpa$mpa == TRUE])

  set.seed(simseed)

  no_mpa <-
    sim_fishery(
      fish = fish,
      fleet = fleet,
      manager = create_manager(year_mpa = year_mpa, mpa_size = 0,
                               mpa_locations = mpa_locations),
      sim_years = sim_years,
      num_patches = num_patches,
      burn_years = burn_years,
      enviro = enviro,
      enviro_strength = enviro_strength,
      rec_driver = rec_driver,
      sprinkler = sprinkler,
      mpa_habfactor = mpa_habfactor,
      tune_costs = FALSE,
      est_msy = FALSE,
      random_mpa = random_mpa,
      min_size = min_size
    ) %>%
    mutate(experiment = 'no-mpa')


  outcomes <- no_mpa %>%
    bind_rows(wi_mpa) %>%
    group_by(year, experiment) %>%
    summarise(
      numbers = sum(numbers),
      biomass = sum(biomass),
      ssb = sum(ssb),
      percent_mpa = mean(mpa),
      catch = sum(biomass_caught),
      profits = sum(profits),
      effort = sum(effort)
    ) %>%
    ungroup()
  raw_outcomes <- no_mpa %>%
    bind_rows(wi_mpa)

  out <- list(outcomes = outcomes,
              raw_outcomes = raw_outcomes)

}