
tune_f <-
  function(params,
           target_depletion,
           linf,
           scientific_name,
           adult_movement,
           larval_movement,
           steepness,
           density_dependence_form,
           target_catch = 100,
           alpha = 0.5,
           sim_years = sim_years,
           burn_years = burn_years,
           num_patches = num_patches) {

    scientific_name <- gsub("(^)([[:alpha:]])", "\\1\\U\\2", scientific_name, perl=TRUE)

    fish <- create_fish(
      scientific_name = scientific_name,
      linf = linf,
      m = NA,
      r0 = exp(params[1]),
      adult_movement = adult_movement,
      larval_movement = larval_movement,
      density_dependence_form = density_dependence_form,
      steepness = steepness,
      query_fishlife = T
    )

    manager <- create_manager(year_mpa = 999)
    unfished <-
      sim_fishery(
        fish = fish,
        fleet = create_fleet(initial_effort = 0, fish = fish, fleet_model = "constant-effort"),
        manager = manager,
        num_patches = num_patches,
        sim_years = sim_years,
        burn_year = burn_years
      )

    ssb0 <- unfished %>%
      filter(year == max(year)) %>%
      summarise(ssb = sum(ssb)) %>% {
        .$ssb
      }

    fished <-
      sim_fishery(
        fish = fish,
        fleet = create_fleet(initial_effort = exp(params[2]), fish = fish, fleet_model = "constant-effort"),
        manager = manager,
        sim_years = sim_years,
        burn_year = burn_years,
        num_patches = num_patches
      )

    ssb <- fished %>%
      filter(year == max(year)) %>%
      summarise(ssb = sum(ssb)) %>% {
        .$ssb
      }

    fish_caught <- fished %>%
      filter(year == max(year)) %>%
      summarise(caught = sum(biomass_caught)) %>% {
        .$caught
      }


    depletion <- ssb / ssb0

    ss <-
      alpha * (log(depletion) - log(target_depletion)) ^ 2 + (1 - alpha) * (log(fish_caught) - log(target_catch)) ^
      2

    print(ss)
    return(ss)
  }
