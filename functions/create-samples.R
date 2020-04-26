create_samples <- function(fishes,
                           f_v_m = 1.5,
                           fleet_q = 0.001,
                           num_patches = 2,
                           year_mpa = 50,
                           sim_years = 100,
                           burn_years = 50,
                           mpa_size = 0.5,
                           samples = 2,
                           time_step = 1,
                           enviro_strength = 1,
                           rec_driver = 'stochastic',
                           sigma_r = 0,
                           cores = 1,
                           cv = .1) {
  simple_fish <- fishes %>%
    mutate(fish = pmap(
      list(
        linf = loo,
        vbk = k,
        m = m,
        length_50_mature = lm
      ),
      create_fish,
      time_step = time_step,
      sigma_r = sigma_r,
      larval_movement = round(num_patches  *.25),
      adult_movement = round(num_patches * .05),
      steepness = 0.8
    )) %>%
    mutate(constant_f = map_dbl(fish, ~ .x$m *  f_v_m, f_v_m = f_v_m)) %>%
    mutate(constant_effort =  ((constant_f / fleet_q) / time_step) * num_patches * targeted) %>%
    mutate(fleet = map2(
      constant_effort,
      fish,
      ~ create_fleet(
        initial_effort = .x,
        fish = .y,
        q = fleet_q,
        length_50_sel = .y$length_50_mature * .9
      ),
      fleet_q = fleet_q
    ))
  
  simple_fish <- simple_fish %>%
    mutate(
      mpa_experiment = pmap(
        list(
          fish = fish,
          fleet = fleet,
          enviro = enviro
        ),
        run_mpa_experiment,
        year_mpa = year_mpa,
        sim_years = sim_years,
        num_patches = num_patches,
        burn_years =  burn_years,
        mpa_size = mpa_size,
        enviro_strength = enviro_strength,
        rec_driver = rec_driver
      )
    )
  
  simple_fish <- simple_fish %>%
    mutate(net_outcomes = map(mpa_experiment, 'outcomes')) %>%
    mutate(raw_outcomes = map(mpa_experiment, 'raw_outcomes')) %>%
    mutate(mpa_plot = map(net_outcomes, ~ .x %>%
                            ungroup() %>%
                            ggplot(
                              aes(year, ssb, color = experiment, linetype = experiment)
                            ) +
                            # geom_point(aes(size = percent_mpa)) +
                            geom_line()))

  calc_mpa_effect <- function(data) {
    data %>%
      group_by(experiment, year) %>%
      summarise(total_biomass = sum(biomass)) %>%
      mutate(log_total_biomass = log(total_biomass)) %>%
      select(-total_biomass) %>%
      spread(experiment, log_total_biomass) %>%
      mutate(mpa_effect = `with-mpa` - `no-mpa`)


  }

  simple_fish <- simple_fish %>%
    mutate(mpa_effect = map(raw_outcomes, calc_mpa_effect))

  return(simple_fish)

}