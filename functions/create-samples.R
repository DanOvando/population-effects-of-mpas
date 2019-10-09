create_samples <- function(fishes,
                           divers,
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
                           cores = 1) {
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
      sigma_r = sigma_r
    )) %>%
    mutate(constant_f = map_dbl(fish, ~ .x$m *  f_v_m, f_v_m = f_v_m)) %>%
    mutate(constant_effort = num_patches / 2 * ((constant_f / fleet_q * targeted) / time_step)) %>%
    mutate(fleet = map2(
      constant_effort,
      fish,
      ~ create_fleet(
        initial_effort = .x,
        fish = .y,
        q = fleet_q
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

  # cool
  # ok now you need to think about how to set up your sampling regime.
  # The idea now is that you can send out a bunch of observers to monitor the population each year...

  go_sample <- function(pop, fish, divers, samples, cores = 1, cv = 0.1)
  {
    annual_data <- pop %>%
      nest(-year, -experiment, .key = 'pop')

    pisco_samples <-
      cross_df(
        list(
          year = unique(annual_data$year),
          experiment = 'with-mpa',
          patches = unique(pop$patch),
          diver = divers$diver,
          sample_event = 1:samples
        )
      ) %>%
      left_join(annual_data, by = c('year', 'experiment')) %>%
      left_join(divers, by = 'diver')

    doParallel::registerDoParallel(cores = cores)
    foreach::getDoParWorkers()

    if (file.exists('scuba-steve-progress.txt')) {
      file.remove('scuba-steve-progress.txt')
    }

    sampled_lengths <-
      foreach::foreach(i = 1:nrow(pisco_samples)) %dopar% {
        write(file = 'scuba-steve-progress.txt', paste0(round(i / nrow(
          pisco_samples
        ) * 100, 2), '% done'))

        sim_scuba_steve(
          pop = pisco_samples$pop[[i]],
          diver = pisco_samples$diver_stats[[i]],
          patches = pisco_samples$patches[i],
          effort = 1,
          cv = cv,
          fish = fish
        )

      }


    pisco_samples <- pisco_samples %>%
      mutate(sampled_lengths = sampled_lengths) %>%
      mutate(density = map_dbl(sampled_lengths, ~ .x$length_samples$weight %>% sum()))

  }

  simple_fish <- simple_fish %>%
    mutate(
      pisco_samples = map2(
        raw_outcomes,
        fish,
        go_sample,
        divers = divers,
        samples = samples,
        cores = cores,
        cv = 0.001
      )
    )

  return(simple_fish)

}