test_performance <-
  function(fishes,
           year_mpa,
           min_year = 75,
           max_year = 100,
           time_step = 1,
           sigma_obs = 0,
           samps_per_patch = 1) {
    
    
    simple_data <- fishes %>%
      select(loo, k, lm, m, targeted, classcode, commonname,raw_outcomes) %>%
      unnest(cols = raw_outcomes) %>%
      group_by(year, patch, experiment, targeted) %>% 
      summarise(biomass_density = list(rlnorm(samps_per_patch,log(sum(biomass) / 100 + 1e-3), sigma_obs))) %>% 
      unnest(cols = biomass_density) %>% 
      ungroup() %>% 
      filter(year > min_year, year < max_year,
             experiment == "with-mpa") %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year)) %>%
      mutate(
        factor_year = as.factor(year),
        log_density = log(biomass_density),
        logical_targeted = targeted > 0,
        post_mpa = year >= year_mpa
      )
    
    years_protected <- unique(simple_data$year) - year_mpa

    bins <-
      c(seq(min(years_protected),-1, by = 5), seq(0, max(years_protected), by = 5))

    year_block <- data_frame(year = unique(simple_data$year)) %>%
      mutate(years_protected = year - year_mpa) %>%
      mutate(protected_block = cut(
        years_protected,
        breaks = bins,
        include.lowest = T
      )) %>%
      select(-years_protected)

    simple_data <- simple_data %>%
      left_join(year_block, by = "year")

    true_effect <- fishes %>%
      select(targeted, classcode, commonname, raw_outcomes) %>%
      unnest(cols = raw_outcomes) %>%
      group_by(year, patch, experiment, targeted) %>%
      summarise(b = sum(biomass)) %>%
      group_by(year, experiment, targeted) %>%
      summarise(mbd = mean(b)) %>%
      ungroup() %>%
      pivot_wider(names_from = experiment,
                  values_from = mbd) %>%
      mutate(mpa_effect = `with-mpa` / `no-mpa` - 1) %>%
      filter(year > min_year, targeted == 1, year < max_year) %>%
      mutate(post_mpa = year > year_mpa) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year)) %>% 
      filter(targeted == 1)
# 
      env <- new.env(parent = parent.frame())
      
      env$simple_data <- simple_data
    
        # did_model <- with(env, {stan_glmer(biomass_density ~ targeted * factor_year + (1|patch),
        #                     data = simple_data,
        #                     iter = 5000,
        #                     chains = 4,
        #                     cores = 4,
        #                     prior = normal(0,2.5, autoscale = TRUE),
        #                     family = Gamma(link = "log"))}
        # )

        
        did_model <- with(env, {stan_glm(biomass_density ~ targeted * factor_year,
                                           data = simple_data,
                                           iter = 3000,
                                           chains = 4,
                                           cores = 4,
                                           prior = normal(0,2.5, autoscale = TRUE),
                                           family = Gamma(link = "log"))}
        )
        # browser()
        # 
        # did_model <- with(env, {lme4::glmer(biomass_density ~ targeted * factor_year + (1|patch),
        #                                    data = simple_data,
        #                                    family = Gamma(link = "log"))}
        # )
        # browser()
    
    # did_model <- stan_glm(log(biomass_density) ~ targeted * factor_year,
    #                         data = simple_data,
    #                         chains = 4,
    #                         cores = 4)


    did_values <-  tidybayes::tidy_draws(did_model) %>%
      select(contains("."), contains("targeted:factor_year")) %>%
      pivot_longer(
        contains("factor_year"),
        names_to = "year",
        values_to = "value",
        names_prefix = "targeted:factor_year",
        names_ptypes = list(year = integer())
      ) %>% 
      mutate(value = exp(value) - 1) %>% 
      left_join(true_effect %>% select(year, mpa_effect), by = "year") %>% 
      mutate(years_protected = year - year_mpa)
    
    did_plot <- did_values %>%
      ggplot() +
      tidybayes::stat_halfeye(aes(years_protected, value))+
      geom_line(
        aes(years_protected, mpa_effect))

    out <- list(
      did_values = did_values
    )
      

   

    return(out)
  }
