test_performance <-
  function(fishes,
           year_mpa,
           min_year = 75,
           max_year = 100,
           time_step = 1,
           sigma_obs = 0) {
    
    
    simple_data <- fishes %>%
      select(loo, k, lm, m, targeted, classcode, commonname,raw_outcomes) %>%
      unnest(cols = raw_outcomes) %>%
      group_by(year, patch, experiment, targeted) %>% 
      summarise(biomass_density = rlnorm(1,log(sum(biomass) + 1e-3), sigma_obs)) %>% 
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
             year = floor(year))
# 
    did_model <- stan_glmer(biomass_density ~ targeted * protected_block + (1|patch),
                            data = simple_data,
                            chains = 4,
                            cores = 4,
                            family = Gamma(link = "log"))
    
    # did_model <- stan_glm(log(biomass_density) ~ targeted * protected_block,
    #                         data = simple_data,
    #                         chains = 4,
    #                         cores = 4)


    did_values <-  tidybayes::tidy_draws(did_model) %>%
      select(contains("."), contains("targeted:")) %>%
      pivot_longer(
        contains("targeted"),
        names_to = "protected_block",
        values_to = "value",
        names_prefix = "targeted:protected_block"
      ) %>% 
      mutate(value = exp(value) - 1)
    
  
    get_range <- function(bin) {
      bin_range <- str_extract(bin, pattern = '(?<=\\().*(?=])')

      mean(str_split(bin_range, ',', simplify = T) %>% as.numeric())

    }
    
    did_values <- did_values %>% 
      mutate(mid_bin = map_dbl(protected_block, get_range))

    mean_effect <- true_effect %>%
      group_by(year) %>%
      summarise(mean_effect = mean(mpa_effect))

    did_plot <- did_values %>% 
      ggplot() +
      tidybayes::stat_halfeye(aes(mid_bin, value))+
      geom_line(
        data = true_effect,
        aes(year - year_mpa, mpa_effect),
        show.legend = F,
        alpha = 0.5
      ) 

    out <- list(
      did_values = did_values,
      mean_effect = mean_effect
    )
      

   

    return(out)
  }
