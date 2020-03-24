estimate_did <-
  function(data,
           site_data,
           life_history_data,
           bin_width = 3,
           min_classcode_years = 15,
           consistent_sites_only = TRUE,
           cdfw_catches,
           data_to_use = "all",
           data_source = "pisco",
           chains = 4,
           cores = 4,
           iter = 2000) {
    
    # data <- pisco_abundance_data
    # 
    # data <- kfm_abundance_data

    if (data_to_use == "mpa_only"){
      
      data <- data %>% 
        filter(eventual_mpa == TRUE)
      
    } else if (data_to_use == "fished_only"){
      
      data <- data %>% 
        filter(eventual_mpa == FALSE)
    }
    
    annual_catches <- cdfw_catches %>%
      filter(classcode %in% unique(data$classcode)) %>%
      group_by(year) %>%
      summarise(var_catch = sum(catch)) %>%
      mutate(var_lag_catch = lag(var_catch, 1))
    
    annual_catches$var_lag_catch[1] <-   annual_catches$var_lag_catch[2] #assume 1999 catches = 2000 catches, not crazy given AC in totals
    
    consistent_did_sites <- data %>%
      select(site_side, year) %>%
      unique() %>%
      group_by(site_side) %>%
      mutate(has_all = all((2003:2003) %in% year)) %>% # results change dramatically if you filter until present time...
      ungroup() %>%
      filter(has_all)
    
    consistent_classcodes <- data %>%
      group_by(classcode) %>%
      filter(any_seen) %>%
      summarise(n_years = n_distinct(year)) %>%
      ungroup() %>%
      filter(n_years >= min_classcode_years)
    
    year_bins <-
      c(1999, seq(2003, max(data$year) + 1, by = bin_width))
    
    did_data <- data %>% {
      if (consistent_sites_only == TRUE) {
        filter(.,
               site_side %in% unique(consistent_did_sites$site_side))
      } else{
        .
      }
    } %>%
    filter(classcode %in% consistent_classcodes$classcode) 
    # rm("data","cdfw_catches","life_history_data")
    if (data_source == "pisco"){
    
# did_data %>% 
#         ggplot(aes(regional_temp_dev)) + 
#         geom_histogram() + 
#         facet_wrap(~region)
      
    did_data <- did_data %>% 
      group_by(year,
               site_side,
               region,
               zone,
               transect,
               eventual_mpa,
               classcode,
               targeted) %>%
      summarise(
        total_classcode_density = sum(exp(log_density)),
        var_tex = sum(cumulative_n_obs),
        var_vis = mean(mean_vis),
        var_temp = mean(mean_temp),
        var_depth = mean(mean_depth),
        var_surge = mean(surge),
        var_kelp = mean(interp_kelp)
      ) %>%  # sum density across all levels of a transect
      group_by(year, site_side, region, eventual_mpa, classcode, targeted) %>%
      summarise(
        md = mean(total_classcode_density),
        var_tex = mean(var_tex),
        var_vis = mean(var_vis),
        var_temp = mean(var_temp),
        var_depth = mean(var_depth),
        var_surge = mean(var_surge),
        var_kelp = mean(var_kelp)
      ) %>% # calculate mean density per year site, side, species, averaging over zone, transect
      group_by(year, site_side, region, eventual_mpa, targeted) %>%
      summarise(
        total_biomass_density = (sum(md) / 1e6) * 10000,
        # calculate total and mean biomass densities across all species per year site side
        mean_biomass_density = (mean(md) / 1e6) * 10000,
        var_tex = mean(var_tex),
        var_vis = mean(var_vis),
        var_temp = mean(var_temp),
        var_depth = mean(var_depth),
        var_surge = mean(var_surge),
        var_kelp = mean(var_kelp)
      ) %>%
      ungroup() %>%
      mutate(fyear = factor(year)) %>%
      mutate(fyear = relevel(fyear, "2003")) %>%
      mutate(year_bins = cut(year, year_bins, include.lowest = TRUE)) %>% 
      left_join(annual_catches, by = "year") %>% 
      group_by(region) %>% 
      mutate(regional_temp_dev = scale(var_temp)) %>% 
      ungroup()

    vars <- which(str_detect(colnames(did_data), "var_"))
    
    nafoo <- function(x){
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    }
    did_data <- purrrlyr::dmap_at(did_data, vars, ~ scale(.x)) %>%
     purrrlyr::dmap_at(vars, nafoo) %>%
      mutate(total_biomass_density = total_biomass_density + 1e-6,
             mean_biomass_density = mean_biomass_density + 1e-6) %>% 
      group_by(site_side, targeted) %>%
      mutate(scaled_total_biomass_density = scale(log(total_biomass_density)))
    # filter(region != "SMI")
    
    # did_data %>%
    #   group_by(region, targeted) %>%
    #   mutate(total_biomass_density = scale(total_biomass_density)) %>%
    #   ggplot(aes(year, total_biomass_density, color = targeted == 1)) +
    #   geom_point() +
    #   geom_smooth() +
    #   facet_wrap( ~ region)
    # 
    # 
    # did_data %>%
    #   # group_by(region, targeted) %>%
    #   # mutate(total_biomass_density =scale(total_biomass_density)) %>%
    #   ggplot(aes(year, total_biomass_density, color = targeted == 1)) +
    #   geom_point() +
    #   geom_smooth()
    
    
    # did_reg <-
    #   stan_glmer(
    #     log(total_biomass_density + 1e-6) ~ targeted * year_bins + (site_side - 1 |
    #                                                                   region) + var_tex + var_vis + var_temp + var_depth + var_surge + var_kelp + var_catch + var_lag_catch,
    #     data = did_data,
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    
    env <- new.env(parent = .GlobalEnv)
    
    env$did_data <- did_data
    
    env$cores <- cores
    
    env$chains <- chains
    
    env$iter <-  iter
    
    # did_reg <- with(env, {
    #   stan_glmer(
    #     log(total_biomass_density + 1e-6) ~ targeted * year_bins + (site_side - 1 |
    #                                                                   region) + var_tex + var_surge + var_kelp + var_catch,
    #     data = did_data,
    #     cores = cores,
    #     chains = chains,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    # })
    did_reg <- with(env, {
      stan_glmer(
        total_biomass_density ~ targeted * year_bins + (site_side - 1 |
                                                          region) + var_tex + var_surge + var_kelp + var_catch + var_temp +
          regional_temp_dev,
        data = did_data,
        cores = cores,
        chains = chains,
        prior_intercept = normal(autoscale = TRUE),
        prior = normal(0, 2, autoscale = FALSE),
        iter = iter,
        family = Gamma(link = "log")
      )
    })
    

    # ln_did_reg <- with(env,{
    #   stan_glmer(
    #     log(total_biomass_density + 1e-6) ~ targeted * year_bins +  (site_side - 1 |
    #                                                                    region)  + var_tex + var_surge + var_kelp + var_catch,
    #     data = did_data,
    #     cores = cores,
    #     chains = chains,
    #     prior_intercept = normal(autoscale = TRUE),
    #     prior = normal(0, 2),
    #     iter = iter
    #   )
    # })
    # browser()
    # 
    # 
    # gr <- broom::augment(gamma_did_reg)
    # 
    # 
    # 
    # lnr <- broom::augment(ln_did_reg)
    # 
    # color_scheme_set("red")
    # 
    # sigma <- sd()
    # 
    # ppc_dens_overlay(y = exp(lnr$log.total_biomass_density...1e.06.),
    #                  yrep = exp(posterior_predict(ln_did_reg, draws = 50)))
    # 
    # ln_ppd <- exp(posterior_predict(ln_did_reg, draws = 200) + 0.8^2/2) %>% 
    #   as.data.frame() %>% 
    #   mutate(i = 1:nrow(.)) %>% 
    #   pivot_longer(-i, names_to = "j", values_to = "value") %>% 
    #   group_by(i) %>% 
    #   summarise(n = sd(value),
    #             m = mean(value))
    # 
    # ln_ppd_ln <- (posterior_predict(ln_did_reg, draws = 200)) %>% 
    #   as.data.frame() %>% 
    #   mutate(i = 1:nrow(.)) %>% 
    #   pivot_longer(-i, names_to = "j", values_to = "value") %>% 
    #   group_by(i) %>% 
    #   summarise(n = sd(value),
    #             m = mean(value))
    # 
    # 
    # ggplot() + 
    #   geom_density(data = ln_ppd_ln, aes(n, fill = "ln"), alpha = 0.5) + 
    #   geom_vline(xintercept = sd((lnr$log.total_biomass_density...1e.06.))) + 
    #   labs(caption = "distributions are posterior predictive standard deviation on the natrual scale. line is empirical SD on natural scale")
    # 
    # 
    # gamma_ppd <- (posterior_predict(gamma_did_reg, draws = 100)) %>% 
    #   as.data.frame() %>% 
    #   mutate(i = 1:nrow(.)) %>% 
    #   pivot_longer(-i, names_to = "j", values_to = "value") %>% 
    #   group_by(i) %>% 
    #   summarise(n = sd(value),
    #             m = mean(value))
    # 
    # ggplot() + 
    #   geom_density(data = ln_ppd, aes(n, fill = "ln"), alpha = 0.5) + 
    #   geom_density(data = gamma_ppd, aes(n, fill = "gamma"), alpha = 0.5) + 
    #   geom_vline(xintercept = sd(exp(lnr$log.total_biomass_density...1e.06.))) + 
    #   labs(caption = "distributions are posterior predictive standard deviation on the natrual scale. line is empirical SD on natural scale")
    # 
    # 
    # ggplot() + 
    #   geom_density(data = ln_ppd, aes(m, fill = "ln"), alpha = 0.5) + 
    #   geom_density(data = gamma_ppd, aes(m, fill = "gamma"), alpha = 0.5) + 
    #   geom_vline(xintercept = mean(exp(lnr$log.total_biomass_density...1e.06.))) + 
    #   labs(caption = "distributions are posterior predictive standard deviation on the natrual scale. line is empirical SD on natural scale")
    # 
    # 
    # color_scheme_set("red")
    # ppc_dens_overlay(y = log(gr$total_biomass_density),
    #                  yrep = log(posterior_predict(gamma_did_reg, draws = 50)))
    # 
    # 
    # gr %>% 
    #   mutate(lresid = log(total_biomass_density) - log(.fitted)) %>% 
    #   ggplot(aes(sample = lresid)) + 
    #   geom_qq() + 
    #   geom_qq_line()
    # 
    # lnr %>% 
    #   ggplot(aes(sample = .resid)) + 
    #   geom_qq() + 
    #   geom_qq_line()
    # 
    # 
    # gr %>% 
    #   ggplot(aes(.fitted, .resid)) + 
    #   geom_point()
    # 
    # gr %>% 
    #   ggplot(aes(var_catch, .resid)) + 
    #   geom_point()
    # 
    # gr %>% 
    #   ggplot(aes(.resid, fill = region)) + 
    #   geom_histogram() + 
    #   facet_wrap(~region, scales = "free")
    # 
    # gr %>% 
    #   ggplot(aes(.resid, fill = year_bins)) + 
    #   geom_histogram() + 
    #   facet_wrap(~year_bins, scales = "free")
    # 
    # lnr %>% 
    #   ggplot(aes(exp(.fitted), .resid)) + 
    #   geom_point()
    # 
    # lnr %>% 
    #   ggplot(aes(.resid, fill = region)) + 
    #   geom_histogram() + 
    #   facet_wrap(~region, scales = "free")
    # 
    # lnr %>% 
    #   ggplot(aes(.resid, fill = year_bins)) + 
    #   geom_histogram() + 
    #   facet_wrap(~year_bins, scales = "free")
    # 
    # lnr %>% 
    #   ggplot(aes(sample = .resid)) + 
    #   geom_qq() + 
    #   geom_qq_line()
    # 
    # gr %>% 
    #   ggplot(aes(sample = .resid)) + 
    #   geom_qq() + 
    #   geom_qq_line()
    
    
    } else if (data_source == "kfm"){
      did_data <- did_data %>% 
        group_by(year,
                 site_side,
                 region,
                 eventual_mpa,
                 classcode,
                 targeted) %>%
        summarise(
          total_classcode_density = sum(exp(log_density)),
          var_temp = mean(mean_temp),
          var_kelp = mean(interp_kelp)
        ) %>%  # sum density across all levels of a transect
        group_by(year, site_side, region, eventual_mpa, classcode, targeted) %>%
        summarise(
          md = mean(total_classcode_density),
          var_temp = mean(var_temp),
          var_kelp = mean(var_kelp)
        ) %>% # calculate mean density per year site, side, species, averaging over zone, transect
        group_by(year, site_side, region, eventual_mpa, targeted) %>%
        summarise(
          total_biomass_density = (sum(md) / 1e6) * 10000,
          # calculate total and mean biomass densities across all species per year site side
          mean_biomass_density = (mean(md) / 1e6) * 10000,
          var_temp = mean(var_temp),
          var_kelp = mean(var_kelp)
        ) %>%
        ungroup() %>%
        mutate(fyear = factor(year)) %>%
        mutate(fyear = relevel(fyear, "2003")) %>%
        mutate(year_bins = cut(year, year_bins,include.lowest = TRUE)) %>% 
        left_join(annual_catches, by = "year") %>% 
        group_by(region) %>% 
        mutate(regional_temp_dev = scale(var_temp)) %>% 
        ungroup()
      
      vars <- which(str_detect(colnames(did_data), "var_"))
      
      nafoo <- function(x){
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
      }
      
      did_data <- purrrlyr::dmap_at(did_data, vars, ~ scale(.x)) %>%
        purrrlyr::dmap_at(vars, nafoo) %>%
        filter(total_biomass_density > 0) %>% 
        group_by(site_side, targeted) %>%
        mutate(scaled_total_biomass_density = scale(log(total_biomass_density))) %>% 
        ungroup()
    
      env <- new.env(parent = .GlobalEnv)
      
      env$did_data <- did_data
      
      env$cores <- cores
      
      env$chains <- chains
      
      env$iter <-  iter
      
      did_reg <- with(env,{
        stan_glmer(
          total_biomass_density ~ targeted * year_bins + (1 |region) + var_kelp + var_catch,
          data = did_data,
          cores = cores,
          chains = chains,
          prior_intercept = normal(0, 2, autoscale = TRUE),
          prior = normal(0, 2),
          family = Gamma(link = "log"),
          iter = iter
        )})
      
      
    }
    
    # med_did_reg <-
    #   stan_glmer(
    #     log(total_biomass_density + 1e-6) ~ targeted * year_bins + (site_side - 1 |
    #                                                                   region),
    #     data = did_data,
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    # 
    # min_did_reg <-
    #   stan_glm(
    #     log(total_biomass_density + 1e-6) ~ targeted * year_bins,
    #     data = did_data,
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    # 
    # 
    # compare <- loo_compare(loo(did_reg), loo(simpler_did_reg))
    
    # 
    # scaled_did_reg <-
    #   stan_glmer(
    #     (scaled_total_biomass_density) ~ targeted * year_bins + (site_side |
    #                                                                region) + var_tex + var_vis + var_temp + var_depth + var_surge + var_kelp,
    #     data = did_data,
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    
    # did_reg <-
    #   stan_glmer(
    #     log(total_biomass_density + 1e-6) ~   year_bins + targeted:year_bins + (region  - 1 |
    #                                                                               targeted) + var_tex + var_vis + var_temp + var_depth + var_surge + var_kelp,
    #     data = did_data,
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    
    # did_reg_gamma <-
    #   stan_glmer(
    #     (total_biomass_density + 1e-6) ~ targeted * year_bins + (site_side - 1 |
    #                                                                region) + var_tex + var_tex_2 + var_vis + var_temp + var_depth + var_surge + var_kelp,
    #     data = did_data,
    #     family = Gamma(link = "log"),
    #     cores = 4,
    #     chains = 4,
    #     prior_intercept = normal(0, 2),
    #     prior = normal(0, 2)
    #   )
    # 
    
    did_results <- tidybayes::tidy_draws(did_reg) %>%
      select(contains("."), contains("targeted:year_bins")) %>%
      pivot_longer(
        contains("targeted"),
        names_to = "year",
        values_to = "did",
        names_prefix = "targeted:year_bins"
      ) %>%
      mutate(did =  exp(did) - 1) %>%
      group_by(year) %>%
      mutate(prank = percent_rank(did)) %>%
      ungroup() 
    
    # 
    # mpa_effect_plot <-  did_results %>%
    #   ggplot(aes(year, did)) +
    #   geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
    #   tidybayes::stat_halfeye(alpha = 0.7,
    #                           .width = c(0.5, 0.95)) +
    #   scale_y_continuous(labels = percent, name = "Estimated MPA Effect") +
      scale_x_discrete(name = "Year Bin")
    # 
    out <- list(did_results = did_results,
                did_reg = did_reg
                )
    
  }