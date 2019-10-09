fit_ahnold <- function(data,
                       non_nested_variables =  c(
                         'targeted',
                         'factor_year',
                         'level',
                         'factor_month',
                         'cumulative_n_obs',
                         'temp_deviation',
                         'surge',
                         'mean_canopy',
                         'mean_depth'
                       ),
                       script_name,
                       seed = 42,
                       run_dir,
                       use_tmb = T,
                       fixed_did = T,
                       include_intercept = T,
                       center_and_scale = T,
                       fixed_regions = F) {


  numeric_species_key <-
    data_frame(classcode = unique(data$classcode)) %>%
    arrange(classcode) %>%
    mutate(numeric_classcode = 1:nrow(.))


  seen_has_important <- data %>%
    filter(any_seen == T) %>%
    select(non_nested_variables) %>%
    mutate(index = 1:nrow(.)) %>%
    na.omit()

  seeing_has_important <- data %>%
    select(non_nested_variables) %>%
    mutate(index = 1:nrow(.)) %>%
    na.omit()

  seen_data <- data %>%
    filter(any_seen == T) %>%
    left_join(numeric_species_key, by = "classcode") %>%
    slice(seen_has_important$index)

  log_density <- seen_data$log_density

  seen_data <- seen_data %>%
    select(-log_density)

  seeing_data <- data %>%
    left_join(numeric_species_key, by = "classcode") %>%
    slice(seeing_has_important$index)

  any_seen <- seeing_data$any_seen

  seeing_data <- seeing_data %>%
    select(-any_seen)


  # prepare data for c++ ----------------------------------------------------


browser()
  seen_cdata <-
    make_c_worthy(seen_data,
                  non_nested_vars = non_nested_variables,
                  fixed_regions = fixed_regions,
                  include_intercept = include_intercept,
                  fixed_did = fixed_did,
                  numeric_species_key = numeric_species_key)

  seeing_cdata <-
    make_c_worthy(seeing_data,
                  non_nested_vars = non_nested_variables,
                  fixed_regions = fixed_regions,
                  include_intercept = include_intercept,
                  fixed_did = fixed_did,
                  numeric_species_key = numeric_species_key)

  # prepare standardized matrices -------------------------------------------

  standard_non_nested <- seeing_cdata$x_non_nested %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_non_nested <-
    standard_non_nested[rep(1, n_distinct(seeing_data$factor_year)), ]

  standard_non_nested <-
    standard_non_nested[colnames(seeing_cdata$x_non_nested)]


  standard_year_species <- seeing_cdata$x_year_species %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_year_species <-
    standard_year_species[rep(1, n_distinct(seeing_data$factor_year)), ]

  standard_year_species <-
    standard_year_species[colnames(seeing_cdata$x_year_species)]


  standard_region_cluster <- seeing_cdata$x_region_cluster %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_region_cluster <-
    standard_region_cluster[rep(1, n_distinct(seeing_data$factor_year)), ]

  standard_region_cluster <-
    standard_region_cluster[colnames(seeing_cdata$x_region_cluster)]

  standard_did_with_mpa <-
    data_frame(year = unique(seeing_data$factor_year)) %>%
    spread_factors(drop_one = fixed_did)

  standard_did_without_mpa <- standard_did_with_mpa

  standard_did_without_mpa[standard_did_without_mpa > 0] = 0



  # fit model ---------------------------------------------------------------

  seen_species_index <- seen_data$numeric_classcode

  n_species <- length(unique(seen_species_index))


  ahnold_data <- list(
    x_seen_non_nested = seen_cdata$x_non_nested,
    x_seen_did = seen_cdata$x_did,
    x_seen_year_species = seen_cdata$x_year_species,
    x_seen_region_cluster = seen_cdata$x_region_cluster,
    x_seeing_non_nested = seeing_cdata$x_non_nested,
    x_seeing_did = seeing_cdata$x_did,
    x_seeing_year_species = seeing_cdata$x_year_species,
    x_seeing_region_cluster = seeing_cdata$x_region_cluster,
    region_cluster_index = seeing_cdata$region_cluster_index,
    log_density = log_density,
    any_seen = any_seen,
    standard_non_nested = standard_non_nested,
    standard_year_species = standard_year_species,
    standard_did_with_mpa = standard_did_with_mpa,
    standard_did_without_mpa = standard_did_without_mpa,
    standard_region_cluster = standard_region_cluster,
    seen_species_index = seen_species_index,
    year_species_index = seen_cdata$year_species_index,
    seen_weights = rep(1, nrow(seen_cdata$x_non_nested)),
    seeing_weights = rep(1, nrow(seeing_cdata$x_non_nested))
  )

  ahnold_data <- map_if(ahnold_data, is.data.frame, ~ as.matrix(.x))

  any_na <- map_lgl(ahnold_data, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in ahnold_data")
  }
  ahnold_params <- list(
    seen_non_nested_betas = rep(0, ncol(ahnold_data$x_seen_non_nested)),
    seen_did_betas = rep(0, ncol(ahnold_data$x_seen_did)),
    seen_year_species_betas = rep(0, ncol(ahnold_data$x_seen_year_species)),
    seen_year_species_sigmas = rep(log(1), n_species),
    seen_region_cluster_betas = rep(0, ncol(ahnold_data$x_seen_region_cluster)),
    seen_region_cluster_sigmas = rep(log(1), n_distinct(ahnold_data$region_cluster_index)),
    seen_density_species_sigma = rep(log(1), n_species),
    seeing_non_nested_betas = rep(0, ncol(ahnold_data$x_seeing_non_nested)),
    seeing_did_betas = rep(0, ncol(ahnold_data$x_seeing_did)),
    seeing_year_species_betas = rep(0, ncol(ahnold_data$x_seen_year_species)),
    seeing_year_species_sigmas = rep(log(1), n_species),
    seeing_region_cluster_betas = rep(0, ncol(ahnold_data$x_seeing_region_cluster)),
    seeing_region_cluster_sigmas = rep(log(1), n_distinct(ahnold_data$region_cluster_index))
  )

  any_na <- map_lgl(ahnold_params, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in ahnold_params")
  }

  if (use_tmb == T) {
    compile(here::here("scripts", paste0(script_name, ".cpp")), "-O0") # what is the -O0?

    dyn.load(dynlib(here::here("scripts", script_name)))

    if (fixed_regions == T){

      randos <- c(
        "seen_year_species_betas",
        "seeing_year_species_betas"
      )
    } else {
      randos = c(
        "seen_year_species_betas",
        "seeing_year_species_betas",
        "seen_region_cluster_betas",
        "seeing_region_cluster_betas"
      )    }

    ahnold_model <-
      MakeADFun(
        ahnold_data,
        ahnold_params,
        DLL = script_name,
        random = randos
      )

    save(file = here::here(run_dir, "ahnold-onestage-tmb-model.Rdata"),
         ahnold_model)


    a <- Sys.time()
    set.seed(seed)
    ahnold_fit <-
      nlminb(
        ahnold_model$par,
        ahnold_model$fn,
        ahnold_model$gr,
        control = list(iter.max = 4000, eval.max = 5000)
      )
    Sys.time() - a

    save(file = here::here(run_dir, "ahnold-onestage-tmb-fit.Rdata"),
         ahnold_fit)

    ahnold_report <- ahnold_model$report()

    ahnold_sd_report <- sdreport(ahnold_model)


    save(file = here::here(run_dir, "ahnold-tmb-onestage-sdreport.Rdata"),
         ahnold_sd_report)

    save(file = here::here(run_dir, "ahnold-tmb-onestage-report.Rdata"),
         ahnold_report)


  }


  ahnold_estimates <-
    summary(ahnold_sd_report) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    set_names(tolower) %>%
    rename(std_error = `std. error`) %>%
    mutate(lower = estimate - 1.96 * std_error,
           upper = estimate + 1.96 * std_error)


  return(
    list(
      seen_cdata = seen_cdata,
      seeing_cdata = seeing_cdata,
      ahnold_fit = ahnold_fit,
      ahmold_model = ahnold_model,
      ahnold_report = ahnold_report,
      ahnold_sd_report = ahnold_sd_report,
      ahnold_estimates = ahnold_estimates
    )
  )
}