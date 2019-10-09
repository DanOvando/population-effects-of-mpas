fit_zissou <- function(data,
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
                       non_nested_did_variables = c(
                         "temp", "kelp", "lag_catch"
                       ),
                       script_name,
                       seed = 42,
                       run_dir,
                       use_tmb = T,
                       fixed_did = T,
                       include_intercept = T,
                       center_scale = T,
                       fixed_regions = F,
                       mpa_only = F,
                       bin_years = FALSE) {

  # data <- data %>%
  #   filter(classcode != "opic")

  if (bin_years == TRUE){

    data$year <- plyr::round_any(data$year, 3)

    data$factor_year <- factor(data$year)

  }


  if (mpa_only == T){

    data <- data %>%
      filter(eventual_mpa == T)
  }


  if (center_scale == T){

    data_recipe <- recipes::recipe(log_density ~ ., data = data) %>%
      recipes::step_center(all_numeric(), -all_outcomes(),-year,-targeted,-geographic_cluster) %>%
      recipes::step_scale(all_numeric(), -all_outcomes(),-year,-targeted,-geographic_cluster)

    prepped_data <- recipes::prep(data_recipe, data, retain = TRUE)

    data <- recipes::juice(prepped_data)

  }

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

  candidate_years <- data$year %>% unique()

  seen_data <- add_mising_years(seen_data, candidate_years)

  log_density <- seen_data$log_density

  seen_data <- seen_data %>%
    select(-log_density)

  seeing_data <- data %>%
    left_join(numeric_species_key, by = "classcode") %>%
    slice(seeing_has_important$index)

  seeing_data <- add_mising_years(seeing_data, candidate_years)

  any_seen <- seeing_data$any_seen

  seeing_data <- seeing_data %>%
    select(-any_seen)

  # prepare data for c++ ----------------------------------------------------

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

  standard_year_species <- expand.grid(year = unique(seeing_data$year), classcode = unique(seeing_data$classcode), stringsAsFactors = F) %>%
    as_data_frame() %>%
    mutate(marker = 1,
           classcode_year = paste(classcode, year, sep = "-")) %>%
    spread(classcode_year, marker, fill = 0) %>%
    arrange(classcode, year) %>%
    left_join(numeric_species_key, by = "classcode")


  annual_data <- data %>%
    group_by(year, classcode) %>%
    summarise(
              temp = mean(temp_deviation),
              targeted = unique(targeted),
              kelp = mean(interp_kelp),
              lag_catch = mean(lag_catch)) %>%
    ungroup()

  did_data <- standard_year_species %>%
    select(year, classcode) %>%
    left_join(annual_data, by = c("year",'classcode')) %>%
    arrange(classcode, year)

  non_nested_did_data <- did_data %>%
    select(non_nested_did_variables) %>%
    # select(enso, pdo, temp) %>%
    mutate(intercept = 1)

  if (any(is.na(non_nested_did_data))){

    non_nested_did_data <- non_nested_did_data %>%
      mutate(index = 1:nrow(.)) %>%
      gather(variable, value,-index) %>%
      group_by(variable) %>%
      mutate(value = zoo::na.approx(value)) %>%
      ungroup() %>%
      spread(variable, value) %>%
      select(-index)


  }

  year_did_data <- did_data %>%
    mutate(targeted_year = paste(targeted,year, sep = '-'),
           marker = 1) %>%
    select(year, classcode, targeted_year, marker) %>%
    spread(targeted_year,marker, fill = 0) %>%
    arrange(classcode, year)  %>%
    select(-year,-classcode)

  targeted <- str_detect(colnames(year_did_data),"1-")

  targeted_year_did_data <- year_did_data[,targeted]

  nontargeted_year_did_data <- year_did_data[,targeted == F]

  species_did_data <- did_data %>%
    select(classcode) %>%
    mutate(marker = 1, index = 1:nrow(.)) %>%
    spread(classcode, marker, fill = 0) %>%
    arrange(index) %>%
    select(-index)


  standard_species <- standard_year_species$numeric_classcode

  standard_year_species <-  standard_year_species %>%
    select(-year,-classcode,-numeric_classcode)

  standard_non_nested <- seeing_cdata$x_non_nested %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_non_nested <-
    standard_non_nested[rep(1, nrow(standard_year_species)), ]

  standard_non_nested <-
    standard_non_nested[colnames(seeing_cdata$x_non_nested)]

  standard_region_cluster <- seeing_cdata$x_region_cluster %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) %>%
    spread(variable, mean_value)

  standard_region_cluster <-
    standard_region_cluster[rep(1, nrow(standard_year_species)), ]

  standard_region_cluster <-
    standard_region_cluster[colnames(seeing_cdata$x_region_cluster)]

  # fit model ---------------------------------------------------------------

  seen_species_index <- seen_data$numeric_classcode

  n_species <- length(unique(seen_species_index))

  zissou_data <- list(
    x_seen_non_nested = seen_cdata$x_non_nested,
    x_seen_year_species = seen_cdata$x_year_species,
    x_seen_region_cluster = seen_cdata$x_region_cluster,
    x_seeing_non_nested = seeing_cdata$x_non_nested,
    x_seeing_year_species = seeing_cdata$x_year_species,
    x_seeing_region_cluster = seeing_cdata$x_region_cluster,
    region_cluster_index = seeing_cdata$region_cluster_index,
    log_density = log_density,
    any_seen = any_seen,
    standard_non_nested = standard_non_nested,
    standard_year_species = standard_year_species,
    standard_region_cluster = standard_region_cluster,
    seen_species_index = seen_species_index,
    year_species_index = seen_cdata$year_species_index,
    seen_weights = rep(1, nrow(seen_cdata$x_non_nested)),
    seeing_weights = rep(1, nrow(seeing_cdata$x_non_nested)),
    non_nested_did_data = non_nested_did_data,
    targeted_year_did_data = targeted_year_did_data,
    nontargeted_year_did_data = nontargeted_year_did_data,
    species_did_data = species_did_data
  )
  zissou_data <- map_if(zissou_data, is.data.frame, ~ as.matrix(.x))

  any_na <- map_lgl(zissou_data, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in zissou_data")
  }
  zissou_params <- list(
    seen_non_nested_betas = rep(0, ncol(zissou_data$x_seen_non_nested)),
    seen_year_species_betas = rep(0, ncol(zissou_data$x_seen_year_species)),
    seen_year_species_sigmas = rep(log(1), n_species),
    seen_region_cluster_betas = rep(0, ncol(zissou_data$x_seen_region_cluster)),
    seen_region_cluster_sigmas = rep(log(1), n_distinct(zissou_data$region_cluster_index)),
    seen_density_species_sigma = rep(log(1), n_species),
    seeing_non_nested_betas = rep(0, ncol(zissou_data$x_seeing_non_nested)),
    seeing_year_species_betas = rep(0, ncol(zissou_data$x_seen_year_species)),
    seeing_year_species_sigmas = rep(log(1), n_species),
    seeing_region_cluster_betas = rep(0, ncol(zissou_data$x_seeing_region_cluster)),
    seeing_region_cluster_sigmas = rep(log(1), n_distinct(zissou_data$region_cluster_index)),
    non_nested_did_betas = rep(0, ncol(non_nested_did_data)),
    targeted_did_betas = rep(0, ncol(targeted_year_did_data)),
    nontargeted_did_betas = rep(0, ncol(nontargeted_year_did_data)),
    species_did_betas = rep(0, ncol(species_did_data)),
    log_targeted_sigma = log(1),
    log_nontargeted_sigma = log(1),
    log_species_sigma = log(1),
    log_did_sigma = log(1)
  )
  any_na <- map_lgl(zissou_params, ~ any(is.na(.x))) %>% any()

  if (any_na) {
    stop("NAs in zissou_params")
  }

  if (use_tmb == T) {

    compile(here::here("src", paste0(script_name, ".cpp")), "-O0") # what is the -O0?

    dyn.load(dynlib(here::here("src", script_name)))

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
        "seeing_region_cluster_betas",
        "targeted_did_betas",
        "nontargeted_did_betas",
        "species_did_betas"
      )    }
    zissou_model <-
      MakeADFun(
        zissou_data,
        zissou_params,
        DLL = script_name,
        random = randos
      )

    # save(file = here::here(run_dir, "zissou-onestage-tmb-model.Rdata"),
    #      zissou_model)
# browser()
#
#
#     fit <- tmbstan::tmbstan(zissou_model, cores = 1,
#                             chains = 2000)
#
    # a <- Sys.time()
    set.seed(seed)
    zissou_fit <-
      nlminb(
        zissou_model$par,
        zissou_model$fn,
        zissou_model$gr,
        control = list(iter.max = 4000, eval.max = 5000)
      )
    # Sys.time() - a

    # save(file = here::here(run_dir, "zissou-onestage-tmb-fit.Rdata"),
    #      zissou_fit)

    zissou_report <- zissou_model$report()

    zissou_sd_report <- sdreport(zissou_model)

    diagnostics = data.frame(
      "name" = names(zissou_model$par),
      "est" = zissou_fit$par,
      "final_gradient" = as.vector(zissou_model$gr(zissou_fit$par))
    )


    # save(file = here::here(run_dir, "zissou-tmb-onestage-sdreport.Rdata"),
    #      zissou_sd_report)
    #
    # save(file = here::here(run_dir, "zissou-tmb-onestage-report.Rdata"),
    #      zissou_report)


  }


  zissou_estimates <-
    summary(zissou_sd_report) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    set_names(tolower) %>%
    rename(std_error = `std. error`) %>%
    mutate(lower = estimate - 1.96 * std_error,
           upper = estimate + 1.96 * std_error)


  return(
    list(
      zissou_fit = zissou_fit,
      zissou_report = zissou_report,
      zissou_sd_report = zissou_sd_report,
      zissou_estimates = zissou_estimates,
      did_data = did_data,
      diagnostics = diagnostics,
      zissou_data = zissou_data,
      seen_data = seen_data
    )
  )
}