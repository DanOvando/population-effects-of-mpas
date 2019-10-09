fit_did <- function(did_data, timing, complexity, dirty_dishes,chains = 4, cores = 4) {
  did_data <- did_data %>%
    nest(-population_structure,
         -data_source,
         -population_filtering,
         -abundance_source)


  # test correlations inspect abundance trends ---------------------------------------
  # did_data <- did_data %>%
  #   mutate(correlation_tests = map(data, test_parallel_trends)) %>%
  #   mutate(
  #     pre_correlation = map_dbl(correlation_tests, ~ .x$overall_correlation_test$estimate),
  #     pre_correlation_signif = map_dbl(correlation_tests, ~ .x$overall_correlation_test$p.value)
  #   )


  # fit DiD estimator on abundance indicies ---------------------------------


  if (timing == 'years'){

    did_term <- 'targeted:factor_year'

  }
  if (timing == 'generations'){

    did_term <- paste(c('generations_protected',
                        'targeted:generations_protected'), collapse = '+')

  }
  if (timing == 'recruits'){

    did_term <- paste(c('recruits_protected',
                        'targeted:recruits_protected'), collapse = '+')

  }

  if (complexity == 'bare_bones'){

    did_reg <-
      paste0('log_abundance_index ~', paste(
        c(
          'targeted',
          'factor_year',
          did_term
        ),
        collapse = '+'
      ))

  }
  if (complexity == 'kitchen_sink'){

    did_reg <-
      paste0('log_abundance_index ~', paste(
        c(
          'targeted',
          'factor_year',
          did_term,
          dirty_dishes),
        collapse = '+'
      ))
  }
  if (str_detect(did_reg,'\\|')){ # if there are random effects

    # fitfoo <- rstanarm::stan_glmer

    did_models <- did_data %>%
      mutate(did_reg = did_reg) %>%
      mutate(did_model = map2(
        data,
        did_reg,
        ~ rstanarm::stan_glmer(
          .y,
          data = .x,
          chains = chains,
          cores = cores
        )
      ))
  } else {

    # fitfoo <- rstanarm::stan_glm

    did_models <- did_data %>%
      mutate(did_reg = did_reg) %>%
      mutate(did_model = map2(
        data,
        did_reg,
        ~  rstanarm::stan_glm(
          .y,
          data = .x,
          chains = chains,
          cores = cores
        )
      ))
  }


} # close fit_did
