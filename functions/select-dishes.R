select_dishes <- function(did_data, run_dir) {

  did_data <- did_data %>%
    filter(population_filtering == 'all',
           data_source == 'length_to_density',
           population_structure == 'one-pop',
           abundance_source == 'glm_abundance_index')

time_enviro_correlations <- did_data %>%
  select(year, mean_enso, mean_pdo, mean_annual_kelp, mean_annual_temp) %>%
  cor()

pdf(glue::glue("{run_dir}/did_covariate_correlations.pdf"))
corrplot::corrplot(time_enviro_correlations)
dev.off()


# dishes <- 'loo + (mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + mean_annual_kelp + temp_deviation |classcode)'

dishes <- c('loo + (mean_enso + mean_annual_kelp + temp_deviation +mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo - 1 |classcode)',
            '(mean_enso + mean_annual_kelp + temp_deviation +mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo - 1 |classcode)',
            'loo + (mean_enso + mean_annual_kelp + temp_deviation - 1|classcode)',
            'loo + (mean_enso + lag1_enso + lag2_enso + lag3_enso + lag4_enso + mean_annual_kelp + temp_deviation - 1 |classcode)',
            'loo + (mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + mean_annual_kelp + temp_deviation - 1 |classcode)',
            'loo + mean_enso + lag1_enso + lag2_enso + lag3_enso + lag4_enso + mean_annual_kelp + temp_deviation',
            'loo + mean_pdo + lag1_pdo + lag2_pdo + lag3_pdo + lag4_pdo + mean_annual_kelp + temp_deviation',
            'loo + mean_enso + mean_annual_kelp + temp_deviation' ,
            'loo + (mean_enso + mean_annual_kelp + temp_deviation - 1|geographic_cluster)')

did_models <-
  cross_df(
    list(
      did_data = list(did_data),
      timing = c('years'),
      complexity = c('kitchen_sink'),
      dirty_dishes = dishes
    )
  )


candidate_models <- list()

candidate_models_kfold <- list()


for (i in 1:nrow(did_models)){

  did_term <- 'targeted:factor_year'


  did_reg <-
    paste0('log_abundance_index ~', paste(
      c(
        'targeted',
        'factor_year',
        did_term,
        did_models$dirty_dishes[i]),
      collapse = '+'
    ))

  if (str_detect(did_reg,'\\|')){

    # fitfoo <- rstanarm::stan_glmer
#
#     candidate_models[[i]] <- rstanarm::stan_glmer(did_reg, data = did_models$data[[i]],
#                                                cores = cores, chains = chains)

        candidate_models[[i]] <- lme4::glmer(did_reg, data = did_models$did_data[[i]])

  } else {


    candidate_models[[i]] <- glm(did_reg, data = did_models$did_data[[i]])

    # candidate_models[[i]] <- rstanarm::stan_glm(did_reg, data = did_models$data[[i]],
    #                                            cores = cores, chains = chains)

  }


  candidate_models_kfold[[i]] <- AIC(candidate_models[[i]]) #  stanarm::kfold(candidate_models[[i]], K = 10)

}

did_models$BIC <- candidate_models_kfold %>% as.numeric()

did_models <- did_models %>%
  ungroup() %>%
  arrange(BIC)

# did_models <- did_models %>%
#   mutate(did_model = pmap(list(did_data = did_data,
#                                timing = timing,
#                                complexity = complexity,
#                                dirty_dishes = dirty_dishes), fit_did,
#                           cores = 1,
#                           chains = 1)) %>%
#   select(-did_data) %>%
#   unnest()
#
#
# did_models <- did_models %>%
#   mutate(did_plot = map(did_model, plot_did),
#          did_diagnostics = map(did_model, diagnostic_plots))
#
# did_models <- did_models %>%
#   mutate(loo_results = map(did_diagnostics,'loo_model')) %>%
#   mutate(elpd_loo = map_dbl(loo_results,'elpd_loo')) %>%
#   arrange(desc(elpd_loo))
#
# rstanarm::compare_models(loos = did_models$loo_results)

best_dish <- did_models$dirty_dishes[[1]]

}


