diagnostic_plots <- function(model){

  # loo_model <- rstanarm::kfold(model, K = 10)
  safe_bayes <- safely(rstanarm::bayes_R2)

  a <- safe_bayes(model)

  if (is.null(a$error)){

  loo_model <- rstanarm::loo(model)

  r2 <- rstanarm::bayes_R2(model)

  r2_hist_plot <- data_frame(r2 = r2) %>%
    ggplot(aes(r2)) +
    geom_histogram(color = 'black', fill = 'grey')

  poster_predictive_plot <- rstanarm::pp_check(model)

  } else{

    loo_model <- NA

    r2 <- NA

    r2_hist_plot <- ggplot()

    poster_predictive_plot <- ggplot()

  }


  augmod <- broom::augment(model)

  normal_qq_plot <- augmod %>%
    ggplot(aes(sample = .resid)) +
    stat_qq() +
    stat_qq_line(color = 'red')

  hist_resid_plot <- augmod %>%
    ggplot(aes(.resid)) +
    geom_histogram() +
    geom_vline(aes(xintercept = 0), color = 'blue', linetype = 2) +
    geom_vline(aes(xintercept = mean(.resid)), color = 'red', linetype = 3)


  resid_v_fitted_plot <- augmod %>%
    ggplot(aes(.fitted, .resid)) +
    geom_point() +
    geom_hline(aes(yintercept = 0), color = 'red')

  obs_v_predicted_plot <- augmod %>%
    ggplot(aes(log_abundance_index, .fitted)) +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
    geom_smooth(method = 'lm') +
    labs(title = glue::glue('Mean Bayesian R2 is {round(mean(r2),2)}'))


  comp_plot <-
    (obs_v_predicted_plot + labs(title = 'A')) + (r2_hist_plot+labs(title = 'B')) +
    plot_layout(ncol = 1,
                nrow = 2,
                heights = c(3, 1))

  if (class(model)[[1]] == 'lmerMod'){

    ind_vars <- model@frame %>% colnames()

  }else{

  ind_vars <- (model$terms %>% as.character())[3]

  ind_vars <- str_split(ind_vars, pattern = '\\*|\\+', simplify = T)

  ind_vars <-  map_chr(ind_vars, ~str_replace_all(.x, ' ',''))
}
  numeric_vars <- colnames(augmod)[map_lgl(augmod, is.numeric)]

  factor_vars <-  colnames(augmod)[map_lgl(augmod, ~!is.numeric(.))]


  numeric_coef_v_resid <- augmod %>%
    select(ind_vars[ind_vars %in% numeric_vars], .resid) %>%
    gather(variable, value, -.resid)

  factor_coef_v_resid <- augmod %>%
    select(ind_vars[ind_vars %in% factor_vars], .resid) %>%
    gather(variable, value, -.resid)

  numeric_coef_v_resid_plot <- numeric_coef_v_resid %>%
    ggplot(aes(value, .resid)) +
    geom_point() +
    facet_wrap(~variable, scales = 'free_x')


  factor_coef_v_resid_plot <- factor_coef_v_resid %>%
    ggplot(aes(value, .resid)) +
    geom_boxplot() +
    facet_wrap(~variable, scales = 'free_y') +
    coord_flip()



out <- list(factor_coef_v_resid_plot = factor_coef_v_resid_plot,
            numeric_coef_v_resid_plot = numeric_coef_v_resid_plot,
            resid_v_fitted_plot = resid_v_fitted_plot,
            hist_resid_plot = hist_resid_plot,
            normal_qq_plot = normal_qq_plot,
            comp_plot = comp_plot,
            poster_predictive_plot = poster_predictive_plot,
            loo_model = loo_model)

}