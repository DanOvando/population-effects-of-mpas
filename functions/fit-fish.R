fit_fish <-
  function(data,
           pop_structure,
           pop_filter,
           ind_covars,
           dep_var,
           family,
           fit_seen = T,
           consistent_sites,
           model_type = 'glm') {
    if (pop_filter == 'mpa-only') {
      data <- data %>%
        filter(eventual_mpa == T)

    }
    if (pop_filter == 'consistent-sites') {
      data <- data %>%
        filter(site %in% consistent_sites$site)
    }

    if (fit_seen == T) {
      data <- data %>%
        filter(any_seen == fit_seen & is.na(any_seen) == F)
    } else{
      data <- data %>%
        filter(is.na(any_seen) == F)
    }



    if (pop_structure == 'one-pop') {
      pop_term <- 'factor_year'

      if (model_type == 'glmer') {
        pop_term <- '(1|factor_year)'
      }

    }
    if (pop_structure == 'regional-pops') {

      num_levels <- length(unique(data$region))

      pop_term <- 'factor_year:region'

      if (num_levels > 1) {
        pop_term <- 'factor_year:region'
      } else {
        pop_term <- 'factor_year'
      }


      if (model_type == 'glmer') {

        if (num_levels > 1) {
          pop_term <- '(factor_year|region)'
        } else {
          pop_term <- '(1|factor_year)'
        }
      }

    }

    if (pop_structure == 'mpa-pops') {

      num_levels <- length(unique(data$eventual_mpa))

      pop_term <- 'factor_year:eventual_mpa'

      if (num_levels > 1) {
        pop_term <- '(factor_year:eventual_mpa)'
      } else {
        pop_term <- 'factor_year'
      }

      if (model_type == 'glmer') {

        if (num_levels > 1) {
          pop_term <- '(factor_year|eventual_mpa)'
        } else {
          pop_term <- '(1|factor_year)'
        }
      }
    }
    ind_vars <- paste(c(pop_term, ind_covars), collapse = '+')
    reg_fmla <- paste(dep_var, ind_vars, sep = '~') %>% as.formula()

    if (model_type == 'gam') {

    }
    if (model_type == 'glm') {
      model <- glm(reg_fmla, data  = data, family = family, y = FALSE, model = FALSE)
    }
    if (model_type == 'glmer') {
      model <-
        lme4::glmer(
          reg_fmla,
          data  = data,
          family = family,
          control = lmerControl(optimizer = 'Nelder_Mead')
        )
    }

return(model)
  }