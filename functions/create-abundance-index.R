create_abundance_index <-
  function(seen_model,
           seeing_model,
           seeing_aug,
           seen_aug,
           pop_structure) {


    # find seen omitted terms
    # omitted_levels <-
    #   data_frame(factor_terms = names(seen_model$xlevels),
    #              omitted_term = NA)
    #
    # seen_coefs <- broom::tidy(seen_model)
    # for (i in seq_along(1:nrow(omitted_levels))) {
    #   possible_levels <- seen_model$xlevels[[i]]
    #
    #   estimated_levels <- seen_coefs %>%
    #     filter(str_detect(term, omitted_levels$factor_terms[i])) %>% {
    #       .$term
    #     }
    #
    #   estimated_levels <-
    #     str_replace(estimated_levels, omitted_levels$factor_terms[i], '')
    #
    #   omitted_level <-
    #     possible_levels[!possible_levels %in% estimated_levels]
    #
    #   omitted_levels$omitted_term[i] <- omitted_level
    #
    # }
    #
    # omitted_levels <- omitted_levels %>%
    #   mutate(omitted_name = map2_chr(factor_terms, omitted_term, ~ paste(.x, .y, sep = '_')))
    #
    # intercept_term <-
    #   paste(omitted_levels$omitted_name, collapse = '-') # a whole ton of work to double check what the intercept is

    # extract and convert intercept and year terms the old fashioned way
    # abundance_index <- seen_coefs %>%
    #   filter(term == "(Intercept)" | str_detect(term, 'factor_year')) %>%
    #   mutate(term = str_replace(term, 'factor_year',''),
    #          abundance_index = exp(estimate + std.error ^ 2 / 2)) %>%
    #   select(term, abundance_index) %>%
    #   filter(term != '(Intercept)') %>%
    #   ungroup() #%>%
    #mutate(abundance_index = center_scale(abundance_index))

    # in order to not drop a year use model to predict abundance over time for the most frequent factor levels held constant

    # seen_aug <-   seen_model %>%
    #   broom::augment()
    if (pop_structure == 'one-pop') {
      seen_series <-
        create_reference_case(seen_aug = seen_aug, seen_model = seen_model) %>%
        mutate(smearing_term = mean(exp(seen_aug$.resid)))
    }
    if (pop_structure == 'regional-pops')
    {

      regional_smear <- seen_aug %>%
        group_by(region) %>%
        summarise(smearing_term = mean(exp(.resid)))

      # regions <- seen_aug$region %>% unique()

      seen_series <- seen_aug %>%
        nest(-region) %>%
        # mutate(data = map2(data,region, ~.x %>% mutate(region = factor(.y, levels = regions)), regions = regions)) %>%
        mutate(seen_series = map(data, create_reference_case, seen_model = seen_model)) %>%
        select(-data) %>%
        unnest() %>%
        left_join(regional_smear, by = 'region') %>%
        ungroup()

    }
    if (pop_structure == 'mpa-pops') {

      smearing_term <- seen_aug %>%
        group_by(eventual_mpa) %>%
        summarise(smearing_term = mean(exp(.resid)))


      seen_series <- seen_aug %>%
        nest(-eventual_mpa) %>%
        mutate(seen_series = map(data, create_reference_case, seen_model = seen_model)) %>%
        select(-data) %>%
        unnest() %>%
        left_join(smearing_term, by = 'eventual_mpa') %>%
        ungroup()


    }

    # smearing_term <- mean(exp(seen_aug$.resid))

    seen_series <- seen_series %>%
      modelr::add_predictions(seen_model) %>%
      mutate(linear_predictions = exp(pred) * smearing_term) %>%
      ungroup() #%>%
    # mutate(linear_predictions = linear_predictions / max(linear_predictions))

    # ggplot() +
    #   geom_point(data = abundance_index, aes(term, abundance_index)) +
    #  geom_point(data = seen_series, aes(factor_year, linear_predictions), color = 'red')


    # create reference frame for probabilities and calculate abundance index
    prob_seen <- predict(seeing_model, newdata = seen_series, type = 'response') # calculate probabilty of seeing
        # predict(seeing_model, newdata = broom::augment(seeing_model), type = 'response') # calculate probabilty of seeing

    if (pop_structure == 'one-pop') {
      seen_series <- seen_series %>%
        mutate(prob_seen = prob_seen) %>%
        mutate(
          abundance_index = linear_predictions * prob_seen)

    }
    if (pop_structure == 'regional-pops')
    {
      seen_series <- seen_series %>%
        mutate(prob_seen = prob_seen) %>%
        group_by(region) %>%
        mutate(
          abundance_index = linear_predictions * prob_seen)


    }
    if (pop_structure == 'mpa-pops') {
      seen_series <- seen_series %>%
        mutate(prob_seen = prob_seen) %>%
        group_by(eventual_mpa) %>%
        mutate(
          abundance_index = linear_predictions * prob_seen)

    }



    # seen_series %>%
    #   ggplot(aes(factor_year, abundance_index)) +
    #   geom_point()
    return(seen_series)

    # predict observation probabilities for reference frame

    # create composite index

  }
