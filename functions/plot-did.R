plot_did <- function(model){


 check <-  model$stan_summary %>%
    as.data.frame() %>%
    mutate(term = rownames(.)) %>%
    filter(str_detect(term, 'targeted'))

 if (any(str_detect(check$term, 'generations'))){

   did_terms <- model$stan_summary %>%
     as.data.frame() %>%
     mutate(term = rownames(.)) %>%
     filter(str_detect(term, 'targeted')) %>%
     mutate(year = map_dbl(term, ~str_extract_all(.x,'(?<=(protected)).*') %>% as.numeric())) %>%
     filter(!is.na(year))

 } else{

   did_terms <- model$stan_summary %>%
     as.data.frame() %>%
     mutate(term = rownames(.)) %>%
     filter(str_detect(term, 'targeted')) %>%
     mutate(year = map_dbl(term, ~str_replace_all(.x,'\\D','') %>% as.numeric())) %>%
     filter(!is.na(year))

 }


  did_terms %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), color = 'red') +
    # geom_vline(aes(xintercept = 2003), color = 'blue') +
    geom_pointrange(aes(
      x = year,
      y = mean,
      ymin = `2.5%`,
      ymax = `97.5%`
    ))

  # unfished_terms <- model$stan_summary %>%
  #   as.data.frame() %>%
  #   mutate(term = rownames(.)) %>%
  #   filter(str_detect(term, 'factor_year') & !str_detect(term,'targeted')) %>%
  #   mutate(year = map_dbl(term, ~str_replace_all(.x,'\\D','') %>% as.numeric())) %>%
  #   filter(!is.na(year))
  #
  # unfished_terms %>%
  #   ggplot() +
  #   geom_hline(aes(yintercept = 0), color = 'red') +
  #   geom_pointrange(aes(
  #     x = year,
  #     y = mean,
  #     ymin = `2.5%`,
  #     ymax = `97.5%`
  #   ))



}