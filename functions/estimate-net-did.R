estimate_net_did <- function(seen_model,
                             seeing_model){


seen_model <- abundance_models$seen_model[[1]]

seeing_model <- abundance_models$seen_model[[1]]


a <- abundance_models %>%
  filter(data_source == 'length_to_density',
         population_structure == 'one-pop',
         population_filtering == 'all'
         ) %>%
  select(classcode,data) %>%
  unnest()

seen_data <- a %>%
  filter(any_seen== T) %>%
  mutate(targeted = (targeted == 'Targeted') %>% as.numeric) %>%
  left_join(enso, by = 'year') %>%
  left_join(pdo, by = 'year') %>%
  left_join(ci_catches, by = c('classcode', 'year')) %>%
  mutate(catch = ifelse(is.na(catch), 0, catch),
         post_mpa = year >=2003)




test <- lm(log_density ~  targeted + post_mpa +targeted:factor_year + region + mean_vis + factor_month +
             trunc_observer + cumulative_n_obs + level + surge + mean_temp +
             mean_enso + mean_pdo + lag1_enso + lag2_enso + lag1_pdo + lag2_pdo, data = seen_data)

dids <- test %>%
  broom::tidy() %>%
  filter(str_detect(term,'targeted:')) %>%
  mutate(year = str_replace_all(term,'\\D','') %>% as.numeric())


dids %>%
  ggplot() +
  geom_pointrange(aes(x = year, y = estimate, ymin =  estimate - 1.96*std.error,
                      ymax = estimate + 1.96*std.error))

}