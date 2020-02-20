process_zissou_fit <- function(fit){

seen_non_nested_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_non_nested_betas")) %>%
  rename(group = variable) %>%
  mutate(variable  = fit$zissou_data$x_seen_non_nested %>% colnames())

seeing_non_nested_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_non_nested_betas")) %>%
  rename(group = variable) %>%
  mutate(variable  = fit$zissou_data$x_seen_non_nested %>% colnames())

seen_year_species_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_year_species_betas")) %>%
  rename(group = variable) %>%
  mutate(variable = fit$zissou_data$x_seen_year_species %>% colnames())

seeing_year_species_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_year_species_betas")) %>%
  # filter(variable == "seeing_year_species_betas") %>%
  rename(group = variable) %>%
  mutate(variable = fit$zissou_data$x_seen_year_species %>% colnames())

seen_region_cluster_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seen_region_cluster_betas")) %>%
  # filter(variable == "seen_region_cluster_betas") %>%
  rename(group = variable) %>%
  mutate(variable = fit$zissou_data$x_seen_region_cluster %>% colnames())

seeing_region_cluster_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "seeing_region_cluster_betas")) %>%
  # filter(variable == "seeing_region_cluster_betas") %>%
  rename(group = variable) %>%
  mutate(variable = fit$zissou_data$x_seen_region_cluster %>% colnames())

did_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "mpa_effect")) %>%
  # filter(variable == "mpa_effect") %>%
  mutate(group = variable) %>%
  mutate(year = fit$did_data$year %>% unique())



non_targeted_abundance_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "nontargeted_did")) %>%
  # filter(variable == "mpa_effect") %>%
  mutate(group = variable) %>%
  mutate(year = fit$did_data$year %>% unique())

targeted_abundance_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "targeted_did") & !stringr::str_detect(variable, "non")) %>%
  # filter(variable == "mpa_effect") %>%
  mutate(group = variable) %>%
  mutate(year = fit$did_data$year %>% unique())


did_betas <- fit$zissou_estimates %>%
  filter(stringr::str_detect(variable, "mpa_effect")) %>%
  # filter(variable == "mpa_effect") %>%
  mutate(group = variable) %>%
  mutate(year = fit$did_data$year %>% unique())


betas <- bind_rows(
  targeted_abundance_betas,
  non_targeted_abundance_betas,
  seen_non_nested_betas,
  seeing_non_nested_betas,
  seen_year_species_betas,
  seeing_year_species_betas,
  seen_region_cluster_betas,
  seeing_region_cluster_betas,
  did_betas %>% select(-year)
) %>%
  as_data_frame()


beta_names <- ls()[str_detect(ls(), "beta")]

out <- purrr::map(beta_names, ~get(.x)) %>% 
  set_names(beta_names)

return(out)

}


