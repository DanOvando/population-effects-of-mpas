process_fits <- function(zissou_fit) {
  seen_non_nested_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seen_non_nested_betas")) %>%
    rename(group = variable) %>%
    mutate(variable  = zissou_fit$zissou_data$x_seen_non_nested %>% colnames())

  seeing_non_nested_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seeing_non_nested_betas")) %>%
    rename(group = variable) %>%
    mutate(variable  = zissou_fit$zissou_data$x_seen_non_nested %>% colnames())

  seen_year_species_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seen_year_species_betas")) %>%
    rename(group = variable) %>%
    mutate(variable = zissou_fit$zissou_data$x_seen_year_species %>% colnames())

  seeing_year_species_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seeing_year_species_betas")) %>%
    # filter(variable == "seeing_year_species_betas") %>%
    rename(group = variable) %>%
    mutate(variable = zissou_fit$zissou_data$x_seen_year_species %>% colnames())

  seen_region_cluster_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seen_region_cluster_betas")) %>%
    # filter(variable == "seen_region_cluster_betas") %>%
    rename(group = variable) %>%
    mutate(variable = zissou_fit$zissou_data$x_seen_region_cluster %>% colnames())

  seeing_region_cluster_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "seeing_region_cluster_betas")) %>%
    # filter(variable == "seeing_region_cluster_betas") %>%
    rename(group = variable) %>%
    mutate(variable = zissou_fit$zissou_data$x_seen_region_cluster %>% colnames())

  did_betas <- zissou_fit$zissou_estimates %>%
    filter(stringr::str_detect(variable, "mpa_effect")) %>%
    # filter(variable == "mpa_effect") %>%
    mutate(group = variable) %>%
    mutate(year = zissou_fit$did_data$year %>% unique())


  betas <- bind_rows(
    seen_non_nested_betas,
    seeing_non_nested_betas,
    seen_year_species_betas,
    seeing_year_species_betas,
    seen_region_cluster_betas,
    seeing_region_cluster_betas,
    did_betas %>% select(-year)
  ) %>%
    as_data_frame()

  non_nested_beta_plot <- betas %>%
    filter(str_detect(group, "non_nested")) %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
    geom_pointrange(aes(
      x = variable,
      y = estimate,
      ymin = lower,
      ymax = upper
    )) +
    facet_wrap( ~ group) +
    coord_flip() +
    theme(axis.text.y = element_text(size = 10))


  year_species_effects_plots <- betas %>%
    filter(str_detect(group, "year_species_betas")) %>%
    mutate(year = str_replace_all(variable, "\\D", "") %>% as.numeric()) %>%
    mutate(classcode = str_split(variable, '-', simplify = T)[, 2]) %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
    geom_ribbon(aes(
      x = year,
      ymin = lower,
      ymax = upper,
      fill = classcode
    ), alpha = 0.25) +
    geom_line(aes(x = year, y = estimate, color = classcode)) +
    facet_wrap( ~ group)

  region_cluster_plots <- betas %>%
    filter(str_detect(group, "region_cluster_betas")) %>%
    mutate(cluster = str_replace_all(variable, "\\D", "")) %>%
    mutate(region = str_split(variable, '-', simplify = T)[, 3]) %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), linetype = 2, color = "red") +
    geom_pointrange(aes(
      x = region,
      y = estimate,
      ymin = lower,
      ymax = upper,
      color = cluster
    )) +
    facet_wrap( ~ group)


  did_plot <- did_betas %>%
    ggplot() +
    geom_vline(
      aes(xintercept = 2003),
      color = 'red',
      linetype = 2,
      size = 2
    ) +
    geom_hline(aes(yintercept = 0)) +
    geom_pointrange(aes(
      year,
      y = estimate,
      ymin = lower,
      ymax = upper
    ),
    size = 1.5) +
    labs(x = "Year", y = "Divergence",
         caption = "'Divergence' refers to difference in mean density of targeted and non-targeted species")

  locals <- ls()

  out <- map(locals, ~ get(.x)) %>%
    set_names(locals)

  return(out)
}
