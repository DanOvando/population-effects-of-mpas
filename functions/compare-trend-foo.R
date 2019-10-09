compare_trend_foo <- function(common_name, data, run_dir) {
  data <- data %>%
    select(
      population_structure,
      data_source,
      population_filtering,
      abundance_index,
      raw_abundance_trend
    ) %>%
    mutate(
      simple_abundance_index = map2(abundance_index, population_structure, simplify_trend),
      simple_raw_index = map2(raw_abundance_trend, population_structure, simplify_trend)
    ) %>%
    select(-abundance_index, -raw_abundance_trend) %>%
    gather(abundance_type, data, contains('_index')) %>%
    unnest()

  length_to_density_pop_trends_plot <- data %>%
    filter(data_source == 'length_to_density') %>%
    group_by(population_level, abundance_type, population_structure,population_filtering) %>%
    mutate(abundance_index = abundance_index / max(abundance_index)) %>%
    ggplot(aes(
      year,
      abundance_index,
      color = population_level,
      linetype = abundance_type
    )) +
    geom_line() +
    facet_grid(population_structure ~ population_filtering) +
    labs(title = common_name)

  raw_density_pop_trends_plot <- data %>%
    filter(data_source == 'supplied_density') %>%
    group_by(population_level, abundance_type, population_structure,population_filtering) %>%
    mutate(abundance_index = abundance_index / max(abundance_index)) %>%
    ggplot(aes(
      year,
      abundance_index,
      color = population_level,
      linetype = abundance_type
    )) +
    geom_line() +
    facet_grid(population_structure ~ population_filtering) +
    labs(title = common_name)

  kfm_density_pop_trends_plot <- data %>%
    filter(data_source == 'kfm_density') %>%
    group_by(population_level, abundance_type, population_structure,population_filtering) %>%
    mutate(abundance_index = abundance_index / max(abundance_index)) %>%
    ggplot(aes(
      year,
      abundance_index,
      color = population_level,
      linetype = abundance_type
    )) +
    geom_line() +
    facet_grid(population_structure ~ population_filtering) +
    labs(title = common_name)

  ggsave(
    file = paste0(
      run_dir,
      '/',
      common_name,
      '-length_to_density_pop_trends.pdf'
    ),
    length_to_density_pop_trends_plot,
    height = 8,
    width = 8
  )

  ggsave(
    file = paste0(run_dir, '/', common_name, '-raw_density_pop_trends.pdf'),
    raw_density_pop_trends_plot,
    height = 8,
    width = 8
  )


  ggsave(
    file = paste0(run_dir, '/', common_name, '-kfm_density_pop_trends.pdf'),
    kfm_density_pop_trends_plot,
    height = 8,
    width = 8
  )


}