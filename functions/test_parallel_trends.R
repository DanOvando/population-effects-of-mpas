test_parallel_trends <- function(data, mpa_year = 2003){

  data <- data %>%
    mutate(targeted = as.numeric(targeted > 0))

  pre_species_correlations <- data %>%
    select(year,population_level,classcode, abundance_index) %>%
    spread(classcode, abundance_index) %>%
    filter(year <= mpa_year) %>%
    select(-year,-population_level) %>%
    corrr::correlate()


  post_species_correlations <- data %>%
    select(year,population_level,classcode, abundance_index) %>%
    spread(classcode, abundance_index) %>%
    filter(year > mpa_year) %>%
    select(-year,-population_level) %>%
    corrr::correlate()

  overall_species_correlations <- data %>%
    select(year,population_level,classcode, abundance_index) %>%
    spread(classcode, abundance_index) %>%
    select(-year,-population_level) %>%
    corrr::correlate()


  trend_data <- data %>%
    mutate(abundance_index = (abundance_index - mean(abundance_index)) / (2 * sd(abundance_index))) %>%
    group_by(year, targeted) %>%
    filter(year <= mpa_year) %>%
    summarise(mean_abundance = mean(abundance_index)) %>%
    spread(targeted, mean_abundance)

  pre_correlation_test <- cor.test(~ `0` + `1`, data = trend_data)


  trend_data <- data %>%
    mutate(abundance_index = (abundance_index - mean(abundance_index)) / (2 * sd(abundance_index))) %>%
    group_by(year, targeted) %>%
    filter(year > mpa_year) %>%
    summarise(mean_abundance = mean(abundance_index)) %>%
    spread(targeted, mean_abundance)

  post_correlation_test <- cor.test(~ `0` + `1`, data = trend_data)


  trend_data <- data %>%
    mutate(abundance_index = (abundance_index - mean(abundance_index)) / (2 * sd(abundance_index))) %>%
    group_by(year, targeted) %>%
    summarise(mean_abundance = mean(abundance_index)) %>%
    spread(targeted, mean_abundance)

  overall_correlation_test <- cor.test(~ `0` + `1`, data = trend_data)

  out <- list(pre_correlation_test = pre_correlation_test,
              post_correlation_test = post_correlation_test,
              overall_correlation_test = overall_correlation_test,
              pre_species_correlations = pre_species_correlations,
              post_species_correlations = post_species_correlations,
              overall_species_correlations = overall_species_correlations)

}
