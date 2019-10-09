add_mising_years <- function(data, candidate_years){


  missing_some <- data %>%
    group_by(classcode) %>%
    summarise(missing_some = n_distinct(year) < length(candidate_years)) %>%
    filter(missing_some)

  fillfoo <- function(data, years_missing){

    blank <- data[1,]

    blank <- map_df(years_missing, ~ mutate(blank,factor_year = as.character(.x)), blank = blank)

    blank$year <- as.numeric(blank$factor_year)

    data <- data %>%
      bind_rows(blank) %>%
      arrange(numeric_classcode, year)

    data$density_g_m2 <- exp(data$log_density)

    data$density_g_m2[data$year %in% years_missing] <-  NA

    density_trend <- data %>%
      group_by(year) %>%
      summarise(mean_density = mean(density_g_m2)) %>%
      ungroup() %>%
      mutate(interped_density = zoo::na.approx(mean_density, rule = 2)) %>%
      select(-mean_density)

    data <- data %>%
      left_join(density_trend, by = "year") %>%
    mutate(log_density = ifelse(is.na(log_density), log(interped_density),log_density )) %>%
      select(-interped_density,-density_g_m2)

    return(data)
}

  needs_filling <- data %>%
    filter(classcode %in% missing_some$classcode) %>%
    nest(-classcode) %>%
    mutate(years_missing = map(data,~candidate_years[!candidate_years %in% unique(.x$year)], candidate_years =candidate_years )) %>%
    mutate(filled_data = map2(data, years_missing, fillfoo)) %>%
    select(classcode, filled_data) %>%
    unnest()

  out <- data %>%
    filter(!classcode %in% missing_some$classcode) %>%
    bind_rows(needs_filling) %>%
    arrange(classcode, year)

  return(out)
}