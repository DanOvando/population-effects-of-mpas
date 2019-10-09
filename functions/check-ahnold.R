check_ahnold <- function(length_to_density_data, abundance_indices, did_models,min_years = 10, min_year = 2000){

  reg_data <- abundance_indices %>%
    filter(data_source == 'length_to_density', population_filtering == 'all',
           population_structure == 'one-pop') %>%
    select(data, classcode) %>%
    unnest()

  did_data <- did_models %>%
    filter(data_source == 'length_to_density', population_filtering == 'all',
           population_structure == 'one-pop',
           timing == 'years', abundance_source == 'glm_abundance_index',
           complexity == 'kitchen_sink') %>%
    select(data) %>%
    unnest()

  check_species_counts <- length_to_density_data %>%
    group_by(year, month, day, site, side, zone, level, transect) %>%
    summarise(num_species = n_distinct(classcode)) %>%
    group_by(site) %>%
    summarise(species_counts_per_event = n_distinct(num_species))

  if (all(check_species_counts$species_counts_per_event == 1)){
    print('pass: adding zeros worked correctly - all events have same # species')

  }else{
    print('fail: different number of species per event - check add_missing_fish')

  }

  check_species_coverage <- reg_data %>%
    group_by(classcode) %>%
    summarise(nyears = n_distinct(year),
              min_year = min(year),
              max_year = max(year))

  if (all(check_species_coverage$nyears >= min_years)){
    print('pass: all species have minimum number of observed years')

  } else{

    print('fail: some species have too few years - check filterfoo')
  }

  if (all(check_species_coverage$min_year <= min_year)){

    print(glue::glue('pass: all species start in {min_year}'))

  } else{
    print('fail: some species do not start in min_year - check filterfoo')
  }

  reg_data %>%
    group_by(targeted) %>%
    summarise(number_of_species = n_distinct(classcode)) %>%
    knitr::kable() %>%
    print()

  catch_check <- did_data %>%
    group_by(targeted) %>%
    summarise(total_reported_catch = sum(catch, na.rm =T))

  print(knitr::kable(catch_check))

  if (catch_check$total_reported_catch[catch_check$targeted == 0] > 0){

    print('Some "non-targeted" species have reported CDFW catches')

    did_data %>%
      filter(targeted == 0) %>%
      group_by(classcode, commonname) %>%
      summarise(total_reported_catch = sum(catch, na.rm = T)) %>%
      filter(total_reported_catch > 0) %>%
      knitr::kable() %>%
      print()

  }


}