vast_prep <- function(data, classcode,site_coords,conditions_data){

  data <- data %>%
    left_join(site_coords, by = c('site', 'side')) %>%
    mutate(
      year = factor_year %>% as.character() %>% as.numeric(),
      areaswept_km2 = (30 * 4) * .001,
      spp = classcode
    )

  if (is.null(data$observer)){
    data$observer = 'unknown'
  }
  if (is.null(data$total_biomass_g)){
    data$total_biomass_g = data$biomass / 1000
  }

  data <- data %>%
    # select(year, latitude, longitude, observer,areaswept_km2, total_biomass_g, spp, mean_vis) %>%
    rename(lat = latitude, lon = longitude, vessel = observer, catch_kg = total_biomass_g) %>%
    filter(is.na(catch_kg) == F)

  data <- data %>%
    group_by(spp,year) %>%
    mutate(min_catch = min(catch_kg, na.rm = T),
           npos = sum(catch_kg > 0)) %>%
    ungroup() %>%
    mutate(catch_kg = ifelse(catch_kg == min_catch, 0, catch_kg)) %>%
    filter(npos > 0)
}