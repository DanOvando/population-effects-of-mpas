#' add missing fish to data
#'
#' @param this_site
#' @param this_side
#' @param this_year
#' @param this_transect
#' @param observations
#' @param species_sightings
#'
#' @return
#' @export
#'
add_missing_fish <-
  function(this_region,
           this_site,
           this_side,
           this_year,
           this_month,
           this_day,
           this_transect,
           this_zone,
           this_level,
           observations,
           species_sightings,
           life_history_vars) {
    sampling_event <- observations %>%
      ungroup() %>%
      filter(
        site == this_site,
        side == this_side,
        year == this_year,
        month == this_month,
        day == this_day,
        transect == this_transect,
        zone == this_zone,
        level == this_level
      )


    if (dplyr::n_distinct(sampling_event$observer) > 1) {
      stop('multiple observers per transect - check add_missing_fish')
    }

    species_seen <- unique(sampling_event$classcode)

    # species_possible <- species_sightings %>%
    #   filter(region == this_region) %>%
    #   select(species_seen) %>%
    #   unlist() %>%
    #   as.character()

    species_possible <- species_sightings %>%
      filter(site == this_site) %>%
      select(species_seen) %>%
      unlist() %>%
      as.character()

    species_missing <-
      species_possible[!species_possible %in% species_seen]

    if (length(species_missing) > 0) {

      blank <- sampling_event[rep(1, length(species_missing)), ] %>%
        mutate(classcode = species_missing)

      blank[, str_detect(colnames(blank),'biomass_g')] <-  0

      # %>%
      #   gather('metric', 'value', dplyr::contains('biomass_g')) %>%
      #   mutate(value = 0) %>%
      #   spread(metric, value)

      missing_fish <- blank[, colnames(sampling_event)]

      missing_fish <- missing_fish %>%
        mutate(count = 0)

      sampling_event <- sampling_event %>%
        bind_rows(missing_fish)

    }
    return(sampling_event)

  }