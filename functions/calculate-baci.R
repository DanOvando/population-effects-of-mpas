#' calculate baci
#'
#' calculates before-after-control-impact metric
#'
#' @param dat the simulation jdata
#' @param distance a matrix of distance to/from MPA border
#' @param type one of 'biased' or "unbiased" where biased is the world with mpa, unbiased without
#' @param before_window the number of years to average over for the 'before' treatment
#'
#' @return a dataframe of BACI estimates over time
#' @export
#'
calculate_baci <- function(dat, distance_grid = distance_grid,type = "biased", before_window = 1){

  # results$mpa_experiment[[1]]$raw_outcomes -> dat

  ex <- ifelse(type == "biased","with-mpa","no-mpa")

  mpa_year <- min(dat$year[dat$mpa == TRUE])

  sum_dat <- dat %>%
    filter(experiment == ex) %>%
    left_join(distance_grid, by = c("patch" = "from")) %>%
    group_by(year, patch, eventual_mpa) %>%
    summarise(biomass = sum(biomass),
              min_fished_dist = unique(min_fished_dist),
              min_dist = unique(min_dist)) %>%
    group_by(year) %>%
    summarise(treated_biomass = weighted.mean(biomass[eventual_mpa == TRUE], min_fished_dist[eventual_mpa == TRUE]),
              control_biomass = weighted.mean(biomass[eventual_mpa == FALSE], min_dist[eventual_mpa == FALSE]))

  before_treatment <-
    mean(sum_dat$treated_biomass[sum_dat$year < mpa_year &
                           sum_dat$year >= mpa_year - before_window])

  before_control <-
    mean(sum_dat$control_biomass[sum_dat$year < mpa_year &
                                   sum_dat$year >= mpa_year - before_window])

  baci <- sum_dat %>%
    group_by(year) %>%
    summarise(baci = (treated_biomass - before_treatment) - (control_biomass - before_control))

}