#' assign tuned variables
#'
#' tunes a control variable to achieve a target depletion level
#'
#' @param control_variable
#' @param target_depletion
#' @param fish
#' @param fleet
#' @param sim_years
#' @param burn_years
#' @param num_patches
#'
#' @return ss for depletion
#' @export
#'
#' @examples
tune_fishery <- function( f_v_m,
                               fish,
                               fleet,
                               sim_years = 50,
                               burn_years = 1,
                               num_patches = 10,
                          sprinkler = sprinkler,
                          mpa_habfactor = mpa_habfactor){
  

  if (fleet$fleet_model == "constant-catch"){

    tol <- 1e-3

    lower <- 0

    upper <- 400

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <-  0

    while(delta_best > tol & counter < 20) {

      counter <- counter + 1

      constant <- (1 - golden) * (upper - lower)

      x1 <- lower + constant

      x2 <- upper - constant

      yield_1 <- estimate_msy(x1, fish = fish, fleet = fleet,        num_patches = 1,
                              sim_years = sim_years,
                              burn_years = burn_years)

      yield_2 <- estimate_msy(x2, fish = fish, fleet = fleet,
                              num_patches = 1,
                              sim_years = sim_years,
                              burn_years = burn_years)

      delta_best <-  (best -  min(yield_1,yield_2))^2

      best <- min(yield_1,yield_2)

      if (yield_1 < yield_2){

        lower <- lower

        upper <- x2
      } else{

        lower <- x1

        upper <- upper

      }

    } # close golden while

    msy_fit <-
      nlminb(
        mean(c(lower, upper)),
        estimate_msy,
        fish = fish,
        fleet = fleet,
        lower = 0,
        num_patches = 1,
        sim_years = sim_years,
        burn_years = burn_years
      )

    fleet$e_msy <- msy_fit$par

    fish$msy <- -msy_fit$objective

    fish$b_msy <- estimate_msy(fleet$e_msy, fish = fish, fleet = fleet, use = "other", num_patches = 1)

    msy <- fish$msy

    fleet$target_catch <- msy * f_v_m
    
  } else if (fleet$fleet_model == "constant-effort"){

    fleet$initial_effort <-   ((fish$m * f_v_m) / fleet$q) * num_patches

  } else if (fleet$fleet_model == "open-access"){

    temp_fleet <- fleet
    
    temp_fleet$initial_effort <-   ((fish$m * f_v_m) / fleet$q) * num_patches
    
    temp_fleet$fleet_model <- "constant-effort"
    
    oa <- sim_fishery(
      fish = fish, 
      fleet = temp_fleet,
      num_patches = num_patches,
      sim_years = sim_years,
      burn_years = burn_years,
      sprinkler = sprinkler,
      mpa_habfactor = mpa_habfactor,
      create_manager(mpa_size = 0)
    )
    
    # plot_spasm(oa, type = "totals")
    
    b_oa <- oa %>% 
      filter(year > (.75 * max(year))) %>% 
      group_by(year) %>% 
      summarise(b = sum(biomass),
                b0 = unique(b0)) %>% 
     mutate(oa_dep = b / b0) %>% 
    summarise(oa_dep = mean(oa_dep))
    
    b_ref_oa <- b_oa$oa_dep # target open access level of depletion
    
    # crs <- seq(0.0001,.9, length.out = 25)
    # 
    # huh <- map_dbl(crs, estimate_costs, fish = fish,
    #                fleet = fleet,
    #                b_ref_oa = b_ref_oa,
    #                lags = 1,
    #                num_patches = num_patches,
    #                sim_years = sim_years,
    #                burn_years = burn_years,
    #                sprinkler = sprinkler,
    #                mpa_habfactor = mpa_habfactor
    # )
    # 
    # browser()
    
    cost_fit <-
      nlminb(
        .5,
        estimate_costs,
        fish = fish,
        fleet = fleet,
        b_ref_oa = b_ref_oa,
        lower = 1e-3,
        upper = 0.9,
        lags = 1,
        num_patches = num_patches,
        sim_years = sim_years,
        burn_years = burn_years,
        sprinkler = sprinkler,
        mpa_habfactor = mpa_habfactor
      )
    # browser()



    # cost_fit <-         nlminb(
    #   mean(c(lower, upper)),
    #   estimate_costs,
    #   fish = fish,
    #   fleet = fleet,
    #   lower = 0,
    #   msy = fish$msy,
    #   e_msy = fleet$e_msy,
    #   b_msy = fish$b_msy$b_msy,
    #   p_response = fleet$max_perc_change_f,
    #   b_v_bmsy_oa =  pmax(.05,2 - f_v_m)
    # )

    fleet$max_cr_ratio <- cost_fit$par

    # fleet$max_perc_change_f <-  20
    # 
    # oa <- sim_fishery(
    #   fish = fish, 
    #   fleet = fleet,
    #   num_patches = num_patches,
    #   sim_years = sim_years,
    #   burn_years = burn_years,
    #   sprinkler = sprinkler,
    #   mpa_habfactor = mpa_habfactor,
    #   create_manager(mpa_size = 0)
    # )
    # 
    # plot_spasm(oa)
    
    # fleet$p_msy <- fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta


  } else if (fleet$fleet_model == "supplied-catch"){


    tol <- 100

    lower <- 0

    upper <- 400

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <-  0

    while(delta_best > tol & counter < 20) {

      counter <- counter + 1

      constant <- (1 - golden) * (upper - lower)

      x1 <- lower + constant

      x2 <- upper - constant

      yield_1 <- estimate_msy(x1, fish = fish, fleet = fleet)

      yield_2 <- estimate_msy(x2, fish = fish, fleet = fleet)

      delta_best <-  (best -  min(yield_1,yield_2))^2

      best <- min(yield_1,yield_2)

      if (yield_1 < yield_2){

        lower <- lower

        upper <- x2
      } else{

        lower <- x1

        upper <- upper

      }

    } # close golden while

    msy_fit <- nlminb(mean(c(lower, upper)), estimate_msy, fish = fish, fleet = fleet, lower = 0)

    fleet$e_msy <- msy_fit$par

    fish$msy <- -msy_fit$objective

    fish$b_msy <- estimate_msy(fleet$e_msy, fish = fish, fleet = fleet, use = "other")

    msy <- fish$msy




  }

  return(list(fish = fish, fleet = fleet))

}




