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

    tol <- 100

    lower <- 0

    upper <- 400

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <-  0

    while(delta_best > tol | counter < 20) {

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

    msy_fit <- nlminb(mean(c(lower, upper)), estimate_msy, fish = fish, fleet = fleet, lower = 0,
                      num_patches = 1)

    fleet$e_msy <- msy_fit$par

    fish$msy <- -msy_fit$objective

    fish$b_msy <- estimate_msy(fleet$e_msy, fish = fish, fleet = fleet, use = "other")

    msy <- fish$msy

    fleet$target_catch <- msy * f_v_m


  } else if (fleet$fleet_model == "constant-effort"){

    fleet$initial_effort <-   ((fish$m * f_v_m) / fleet$q) * num_patches

  } else if (fleet$fleet_model == "open-access"){


    # estimate MSY

    tol <- 100

    lower <- 0

    upper <- 400

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <-  0

    while(delta_best > tol | counter < 20) {

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

    fleet$initial_effort <- fleet$e_msy

    fish$msy <- -msy_fit$objective

    fish$b_msy <- estimate_msy(fleet$e_msy, fish = fish, fleet = fleet, use = "other")

    msy <- fish$msy

    # tune costs

    tol <- .01

    lower <- 0

    upper <- 1

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <- 0

    set.seed(24)

    while(delta_best > tol) {

      counter <- counter + 1

      constant <- (1 - golden) * (upper - lower)

      x1 <- lower + constant

      x2 <- upper - constant

      # ss_1 <- estimate_costs(cost = x1,
      #                        fish = fish,
      #                        fleet = fleet,
      #                        msy = fish$msy,
      #                        e_msy = fleet$e_msy,
      #                        b_msy = fish$b_msy$b_msy,
      #                        p_response = fleet$max_perc_change_f,
      #                        b_v_bmsy_oa = pmax(.05,2 - f_v_m),
      #                        sim_years = 100)

      ss_1 <-
          estimate_costs(
          max_cr_ratio = x1,
          fish = fish,
          fleet = fleet,
          b_ref_oa = pmax(.05,(2 - f_v_m) / 2),
          lags = 10,
          num_patches = num_patches,
          sim_years = sim_years,
          burn_years = burn_years,
        )


      # ss_2 <- estimate_costs(cost = x2,
      #                        fish = fish,
      #                        fleet = fleet,
      #                        msy = fish$msy,
      #                        e_msy = fleet$e_msy,
      #                        b_msy = fish$b_msy$b_msy,
      #                        p_response = fleet$max_perc_change_f,
      #                        b_v_bmsy_oa = pmax(.05,2 - f_v_m),
      #                        sim_years = 100)

      ss_2 <-
        estimate_costs(
          max_cr_ratio = x2,
          fish = fish,
          fleet = fleet,
          b_ref_oa = pmax(.05,(2 - f_v_m) / 2),
          lags = 10,
          num_patches = num_patches,
          sim_years = sim_years,
          burn_years = burn_years,
        )


      delta_best <-  (0 -  min(ss_1,ss_2))^2

      best <- min(ss_1,ss_2)

      if (ss_1 < ss_2){

        lower <- lower

        upper <- x2
      } else{

        lower <- x1

        upper <- upper

      }

      if (counter > 20){
        delta_best <- 0
      }

    } # close golden while


    crs <- seq(0.01,0.8, by = 0.05)

    huh <- map_dbl(crs, estimate_costs, fish = fish,
                   fleet = fleet,
                   b_ref_oa = pmax(.05,(2 - f_v_m) / 2),
                   lags = 1,
                   num_patches = num_patches,
                   sim_years = sim_years,
                   burn_years = burn_years
    )

    cost_fit <-
      nlminb(
        crs[huh == min(huh)],
        estimate_costs,
        fish = fish,
        fleet = fleet,
        b_ref_oa = pmax(.05,(2 - f_v_m) / 2),
        lower = 1e-3,
        upper = 0.8,
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

    # fleet$p_msy <- fish$price * fish$msy - fleet$cost * fleet$e_msy ^ fleet$beta


  } else if (fleet$fleet_model == "supplied-catch"){


    tol <- 100

    lower <- 0

    upper <- 400

    golden <- (sqrt(5) -1)/2

    best <- 1000

    delta_best <- 100

    counter <-  0

    while(delta_best > tol | counter < 20) {

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




