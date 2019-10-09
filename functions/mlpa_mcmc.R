#' Run MCMC on mosquitoes
#'
#' \code{Mozzy_MCMC} runs an MCMC on
#' mosquito count data
#' @param par_init initival vector of parameter guesses
#' @param varnames names of the parameters
#' @param dat are the raw stream-count-distance data
#' @param covar is the variance-covariance matrix of pars
#' @param n_sim number of simulations including burn in
#' @param n_burn number of simulations to burn in
#' @param n_thin the number of iterations to thin by
#' @param vcov_augment factor to multiply vcov by in the burn in period
#' @param prog_bar T or F to output progress bar
#' @param seed if you want to set a specific seed for this run
#' @param targ_accept_rate the target acceptance rate to tune the model
#' @param jumpyness parameter to modify acceptance ratio

mlpa_mcmc <- function(par_init,parm.names,dat,vcov,n_sim=1000,n_burn=0,n_thin=1,
                     vcov_augment = 1,prog_bar = T,jumpyness = 1, seed = NA,targ_accept_rate = 0.2)
{

  if (is.na(seed) == F)
  {
    set.seed(seed)
  }

  # Set things up ----

  par_curr <- par_init

  if (class(vcov) == 'data.frame'){ vcov <- as.matrix(vcov) }

  likelihoods <- mlpa_delta_likelihood(parm = par_curr, Data = dat)

  ll_curr <- likelihoods$LP

  deviance_curr <- likelihoods$Dev

  outvars <- c(tolower(parm.names),'ll','deviance')

  par_out <- matrix(NA, nrow = (n_sim - n_burn), ncol = length(outvars))

  i_count <- 0

  if(prog_bar == T)
  {
    pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  }

  #   a <- proc.time()


  sel_count <- 0

  orig_vcov <- vcov

  vcov_scalar <- vcov_augment

  track_rate <- rep(NA,n_sim)

  ### run MCMC----
  for (i in 1 : n_sim)
  {

    vcov <- orig_vcov * vcov_scalar

    par_next <- rmvnorm(n = 1, mean= par_curr)

    colnames(par_next) <- parm.names

    likelihoods <- mlpa_delta_likelihood(parm = par_next, Data = dat)

    deviance_next <- likelihoods$Dev

    ll_next <- likelihoods$LP

    if(prog_bar == T)
    {

      setTxtProgressBar(pb, i)

    }

    rand_accept <- log(runif(1,0,jumpyness))

    if (ll_next > ll_curr + rand_accept)
    {
      par_curr <- par_next
      ll_curr <- ll_next
      deviance_curr <- deviance_next
      sel_count <- sel_count + 1
    }

    track_rate[i] <- sel_count/i
#     a <- proc.time()
    if (i > n_burn & (i %% n_thin)==0) # store results
    {
      i_count <- i_count + 1
#       par_out[i_count,] <- c(par_curr,ll_curr,deviance_curr)
    } #else{
#     proc.time() - a
    if (i %% round(.2*n_burn,0) == 0 & i <= n_burn) { #Tune variance covariance

      acceptance_rate <- sel_count/i

      perc_off <- (acceptance_rate/targ_accept_rate) * 1.1
      vcov_scalar <- vcov_augment * perc_off

    }

  }

  # proc.time() - a
  par_out <- as.data.frame(par_out)

  colnames(par_out) <- outvars

  accepted_runs <- data.frame(num_selected = sel_count, num_simulated = n_sim, perc_selected = sel_count/n_sim,
                              perc_rejected = 1 - sel_count/n_sim, num_kept = n_sim - n_burn)

  return(list(posteriors = par_out[1:i_count,], initial_par = par_init,
              accepted_runs = accepted_runs))
}