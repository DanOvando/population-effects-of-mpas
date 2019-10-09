thin_mcmc <- function(chains, thin_every =  1)
{

  thinned_chains <- chains[(1:dim(chains)[1]) %% thin_every == 0,] #only select chains divisible by the thinning interval

  return(thinned_chains)
}