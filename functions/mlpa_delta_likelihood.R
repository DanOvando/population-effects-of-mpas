mlpa_delta_likelihood <- function(parm, Data,reg_model = 'tobit')
{
  ### Separate out parameters for observed density  ----
  beta <- parm[Data$pos_den_beta]

  parm[Data$pos_den_sigma] <- interval(parm[Data$pos_den_sigma], 1e-100, Inf)

  #   parm[Data$pos_den_sigma] <- exp(parm[Data$pos_den_sigma])

  sigma_year <- parm[Data$pos_sigma_year]

  sigma_bi_year <- parm[Data$pos_sigma_bi_year]

  sigma_region <- parm[Data$pos_sigma_region]

  sigma_density <- parm[Data$pos_sigma_density]

  ### Log-Priors for observed ----

  year_priors <- sum(dnorm(parm[Data$pos_den_time_terms],0,sigma_year, log = T))

  bi_year_priors <- sum(dnorm(parm[Data$pos_bi_time_terms],0,sigma_bi_year, log = T))


  sigma_year_prior <- dgamma(sigma_year, 2,0.5, log = T)

  sigma_bi_year_prior <- dgamma(sigma_bi_year, 2,0.5, log = T)

  region_priors <- sum(dnorm(parm[Data$pos_den_region_terms],
                             0,sigma_region, log = T))

  sigma_region_prior <- dgamma(sigma_region, 2,.5, log = T)

  sigma_density_prior <- dgamma(sigma_density, 2,.5, log = T)


  ### Hurdle Log-Likelihood ----

  bi_beta <- parm[Data$beta_to_use_binom]

  bi_dat <- Data$bi_reg_mat

  bi_hat <- pmin(10,bi_dat %*% bi_beta)
# browser()
#   bi_hat <- pmin(10,rowSums(bi_dat * bi_beta))
#

  prob_hat <- exp(bi_hat)/(1 + exp(bi_hat))

  bi_loglike <-  sum(dbinom(Data$binom_dep_var,Data$sites_checked,pmax(1e-15,prob_hat), log = T)) # actual values of 1 return -inf if not fulfilled

#   bi_loglike <-  sum(dbinom(Data$binom_dep_var,1,pmax(1e-15,prob_hat), log = T)) # actual values of 1 return -inf if not fulfilled

#   wtf <- bi_hat[is.na(bi_loglike)]

  ### Density Log-likelihood ----

  mu <- Data$den_reg_mat %*% beta

#   observed_density <- Data$dep_var
  density_loglike <- sum(dnorm(Data$dep_var, mu,sigma_density, log=TRUE)^(Data$any_seen))

  ### Log-Posterior

  LP <- density_loglike + bi_loglike  + year_priors + bi_year_priors +
    sigma_year_prior + sigma_bi_year_prior + region_priors + sigma_region_prior + sigma_density_prior

#   if (is.finite(LP) == F){browser()}

  LL <- density_loglike + bi_loglike

  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=LP,yhat = 1,
                   parm=parm)

  return(Modelout)
}
