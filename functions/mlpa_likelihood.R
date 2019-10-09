mlpa_likelihood <- function(parm, Data,reg_model = 'tobit')
{
  ### Parameters
  #         a <- proc.time()
  beta <- parm[Data$pos.beta]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)

  parm[Data$pos.sigma] <- sigma
  ### Log-Priors

  year_priors <- sum(dnorm(parm[Data$pos.time_terms],0,parm[Data$parm.names == 'sigma_year'], log = T))

  sigma_year_prior <- dnorm(parm[Data$parm.names == 'sigma_year'], 1,.2, log = T)

#   sigma_year_prior <- dhalfcauchy(parm[Data$parm.names == 'sigma_year'],25,log = T) #think through this

  region_priors <- sum(dnorm(parm[Data$parm.names %in% Data$site_vars ],
                             0,parm[Data$parm.names == 'sigma_region'], log = T))

  sigma_region_prior <- dnorm(parm[Data$parm.names == 'sigma_region'], .1,.2, log = T)

  sigma_density_prior <- dnorm(parm[Data$parm.names == 'sigma_density'], .1,.2, log = T)


  ### Log-Likelihood
  mu <- Data$reg_dat %*% beta

  y_star <- Data$dep_var

  #     if (reg_model == 'tobit')
  #     {

  min_val <- min(Data$dep_var)

  censored <- y_star <= min_val


#   thing <- NULL
#   for (i in 1:sum(censored))
#   {
#     thing[i] <- rtnorm(1,mean = mu[censored][i],sd = parm[Data$parm.names == 'sigma_density'], upper = min_val)
#   }
#   browser()

  y_star[censored] <- rtnorm(sum(censored),mean = mu[censored],sd = parm[Data$parm.names == 'sigma_density'], upper = min_val)
  # }

  LL <- sum(dnorm(y_star, mu, parm[Data$parm.names == 'sigma_density'], log=TRUE))
  ### Log-Posterior
  LP <- LL + year_priors + sigma_year_prior +  region_priors + sigma_region_prior + sigma_density_prior

  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)

  #         show(proc.time() - a)
  return(Modelout)
}
