#' Function to run MLPA demon using a
#' delta/hurdle method
#'
#' \code{run_delta_demon} uses Laplaces Demon to
#' fit a bayesian hierarchechal model to the MLPA
#' data
#' @param dat dataframe of regression variables
#' @param dep_var the name of the dependent variable
#' in the regression
#' @param pos_vars vector of variables to be potentially included
#' in the regression
#' @param delta_vars vector of variables to be included in the
#' binomial part of the regression
#' @param iterations the number of runs for the MCMC
#' @param status the fraction of iterations at which to display status
#' @param thin the thinning rate for the MCMC
#' @param burn the amount of the chain to burn off
#' @param scale_numerics T or F to center and scale numeric covariates
#' @param runpath the folder to store results
#' @param acceptance_rate the target acceptance rate
#' @param method Summon Demon for standard MCMC, Summon Reversible Demon for
#' reversible MCMC
#'

run_delta_demon <- function(dat,dep_var,pos_vars,delta_vars,iterations = 1000,status = .05,thin = 1,burn = .5,
                            scale_numerics = F,
                            runpath, acceptance_rate  = 0.234, method = 'Summon Demon',
                            num_chains = 2) {


#   dat = reg_data #[1:6000,]
#   method = 'Jagged Demon'
#   dep_var = dep_var
#   pos_vars = pos_vars
#   delta_vars = delta_vars
#   runpath = runpath
#   scale_numerics = scale_numerics
#   iterations = its
#   status = .01
#   acceptance_rate = 0.43
#   thin = its/1e4


  # Convert data to regression format----
  observed_dat <- prep_demon(dat,pos_vars = pos_vars, scale_numerics = scale_numerics)

  binom_dat <- prep_demon(dat,pos_vars = delta_vars, scale_numerics = scale_numerics)

  vars_to_use_binom <- colnames(binom_dat)

  colnames(binom_dat) <- paste('bi',colnames(binom_dat), sep = '.')

  reg_dat <- cbind(observed_dat, binom_dat)

  J <- dim(reg_dat)[2]

  mon.names <- "LP"

  names <- colnames(reg_dat)

  beta_to_use_binom <- grepl('bi.', names, fixed = T)

  sigmas <- c('sigma_density','sigma_year','sigma_region','sigma_bi_year')

  time_vars <- names[grepl('_year.',names, fixed = T)]

  site_vars <- names[grepl('region.',names, fixed = T)]

  species_vars <- names(grepl('trophic.',names, fixed = T))

  parm.names <- c(colnames(reg_dat),sigmas)

  # Prepare reversible jump parameteres ----

  off_the_table <- (grepl('fished', parm.names, fixed = T) | grepl('mpa_applied', parm.names,fixed = T) |
                      grepl('fished_x_mpa', parm.names,fixed = T) | grepl('constant', parm.names,fixed = T) |
                      grepl('.factor', parm.names,fixed = T) | grepl('sigma', parm.names,fixed = T))

  selectable <- rep(1,length(parm.names))

  selectable[off_the_table] <- 0

  selected <- selectable

  bin.n <- J #Maximum allowable model size

  bin.p <- 0.9 #Most probable size:  bin.p x bin.n is binomial mean and median

  parm.p <- rep(1/J,J+1)

  # Prepare indices for MCMC ----

  vars_for_binom <- grepl('bi.', parm.names, fixed = T)

  pos_vars_for_binom <- which(vars_for_binom)

  parm <- rep(0,length(parm.names))

  pos_den_beta <- which(parm.names %in% colnames(reg_dat) & vars_for_binom == F)

  pos_bi_beta <- which(parm.names %in% colnames(reg_dat) & vars_for_binom == T)

  pos_den_sigma <- which(parm.names %in% sigmas & vars_for_binom == F)

  pos_den_time_terms <- which(parm.names %in% time_vars & vars_for_binom == F)

  pos_den_region_terms <- which(parm.names %in% site_vars & vars_for_binom == F)

  pos_bi_time_terms <- which(parm.names %in% time_vars & vars_for_binom == T)

  pos_any_betas <- which(grepl('sigma_',parm.names, fixed = T) == F)

  pos_any_sigmas <-  which(grepl('sigma_',parm.names, fixed = T))

  pos_sigma_year <- which(grepl('sigma_year',parm.names, fixed = T))

  pos_sigma_region <- which(grepl('sigma_region',parm.names, fixed = T))

  pos_sigma_density <- which(grepl('sigma_density',parm.names, fixed = T))

  pos_sigma_bi_year <- which(grepl('sigma_bi_year',parm.names, fixed = T))

  PGF <- function(Data) {
    beta <- rnorm(length(Data$pos_any_betas))
    sigma <- runif(length(Data$pos_any_sigmas))
    return(c(beta, sigma))
  }

  # Subset data to possible for regression ----

  possible <- matrix(NA,nrow = dim(reg_dat)[1], ncol = 2)

  possible[,1] <- is.na(as.matrix(reg_dat[,beta_to_use_binom == F]) %*% parm[pos_den_beta]) == F

  possible[,2] <- is.na(as.matrix(reg_dat[,beta_to_use_binom == T]) %*% parm[pos_bi_beta]) == F

  which_contrains <- which(colSums(possible) == min(colSums(possible)))[1]

  reg_dat <- reg_dat[possible[,which_contrains],]

  dat <- dat[possible[,which_contrains],]

  binom_dep_var <- dat$sites_seenat

  #   any_seen <- as.numeric(dat$sites_seenat > 0)

  any_seen <- as.numeric(dat[,dep_var] > min(dat[,dep_var]))

  N <- dim(dat)[1]

  # Prepare priors ----

  dep_sd <- sd(dat$log_density)

  # Prepare the Demon's snacks ----
  #

  pos_vars_for_binom <- which(vars_for_binom)

  pos_beta_to_use_binom <- which(beta_to_use_binom)

  bi_reg_mat = as.matrix(reg_dat[,beta_to_use_binom])

  den_reg_mat = as.matrix(reg_dat[,beta_to_use_binom == F])


  Data <- list(N = N,
               J=J,
               PGF=PGF,
               bi_reg_mat = bi_reg_mat,
               den_reg_mat = den_reg_mat,
               reg_dat = as.matrix(reg_dat),
               mon.names = mon.names,
               binom_dep_var = binom_dep_var,
               parm.names = parm.names,
               beta_to_use_binom = beta_to_use_binom,
               vars_for_binom = vars_for_binom,
               pos_beta_to_use_binom = pos_beta_to_use_binom,
               pos_den_beta = pos_den_beta,
               pos_den_sigma = pos_den_sigma,
               pos_den_time_terms = pos_den_time_terms,
               pos_den_region_terms = pos_den_region_terms,
               pos_any_betas = pos_any_betas,
               pos_any_sigmas = pos_any_sigmas,
               pos_sigma_year = pos_sigma_year,
               pos_sigma_density = pos_sigma_density,
               pos_sigma_region = pos_sigma_region,
               pos_sigma_bi_year = pos_sigma_bi_year,
               pos_bi_time_terms = pos_bi_time_terms,
               dep_var = as.matrix(dat[,dep_var]),
               time_vars = time_vars,
               site_vars = site_vars,
               species_vars = species_vars,
               sites_checked = dat$sites_checked,
               any_seen = any_seen)

  Initial.Values <- GIV(mlpa_delta_likelihood, Data, PGF=TRUE)

  dense_vars <- colnames(Data$reg_dat[,Data$pos_den_beta])

  binary_vars <- colnames(Data$reg_dat[,Data$pos_beta_to_use_binom])

  quickreg <- as.data.frame(Data$reg_dat)

  quickreg$mean_density <- exp(Data$dep_var)

  quickreg$mean_density[quickreg$mean_density == min(quickreg$mean_density)] <- 0

  #Run a quick delta lognormal glm for starting guesses
  density_fmla <- formula(paste('mean_density ~ -1 +',paste(dense_vars, collapse = '+'), sep = ''))

  logit_fmla <- formula(paste('~ -1 +',paste(binary_vars, collapse = '+'), sep = ''))

  delta_glm <- deltaLN(ln.form = density_fmla,binary.form = logit_fmla, data =quickreg)

  Initial.Values[Data$pos_den_beta] <- as.numeric(delta_glm$coefs$ln)

  Initial.Values[Data$pos_beta_to_use_binom] <- as.numeric(delta_glm$coefs$binary)

  vcov <- cov(reg_dat)

  full_vcov <- diag(x = .1, nrow = length(parm.names), ncol = length(parm.names))

  rownames(full_vcov) <- parm.names

  colnames(full_vcov) <- parm.names

  full_vcov[rownames(full_vcov) %in% colnames(reg_dat), colnames(full_vcov) %in% colnames(reg_dat) ] <- vcov


  #   jags_demon <-
  # Run Demon ----

  if (method == 'Jagged Demon')
  {

    den_dat <- Data$den_reg_mat #[Data$dep_var > min(Data$dep_var),]

    model_dat <- list(den_dat = den_dat, binom_dat = Data$bi_reg_mat,
                      n = dim(den_dat)[1], d = dim(Data$den_reg_mat)[2],
                      b = dim(Data$bi_reg_mat)[2],any_seen = any_seen,
                      den_dep_var = Data$dep_var, binom_dep_var = Data$binom_dep_var,
                      sites_checked = Data$sites_checked)
    #
    #     model_dat <- list(den_dat = den_dat, binom_dat = Data$bi_reg_mat,
    #                       n = dim(den_dat)[1], d = dim(Data$den_reg_mat)[2],
    #                       b = dim(Data$bi_reg_mat)[2],z = dim(Data$bi_reg_mat)[1], any_seen = any_seen,
    #                       den_dep_var = Data$dep_var, binom_dep_var = Data$binom_dep_var,
    #                       sites_checked = Data$sites_checked)
    #

    inits <- list(den_beta = as.numeric(delta_glm$coefs$ln),
                  binom_beta = as.numeric(delta_glm$coefs$binary),
                  sigma_density = 0.05, .RNG.name = "base::Super-Duper", .RNG.seed=1)

    model <-  'model{
    for (i in 1:n) {

    den_hat[i] <- (den_dat[i,] %*% den_beta)
    den_dep_var[i,] ~ dnorm(den_hat[i],sigma_density)
#     den_dep_var[i,] ~ ifelse( den_dep_var[i,] > 0,dnorm(den_hat[i],sigma_density),1)
#     den_dep_var[i,] <- den_dep_var[i,]^(any_seen[i])
#     }
# for (l in 1:z){
    binom_hat[i] <- min(10, binom_dat[i,] %*% binom_beta)
    prob_hat[i] <- exp(binom_hat[i])/(1 + exp(binom_hat[i]))
    binom_dep_var[i] ~ dbinom(prob_hat[i], sites_checked[i])
    }

sigma_density ~ dunif(1e-5,5)

for (j in 1:d){
den_beta[j] ~ dunif(-1000,1000)
}
for (k in 1:b){
binom_beta[k] ~dunif(-1000,1000)
}

} #close model'

# jagged_demon <- run.jags(model=model, monitor=c("den_beta",'binom_beta','sigma_density'),
#                          data=model_dat, n.chains=1, method="rjags", inits=inits,plots=F,monitor.deviance=F,
#                          silent.jag=F,burnin = 1000,thin = 10)
jagged_demon <- run.jags(model=model, monitor=c("den_beta",'binom_beta','sigma_density'),
                         data=model_dat, n.chains=1, method="rjags", inits=inits,plots=F,monitor.deviance=F,
                         silent.jag=F,modules=c("dic","glm","bugs"),
                         sample=1000,adapt=15000,burnin=150000,thin=500)

  }


if (method == 'Banish Demon')
{
  a <- proc.time()
  Fit <- mlpa_mcmc(par_init = Initial.Values,parm.names = parm.names,
                   dat = Data,vcov = full_vcov,n_sim  = iterations,
                   n_burn =  burn*iterations, targ_accept_rate = 0.25,
                   vcov_augment = (2.4/sqrt(length(parm.names)))^2, jumpyness = 1)
  show(proc.time() - a)
}

if (method == 'Summon Demon')
{
  #     a <- proc.time()
  Fit <- LaplacesDemon(mlpa_delta_likelihood, Data=Data, Initial.Values = Initial.Values,
                       Covar=NULL, Iterations=iterations, Status = iterations*status, Thinning=thin,
                       Algorithm = 'HARM', Specs=list(alpha.star = acceptance_rate, B = NULL),
                       parm.names = parm.names)
  #     show(proc.time() - a)
  #     browser()

}
if (method == 'Summon Parallel Demon')
{

  jitter_inits <- c(rep(Initial.Values,num_chains))

  disperse_mat <- matrix(runif(length(jitter_inits),-2,2),nrow = num_chains, ncol = length(Initial.Values))

  par_initial_values <- matrix(jitter_inits,nrow = num_chains,
                               ncol = length(Initial.Values), byrow = T) + disperse_mat

  Fit <- LaplacesDemon.hpc(mlpa_delta_likelihood, Data=Data, Initial.Values = par_initial_values,
                           Covar=NULL, Iterations=100, Status=iterations*status, Thinning=thin,
                           Algorithm = 'HARM', Specs=list(alpha.star=acceptance_rate, B = NULL),
                           Chains = num_chains, CPUs = num_chains)

}
if (method == 'Summon Reversible Demon'){

  Fit <- LaplacesDemon(mlpa_delta_likelihood, Data=Data, Initial.Values = Initial.Values,
                       Covar=NULL, Iterations=iterations, Status=iterations*status, Thinning=thin,
                       Algorithm = 'RJ', Specs=list(bin.n=bin.n, bin.p=bin.p,
                                                    parm.p=parm.p, selectable=selectable,
                                                    selected=selected), parm.names = parm.names)

}

return(list(demon_fit = Fit,Data = Data))
}
