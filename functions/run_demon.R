#' Function to run MLPA demon
#'
#' \code{run_mlpa_demon} uses Laplaces Demon to
#' fit a bayesian hierarchechal model to the MLPA
#' data
#'

run_mlpa_demon <- function(dat,dep_var,pos_vars,iterations = 1000,status = .05,thin = 1,burn = .5,
                           scale_numerics = F,
                           runpath, acceptance_rate  = 0.234, method = 'Summon Demon') {

  observed_dep_var <- dat[,dep_var] > min(dat[,dep_var])

  reg_dat <- prep_demon(dat,pos_vars = pos_vars, scale_numerics = scale_numerics)

  binom_dat <- prep_demon(dat,pos_vars = pos_vars, scale_numerics = scale_numerics)

  reg_frame <- data.frame(dat[,dep_var],reg_dat)

  fmla <- as.formula(paste("log_density ~ ", paste(colnames(reg_dat), collapse= "+"),'- 1'))

  tobit_reg <- tobit(fmla,
                     data = reg_frame,
                     left = min(dat[,dep_var], na.rm = T))

  vcov <- tobit_reg$var

  vcov <- vcov[rownames(vcov) != 'Log(scale)', colnames(vcov) != 'Log(scale)']

  J <- dim(reg_dat)[2]

  names <- colnames(reg_dat)

  sigmas <- c('sigma_density','sigma_year','sigma_region')

  #   mus <- c('mu_year')

  time_vars <- names[grepl('_year.2',names)]

  site_vars <- names[grepl('region.',names)]

  species_vars <- names(grepl('trophic.',names))

  mon.names <- "LP"

  parm.names <- c(colnames(reg_dat),sigmas)

  full_vcov <- diag(x = .1, nrow = length(parm.names), ncol = length(parm.names))

  rownames(full_vcov) <- parm.names

  colnames(full_vcov) <- parm.names

  full_vcov[rownames(full_vcov) %in% colnames(reg_dat), colnames(full_vcov) %in% colnames(reg_dat) ] <- vcov

  parm <- rep(0,length(parm.names))

  pos.beta <- which(parm.names %in% colnames(reg_dat))

  pos.sigma <- which(parm.names %in% sigmas)

  pos.time_terms <- which(parm.names %in% time_vars)

  PGF <- function(Data) {
    beta <- rnorm(length(Data$pos.beta))
    sigma <- runif(length(Data$pos.sigma))
    return(c(beta, sigma))
  }

  possible <- is.na(as.matrix(reg_dat) %*% parm[pos.beta]) == F

  reg_dat <- reg_dat[possible,]

  dat <- dat[possible,]

  N <- dim(dat)[1]

  Data <- list(N = N,J=J, PGF=PGF, reg_dat = as.matrix(reg_dat), mon.names=mon.names,
               parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma,pos.time_terms = pos.time_terms,
               dep_var = as.matrix(dat[,dep_var]),time_vars = time_vars,site_vars = site_vars,
               species_vars = species_vars)
  Initial.Values <- GIV(mlpa_likelihood, Data, PGF=TRUE)

  if (method != 'Summon Demon')
  {
    Fit <- mlpa_mcmc(par_init = Initial.Values,parm.names = parm.names,
                     dat = Data,vcov = full_vcov,n_sim  = iterations,
                     n_burn =  burn*iterations, targ_accept_rate = 0.25,
                     vcov_augment = (2.4/sqrt(45))^2, jumpyness = 1)
  }else{

    Fit <- LaplacesDemon(mlpa_likelihood, Data=Data, Initial.Values = Initial.Values,
                         Covar=NULL, Iterations=iterations, Status=iterations*status, Thinning=1,
                         Algorithm = 'HARM', Specs=list(alpha.star=acceptance_rate, B = NULL))
  }

  return(Fit)
}
