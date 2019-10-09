# mlpa_delta_likelihood <- function(parm, Data,reg_model = 'tobit')
# {

'model{
  for (i in 1:n){

    den_hat[n] <- den_dat[n,1:d] %*% den_beta

    dep_var ~ dnorm(den_hat[n],sigma_density)^(any_seen)

    binom_hat[n] <- min(10, binom_dat[n,1:b] %*% binom_beta)

    prob_hat[n] <- exp(binom_hat[n])/(1 + exp(binom_hat[n]))

    binom_dep_var ~ dbinom(sites_checked,max(1e-15,prob_hat))

  }
}'
