#' Function to estimate predicted values from
#'
#' \code{apply_demon} applies the posterior from
#' a Laplaces Demon call to a supplied and appropriate
#' dataframe to obtain predictions and posterior predictions
#' @param demonpost the burned and thinned posterior
#' @param demon an object of class demonoid produced by
#' LaplacesDemon
#' @param dat the dataframe associated with the
#' demonoid object

apply_demon <- function(demonpost,dat, raw_data)
{

#     demonpost = thinned_post
#
#     demon <- bayes_reg$demon_fit
#
#     dat <- bayes_reg$Data

    rows_used <- as.numeric(rownames(dat$reg_dat))

    odat <- data.frame(observation = rows_used,
                       raw_data[rows_used,], stringsAsFactors = F)

  # Obtain predicted values ----

  betas <- demonpost[,dat$pos_den_beta]

  ind_vars <- dat$reg_dat[,dat$beta_to_use_binom == F]

  predicted_density <- as.data.frame(t(betas %*% t(ind_vars)))
  predicted_density$obs_log_density <- dat$dep_var

  sigma_density <- as.data.frame(demonpost[,'sigma_density'])

  colnames(sigma_density) <- 'sigma_density'

  sigma_density$chain <- 1:dim(sigma_density)[1]

  predicted_density$observation <- rows_used

  posterior <- predicted_density %>%
    gather('chain','pred_log_density',which(grepl('V', colnames(.), fixed = T))) %>%
    mutate(chain = as.numeric(gsub('V', '',chain))) %>%
    left_join(sigma_density, by = 'chain') %>%
    mutate(post_predict = rnorm(length(pred_log_density),pred_log_density, sigma_density),
           resid = pred_log_density - obs_log_density)

#   subset(posterior, observation <40) %>%
#     ggplot(aes(post_predict)) +
#     geom_histogram() +
#     geom_vline(aes(xintercept = obs_log_density)) +
#     facet_wrap(~observation)

  # Obtain predicted zero probs ---------------------------------------------

  betas <- demonpost[,dat$beta_to_use_binom]

  ind_vars <- dat$reg_dat[,dat$beta_to_use_binom == T]

  prob_zero <- as.data.frame(t(betas %*% t(ind_vars)))

  prob_zero <- exp(prob_zero) / (1+exp(prob_zero))

  prob_zero$log_density <- dat$dep_var

  prob_zero$is_zero <- prob_zero$log_density > min(prob_zero$log_density)

  prob_zero$observation <- rows_used

  post_prob_zero <- prob_zero %>%
    gather('chain','prob_zero',which(grepl('V', colnames(.), fixed = T))) %>%
    mutate(chain = as.numeric(gsub('V', '',chain)))

  mean_prob_zero <- post_prob_zero %>%
    group_by(observation) %>%
    summarise(mean_prob_zero = mean(prob_zero))


  # Place mean observations in data -----------------------------------------------------------

  post_dat <- posterior %>%
    group_by(observation) %>%
    summarise(mean_pred_log_den = mean(pred_log_density),
              mean_sigma_den = mean(sigma_density),mean_post_pred_log_den = mean(post_predict),
              mean_resid = mean(resid)) %>%
    right_join(odat, by = 'observation') %>%
    left_join(mean_prob_zero, by = 'observation')


  return(list(posterior = posterior, post_prob_zero = post_prob_zero, post_dat = post_dat))
}