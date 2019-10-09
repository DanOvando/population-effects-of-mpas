#' analyze demon object
#'
#' \code{process_demon} loads saved data
#' from run_demon and produces diagnostics
#' fits etc
#' @param runfolder the folder where results are stored
#' @param fontsize numeric fontsize to use
#' @param post_sample_size the desired final number of samples
#' from the MCMC after burning and thinning
#' @param burn proportion 0 to 1 of iterations to burn

process_demon <- function(runfolder,fontsize = 14,post_sample_size = 1000, burn = 0.6)
{

  runpath <- paste('results/',runfolder,'/', sep = '')

  load(paste(runpath,'reg_data.Rdata', sep = ''))

  load(paste(runpath,'MCMC results.Rdata', sep = ''))

  plot_theme <- theme_classic() + theme(text = element_text(size = fontsize, family = 'Helvetica'))

  post <- reg_results$post

  #   post <- bayes_reg$demon_fit$Posterior1

  its <- dim(post)[1]
  thinned_post <- thin_mcmc(post[(burn * its):its, ], thin_every = (its*(1-burn))/post_sample_size)
  ggmcmc(ggs(mcmc(thinned_post)), file = paste(runpath,'ggMCMC Diagnostics.pdf', sep = ''))

  predictions <- apply_demon(demonpost = thinned_post,
                             dat = reg_results$Data, raw_data = reg_data)


  # Make three/4 plots:
  # map/variogram of residuals by site/region
  # residuals by time
  # overlayed distributions of residuals by species group

  # Diagnostic plots ----

  resids_by_time_plot <- predictions$post_dat %>%
    subset(mean_density >0) %>%
    group_by(year) %>%
    mutate(mean_mean_resid = mean(mean_resid)) %>%
    ggplot(aes(factor(year),mean_resid)) +
    geom_boxplot() +
    geom_point(aes(factor(year),mean_mean_resid),
               shape = 17) +
    geom_hline(aes(yintercept = 0), color = 'red') +
    xlab('Year') +
    ylab('Residuals') +
    plot_theme

  resids_by_site_plot <- predictions$post_dat %>%
    subset(mean_density >0) %>%
    group_by(site) %>%
    mutate(mean_mean_resid = mean(mean_resid)) %>%
    ggplot(aes(factor(site),mean_resid, fill = factor(region))) +
    geom_boxplot(alpha = 0.85, notch = T) +
    geom_point(aes(factor(site),mean_mean_resid, fill = factor(region)),
               shape = 17) +
    scale_fill_discrete(name = 'Region') +
    geom_hline(aes(yintercept = 0), color = 'red') +
    plot_theme +
    ylab('Residuals') +
    xlab('') +
    plot_theme +
    coord_flip()

  resids_by_species_plot <- predictions$post_dat %>%
    subset(mean_density >0) %>%
    group_by(species) %>%
    mutate(mean_mean_resid = mean(mean_resid)) %>%
    ggplot(aes(factor(species),mean_resid, fill = factor(trophic.group))) +
    geom_boxplot(alpha = 0.85, notch = T) +
    geom_point(aes(factor(species),mean_mean_resid, fill = factor(trophic.group)),
               shape = 17) +
    scale_fill_discrete(name = 'Trophic Group') +
    geom_hline(aes(yintercept = 0), color = 'red') +
    plot_theme +
    ylab('Residuals') +
    xlab('') +
    plot_theme +
    coord_flip()

  resids_by_trophic_plot <- predictions$post_dat %>%
    subset(mean_density >0) %>%
    group_by(trophic.group) %>%
    mutate(mean_mean_resid = mean(mean_resid)) %>%
    ggplot(aes(factor(trophic.group),mean_resid)) +
    geom_boxplot(alpha = 0.85, notch = T) +
    geom_point(aes(factor(trophic.group),mean_mean_resid),
               shape = 17) +
    geom_hline(aes(yintercept = 0), color = 'red') +
    plot_theme +
    ylab('Residuals') +
    xlab('') +
    plot_theme +
    coord_flip()

  resids <- predictions$posterior %>%
    subset(obs_log_density > min(obs_log_density)) %>%
    group_by(observation) %>%
    summarise(mean_resid = mean(resid), mean_obs = mean(obs_log_density),mean_pred = mean(pred_log_density), mean_post_pred = mean(post_predict),
              median_post_pred = median(post_predict))

  resid_plot <- resids %>%
    ggplot(aes(mean_resid)) +
    geom_histogram(fill = 'steelblue4',color = 'black') +
    plot_theme +
    xlab('Residual') +
    ylab('Count')


  density_trend <- predictions$post_dat %>%
    subset(log_density > min(log_density)) %>%
    group_by(year,fished) %>%
    summarize(mean_density = mean(log_density), mean_predicted_density = mean(mean_pred_log_den)) %>%
    gather('data_source','density',mean_density:mean_predicted_density)


  obs_pre_trend_plot <- (ggplot() +
                           geom_point(data = filter(density_trend,data_source == 'mean_density'),aes(year,density,fill = factor(fished)), shape = 21) +
                           geom_line(data = filter(density_trend,data_source == 'mean_predicted_density'), aes(year,density,color = factor(fished))))


  density_trend <- reg_data %>%
    subset(log_density > min(log_density)) %>%
    group_by(year,fished) %>%
    summarize(mean_density = mean(log_density)) #%>%
#     gather('data_source','density',mean_density:mean_predicted_density)


  obs_pre_trend_plot <- (ggplot(density_trend) +
#                            geom_point(data = filter(density_trend,data_source == 'mean_density'),aes(year,density,fill = factor(fished)), shape = 21) +
                           geom_line(aes(year,mean_density,color = factor(fished))))


#   density_trend <- reg_data %>%
#     filter(log_density > min(log_density)) %>%
#     group_by(year,fished) %>%
#     summarize(mean_density = mean(log_density),
#               upper95 = quantile(log_density,.75),
#               lower95 = quantile(log_density,0.25) ) %>%
#     mutate(`Fished?` = 'Yes')
#
#   density_trend$`Fished?`[density_trend$fished == 0] <- 'No'


#   reg_data$`Fished?` <- 'Yes'
#
#   reg_data$`Fished?`[reg_data$fished == 0] <- 'No'
#

#   density_plot <- ggplot(density_trend,aes(year,mean_density,color = `Fished?`, fill = `Fished?`)) +
#     #     geom_smooth(span = .6,size = 2) +
#     geom_line(size = 2) +
#     #     geom_ribbon(aes(x = year, ymax = upper95, ymin = lower95, fill = `Fished?`)) +
#     geom_vline(aes(xintercept = 2003), linetype = 'longdash', color = 'black', size = 1.5) +
#     xlab('Year') +
#     ylab('Mean Log Density') +
#     plot_theme

  qq_plot <- resids %>%
    ggplot(aes(sample = mean_resid)) +
    stat_qq(shape = 21, size = 4, alpha = 0.6, fill = 'steelblue4') +
    plot_theme +
    xlab('Theoretical') +
    ylab('Sample')

  misspec_plot <- resids %>%
    ggplot(aes(mean_pred,mean_resid)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +   geom_hline(aes(yintercept = 0), linetype = 'longdash') +
    geom_smooth(method = 'lm', color = 'red', se = F) +
    plot_theme +
    xlab('Predicted Density') +
    ylab('Residual')


  lm_eqn = function(m) {

    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              r2 = format(summary(m)$r.squared, digits = 3));

    if (coef(m)[2] >= 0)  {
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    } else {
      eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    }

    as.character(as.expression(eq));
  }

  fit_plot <- resids %>%
    ggplot(aes(mean_obs,mean_pred)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_smooth(method = 'lm', color = 'red',linetype = 'longdash') +
    plot_theme +
    xlab('Observed') +
    ylab('Predicted')
  #+ geom_text(x = -8,y= -3, label = lm_eqn(lm(mean_pred ~ mean_obs,resids )), parse = T)

  post_pred_plot <- resids %>%
    ggplot(aes(mean_obs,median_post_pred)) +
    geom_point(shape = 21, fill = 'steelblue4', size = 3, alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_smooth(method = 'lm', color = 'red',linetype = 'longdash') +
    plot_theme +
    xlab('Observed') +
    ylab('Median Posterior Predicted')

  r2 <- round(summary(lm(median_post_pred ~ median_post_pred,resids ))$r.squared, digits = 2)

  gweke_plot <- ggs_geweke(ggs(mcmc(thinned_post))) +
    plot_theme

  eff_sample_plot <- data.frame(effective.ss = effectiveSize(mcmc(thinned_post))) %>%
    mutate(parameter = rownames(.)) %>%
    subset(!parameter %in% c('ll','deviance')) %>%
    ggplot(aes(factor(parameter),effective.ss)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    coord_flip() +
    xlab('Parameter') +
    ylab('Effective Sample Size') +
    plot_theme


  lag_1 <- rep(NA,dim(thinned_post)[2])
  for (j in 1:dim(thinned_post)[2])
  {
    lag_1[j] <- acf(thinned_post[,j], plot = F)$acf[2]
  }


  bino_plot <- predictions$post_prob_zero %>%
    group_by(observation) %>%
    summarize(is_zero = mean(is_zero), prob_zero = mean(prob_zero)) %>%
    mutate(bin = cut(prob_zero,5, dig.lab = 2)) %>%
    group_by(bin) %>%
    summarize(seen = mean(is_zero)) %>%
    ggplot(aes(bin,seen)) +
    geom_bar(stat = 'identity', position = 'dodge',
             fill = 'steelblue4', color = 'black') +
    plot_theme +
    xlab('Proportion Observed > 0') +
    ylab('Mean P(D)')


  acf_hist_plot <- ggplot(data.frame(acf = lag_1), aes(acf)) +
    geom_histogram(binwidth = .05, fill = 'steelblue4', color = 'black') +
    xlab('Lag 1 ACF') +
    ylab('Count') +
    plot_theme

  crosscor_plot <- ggs_crosscorrelation(ggs(mcmc(thinned_post))) +
    plot_theme

  temperature_effect_plot <- as.data.frame(thinned_post) %>%
    select(contains('temp')) %>%
    select(-contains('bi'))  %>%
    #     rename('mean_temp_lag0' = na_temp) %>%
    gather('lag_name','coefficient') %>%
    mutate( lag = str_match(lag_name, paste('lag', '(.+)', sep=''))[,2]) %>%
    group_by(lag) %>%
    ggplot(aes(lag,coefficient)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = 0)) +
    plot_theme +
    xlab('Temperature Lag') +
    ylab('Coefficient')

  mpa_time_effect_plot <- as.data.frame(thinned_post) %>%
    select(contains('fishedeffect')) %>%
    select(-contains('bi'))  %>%
    gather('year_term','coefficient') %>%
    mutate( year = str_match(year_term, paste('_', '(.+)', '_', sep=''))[,2]) %>%
    group_by(year) %>%
    ggplot(aes(year,coefficient)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = 0)) +
    plot_theme +
    xlab('Year') +
    ylab('Coefficient')

  outs_of_interest_plot <- as.data.frame(thinned_post) %>%
    dplyr::select(fished,mpa_applied) %>%
    rename('Species Fished' = fished,'MPAs Applied' = mpa_applied) %>%
    gather('Variable','Coefficient') %>%
    group_by(Variable) %>%
    mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
    ggplot(aes(Coefficient, fill = Variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_histogram(color = 'black') +
    geom_vline(aes(xintercept = lower95),alpha= 0.75) +
    geom_vline(aes(xintercept = upper95),alpha= 0.75)+
    geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
    facet_grid(.~Variable, scales = 'free') +
    theme_light() +
    theme(text = element_text(size = 12),strip.text.y = element_text(size = 11)) +
    ylab('Density')

    density_summary_plot <-  arrangeGrob(grid.arrange(mpa_time_effect_plot,outs_of_interest_plot, newpage = F,nrow = 2,ncol = 1))

  outs_of_interest_sigma_plot <- as.data.frame(thinned_post) %>%
    dplyr::select(contains('sigma')) %>%
    #     rename(Fished = fished,'Years MLPA' = years_mlpa_mpas, 'Fished X Years MLPA' = fished_x_yearsmlpa) %>%
    gather('Variable','Coefficient') %>%
    group_by(Variable) %>%
    mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
    ggplot(aes(Coefficient, fill = Variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_density() +
    geom_vline(aes(xintercept = lower95),alpha = 0.75) +
    geom_vline(aes(xintercept = upper95),alpha = 0.75) +
    geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
    facet_grid(Variable~., scales = 'free') +
    theme_light() +
    theme(text = element_text(size = 12), strip.text.y = element_text(size = 11)) +
    ylab('Density')

  outs_of_interest_binomial_plot <- as.data.frame(thinned_post) %>%
    dplyr::select(bi.fished,bi.mpa_applied) %>%
    rename('Species Fished' = bi.fished,'MPA Applied' = bi.mpa_applied) %>%
    gather('Variable','Coefficient') %>%
    group_by(Variable) %>%
    mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
    ggplot(aes(Coefficient, fill = Variable)) +
    scale_fill_brewer(guide = F, palette = 'Spectral') +
    geom_density() +
    geom_vline(aes(xintercept = lower95),alpha = 0.75) +
    geom_vline(aes(xintercept = upper95),alpha = 0.75) +
    geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
    facet_grid(Variable~.) +
    theme_light() +
    theme(text = element_text(size = 12),strip.text.y = element_text(size = 11)) +
    ylab('Density')




  #   outs_of_interest_plot <- as.data.frame(thinned_post) %>%
  #     dplyr::select(fished,years_mlpa_mpas,fished_x_yearsmlpa) %>%
  #     rename(Fished = fished,'Years MLPA' = years_mlpa_mpas, 'Fished X Years MLPA' = fished_x_yearsmlpa) %>%
  #     gather('Variable','Coefficient') %>%
  #     group_by(Variable) %>%
  #     mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
  #     ggplot(aes(Coefficient, fill = Variable)) +
  #     scale_fill_brewer(guide = F, palette = 'Spectral') +
  #     geom_density() +
  #     geom_vline(aes(xintercept = lower95),alpha= 0.75) +
  #     geom_vline(aes(xintercept = upper95),alpha= 0.75)+
  #     geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
  #     facet_grid(Variable~., scales = 'free') +
  #     theme_light() +
  #     theme(text = element_text(size = 12),strip.text.y = element_text(size = 11)) +
  #     ylab('Density')
  #
  #   outs_of_interest_sigma_plot <- as.data.frame(thinned_post) %>%
  #     dplyr::select(contains('sigma')) %>%
  # #     rename(Fished = fished,'Years MLPA' = years_mlpa_mpas, 'Fished X Years MLPA' = fished_x_yearsmlpa) %>%
  #     gather('Variable','Coefficient') %>%
  #     group_by(Variable) %>%
  #     mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
  #     ggplot(aes(Coefficient, fill = Variable)) +
  #     scale_fill_brewer(guide = F, palette = 'Spectral') +
  #     geom_density() +
  #     geom_vline(aes(xintercept = lower95),alpha = 0.75) +
  #     geom_vline(aes(xintercept = upper95),alpha = 0.75) +
  #     geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
  #     facet_grid(Variable~., scales = 'free') +
  #     theme_light() +
  #     theme(text = element_text(size = 12), strip.text.y = element_text(size = 11)) +
  #     ylab('Density')
  #
  #   outs_of_interest_binomial_plot <- as.data.frame(thinned_post) %>%
  #     dplyr::select(bi.fished,bi.years_mlpa_mpas,bi.fished_x_yearsmlpa) %>%
  #     rename(Fished = bi.fished,'Years MLPA' = bi.years_mlpa_mpas, 'Fished X Years MLPA' = bi.fished_x_yearsmlpa) %>%
  #     gather('Variable','Coefficient') %>%
  #     group_by(Variable) %>%
  #     mutate(lower95 = quantile(Coefficient, 0.025), upper95 = quantile(Coefficient, 0.975)) %>%
  #     ggplot(aes(Coefficient, fill = Variable)) +
  #     scale_fill_brewer(guide = F, palette = 'Spectral') +
  #     geom_density() +
  #     geom_vline(aes(xintercept = lower95),alpha = 0.75) +
  #     geom_vline(aes(xintercept = upper95),alpha = 0.75) +
  #     geom_vline(aes(xintercept = 0), linetype = 'longdash', color = 'red')+
  #     facet_grid(Variable~.) +
  #     theme_light() +
  #     theme(text = element_text(size = 12),strip.text.y = element_text(size = 11)) +
  #     ylab('Density')

  #   year_trend_plot <- as.data.frame(thinned_post) %>%
  #     select(contains('factor_year')) %>%
  #     select(which(!grepl('bi.',colnames(.), fixed = T))) %>%
  #     gather('factor_year','Coef') %>%
  #     mutate(Year = as.numeric(factor_year) + 2001 -1) %>%
  #     ggplot(aes(factor(Year),Coef)) +
  #     geom_boxplot(fill = 'steelblue4', color = 'black') +
  #     geom_hline(aes(yintercept = 0))+
  #     plot_theme +
  #     theme(axis.text.x = element_text(size = 12)) +
  #     xlab('Year') +
  #     ylab('Effect Relative to 2000')

  local_files <- ls()

  plot_files <- local_files[grep('_plot', local_files, fixed = T)]

  not_plot_files <- local_files[!grepl('_plot', local_files, fixed = T)]


  plot_list <- list()

  for (i in 1:length(plot_files))
  {
    eval(parse( text = paste('plot_list$',plot_files[i],' <- ',plot_files[i], sep = '')))
  }

  predicted_data <- predictions$post_dat

  drop <- local_files[!local_files %in% c('plot_list','resids','thinned_post','predicted_data','predictions')]

  rm(list = drop)

  return(list(plot_list = plot_list, resids = resids, thinned_post = thinned_post,predictions = predicted_data, post_pred = predictions))
}