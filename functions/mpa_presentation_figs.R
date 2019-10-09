mpa_presentation_figs <- function(run, out,font_size = 18, font_family = 'Helvetica',
                                  base_theme = theme_grey(), fig_width = 10, fig_height = 8)
{


  plot_theme <- base_theme +
    theme(text = element_text(size = font_size, family = font_family),
          axis.title.y = element_text(angle = 0, hjust = 0))

  runfolder <- run

  runpath <- paste('results/',runfolder,'/', sep = '')

  out_folder <- paste(runpath,out,'/',sep = '')


  if (dir.exists(out_folder) == F)
  {
    dir.create(out_folder, recursive = T)
  }

  load(paste(runpath,'reg_data.Rdata', sep = ''))

  load(paste(runpath,'MCMC results.Rdata', sep = ''))

  load(paste(runpath,'Processed MCMC results.Rdata', sep = ''))

  .simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
  }

  reg_data$species <- sapply(reg_data$species, .simpleCap, USE.NAMES = F)

  species_summary <- reg_data %>%
    group_by(species) %>%
    summarize(observations = length(species),fished = unique(fished == 1)) %>%
    mutate(obs_order = order(observations),
           species = factor(species, levels = species[obs_order]))

  species_summary$fished[species_summary$fished == TRUE] <- 'Yes'

  species_summary$fished[species_summary$fished == FALSE] <- 'No'


  density_plot <- ggplot(reg_data, aes(log_density)) +
    geom_histogram(fill = 'lightseagreen', color = 'black') +
    ylab('Log Density (tons/hectare') +
    xlab('# of Observations') +
    plot_theme

  ggsave(filename = paste(out_folder,'Density Summary Plot.pdf', sep = ''), density_plot, width = fig_width,
         height = fig_height)

  species_plot <- ggplot(species_summary,aes(species,observations, fill = fished)) +
    geom_bar(stat = 'identity', color = 'black') +
    coord_flip() +
    ylab('Observations') +
    xlab('') +
  scale_fill_discrete(name = 'Fished?') +
  plot_theme

  ggsave(filename = paste(out_folder,'Species Summary Plot.pdf', sep = ''), species_plot, width = fig_width,
         height = fig_height)

  density_trend <- reg_data %>%
    filter(log_density > min(log_density)) %>%
    group_by(year,fished) %>%
    summarize(mean_density = mean(log_density),
              upper95 = quantile(log_density,.75),
              lower95 = quantile(log_density,0.25) ) %>%
    mutate(`Fished?` = 'Yes')

  density_trend$`Fished?`[density_trend$fished == 0] <- 'No'


  reg_data$`Fished?` <- 'Yes'

  reg_data$`Fished?`[reg_data$fished == 0] <- 'No'


  density_plot <- ggplot(density_trend,aes(year,mean_density,color = `Fished?`, fill = `Fished?`)) +
#     geom_smooth(span = .6,size = 2) +
    geom_line(size = 2) +
#     geom_ribbon(aes(x = year, ymax = upper95, ymin = lower95, fill = `Fished?`)) +
    geom_vline(aes(xintercept = 2003), linetype = 'longdash', color = 'black', size = 1.5) +
    xlab('Year') +
    ylab('Density (log)') +
    plot_theme


  ggsave(filename = paste(out_folder,'Density Trend Plot.pdf', sep = ''), density_plot, width = fig_width,
         height = fig_height)

  density_trend2 <- reg_data %>%
#     filter(log_density > min(log_density)) %>%
    group_by(year,fished) %>%
    summarize(mean_density = mean(log_density),
              upper95 = quantile(log_density,.75),
              lower95 = quantile(log_density,0.25) ) %>%
    mutate(`Fished?` = 'Yes')

  density_trend2$`Fished?`[density_trend2$fished == 0] <- 'No'


  density_plot2 <- ggplot(density_trend2,aes(year,mean_density,color = `Fished?`, fill = `Fished?`)) +
    #     geom_smooth(span = .6,size = 2) +
    geom_line(size = 2) +
    #     geom_ribbon(aes(x = year, ymax = upper95, ymin = lower95, fill = `Fished?`)) +
    geom_vline(aes(xintercept = 2003), linetype = 'longdash', color = 'black', size = 1.5) +
    xlab('Year') +
    ylab('Density (log)') +
    plot_theme


  ggsave(filename = paste(out_folder,'Density Trend Plot 2.pdf', sep = ''), density_plot2, width = fig_width,
         height = fig_height)



  temp_trend_plot <- reg_data %>%
  group_by(year) %>%
  mutate(mmt = mean(mean_temp)) %>%
  ggplot(aes(factor(year),mean_temp)) +
  geom_point(alpha = 0.75,size = 4,shape = 21, aes(fill = mean_temp)) +
  geom_smooth(method = 'lm', aes(group = 1), color = 'red') +
  scale_fill_gradient(low = 'steelblue1', high = 'orangered', guide = F) +
    ylab(expression(paste(degree ~ C, sep = ''))) +
  xlab('Year') +
#   geom_hline(aes(yintercept = mean(mean_temp)), linetype = 'longdash')
  plot_theme

ggsave(filename = paste(out_folder,'Temp Trend Plot.pdf', sep = ''), temp_trend_plot, width = fig_width,
       height = fig_height)

  temperature_effect_plot <- as.data.frame(processed_demon$thinned_post) %>%
    select(contains('temp')) %>%
    select(-contains('bi'))  %>%
#     rename('mean_temp_lag0' = na_temp) %>%
    gather('lag_name','coefficient') %>%
    mutate( lag = str_match(lag_name, paste('lag', '(.+)', sep=''))[,2]) %>%
    group_by(lag) %>%
    mutate(mean_effect = mean(coefficient)) %>%
    ggplot(aes(lag,coefficient)) +
    geom_boxplot(aes(fill = mean_effect))+
    geom_hline(aes(yintercept = 0)) +
    plot_theme +
    xlab('Temperature Lag') +
    ylab('Effect on Densities') +
    scale_fill_gradient2(guide = F,low = 'red',mid = 'white',high = 'green',midpoint = 0) +
    geom_hline(yintercept = 0, size = 2, linetype = 'longdash')

  ggsave(filename = paste(out_folder,'Temperature Effect.pdf', sep = ''), temperature_effect_plot, width = fig_width,
         height = fig_height)

  mpa_time_effect_plot <- as.data.frame(processed_demon$thinned_post) %>%
    select(contains('fishedeffect')) %>%
    select(-contains('bi'))  %>%
    gather('year_term','coefficient') %>%
    mutate( year = str_match(year_term, paste('_', '(.+)', '_', sep=''))[,2]) %>%
    group_by(year) %>%
    mutate(mean_effect = mean(coefficient)) %>%
    ggplot(aes(year,coefficient)) +
    geom_boxplot(aes(fill = mean_effect)) +
    plot_theme +
    xlab('Year') +
    ylab('Effect on Densities') +
    scale_fill_gradient2(guide = F,low = 'red',mid = 'white',high = 'green',midpoint = 0) +
    geom_hline(yintercept = 0, size = 2, linetype = 'longdash')

  mpa_time_effect2_plot <- as.data.frame(processed_demon$thinned_post) %>%
    select(contains('fishedeffect')) %>%
    select(-contains('bi'))  %>%
    gather('year_term','coefficient') %>%
    mutate( year = str_match(year_term, paste('_', '(.+)', '_', sep=''))[,2]) %>%
    group_by(year) %>%
    mutate(mean_effect = mean(coefficient)) %>%
    ggplot(aes(year,coefficient)) +
    geom_violin(aes(fill = mean_effect)) +
    plot_theme +
    xlab('Year') +
    ylab('Effect on Densities') +
    scale_fill_gradient2(guide = F,low = 'red',mid = 'white',high = 'green',midpoint = 0) +
    geom_hline(yintercept = 0, size = 2, linetype = 'longdash')

  ggsave(filename = paste(out_folder,'MPA Effect.pdf', sep = ''), mpa_time_effect_plot, width = fig_width,
         height = fig_height)

  ggsave(filename = paste(out_folder,'MPA Effect 2.pdf', sep = ''), mpa_time_effect2_plot, width = fig_width,
         height = fig_height)
  bimpa_time_effect_plot <- as.data.frame(processed_demon$thinned_post) %>%
    select(contains('fishedeffect')) %>%
    select(contains('bi'))  %>%
    gather('year_term','coefficient') %>%
    mutate( year = str_match(year_term, paste('_', '(.+)', '_', sep=''))[,2]) %>%
    group_by(year) %>%
    mutate(mean_effect = mean(coefficient)) %>%
    ggplot(aes(year,coefficient)) +
    geom_boxplot(aes(fill = mean_effect)) +
    plot_theme +
    xlab('Year') +
    ylab('Effect on P(D>0)') +
    scale_fill_gradient2(guide = F,low = 'red',mid = 'white',high = 'green',midpoint = 0) +
    geom_hline(yintercept = 0, size = 2, linetype = 'longdash')

  ggsave(filename = paste(out_folder,'Logit MPA Effect.pdf', sep = ''), bimpa_time_effect_plot, width = fig_width,
         height = fig_height)


  den_post_preds <- processed_demon$post_pred$posterior %>%
    subset(obs_log_density > min(obs_log_density))

  post_summary <- den_post_preds %>%
    group_by(chain) %>%
    summarize(min_pp = min(post_predict), max_pp = max(post_predict), mean_pp = mean(post_predict),
              sd_pp = sd(post_predict), min_real = min(obs_log_density), max_real = max(obs_log_density), mean_real = mean(obs_log_density),
              sd_real = sd(obs_log_density))

  reals <- gather(post_summary,'statistic','real',min_real:sd_real) %>%
    mutate(stat_name = gsub('\\_.*', '', statistic))

  post_summary_plot <- post_summary %>%
    gather('statistic','value',min_pp:sd_pp) %>%
    group_by(statistic) %>%
    mutate(lower95 = quantile(value,0.025), upper95 = quantile(value,0.975)) %>%
    ungroup() %>%
    mutate(stat_name = gsub('\\_.*', '', statistic))%>%
    left_join(select(reals,stat_name,real), by = 'stat_name') %>%
    #group_by(statistic) %>%
    ggplot(aes(value)) +
    geom_histogram() +
    geom_vline(aes(xintercept = lower95)) +
    geom_vline(aes(xintercept = upper95)) +
    geom_vline(aes(xintercept = real), linetype = 'longdash', color = 'red') +
    facet_wrap(~stat_name, scales = 'free_x') +
    plot_theme +
    xlab('Value') +
    ylab('Count')


  ggsave(filename = paste(out_folder,'Post Predict Diagnostics.pdf', sep = ''), post_summary_plot, width = fig_width,
         height = fig_height)



  lm_eqn = function(m) {

    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              r2 = format(summary(m)$r.squared, digits = 3));

    if (coef(m)[2] >= 0)  {
      eq <- substitute(italic(r)^2~"="~r2,l)
    } else {
#       eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    }

    as.character(as.expression(eq));
  }


  r2 = format(summary(lm(median_post_pred ~ mean_obs, data = processed_demon$resids))$r.squared, digits = 2)

  post_pred_plot <- processed_demon$resids %>%
    ggplot(aes(mean_obs,median_post_pred)) +
    geom_point(shape = 21, fill = 'lightseagreen', size = 3, alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0), size = 1.5) +
    geom_smooth(method = 'lm', color = 'red',linetype = 'longdash', size = 2, se = F) +
    plot_theme +
    xlab('Observed') +
    ylab('Median Posterior Predicted') +
    geom_text(x = -12,y = -2.5, label = paste('R^2 ==', r2 ), parse = T, size = 6)

  ggsave(filename = paste(out_folder,'Obs vs predicted.pdf', sep = ''), post_pred_plot, width = fig_width,
         height = fig_height)

#   reorder_size <- function(x) {
#     factor(x, levels = names(sort(table(x))))
#   }
#
# ggplot(reg_data,aes(reorder_size(species), fill = reorder_size(species))) +
#     geom_bar() +
#   coord_flip() +
#   ylab('Observations') +
#   xlab('') +
#   scale_fill_brewer(guide = F,palette = 'Spectral')

  local_files <- ls()

  plot_files <- local_files[grep('_plot', local_files, fixed = T)]

  plot_list <- list()

  for (i in 1:length(plot_files))
  {
    eval(parse( text = paste('plot_list$',plot_files[i],' <- ',plot_files[i], sep = '')))
  }



return(pres_plots = plot_list)
}