library(tidyverse)
library(TMB)
#library(VAST)
library(sdmTMB)


# assume loading in data already

pisco_abundance_data <-
  abundance_data$data[abundance_data$data_source == "pisco"][[1]]



empirical_biomass_densities  <- pisco_abundance_data %>%
  filter(classcode %in% unique(fitted_data$classcode)) %>%
  group_by(year, site_side, region, zone, eventual_mpa, classcode, targeted) %>%
  summarise(total_classcode_density = sum(exp(log_density)),
            mean_experience  = mean(cumulative_n_obs),
            mean_vis = mean(mean_vis),
            mean_depth = mean(mean_depth),
            mean_canopy = mean(mean_canopy)) %>%
  group_by(year, site_side, classcode) %>%
  summarise(
    total_biomass_density = sum(total_classcode_density),
    mean_biomass_density = mean(total_classcode_density),
    mean_experience  = mean(mean_experience),
    mean_vis = mean(mean_vis),
    mean_depth = mean(mean_depth),
    mean_canopy = mean(mean_canopy)
  ) %>%
  ungroup()

empirical_biomass_densities %>% ggplot(aes(mean_biomass_density)) + geom_histogram() +
  facet_wrap( ~ classcode, scales = "free_x") + 
  theme_minimal()


classcodes <- unique(empirical_biomass_densities$classcode)


standardizer <- function(temp_classscode, pknots = 0.5){

temp <- empirical_biomass_densities %>%
  filter(classcode == temp_classscode) %>%
  # select(year, site_side, mean_biomass_density) %>%
  left_join(
    site_data %>% select(site, side, lon_wgs84, lat_wgs84) %>% unite(col = site_side, site:side, sep = '-'),
    by = c("site_side")
  ) %>%
  rename(X = lon_wgs84,
         Y = lat_wgs84) %>% 
  mutate(mean_experience = scale(mean_experience)) %>% 
  mutate(mean_experience_2 = mean_experience^2,
         mean_vis = scale(mean_vis))

temp %>%
  ggplot(aes(X, Y, color = mean_biomass_density)) +
  geom_point() +
  facet_wrap( ~ year) +
  theme_minimal() +
  scale_color_viridis()


# test out sdmTMB

temp_spde <-
  make_spde(x = temp$X,
            y = temp$Y,
            n_knots = round(n_distinct(temp$site_side) * pknots))

# plot_spde(temp_spde)

m <- sdmTMB(
  data = temp,
  formula = mean_biomass_density ~ 0 + as.factor(year) + mean_experience + mean_experience_2 + mean_vis,
  time = "year",
  spde = temp_spde,
  family = tweedie(link = "log")
)

ogtemp <- temp
# temp$resids <- residuals(m) # randomized quantile residuals
# hist(temp$resids)
# qqnorm(temp$resids)
# abline(a = 0, b = 1)
# ggplot(temp, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#   geom_point() + facet_wrap( ~ year) + coord_fixed() +
#   theme_minimal()


pred_grid <-
  expand_grid(X = seq(min(temp$X), max(temp$X), by = 0.025),
              Y = seq(min(temp$Y), max(temp$Y), by = 0.025))

pred_grid <-
  tibble(year = unique(ogtemp$year), x = list(pred_grid)) %>%
  unnest(cols = x) %>% 
  mutate(mean_experience = 0,
         mean_experience_2 = 0,
         mean_vis = 0)

temp_hat <-
  predict(m, return_tmb_object = TRUE, newdata = pred_grid)

plot_map <- function(dat, column) {
  ggplot(dat) +
    geom_tile(aes_string("X", "Y", fill = column)) +
    facet_wrap( ~ year) +
    coord_fixed()
}


# points <- temp %>%
#   select(X, Y) %>%
#   unique()

# plot_map(temp_hat$data %>% filter(year == 2014), "exp(est)") +
#   scale_fill_viridis_c(trans = "sqrt") +
#   ggtitle("Prediction (fixed effects + all random effects)") +
#   geom_point(data = points, aes(x = X, y = Y), color = "red")

ind <- get_index(temp_hat, bias_correct = FALSE)

out <- list(fit = m,
            prediction = temp_hat,
            index = ind)

return(out)
# scale <-
#   1 #2 * 2 / 1000 # 2 x 2 km grid and converted from kg to tonnes
# ggplot(ind, aes(year, est * scale)) + geom_line() +
#   geom_ribbon(aes(ymin = lwr * scale, ymax = upr * scale), alpha = 0.4) +
#   xlab('Year') + ylab('Biomass estimate (metric tonnes)')
# 
# a = temp %>%
#   group_by(year) %>%
#   summarise(m = mean(mean_biomass_density))
# 
# plot(scale(ind$est))
# lines(scale(a$m))

}

test <- map(classcodes, standardizer)


standardized_index <- tibble(classcode = classcodes, index = map(test,"index")) %>% 
  unnest(cols = "index") %>% 
  group_by(classcode) %>% 
  mutate(scaled_est = scale(est)) %>% 
  left_join(life_history_data %>% select(classcode, targeted), by = "classcode") %>% 
  ungroup() %>% 
  mutate(fyear = factor(year)) %>%
  mutate(fyear = relevel(fyear, ref = "2002"))


standardized_index %>% 
  group_by(classcode) %>% 
  mutate(est = (est)) %>% 
  ggplot(aes(year, est)) + 
  geom_line() + 
  facet_wrap(~classcode, scales = "free_y") + 
  theme_minimal()

zissou_index <-
  data_frame(
    abundance_hat = zissou_fit$zissou_report$abundance_hat,
    classcode = rep(unique(fitted_data$classcode), each = length(years))
  ) %>%
  mutate(log_abundance_hat = log(abundance_hat)) %>%
  group_by(classcode) %>%
  mutate(year = 1999 + 1:length(abundance_hat)) %>%
  mutate(scaled_abundance_hat = (abundance_hat - mean(abundance_hat)) / sd(abundance_hat)) %>%
  ungroup()


raw_index <- empirical_biomass_densities %>%
  group_by(year, classcode) %>%
  summarise(density = mean(mean_biomass_density)) %>% 
group_by(classcode) %>%
  mutate(scaled_density = scale(density)) %>%
  ungroup()

compare <- standardized_index %>% 
  select(year, classcode, est)



ggplot() +
  geom_line(data = standardized_index, aes(year, scaled_est)) +
  geom_line(data = zissou_index, aes(year,scaled_abundance_hat, color = "zissou")) + 
  geom_point(data = raw_index, aes(year, scaled_density, color = "empirical")) + 
  facet_wrap( ~ classcode) + 
  theme_minimal()


standardized_did <- stan_glm(scaled_est ~ targeted + fyear + targeted:fyear, data = standardized_index)

standardized_did_plot <-
  bayesplot::mcmc_areas(as.matrix(standardized_did),
                        regex_pars = ":fyear") +
  geom_vline(aes(xintercept = 0), linetype = 2, color = "red") +
  scale_y_discrete(labels = c(2000:2001, 2003:2017)) +
  coord_flip() +
  scale_x_percent() 
