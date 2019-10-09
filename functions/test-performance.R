test_performance <-
  function(fishes,
           year_mpa,
           min_year = 75,
           max_year = 100,
           time_step = 1) {
    simple_data <- fishes %>%
      select(loo, k, lm, m, targeted, classcode, commonname, pisco_samples) %>%
      unnest() %>%
      filter(year > min_year, year < max_year) %>%
      select(-pop,-sampled_lengths,-diver_stats) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year)) %>%
      mutate(
        factor_year = as.factor(year),
        log_density = log(density),
        logical_targeted = targeted > 0,
        any_seen = density > 0,
        region = patches,
        post_mpa = year >= year_mpa
      )

    years_protected <- unique(simple_data$year) - year_mpa

    bins <-
      c(seq(min(years_protected),-1, by = 5), seq(0, max(years_protected), by = 5))

    year_block <- data_frame(year = unique(simple_data$year)) %>%
      mutate(years_protected = year - year_mpa) %>%
      mutate(protected_block = cut(
        years_protected,
        breaks = bins,
        include.lowest = T
      )) %>%
      select(-years_protected)

    simple_data <- simple_data %>%
      left_join(year_block, by = "year")

    true_effect <- fishes %>%
      select(classcode, targeted, mpa_effect) %>%
      unnest() %>%
      # mutate(mpa_effect = log(`with-mpa`) - log(`no-mpa`)) %>%
      filter(year > min_year, targeted == 1) %>%
      mutate(post_mpa = year > year_mpa) %>%
      mutate(year = year * time_step) %>%
      mutate(subyear = year - floor(year),
             year = floor(year))


    bare_bones_model <-
      lm(log_density ~ targeted + protected_block + targeted:protected_block,
         data = simple_data)
    if (n_distinct(simple_data$diver) > 1) {
      # mixed_effect_model <-
      #   lme4::lmer(
      #     log_density ~ (1 + enso|classcode) + diver + targeted*protected_block ,
      #     data = simple_data,
      #     verbose = 1
      #   )

      mixed_effect_model <-
        lme4::lmer(
          log_density ~ (1 | classcode) + diver + targeted * protected_block ,
          data = simple_data,
          verbose = 1
        )

    } else
    {
      mixed_effect_model <-
        lme4::lmer(
          log_density ~ (1 | classcode) + targeted * protected_block ,
          data = simple_data,
          verbose = 1
        )

    }

    pre_post_model <-
      lm(log_density ~ targeted + post_mpa + targeted:post_mpa, data = simple_data)

    get_range <- function(bin) {
      bin_range <- str_extract(bin, pattern = '(?<=\\().*(?=])')

      mean(str_split(bin_range, ',', simplify = T) %>% as.numeric())

    }


    mean_effect <- true_effect %>%
      group_by(year) %>%
      summarise(mean_effect = mean(mpa_effect))

    bare_bones_did <- broom::tidy(bare_bones_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(year = map_dbl(term, get_range)) %>%
      ggplot() +
      geom_pointrange(
        aes(
          x = year,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      ) +
      geom_line(
        data = true_effect,
        aes(year - year_mpa, mpa_effect, color = classcode),
        show.legend = F,
        alpha = 0.5
      ) +
      labs(title = 'bare bones') +
      geom_line(
        data = mean_effect,
        aes(year - year_mpa, mean_effect),
        color = "red",
        size = 1.5,
        linetype = 2
      )

    mixed_effect_did <- broom::tidy(mixed_effect_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(year = map_dbl(term, get_range)) %>%
      ggplot() +
      geom_pointrange(
        aes(
          x = year,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      ) +
      geom_line(
        data = true_effect,
        aes(year - year_mpa, mpa_effect, color = classcode),
        show.legend = F
      ) +
      labs(title = 'mixed effects') +
      geom_line(
        data = mean_effect,
        aes(year - year_mpa, mean_effect),
        color = "red",
        size = 1.5,
        linetype = 2
      )


    check_block_did <- broom::tidy(pre_post_model) %>%
      filter(str_detect(term, 'targeted:')) %>%
      mutate(post_mpa = TRUE) %>%
      ggplot() +
      geom_boxplot(data = true_effect,
                   aes(post_mpa, mpa_effect),
                   color = 'red') +
      geom_pointrange(
        aes(
          x = post_mpa,
          y = estimate,
          ymin = estimate - 1.96 * std.error,
          ymax = estimate + 1.96 * std.error
        )
      )

    out_plot <-
      {
        bare_bones_did + mixed_effect_did + plot_layout(ncol = 1)
      } + check_block_did + plot_layout(ncol = 2)

    out <- list(
      bare_bones_did = bare_bones_did,
      mixed_effect_did = mixed_effect_did,
      out_plot = out_plot,
      bare_bones_model = bare_bones_model,
      mixed_effect_model = mixed_effect_model,
      pre_post_model = pre_post_model
    )

    return(out)
  }
