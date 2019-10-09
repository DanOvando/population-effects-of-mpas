sim_scuba_steve <- function(pop,diver, fish, patches, effort = 1, cv = 0){

pop_at_age <- pop %>%
  group_by(age) %>%
  summarise(number = sum(numbers))

pop_at_age <- pop_at_age %>%
  mutate(mean_length = fish$linf * (1 - exp(-fish$vbk * (age - fish$t0)))) %>%
  mutate(selectivity = (1 / (1 + exp(-log(
  19
) * ((mean_length - diver$sel_size_50) / (diver$sel_size_delta)
)))))

pop_at_age <- pop_at_age %>%
  mutate(potential_seen = number * (1 - exp(-diver$q * effort * selectivity)))

pop_at_age <- pop_at_age %>%
  mutate(numbers_seen = rmultinom(1, sum(potential_seen), potential_seen) %>% as.numeric())


age_samples <- pop_at_age


length_samples <- spasm::sample_lengths(age_samples, sample_col = quo(numbers_seen),sample_type = 'other',
                                        cv = cv,
                                        k = fish$vbk,
                                        linf = fish$linf,
                                        t0 = fish$t0,
                                        percent_sampled = 1,
                                        linf_buffer = 1.25
                                        ) %>%
  mutate(weight = numbers * (fish$weight_a * length_bin ^fish$weight_b))




return(list(age_samples = age_samples,
            length_samples = length_samples))



}