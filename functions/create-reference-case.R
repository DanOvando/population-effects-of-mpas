create_reference_case <- function(seen_aug, seen_model){

seen_factors <- colnames(seen_aug)[map_lgl(seen_aug, ~class(.x) == 'factor' |class(.x) == 'character')]

top_factors <- seen_aug %>%
  group_by_at(vars(one_of(seen_factors[seen_factors != 'factor_year']))) %>%
  count() %>%
  arrange(desc(n)) %>%
  ungroup() %>%
  slice(1) %>%
  select(-n)

seen_series <- seen_aug %>%
  purrrlyr::dmap_if(is.character, as.factor)

factor_names <- colnames(top_factors)

for (i in factor_names){ # filter down to the most common thing

  where_seen <- (seen_series[,i] == (top_factors[,i][[1]]))

  seen_series <- seen_series[where_seen,]

}

seen_series <- seen_series %>%
  slice(1)

if (class(seen_model) == 'lmerMod'){

  num_original_years <- seen_model@flist$factor_year %>% unique()

  seen_series <- seen_series[rep(1, length(num_original_years)),] #replicate the number of years

  seen_series$factor_year <- levels(num_original_years) %>% as.factor()

} else {
num_original_years <- seen_model$xlevels$factor_year

seen_series <- seen_series[rep(1, length(num_original_years)),] #replicate the number of years

seen_series$factor_year <- seen_model$xlevels$factor_year %>% as.factor()
}
return(seen_series)
}
