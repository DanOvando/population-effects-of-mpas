spread_factors <- function(data, drop_one = T) {
  factors <-
    colnames(data)[map_lgl(data, ~ class(.x) %in% c('character', 'factor'))]
  for (i in seq_along(factors)) {
    var <- factors[i]

    data <- data %>%
      ungroup() %>%
      mutate(dummy = 1, index = 1:nrow(.)) %>%
      spread_(var, 'dummy', fill = 0, sep = 'dummy') %>%
      arrange(index) %>% {
        if (drop_one == T) {
          dplyr::select(., -index, -contains('dummy')[1]) #drop the first dummy
        } else{
          dplyr::select(., -index)
        }
      } %>%
      set_names(str_replace(colnames(.), 'dummy', '-'))

  }

  return(data)
}