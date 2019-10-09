center_scale <- function(x, xname = '', omit_names = 'none') {
  if (is.numeric(x) & !all(unique(x) %in% c(1, 0)) &
      !xname %in% omit_names) {
    x <- (x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T))

  }

  return(x)

}