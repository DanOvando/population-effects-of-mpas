refit_lw <- function(linf, a, b, linf_factor) {
   old_lengths <- seq(1, linf, length.out = 100)

  pesos <- a * old_lengths ^ b

  new_lengths <- old_lengths * linf_factor

  fitfoo <- function(ws, sizes, peso) {
    new_weight <- ws[1] * sizes ^ ws[2]

    ss <- sum((new_weight - peso) ^ 2)

    return(ss)

  }

  news <-
    nlminb(
      start = c(a, b),
      fitfoo,
      sizes = new_lengths,
      peso = pesos,
      lower = c(0, 2),
      upper = c(1e3, 4)
    )

  out <- list(a = news$par[1], b = news$par[2])

  return(out)


}