#' Prepare objects required for Laplaces Demon
#'
#' \code{prep_demon} takes data and options and
#' returns a list with components required by Laplaces Demon
#' @param demondat data to be used
#' @param pos_vars possible variables to include in the model
#' @param  scale_numerics scale T or F to center and scale variables


prep_demon <- function(demondat, pos_vars,scale_numerics = F,constant = T)
{
  ind_vars <- demondat[,pos_vars] #pull out things you need

  var_types <- sapply(ind_vars,class)

  factors <- pos_vars[which(var_types == 'character' | var_types == 'factor') ]

  binaries <- ind_vars %>%
    gather('var','value',convert = T) %>%
    subset(!var %in% factors) %>%
    mutate(value = as.numeric(value)) %>%
    group_by(var) %>%
    summarise(is_binary = sum(value == 0 | value == 1) == length(value)) %>%
    subset(is_binary == T)

  numerics <- pos_vars[which(var_types == 'numeric')]

  leave_alone <- numerics[grepl('fished',numerics) | grepl('year',numerics) | grepl('temp',numerics)| grepl('lag',numerics)]

  numerics <- numerics[!numerics %in% leave_alone]

  if (scale_numerics == T)
  {
    for (j in 1:length(numerics)){
      ind_vars[,numerics[j]] <- CenterScale(as.matrix(ind_vars[,numerics[j]]))
    }
  }
  for (f in 1:length(factors))
  {
    ind_vars <- spread_factor(ind_vars,var = factors[f])
  }

  if (constant == T){
    ind_vars$constant <- 1
  }

  return(ind_vars)
}