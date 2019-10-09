get_fish_life <- function(genus, species) {
  Predict = Plot_taxa(
    Search_species(Genus = genus, Species = species)$match_taxonomy,
    mfrow = c(2, 2),
    partial_match = T,
    verbose = F
  )
  out <- Predict[[1]]$Mean_pred %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()

  out[colnames(out) != 'Temperature'] <-
    exp(out[colnames(out) != 'Temperature'])

  return(out)

}
