#' length_to_weight
#' Esimate weight from one row of observations in pisco data, where one
#' observation is a defined by an integer count of fish of one species
#' of a given length or within a range of lengths defined by a min, max
#' and mean length
#'
#' @param mean_length the mean length of observations
#' @param min_length min observed length
#' @param max_length max observed length
#' @param count total count for that observation
#' @param weight_a a coefficient for weight
#' @param weight_b b coefficient for weight
#' @param length_units units of length, e.g. cm, mm
#' @param weight_units output units for weight, e.g. kg
#' @param length_for_weight_units units from length to weight relationship
#' @param length_type_for_weight type of length, e.g. SL, TL
#' @param tl_sl_a params for conversion to/from SL to TL
#' @param tl_sl_b params for conversion to/from SL to TL
#' @param tl_sl_type params for conversion to/from SL to TL
#' @param tl_sl_formula params for conversion to/from SL to TL
#' @param length_prob type of distirbution to use to generate multiple length comps
#'
#' @return a real number of total estiamted weight for that observation
#' @export
#'
length_to_weight <-
  function(mean_length,
           min_length,
           max_length,
           count,
           weight_a,
           weight_b,
           length_units = 'cm',
           weight_units = 'g',
           length_for_weight_units = 'mm',
           length_type_for_weight,
           tl_sl_a,
           tl_sl_b,
           tl_sl_type,
           tl_sl_formula,
           length_prob = "multinomial") {
    #
    # generate_lengths <- function(count,mean_length, min_length, max_length){

    if (is.na(count) | count == 0){

      outweight <-  0
    }  else {

      if (is.na(min_length) |
          is.na(max_length)) {
        #generate distribution of lengths

        lengths <-  rep(mean_length, count)

      } else{

        if(length_prob == "unif"){

        lengths <- runif(count, min = min_length, max = max_length)
        } else {

          probs <- rep(1, round(max_length) - round(min_length) + 1)

        lengths <- as.numeric(rmultinom(1, count, prob = probs))

        lengths <- data_frame(length = round(min_length):round(max_length), counts = lengths) %>%
          mutate(lengths = map2(length, counts, ~rep(.x,.y)))

        lengths <- unlist(lengths$lengths)

        }
      }

      if (length_type_for_weight == 'SL') {
        if (tl_sl_type  == 'TYPICAL') {

          weight_lengths <-  lengths * tl_sl_a + tl_sl_b
        } else{
          weight_lengths <- (lengths - tl_sl_b) / tl_sl_a

        } # convert from observed total lengths to standard lengths

      } else {
        weight_lengths <-  lengths
      } # convert from standard length to total length if needed

      if (length_for_weight_units == 'mm'){

       weight_a <-  weight_a * 10^weight_b
      } # get a units correct for dealing with length observations in cm instead of mm

      if (weight_units == 'kg'){ weight_a <-  weight_a / 1000} # convert from kg to grams

      weight <-  weight_a * weight_lengths ^ weight_b

      outweight = sum(weight)
    }
    return(outweight)
  } #close function
