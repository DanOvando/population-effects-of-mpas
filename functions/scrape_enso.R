#' \code{scrape_enso} extracts El Niño data from
# [here]('https://www.esrl.noaa.gov/psd/enso/mei/table.html)
#'
#' @return a saved csv of ENSO idices
#' @export
#' @import readr
#'
#' @examples
#' scrape_enso(outdir = '../data/')
scrape_enso <- function(outdir = '../data/'){


    # enso <- readr::read_lines('https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.anom.data')

    enso <- readr::read_table('https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.anom.data',
                              col_names = F, skip = 1,n_max = 2017-1870 + 1) %>%
      set_names(c('year',1:12)) %>%
      gather(month, enso,-year) %>%
      mutate(month =  as.numeric(month)) %>%
      mutate(date = lubridate::ymd(paste(year,month,1, sep = '-'))) %>%
      mutate(enso = ifelse(enso == -99.99, NA, enso))


readr::write_csv(enso, path = paste0(outdir,'enso.csv'))


}