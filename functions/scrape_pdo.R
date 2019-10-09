#' \code{scrape_pdo} extracts pacific decadal oscilation data from
# [here]('https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/pdo.long.data')
#'
#' @return a saved csv of PDO idices
#' @export
#'
#' @examples
#' scrape_pdo(outdir = '../data/')
scrape_pdo <- function(outdir = '../data/'){

  pdo <- read_table("https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/pdo.long.data",
                    na = c("-99.99", "99.99",'-99'), skip = 1, n_max = lubridate::year(Sys.time()) - 1900,
                    col_names = c("year", 1:12)) %>%
    gather(month, pdo, -year) %>%
    mutate(month = as.double(month),
           date = lubridate::ymd(paste(year,month,'01', sep = '-')))
  readr::write_csv(pdo, path = paste0(outdir,'pdo.csv'))


}