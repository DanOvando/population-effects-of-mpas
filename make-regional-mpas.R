# Run Complete Zissou Analysis --------------------------------------------


# libraries ---------------------------------------------------------------
set.seed(42)
library(scales)
library(viridis)
library(ggmap)
library(forcats)
library(stringr)
library(lubridate)
library(lme4)
library(TMB)
library(FishLife)
library(patchwork)
library(rstan)
library(rstanarm)
library(extrafont)
library(hrbrthemes)
library(here)
library(doParallel)
library(spasm)
library(furrr)
library(ggsci)
library(recipes)
library(sf)
library(ggmap)
library(rEDM)
library(bayesplot)
library(tabulizer)
library(tidyverse)
extrafont::loadfonts()
rstan_options(auto_write = TRUE)


functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions


# options -----------------------------------------------------------------

run_name <- 'v6.0'

run_description <- "PNAS R and R with simplified DiD"

# the following analysis run the complete contents of "regional-effects-of-mpas". Each section depends on the out
# outcomes of the prior section, but will load relevant saved files.
# So, once you've run simulate_mpas, you can set it to FALSE and validata_mpas will work


run_did <- TRUE # run difference in difference on data from the CINMS

run_tmb <- FALSE

simulate_mpas <- TRUE # simulate MPA outcomes

validate_mpas <- TRUE

process_results <- TRUE

get_cdfw_catches <-  FALSE

knit_paper <- FALSE

sim_years <- 50

num_patches <- 50

n_cores <- 8

# prepare run -------------------------------------------------------------

run_dir <- here::here("results", run_name)

if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

data_dir <- "data"

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))

plot_theme <- hrbrthemes::theme_ipsum(base_size = 14,
                                      axis_title_size = 16,
                                      base_family = "Fira Sans")

theme_set(plot_theme)

fig_name <- "presentations"

fig_width <- 12

fig_height <- fig_width / 1.333

device <- "pdf"

fig_dir <- file.path(run_dir,"figures",fig_name)

if (!dir.exists(fig_dir)){
  dir.create(fig_dir, recursive = TRUE)
}

# run MLPA difference-in-difference ---------------------------------------
if (run_did == TRUE){


run_length_to_density <-  FALSE

tmb_to_stan <- FALSE # fit the model in stan instead of TMB

max_generations <- 4

run_vast <- FALSE # run VAST, best to leave off for now

num_knots <-  10

channel_islands_only <- T # only include channel islands, leave T

min_year <- 1999 # included years must be greater than this

occurance_ranking_cutoff <- 0.5

small_num <-  0 # no idea

use_mpa_site_effects <- F # no idea

use_cfdw_fished <-  F

year_mpa <- 2003

rank_targeting <- F

max_generations <- 5

max_year <- 2017

bin_years <- FALSE # whether or not to group did terms into bins

length_data <- read.csv(here::here(data_dir,"UCSB_FISH.csv")) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  mutate(classcode = tolower(classcode)) %>%
  mutate(observer = as.character(observer)) %>%
  mutate(observer = ifelse(observer == '', 'unknown', observer))

# Filter data per operations in Fish size biomass processing CIMPA.sas file

length_data <- length_data %>%
  filter(
    level != 'CAN',
    campus == 'UCSB',
    method %in%  c(
      'SBTL_FISH'
    ),!(site == 'SCI_PELICAN' & side == 'FAR WEST'),!(toupper(classcode) %in% c('NO_ORG', 'LDAL', 'CNIC'))
  )

yoy_foo <- function(classcode, fish_tl) {
  new_classcode <- classcode
  if (fish_tl <= 5 &
      is.na(fish_tl) == F &
      classcode %in% c('cpun', 'bfre', 'ocal', 'hsem')) {
    new_classcode <- glue::glue('{classcode}_yoy')

  } # close less than 5cm

  if (fish_tl <= 10 & is.na(fish_tl) == F) {
    if (toupper(classcode) %in% c('PCLA', 'SMYS', 'SPAU', 'SPIN', 'SPUL')) {
      new_classcode <- glue::glue('{classcode}_yoy')
    }

    if (toupper(classcode)  %in% c('SATR', 'SCAR', 'SCAU', 'SCHR', 'GBY')) {
      new_clascode = 'kgb'
    }

    if (toupper(classcode) %in% c('SFLA', 'SSER', 'SMEL', 'OYT')) {
      new_classcode = 'oyb'
    }

    if (toupper(classcode) %in% c('SMIN',
                                  'SNEB',
                                  'SRAS',
                                  'STRE',
                                  'SSAX',
                                  'SDAL',
                                  'SDIP',
                                  'SEBSPP')) {
      new_classcode = 'r_yoy'
    }


  } #close less than than 10cm

  return(new_classcode)

} # close yoy_foo


length_data <- length_data %>%
  mutate(classcode = map2_chr(classcode, fish_tl, yoy_foo))

# flip situations where min_tl is greater than max_tl

flipped <- length_data$max_tl < length_data$min_tl & !is.na(length_data$max_tl) & !is.na(length_data$min_tl)

length_data[flipped,c("min_tl","max_tl")] <-  length_data[flipped,c("max_tl","min_tl")]


# add in covariates -------------------------------------------------------

life_history_data <-
  read_csv(here::here(
    data_dir,
    'VRG Fish Life History in MPA_04_08_11_12 11-Mar-2014.csv'
  )) %>%
  rename(classcode = pisco_classcode) %>%
  mutate(classcode = tolower(classcode)) %>%
  magrittr::set_colnames(., tolower(colnames(.)))

fish_life <-
  life_history_data$taxa %>% str_split(' ', simplify = T) %>%
  as_data_frame() %>%
  select(1:2) %>%
  set_names(c('genus', 'species'))


fish_life <- fish_life %>%
  mutate(life_traits = map2(genus, species, safely(get_fish_life)))

fish_life <- fish_life %>%
  mutate(fish_life_worked = map(life_traits, 'error') %>% map_lgl(is.null)) %>%
  filter(fish_life_worked) %>%
  mutate(life_traits = map(life_traits, 'result')) %>%
  unnest() %>%
  mutate(taxa = glue::glue('{genus} {species}')) %>%
  set_names(tolower)

life_history_data <- life_history_data %>%
  left_join(fish_life, by = 'taxa')

# check age at recruitment to survey

smallest_length_seen <- length_data %>%
  group_by(classcode) %>%
  summarise(smallest_length_seen = min(fish_tl, na.rm = T))

life_history_data <- life_history_data %>%
  left_join(smallest_length_seen, by = 'classcode') %>%
  mutate(first_age_seen = log(1 - smallest_length_seen / loo) / -k)

#Convert to fished species things with CDFW catches

caselle_fished_species <-
  read_csv(file = here::here(data_dir, 'caselle-2015-fished-species-list.csv')) %>%
  mutate(commonname = tolower(commonname),
         mod_commonname = tolower(mod_commonname)) %>%
  select(mod_commonname, caselle_targeted)

# make sure targeted list matches that from Caselle 2015
life_history_data <- life_history_data %>%
  left_join(caselle_fished_species, by = c('commonname' = 'mod_commonname')) %>%
  mutate(targeted = ifelse(is.na(caselle_targeted), targeted, caselle_targeted))


if (get_cdfw_catches){
  
  file_names <- list.files(here('data','cdfw-data'))
  
  
  process_cdfw <- function(file_name) {
    dat <-
      tabulizer::extract_tables(here("data","cdfw-data", file_name))
    
    year <- str_extract(file_name, pattern = '(?<=landings).*(?=_)')
    
    year <- as.numeric(paste0('20', year))
    
    collapse_pages <- function(x) {
      missing_all <-
        map_lgl(x %>% as_tibble, ~ all(str_count(.x) == 0)) == F
      
      x <- x[, missing_all]
      
      # if (dim(x)[2] == 15){ #super hacky fix for random blank column in some years
      #
      #   x <- x[,-14]
      #
      # }
      
      x <- as_tibble(x)
      
      x <- slice(x,-(1:2))
      
      species <- x[, 1]
      species <-
        str_replace_all(species$V1, "(\\.)|^[ \t]|[ \t]$", '')
      
      species_mat <-
        str_split(species, pattern = ',|:', simplify = T)[, c(2, 1)]
      
      species <-
        unite(species_mat %>% as_data_frame(),
              col = species,
              V1,
              V2,
              sep = ' ') %>%
        mutate(species = str_replace(species, "^[ \t]|[ \t]$", '') %>% tolower())
      
      x[, 1] <- species$species
      numfoo <- function(z) {
        z <- str_replace_all(z, '\\.|\\,', '')
        
        if (any(is.na(as.numeric(z)) == F)) {
          z = as.numeric(z)
        } else{
          z = z
        }
        
      }
      
      x <- map_df(x, numfoo)
      
      
      return(x)
    }
    
    flat_dat <- map_df(dat, collapse_pages)
    
    
    flat_dat <- flat_dat %>%
      set_names(
        c(
          'species',
          'january',
          'february',
          'march',
          'april',
          'may',
          'june',
          'july',
          'august',
          'september',
          'october',
          'november',
          'december',
          'total'
        )
      ) %>%
      select(-total) %>%
      gather('month', 'pounds_caught',-species) %>%
      mutate(year = year) %>%
      filter(str_detect(species, '(total)|(Total)|(ishes)|(aters)') == F)
  }
  
  flat_cdfw <- map_df(file_names, process_cdfw)
  
  write_csv(flat_cdfw, path = here("data","clean_flat_cdfw.csv"))
  
  # flat_cdfw <- read_csv( "processed_data/clean_flat_cdfw.csv")
  
  # out <- bind_rows(flat_cdfw)
  
  commnames <- unique(flat_cdfw$species)
  
  sci_names <- taxize::comm2sci(commnames, db = "worms")
  
  names <- tibble(com_name =commnames, sci_name =  map_chr(sci_names,~.x[1])) %>%
    mutate(no_space =  str_remove_all(com_name," ")) %>%
    filter(!is.na(sci_name))
  
  cdfw_catches <- flat_cdfw %>%
    mutate(no_space =  str_remove_all(species," ")) %>%
    left_join(names, by = 'no_space') %>%
    mutate(sci_name = tolower(sci_name)) %>%
    group_by(no_space, sci_name, year) %>%
    summarise(pounds_caught = sum(pounds_caught))
  
  # cdfw_catches %>%
  #   filter(!is.na(sci_name)) %>%
  #   ggplot(aes(year, pounds_caught, fill = sci_name)) +
  #   geom_area(show.legend = FALSE)
  
  write_csv(cdfw_catches, path = here('data','cdfw-catches.csv'))
  
} else {
  
  cdfw_catches <-
    read_csv(file = here::here("data", 'cdfw-catches.csv'))
}

cdfw_catches <-  cdfw_catches %>%
  group_by(sci_name, year) %>%
  summarise(catch = sum(pounds_caught, na.rm = T)) %>%
  left_join(
    life_history_data %>% select(taxa, classcode) %>% mutate(taxa = tolower(taxa)),
    by = c("sci_name" = "taxa")
  ) %>%
  filter(!is.na(classcode))




targeting_rank <- cdfw_catches %>%
  group_by(sci_name) %>%
  summarise(total_catch = sum(catch)) %>%
  arrange(desc(total_catch)) %>%
  ungroup() %>%
  mutate(catch_rank = percent_rank(total_catch))

fished_species <-
  data_frame(classcode = unique(cdfw_catches$classcode),
             fished = 1)

c99 <- cdfw_catches %>%
  ungroup() %>%
  filter(year == min(year)) %>%
  mutate(year = 1999)

lag_catches <-  cdfw_catches %>%
  bind_rows(c99) %>%
  group_by(classcode) %>%
  arrange(year) %>%
  mutate(year = year + 1,
         lag_catch = catch) %>%
  select(-catch)

if (use_cfdw_fished == T) {
  life_history_data <- life_history_data %>%
    left_join(fished_species, by = 'classcode') %>%
    mutate(targeted = ifelse(fished == 1 &
                               is.na(fished) == F,
                             'Targeted',
                             targeted)) %>%
    select(-fished)
}

life_history_data$targeted <-
  as.numeric(life_history_data$targeted == 'Targeted')

if (rank_targeting == T) {
  life_history_data <- life_history_data %>%
    left_join(targeting_rank, by = 'classcode') %>%
    mutate(targeted = ifelse(targeted == 1 &
                               (is.na(catch_rank) == F), catch_rank, targeted)) %>%
    mutate(targeted = ifelse(is.na(catch_rank) &
                               targeted == 1, 0.5, targeted))

}


site_data <- read_csv(here::here(data_dir,'Final_Site_Table_UCSB.csv')) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  select(
    site,
    side,
    mpagroup,
    mpa_status,
    reserve,
    region,
    year_mpa,
    mpaareanm2,
    lat_wgs84,
    lon_wgs84
  ) %>%
  rename(lat_wgs84 = lon_wgs84,
         lon_wgs84 = lat_wgs84) %>%
  unique() %>%
  mutate(eventual_mpa = (year_mpa > 0))

# create geographic clustering of species

sightings <- length_data %>%
  filter(is.na(count) == F) %>%
  select(classcode, site, side) %>%
  left_join(site_data, by = c('site', 'side')) %>%
  filter(is.na(lat_wgs84) == F & is.na(lon_wgs84) == F,
         region %in% c('ANA', 'SCI', 'SMI', 'SRI')) %>%
  group_by(classcode) %>%
  summarise(
    mean_lat = mean(lat_wgs84),
    mean_long = mean(lon_wgs84),
    min_lat = min(lat_wgs84),
    max_lat = max(lat_wgs84),
    min_long = min(lon_wgs84),
    max_long = max(lon_wgs84)
  ) %>%
  ungroup()


num_clusters <- data_frame(clusters = 1:20) %>%
  mutate(within_ss = map_dbl(clusters, ~ sum(
    kmeans(
      sightings %>% select(contains('_')),
      centers = .x,
      nstart = 25,
      iter.max = 1000
    )$withinss
  )))

cluster_classcodes <-
  kmeans(
    sightings %>% select(contains('_')),
    centers = 5,
    nstart = 25,
    iter.max = 1000
  )

sightings <- sightings %>%
  mutate(geographic_cluster = cluster_classcodes$cluster)

# geographic_cluster_plot <-
#   ggmap::qmplot(mean_long,
#                 mean_lat ,
#                 color = factor(geographic_cluster),
#                 data = sightings) + theme_classic()

life_history_data <- life_history_data %>%
  left_join(sightings %>% select(classcode, geographic_cluster),
            life_history_data,
            by = 'classcode')

# add life history data into length data
length_data <- length_data %>%
  left_join(life_history_data, by = 'classcode')


# load kelp data

kelp_conn <- ncdf4::nc_open(here::here("data","LandsatKelpBiomass_2017.nc"))

kelp_vars <- c('lat','lon','year','quarter', 'biomass')

kelp <- ncdf4::ncvar_get(kelp_conn, varid = c("biomass")) %>%
  as_data_frame() %>%
  mutate(lat = ncdf4::ncvar_get(kelp_conn, varid = c("lat")),
         lon = ncdf4::ncvar_get(kelp_conn, varid = c("lon"))) %>%
  gather(tdex, kelp_biomass, -lat,-lon) %>%
  mutate(tdex = str_replace_all(tdex,"\\D","") %>% as.numeric()) %>%
  na.omit()

time <-
  data_frame(
    year = ncdf4::ncvar_get(
      kelp_conn,
      varid = c("year")),
    quarter = ncdf4::ncvar_get(kelp_conn, varid = c("quarter"))) %>%
  mutate(tdex = 1:nrow(.))

kelp <- kelp %>%
  left_join(time, by = "tdex") %>%
  select(-tdex) %>%
  mutate(rounded_lat = plyr::round_any(lat, 0.1),
         rounded_lon = plyr::round_any(lon, 0.1)) %>%
  group_by(rounded_lat, rounded_lon, year, quarter) %>%
  summarise(mean_kelp = mean(kelp_biomass, na.rm = T))


kelp_recipe <- recipes::recipe(mean_kelp ~ ., data = kelp)


## process kfm dat
if (!file.exists(here::here("data","SBCMBON_integrated_fish_20181129.csv"))){
  
  message("KFM data not loaded, downloading, requires internet connection")
  
  download.file("https://pasta.lternet.edu/package/data/eml/edi/5/3/4b1aed53da2cdbf25f9f9001bd74ace8",
                destfile = here::here("data","SBCMBON_integrated_fish_20181129.csv"))
  
}

if (!file.exists(here::here("data","SBCMBON_site_geolocation_20170523.csv"))){
  
  message("KFM data not loaded, downloading, requires internet connection")
  
  download.file("https://pasta.lternet.edu/package/data/eml/edi/5/3/f9eee1db42cd598cfb7b2578a35662fa",
                destfile = here::here("data","SBCMBON_site_geolocation_20170523.csv"))
  
}



kfm_data <-
  read_csv(here::here("data",'SBCMBON_integrated_fish_20181129.csv'))

kfm_locations <-
  read_csv(here::here('data','SBCMBON_site_geolocation_20170523.csv'))

conditions_data <- length_data %>%
  group_by(site, side, classcode, year) %>%
  summarise(
    mean_temp = mean(temp, na.rm = T),
    mean_kelp = mean(pctcnpy, na.rm = T),
    mean_vis = mean(vis, na.rm = T)
  ) %>%
  group_by(year) %>%
  mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  ) %>%
  ungroup() %>%
  mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp),
    mean_kelp = ifelse(is.na(mean_kelp), mean(mean_kelp, na.rm = T), mean_kelp)
  )

observer_experience <- length_data %>%
  group_by(year, month, observer) %>%
  summarise(n_obs = length(side)) %>%
  arrange(observer, year) %>%
  group_by(observer) %>%
  mutate(cumulative_n_obs = cumsum(n_obs),
         lifetime_obs = sum(n_obs)) %>%
  ungroup() %>%
  mutate(observer_ranking = percent_rank(lifetime_obs)) %>%
  mutate(
    trunc_observer = ifelse(
      observer_ranking > 0.5 | observer == 'unknown',
      observer,
      'infrequent'
    )
  ) %>%
  arrange(desc(observer_ranking)) %>%
  mutate(
    cumulative_n_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      cumulative_n_obs
    ),
    n_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      n_obs
    ),
    lifetime_obs = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      lifetime_obs
    ),
    observer_ranking = ifelse(
      trunc_observer == 'infrequent' |
        trunc_observer == 'unknown',
      0,
      observer_ranking
    )
  ) %>%
  mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2)


length_data <- length_data %>%
  left_join(observer_experience, by = c('year', 'month', 'observer'))


density_data <- read_csv(here::here(data_dir,"ci_reserve_data_final3 txt.csv")) %>%
  magrittr::set_colnames(., tolower(colnames(.))) %>%
  gather('concat.name', 'value', grep('_', colnames(.)), convert = T) %>%
  mutate(
    data.type = gsub('\\_.*', '', concat.name),
    classcode = gsub('.*\\_', '', concat.name)
  ) %>%
  mutate(value = as.numeric(value),
         index = 1:nrow(.)) %>%
  spread(data.type, value) %>%
  rename(site_side = site.side) %>%
  arrange(index) %>%
  select(-index)

site_coords <- density_data %>%
  group_by(site, side) %>%
  summarise(latitude = mean(lon.wgs84, na.rm = T),
            longitude = mean(lat.wgs84, na.rm = T))

species_distributions <- length_data %>%
  left_join(site_coords, by = 'site') %>%
  ungroup() %>%
  filter(is.na(latitude) == F & is.na(longitude) == F) %>%
  group_by(classcode) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  summarise(
    mean_latitude = sum(count * latitude, na.rm = T) / sum(count, na.rm = T),
    mean_longitude = sum(count * longitude, na.rm = T) / sum(count, na.rm = T),
    max_latitude = max(latitude, na.rm = T),
    min_latitude = min(latitude, na.rm = T),
    max_longitude = max(longitude, na.rm = T),
    min_longitude = min(longitude, na.rm = T)
  )

# ggmap::qmplot(min_longitude,min_latitude, data = species_distributions)


# ggmap::qmplot(longitude,latitude, color = side, data = site_coords)


if (file.exists(here::here(data_dir,'enso.csv'))) {
  enso <- read_csv(here(data_dir,'enso.csv'))
} else {

  scrape_enso(outdir = paste0(data_dir,'/'))

  enso <- read_csv(here(data_dir,'enso.csv'))
}

if (file.exists(here(data_dir,'pdo.csv'))) {
  pdo <- read_csv(here(data_dir,'pdo.csv'))

} else {
  scrape_pdo(outdir = paste0(data_dir,'/'))

  pdo <- read_csv(here(data_dir,'pdo.csv'))
}



# convert transect data to density estimates ------------------------------


if (file.exists(file.path(run_dir,"pisco-data.Rdata")) == F |
    run_length_to_density == T) {


  length_example <-   length_data %>%
    filter(is.na(commonname) == F) %>%
    mutate(biomass_g = pmap_dbl(
      list(
        mean_length = fish_tl,
        min_length = min_tl,
        max_length = max_tl,
        count = count,
        weight_a = wl_a,
        weight_b = wl_b,
        length_type_for_weight = wl_input_length,
        length_for_weight_units = wl_l_units,
        weight_units = wl_w_units,
        tl_sl_a = lc.a._for_wl,
        tl_sl_b = lc.b._for_wl,
        tl_sl_type = lc_type_for_wl,
        tl_sl_formula = ll_equation_for_wl
      ),
      length_to_weight
    ))

  pisco_data <- length_example %>%
    mutate(
      observer = ifelse(is.na(observer), 'unknown', observer),
      surge = ifelse(is.na(surge), 'unknown', surge)
    ) %>%
    rename(
      total_biomass_g = biomass_g,
      mean_temp = temp,
      mean_vis = vis,
      mean_depth = depth,
      mean_canopy = pctcnpy
    )

  species_sightings <- pisco_data %>%
    left_join(site_data, by = 'site') %>%
    group_by(site) %>%
    summarise(species_seen = list(unique(classcode)))

  pisco_data <- pisco_data %>%
    ungroup() %>%
    left_join(site_data %>% select(site, region), by = 'site') %>%
    select(region, site, side, year, month, day, zone, level, transect) %>%
    unique() %>%  {
      pmap(
        list(
          this_region = .$region,
          this_site = .$site,
          this_side = .$side,
          this_year = .$year,
          this_month = .$month,
          this_day = .$day,
          this_transect = .$transect,
          this_zone = .$zone,
          this_level = .$level
        ),
        add_missing_fish,
        observations = pisco_data,
        species_sightings = species_sightings,
        life_history_vars = colnames(life_history_data)
      )
    } %>%
    bind_rows()


  classcodes <- unique(pisco_data$classcode)

  life_history_vars <-
    which(colnames(pisco_data) %in% colnames(life_history_data))

  for (i in 1:length(classcodes)) {
    where_class <-
      pisco_data$classcode == classcodes[i] &
      pisco_data$count == 0

    temp_life_history <- life_history_data %>%
      filter(classcode == classcodes[i])

    temp_life_history <-
      temp_life_history[, colnames(pisco_data)[life_history_vars]]

    if (nrow(temp_life_history) > 1) {
      stop('multiple classcodes')
    }

    pisco_data[where_class, life_history_vars] <-
      temp_life_history

    if (any((
      colnames(temp_life_history) == colnames(pisco_data[1, life_history_vars])
    ) == F)) {
      stop()
    }

    if (any(is.na(pisco_data$classcode))) {
      stop()
    }


  }

  check_life_history_fill <- pisco_data %>%
    group_by(classcode) %>%
    summarise(a = n_distinct(commonname)) %>%
    arrange(desc(a))

  if (any(check_life_history_fill$a > 1)) {
    stop('multiple species per classcode')
  }

  save(file = file.path(run_dir,'pisco-data.Rdata'),
       pisco_data)

} else {
  load(file.path(run_dir,'pisco-data.Rdata'))

}


# sum biomass across all observed sizes ---------------------------------


transect_covariates <-
  c(
    'mean_depth',
    'mean_vis',
    'mean_temp',
    'surge',
    'mean_canopy',
    'observer',
    'n_obs',
    'cumulative_n_obs'
  )

mean_foo <- function(var) {
  if (class(var) != 'character') {
    out <- mean(var, na.rm = T)
  } else{
    out <- paste(unique(var), collapse = '-')
  }

}


channel_sides <-
  data_frame(region = c("ANA", "SCI", "SRI", "SMI"),
             midline = c(34, 34, 34, 34.03))

pisco_data <- pisco_data %>%
  mutate(year_month = year + (month / 12 - .1))

rolling_mean_temperatures <- pisco_data %>%
  mutate(year_month = year + (month / 12 - .1)) %>%
  left_join(site_data, by = c('site','side')) %>%
  filter(region %in% c("ANA", "SCI","SRI","SMI")) %>%
  group_by(year) %>%
  summarise(mean_temperature = mean(mean_temp)) %>%
  mutate(rolling_mean = RcppRoll::roll_mean(mean_temperature,4, align = "right",
                                            fill = NA)) %>%
  mutate(rolling_mean3 = RcppRoll::roll_mean(mean_temperature,3, align = "right",
                                             fill = NA))  %>%
  mutate(rolling_mean2 = RcppRoll::roll_mean(mean_temperature,2, align = "right",
                                             fill = NA)) %>%
  mutate(rolling_mean1 = RcppRoll::roll_mean(mean_temperature,1, align = "right",
                                             fill = NA)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean3, rolling_mean)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean2, rolling_mean)) %>%
  mutate(rolling_mean = ifelse(is.na(rolling_mean), rolling_mean1, rolling_mean)) %>%
  select(year, rolling_mean) %>%
  mutate(rolling_mean = zoo::na.approx(rolling_mean)) %>%
  ungroup()


# sum biomass across all lengths
pisco_data <- pisco_data %>%
  mutate(sampling_event = paste(campus, method, year, month, day, site, side, zone, transect, level) %>%
           as.factor() %>% as.numeric()) %>%
  group_by(sampling_event,
           classcode) %>%
  mutate(total_biomass_g = sum(total_biomass_g)) %>% 
  mutate(
    density_g_m2 = total_biomass_g / 60,
    total_count = sum(count),
    mean_length = mean(fish_tl, na.rm = T),
    counter = 1:length(total_biomass_g)
  ) %>%
  ungroup() %>%
  filter(counter == 1) %>%
  select(-counter)

pisco_data <- pisco_data %>%
  mutate(
    any_seen = total_biomass_g > 0,
    factor_year = factor(year),
    log_density = log(density_g_m2),
    factor_month = factor(month),
    site_side = paste(site,side, sep = '-')
  ) %>%
  group_by(site, side, month, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(site, side, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(site, year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  group_by(year) %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>%
  ungroup() %>%
  mutate(
    mean_vis = ifelse(is.na(mean_vis), mean(mean_vis, na.rm = T), mean_vis),
    mean_canopy = ifelse(is.na(mean_canopy), mean(mean_canopy, na.rm = T), mean_canopy)
  ) %>% mutate(
    mean_temp = ifelse(is.na(mean_temp), mean(mean_temp, na.rm = T), mean_temp)) %>%
  mutate(temp_deviation = abs(mean_temp - temperature)) %>%
  mutate(generations_protected = pmin(round((year - year_mpa - 1) / tm), max_generations))


# Process kfm count data

kfm_sites <- kfm_locations %>%
  select(latitude, longitude) %>%
  as.matrix() %>%
  na.omit()

pisco_sites <- site_coords %>%
  ungroup() %>%
  select(latitude, longitude) %>%
  as.matrix()

nearest_pisco_site <- RANN::nn2(pisco_sites, kfm_sites)

nearest_pisco_site <- site_coords[nearest_pisco_site$nn.idx[, 1], ]

kfm_pisco_locations <- kfm_locations %>%
  na.omit() %>%
  bind_cols(nearest_pisco_site)


kfm_data <- kfm_data %>%
  filter(sample_method != "visualfish" &
           sample_method != "crypticfish",
         data_source == 'kfm') %>%
  left_join(life_history_data, by = c('taxon_name' = 'taxa')) %>%
  left_join(kfm_pisco_locations %>% select(-data_source), by = 'site_id') %>%
  filter(str_detect(geolocation, '_island'),
         data_source != 'lter',
         is.na(geolocation) == F) %>%
  mutate(region = str_extract(geolocation, '.*(?=_island)')) %>%
  filter(region %in% c('anacapa', 'santa_cruz', 'santa_rosa', 'san_miguel'))

kfm_data <- kfm_data %>%
  mutate(
    count = as.numeric(count),
    density = count / area,
    year = lubridate::year(date),
    month = lubridate::month(date)
  )

region_table <- tribble(
  ~ region,
  ~ short_region,
  'anacapa',
  'ANA',
  'santa_cruz',
  'SCI',
  'santa_rosa',
  'SRI',
  'san_miguel',
  'SMI'
)

kfm_data <- kfm_data %>%
  left_join(region_table, by = 'region') %>%
  select(-region) %>%
  rename(region = short_region)

kfm_data <- kfm_data %>%
  mutate(
    log_density = log(density),
    any_seen = density > 0,
    factor_year = as.factor(year),
    factor_month = as.factor(month),
    total_biomass_g = density
  ) %>%
  select(-region) #for compatibility with things later on


# save raw-isa data -------------------------------------------------------

save(
  file = file.path(run_dir,'rawish_zissou_data.Rdata'),
  life_history_data,
  pisco_data,
  kfm_data
)

# filter data -------------------------------------------------------------

consistent_sites <- pisco_data %>%
  group_by(site) %>%
  summarise(
    num_years = length(unique(year)),
    min_year = min(year),
    max_year = max(year)
  ) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

consistent_regions <- pisco_data %>%
  left_join(site_data, by = 'site') %>%
  group_by(region) %>%
  summarise(
    num_years = length(unique(year)),
    min_year = min(year),
    max_year = max(year)
  ) %>%
  arrange(desc(num_years)) %>%
  filter(min_year <= 2000,
         num_years > 10)

nobs_quantiles <- pisco_data %>%
  filter(any_seen == T) %>%
  group_by(year, classcode) %>%
  summarise(nobs = sum(any_seen)) %>%
  ungroup() %>%
  {
    quantile(.$nobs)
  }

well_observed_species <- pisco_data %>%
  filter(year > 1999) %>%
  group_by(year, classcode, commonname, targeted) %>%
  summarise(nseen = sum(any_seen, na.rm = T)) %>%
  group_by(commonname, classcode, targeted) %>%
  summarise(min_seen = min(nseen, na.rm = T)) %>%
  arrange(desc(min_seen)) %>%
  ungroup() %>%
  filter(min_seen > 2) %>%
  # filter(min_seen > nobs_quantiles[3]) %>%
  mutate(classcode = (classcode))


filterfoo <-
  function(x,
           y,
           min_seen_years = 14,
           mpa_start_year = 2003,
           min_year,
           filter_level) {
    # only years above 1999 at the main channel islands
    x <- x %>%
      mutate(classcode = y) %>%
      filter(is.na(log_density) == F) %>%
      filter(year > min_year &
               region %in% c('ANA', 'SCI', 'SMI', 'SRI')) %>%
      group_by(!!filter_level) %>%
      mutate(
        num_years_observed = length(unique(year[any_seen == T])),
        min_year = min(year[any_seen == T], na.rm = T),
        max_year = max(year[any_seen == T], na.rm = T)
      ) %>%
      ungroup() %>%
      filter(
        num_years_observed >= min_seen_years,
        min_year <= mpa_start_year,
        max_year >= mpa_start_year
      ) %>%
      select(-classcode)# actually observed for at least 10 years
  }


# create dropped pairs


site_data$rounded_lat <- site_data$lat_wgs84

site_data$rounded_lon <- site_data$lon_wgs84

abundance_data <- pisco_data %>%
  left_join(site_data, by = c("site","side")) %>%
  left_join(lag_catches %>% select(year, classcode, lag_catch), by = c('year', 'classcode')) %>%
  mutate(data_source = 'pisco') %>%
  group_by(data_source, classcode) %>%
  nest() %>%
  bind_rows(kfm_data %>%  mutate(data_source = 'kfm') %>%
              left_join(site_data, by = c("site","side")) %>%
              left_join(lag_catches %>% select(year, classcode, lag_catch), by = c('year', 'classcode')) %>%
              mutate(year_month = year + (month / 12 - .1),
                     site_side = site_id) %>%
              group_by(data_source, classcode) %>%
              nest()
  ) %>%
  filter(classcode %in% well_observed_species$classcode) %>%
  # mutate(data = map(data, ~ left_join(.x, site_data, by = c('site', 'side')))) %>%
  # mutate(data = map(data, ~ left_join(.x, channel_sides, by = c("region")))) %>%
  # mutate(data = map(data, ~ mutate(.x, channel_zone = ifelse(lat_wgs84 > midline, "inner", "outer")))) %>%
  mutate(data = map(data, ~ left_join(.x, rolling_mean_temperatures, by = c("year")))) %>%
  mutate(data = map(data, ~ mutate(.x, temp_deviation = abs(rolling_mean - temperature)))) %>%
  mutate(data = map(data, ~ mutate(.x, month = as.numeric(as.character(factor_month))))) %>%
  mutate(data = map(data, ~ left_join(.x, enso, by = c('year', 'month')))) %>%
  mutate(data = map(data, ~ left_join(.x, pdo, by = c('year', 'month')))) %>%
  mutate(data = map(data, ~ mutate(.x, lag_catch = ifelse(is.na(lag_catch), 0, lag_catch)))) %>%
  mutate(
    data = map2(
      data,
      classcode,
      filterfoo,
      min_year = min_year,
      min_seen_years = 14,
      filter_level = quo(classcode)
    )
  ) %>% # filter out things
  mutate(dim_data = map_dbl(data, nrow)) %>%
  filter(dim_data > 0) %>%
  left_join(life_history_data %>% select(classcode, commonname, targeted),
            by = 'classcode') %>%
  filter(
    str_detect(commonname, 'YOY') == F,
    is.na(targeted) == F,
    str_detect(classcode, '_yoy') == F
  )

pisco_cols <- abundance_data$data[abundance_data$data_source == "pisco"][1][[1]] %>% colnames()

abundance_data <- abundance_data %>%
  mutate(data = map(data,~.x[,map_lgl(.x, ~!all(is.na(.x)))])) %>% # unnest having tough time with unnesting all NA columns
  ungroup() %>% 
  select(classcode, data_source, data) %>%
  unnest(col = data, keep_empty = TRUE) %>%
  group_by(data_source) %>%
  nest()


# add kelp data

kelp_model <- caret::knnreg(mean_kelp ~ ., data = kelp, k = 5)

abundance_data <- abundance_data %>%
  mutate(data = map2(data,list(kelp_model),  ~ mutate(.x, interp_kelp = predict(
    .y,
    newdata = .x %>% mutate(quarter = (month/3) %>% ceiling()) %>%  select(year, quarter, rounded_lat, rounded_lon)
  ))))

save(
  file = file.path(run_dir,"abundance_data.Rdata"),
  abundance_data
)

# fit model ---------------------------------------------------------------


tmb_abundance_data <- abundance_data %>%
  mutate(data = map(data, ~select(.,log_density,
                                  targeted,
                                  site_side,
                                  region,
                                  geographic_cluster,
                                  level,
                                  factor_month,
                                  cumulative_n_obs,
                                  surge,
                                  mean_canopy,
                                  mean_vis,
                                  mean_depth,
                                  lag_catch,
                                  factor_year,
                                  classcode,
                                  interp_kelp,
                                  eventual_mpa,
                                  temp_deviation,
                                  any_seen,
                                  year
  )))


script_name <- "fit_zissou"



var_options <-  data_frame(var_names = c("pisco_a","kfm"),
                           non_nested_variables =
                             list(c(
                               "site_side",
                               'level',
                               'factor_month',
                               'cumulative_n_obs',
                               'surge',
                               'mean_depth',
                               'mean_vis'
                             ),
                             c(
                               "site_side",
                               'factor_month'
                             )
                             ))


stupid_filter <- function(x){

  missing = x %>%
    group_by(classcode) %>%
    summarise(nyd = n_distinct(year)) %>%
    ungroup() %>%
    filter(nyd < 17)

  out <- x %>%
    filter(!(classcode %in% missing$classcode))

}

tmb_abundance_data <- tmb_abundance_data %>%
  mutate(data = map(data, stupid_filter))

model_runs <- cross_df(
  list(
    data_source = abundance_data$data_source,
    var_names = var_options$var_names,
    data_to_use = c("all","mpa_only", "fished_only"),
    center_scale = c(TRUE)
  )
) %>%
  left_join(var_options, by = "var_names") %>%
  left_join(abundance_data, by = "data_source") %>%
  filter(!(data_source == "kfm" & str_detect(var_names,"pisco"))) %>%
  filter(!(data_source == "kfm" & (data_to_use %in% c("mpa_only","fished_only"))))


if (run_tmb == T){

  if (file.exists("fit-progress.txt")){
    file.remove("fit-progress.txt")
  }

  # browser()
  # future::plan(future::multiprocess, workers = 4)

  doParallel::registerDoParallel(cores = n_cores)
  #
  # model_runs <- model_runs %>%
  #   slice(3)

  compile(here::here("src", paste0(script_name, ".cpp")), "-O0") # what is the -O0?
  
  dyn.load(dynlib(here::here("src", script_name)))
  
  fits <- foreach::foreach(i = 1:nrow(model_runs)) %dopar% {

    # fits <- list()

    # for (i in 1:nrow(model_runs)){

    sfz <- safely(fit_zissou)

    fits <- sfz(data = model_runs$data[[i]],
                non_nested_variables = model_runs$non_nested_variables[[i]],
                data_to_use = model_runs$data_to_use[[i]],
                center_scale = model_runs$center_scale[[i]],
                run_dir = run_dir,
                script_name = script_name,
                fixed_regions = FALSE,
                include_intercept = TRUE,
                fixed_did = FALSE,
                non_nested_did_variables = c(
                  "temp", "kelp", "lag_catch"
                ),
                bin_years = bin_years
    )

    write(glue::glue("{round(100*i/nrow(model_runs),2)}% done with model fits"), file = "fit-progress.txt",
          append = T)

    write(glue::glue("error message:{fits$error}"), file = "fit-errors.txt",
          append = T)

    out <- fits

  } # close dopar

  save(file = file.path(run_dir, 'model_fits.Rdata'),
       fits)

  model_runs$tmb_fit <- fits

  save(file = file.path(run_dir, 'model_runs.Rdata'),
       model_runs)


} else {
  
  

  # load(file = file.path(run_dir, 'model_runs.Rdata'))
}
model_runs <- model_runs %>%
  mutate(
    did_fit = pmap(
      list(
        data = data,
        data_to_use = data_to_use,
        data_source = data_source
      ),
      estimate_did,
      site_data = site_data,
      cdfw_catches = cdfw_catches,
      life_history_data = life_history_data
    )
  )


# did_worked <- map(model_runs$did_fit, "error") %>% map_lgl(is.null)

# model_runs <- model_runs %>% 
#   filter(did_worked)

# model_runs <- model_runs %>% 
#   filter(did_worked) %>% 
#   mutate(did_fit = map(did_fit, "result"))

  
did_fits <- model_runs %>% 
  select(-data)
  
# print(object.size(did_fits), units = "Mb")


write_rds(did_fits, path = file.path(run_dir,"did_fits.rds"))

test <- read_rds( file.path(run_dir,"did_fits.rds"))
}
# simulate mpa outcomes ---------------------------------------------------
 if (simulate_mpas == TRUE){

   sim_years <- 50

   burn_years <- 25

   num_patches <- 50

   run_experiments <- TRUE

   save_experiment <- TRUE

   create_grid <- TRUE

   samps <- 20

   grid_search <-  FALSE


   # prepare data -----------------------------------------------------

   experiment_dir <- file.path(run_dir, "experiments")

   load(file = file.path(run_dir, "rawish_zissou_data.Rdata"))

   model_runs <- readr::read_rds(file.path(run_dir, "did_fits.rds"))
   
   load(file = file.path(run_dir, "abundance_data.Rdata"))

   fitted_data <- abundance_data$data[abundance_data$data_source == "pisco"][[1]]

   seen_species <- life_history_data %>%
     filter(classcode %in% (fitted_data$classcode %>% unique())) %>%
     rename(sci_name = taxa,
            linf = vbgf.linf,
            common_name = commonname) %>%
     mutate(sci_name = tolower(sci_name)) %>%
     filter(is.na(linf) == F)

   fitted_data <- fitted_data %>%
     left_join(life_history_data %>% select(classcode, taxa) %>% unique(), by = "classcode")

   # prepare experiments -----------------------------------------------------


   # run experiments ---------------------------------------------------------
   if (run_experiments == T) {
     if (create_grid == TRUE) {

       # switch(
       #   utils::menu(c("y", "n"), title = "STOP!!! Are you sure you want to overwrite the experiment grid? y/n", graphics = TRUE) + 1,
       #   cat("Nothing done\n"),
       #   cat("OK, sit back, this might take a while"),
       #   stop("Got it, canceling run to avoid overwriting"))
       #
       #
       # stop("force stop")

       if (grid_search == T) {
         sim_grid <- expand.grid(
           scientific_name = unique(fitted_data$taxa),
           steepness = seq(0.6, 1, by = .2),
           adult_movement = seq(1, 20, length.out = 3),
           larval_movement = seq(1, 20, length.out = 3),
           density_movement_modifier = c(0.25, 1),
           density_dependence_form = 1:5,
           mpa_size = c(.1, .3, .75),
           f_v_m = seq(.01, 1.25, by = 0.5),
           fleet_model = c("constant-catch"),
           effort_allocation = c("profit-gravity", "simple"),
           stringsAsFactors = F
         )
       } else{
         sim_grid <-
           tibble(
             scientific_name = sample(rfishbase::taxonomy(), samps, replace = T),
             steepness = runif(samps, min = 0.6, max = 0.95),
             adult_movement = sample(0:round(0.25 * num_patches), samps, replace = T),
             larval_movement = sample(0:round(0.25 * num_patches), samps, replace = T),
             density_movement_modifier = sample(c(0.25, 1), samps, replace = T),
             density_dependence_form = sample(1:3, samps, replace = T),
             mpa_size = runif(samps, min = 0.01, max = 1),
             f_v_m = runif(samps, min = 0.01, max = 4),
             fleet_model = sample(
               c("open-access", "constant-effort", "constant-catch"),
               samps,
               replace = T
             ),
             effort_allocation = sample(
               c("profit-gravity", "simple", "gravity"),
               samps,
               replace = T
             ),
             year_mpa = sample(5:(sim_years / 1.5), samps, replace = T),
             sprinkler = sample(c(TRUE, FALSE), samps, replace = TRUE),
             mpa_reaction   =  sample(c("stay", "leave"), samps, replace = TRUE),
             min_size = runif(samps, min = 0.01, max = 0.75),
             mpa_habfactor = sample(c(1, 4), samps, replace = TRUE),
             size_limit = runif(samps, 0.1, 1.25),
             random_mpa = sample(c(TRUE, FALSE), samps, replace = TRUE),
             sigma_r = sample(c(0,.05,.1,.2), samps, replace = TRUE),
             rec_ac =  sample(c(0,.05,.1,.2), samps, replace = TRUE)
           )



       }

       # create fish objects
       sim_grid <- sim_grid %>%
         mutate(fish = pmap(
           list(
             scientific_name = scientific_name,
             steepness = steepness,
             adult_movement = adult_movement,
             larval_movement = larval_movement,
             density_dependence_form = density_dependence_form,
             density_movement_modifier = density_movement_modifier,
             sigma_r = sigma_r,
             rec_ac = rec_ac
           ),
           safely(create_fish),
           price = 10
         ))

       fish_worked <- map(sim_grid$fish,"error") %>% map_lgl(is_null)

       sim_grid <- sim_grid %>%
         filter(fish_worked) %>%
         mutate(fish = map(fish, "result"))

       # create fleet objects
       sim_grid <- sim_grid %>%
         mutate(fleet = pmap(
           list(
             fish = fish,
             fleet_model = fleet_model,
             effort_allocation = effort_allocation,
             mpa_reaction = mpa_reaction,
             length_50_sel = size_limit * map_dbl(sim_grid$fish, "length_50_mature")
           ),
           create_fleet,
           q = .1
         ))

       # tune fleet objects

       # future::plan(future::multiprocess, workers = n_cores)
       #
       # message("starting fishery tuning")
       # sim_grid <- sim_grid %>%
       #   mutate(
       #     tuned_fishery = future_pmap(
       #       list(
       #         f_v_m = f_v_m,
       #         fish = fish,
       #         fleet = fleet,
       #         sprinkler = sprinkler,
       #         mpa_habfactor = mpa_habfactor             ),
       #       safely(tune_fishery),
       #       num_patches = num_patches,
       #       sim_years = sim_years,
       #       burn_years = burn_years,
       #       .progress = T
       #     )
       #   )

       doParallel::registerDoParallel(cores = n_cores)

       if (dir.exists(experiment_dir) == F) {
         dir.create(experiment_dir, recursive = T)
       }

       sft = purrr::safely(tune_fishery)

       tuned_fisheries <-
         foreach::foreach(i = 1:nrow(sim_grid)) %dopar% {


           tuned_fishery = sft(
                     f_v_m = sim_grid$f_v_m[i],
                     fish = sim_grid$fish[[i]],
                     fleet = sim_grid$fleet[[i]],
                     sprinkler = sim_grid$sprinkler[i],
                     mpa_habfactor = sim_grid$mpa_habfactor[i],
                   num_patches = num_patches,
                   sim_years = sim_years,
                   burn_years = burn_years
                 )

           filename <- glue::glue("fishery_{i}.rds")


           if (save_experiment == TRUE) {
             saveRDS(tuned_fishery, file = file.path(experiment_dir,filename))

           } else {
             out <- tuned_fishery

           }
         } # close tuning


       loadfoo <-
         function(fishery, experiment_dir) {
           ex <-
             readRDS(file.path(experiment_dir,glue::glue("fishery_{fishery}.rds")))
         }




      tuned_results = map(1:nrow(sim_grid), safely(loadfoo), experiment_dir = experiment_dir)

      tuning_worked <- map(tuned_results,"error") %>% map_lgl(is.null)

      # annoying process to catch things that segfaulted in the parallel process
      sim_grid <- sim_grid %>%
        filter(tuning_worked) %>%
        mutate(tuned_fishery = map(tuned_results[tuning_worked],"result"))


       save(sim_grid, file = file.path(run_dir, "sim_grid.Rdata"))
       message("finished fishery tuning")

     } else{
       load(file = file.path(run_dir, "sim_grid.Rdata"))

     }

     doParallel::stopImplicitCluster()

     tuning_worked <-
       map(sim_grid$tuned_fishery, "error") %>% map_lgl(is.null)

     sim_grid <- sim_grid %>%
       filter(tuning_worked) %>%
       mutate(fish = map(tuned_fishery, c("result", "fish")),
              fleet = map(tuned_fishery, c("result", "fleet")))

     sim_grid$tuned_fishery <- map(sim_grid$tuned_fishery, "result")

     sim_grid <- sim_grid %>%
       select(-tuned_fishery)

     sim_grid$experiment <- 1:nrow(sim_grid)

     if (dir.exists(experiment_dir) == F) {
       dir.create(experiment_dir, recursive = T)
     }

     # library(progress)

     message("starting mpa experiments")

     # sim_grid <- sample_n(sim_grid, 100)
     obs = ls()
     size_grid <- tibble(obs = obs) %>%
       mutate(size = map_dbl(obs, ~as.numeric(object.size(get(.x))))) %>%
       arrange(desc(size))

     doParallel::registerDoParallel(cores = n_cores)


     mpa_experiments <-
       foreach::foreach(i = 1:nrow(sim_grid), .noexport = c("model_runs","pisco_data",'abundance_data',"kfm_data", 'fitted_data')) %dopar% {
         # pb$tick()
         results <- sim_grid %>%
           slice(i) %>%
           mutate(
             mpa_experiment = pmap(
               list(
                 fish = fish,
                 fleet = fleet,
                 mpa_size = mpa_size,
                 year_mpa = year_mpa,
                 sprinkler = sprinkler,
                 mpa_habfactor = mpa_habfactor,
                 min_size = min_size,
                 random_mpa = random_mpa
               ),
               run_mpa_experiment,
               sim_years = sim_years,
               burn_years = burn_years,
               num_patches = num_patches
             )
           )

         # results$mpa_experiment[[1]]$raw_outcomes %>% filter(year == max(year)) -> a
         #
         # a %>% filter(experiment == "with-mpa") %>% group_by(patch) %>% summarise(mpa = unique(mpa)) %>% ggplot(aes(patch, mpa)) + geom_col()

         filename <- glue::glue("experiment_{i}.rds")
         #
         # out <- results

         if (save_experiment == TRUE) {
           saveRDS(results, file = file.path(experiment_dir,filename))

         } else {
           out <- results

         }


         rm(list= ls())
         gc()
       } # close dopar

     doParallel::stopImplicitCluster()

     # mpa_experiments[[1]]$mpa_experiment[[1]]$raw_outcomes %>% filter(year == max(year)) -> a
     #
     # a %>% filter(experiment == "with-mpa") %>% group_by(patch) %>% summarise(mpa = unique(mpa)) %>% ggplot(aes(patch, mpa)) + geom_col()
     #
     # filename <- glue::glue("experiment_{i}.rds")

     # out <- results


   } # close run experiments

   message("finished mpa experiments")

   # process outcomes --------------------------------------------------------




   loadfoo <-
     function(experiment, experiment_dir, output = "mpa-effect") {
       ex <-
         readRDS(file.path(experiment_dir, glue::glue("experiment_{experiment}.rds")))

       # if (output == "mpa-effect") {
       ex$msy <- ex$fish[[1]]$msy

       ex$b_msy <- ex$fish[[1]]$b_msy$b_msy

       ex <- purrr::map_df(list(ex), study_mpa)

       ex <- ex %>%
         select(-mpa_experiment, -fish,-fleet)

       return(ex)
     }


   # if (logged == TRUE) {
   future::plan(future::multiprocess, workers = n_cores)

   processed_grid <-
     future_map(1:nrow(sim_grid),
                safely(loadfoo),
                experiment_dir = experiment_dir,
                .progress = T)

   grid_worked <- map(processed_grid, "error") %>% map_lgl(is_null)

   processed_grid <- processed_grid %>%
     keep(grid_worked)

   processed_grid <- map(processed_grid, "result") %>%
     bind_rows()


   save(processed_grid, file = file.path(run_dir, "processed_grid.Rdata"))

   outcomes <- processed_grid %>%
     unnest()

   # outcomes %>%
   #   ggplot(aes(year,pmin(4,mpa_effect), group = experiment)) +
   #   geom_path() +
   #   facet_wrap(~density_movement_modifier)



 }

# validate estimation strategy --------------------------------------------
if (validate_mpas == TRUE){

  simulate_samples <- TRUE

  burn_years <- 1

  sim_years <- 75

  year_mpa <- 50

  num_patches <-  2

  mpa_size <- 0.5

  n_samples <- 5

  time_step <-  1


  plot_theme <- hrbrthemes::theme_ipsum(base_size = 14,
                                        axis_title_size = 16)


  theme_set(plot_theme)

  # prepare data ------------------------------------------------------------

  load(file.path(run_dir, 'abundance_data.Rdata'))

  load(file.path(run_dir, 'rawish_zissou_data.Rdata'))

  enso <- read_csv(here::here('data','enso.csv'))

  pisco <- abundance_data$data[[1]]

  pisco_fish <- life_history_data %>%
    filter(classcode %in% unique(pisco$classcode)) %>%
    mutate(enviro_effect = ifelse(geographic_cluster > 1,-1, 1)) %>%
    mutate(enviro = list(enso$enso[(nrow(enso) - (sim_years + burn_years)):(nrow(enso) - 1)])) %>%
    mutate(enviro = map2(enviro, enviro_effect, ~ .x * .y))

  n_groups <- 5

  simple_fish <-
    data_frame(
      loo = c(rnorm(n_groups, 120, 0.001), rnorm(n_groups, 100, 0.001)),
      k = 0.4,
      lm = .75 * loo,
      m = 0.2,
      targeted = rep(c(1, 0), each = n_groups),
      classcode = fruit[1:(n_groups * 2)],
      commonname = colors()[1:(n_groups * 2)],
      enviro = NA
    )


  diver <- list(q = .1,
                sel_size_50 = 2 ,
                sel_size_delta = 2)



  pisco_divers <- data_frame(diver = fruit[1:3],
                             diver_stats = list(
                               list(
                                 q = .01,
                                 sel_size_50 = 2 ,
                                 sel_size_delta = 2
                               ),
                               b = list(
                                 q = .067,
                                 sel_size_50 = 10 ,
                                 sel_size_delta = 2
                               ),
                               d = list(
                                 q = .1,
                                 sel_size_50 = 4 ,
                                 sel_size_delta = 2
                               )
                             ))

  if (simulate_samples == T) {
    simple_fish <- create_samples(
      fishes = simple_fish,
      divers = pisco_divers %>% slice(1),
      mpa_size = mpa_size,
      burn_years = burn_years,
      sim_years = sim_years,
      year_mpa = year_mpa,
      num_patches = num_patches,
      samples = n_samples,
      rec_driver = 'stochastic',
      enviro_strength = 1,
      sigma_r = 0,
      cores = n_cores,
      time_step = time_step
    )

    pisco_fish <- create_samples(
      fishes = pisco_fish,
      divers = pisco_divers,
      mpa_size = mpa_size,
      burn_years = burn_years,
      sim_years = sim_years,
      year_mpa = year_mpa,
      num_patches = num_patches,
      samples = n_samples,
      rec_driver = 'environment',
      enviro_strength = 1,
      sigma_r = 0.1,
      cores = n_cores,
      time_step = time_step
    )

    save(file = file.path(run_dir, 'simulated-data.Rdata'),
         simple_fish,
         pisco_fish)

  } else {
    load(file = file.path(run_dir, 'simulated-data.Rdata'))
  }


  # fit simple model --------------------------------------------------------

  a <- pisco_fish %>%
    select(classcode, targeted, pisco_samples) %>%
    unnest() %>%
    select(-pop, -diver_stats, -sampled_lengths) %>%
    group_by(classcode) %>%
    mutate(density = scale(density))

  a %>%
    ggplot(aes(year, density, color = classcode)) +
    geom_smooth() +
    facet_wrap(~targeted, scales = "free_y")


  simple_performance <-
    test_performance(
      simple_fish,
      year_mpa = year_mpa + burn_years,
      min_year = 45,
      time_step = time_step
    )

  pisco_performance <-
    test_performance(pisco_fish,
                     year_mpa + burn_years,
                     min_year = 45,
                     time_step = time_step)

  save(file = file.path(run_dir,'simulated_did.Rdata'), simple_performance, pisco_performance)


}

#process results ------------------------------------------------------------

if (process_results == TRUE){

  short_frame <- 15

  if (!dir.exists(fig_dir)){
    dir.create(fig_dir, recursive = TRUE)
  }

  write(glue::glue("Figures generated on {Sys.Date()} using results {run_name}"),
        file = file.path(fig_dir,"README.txt"))

  experiment_dir <- here::here("results",run_name, "experiments")

  load(file = here::here("results",run_name, "rawish_zissou_data.Rdata"))

  # load(file = here::here("results",run_name, "model_runs.Rdata"))
  
  # load(file = here::here("results",run_name, "model_runs.Rdata"))
  
  model_runs <- read_rds( file.path(run_dir,"did_fits.rds"))
  
  load(file = here::here("results",run_name,"processed_grid.Rdata"))

  load(file = here::here("results",run_name,"simulated_did.Rdata"))

  load(file = file.path(run_dir, "abundance_data.Rdata"))

  channel_islands <- readRDS(here::here("data","channel_islands_map.rds"))

  ca_mpas <- sf::st_read(here::here("data","MPA_CA_Existing_160301")) %>%
    # rmapshaper::ms_simplify() %>%
    sf::st_transform(crs = 4326)
  
  pisco_abundance_data <- abundance_data$data[abundance_data$data_source == "pisco"][[1]]
  

#
#   zissou_theme <-
#     theme_ipsum(
#       base_size = 22,
#       axis_title_size = 26,
#       strip_text_size = 26
#     )
#
#   theme_set(zissou_theme)

  gc <- guide_colorbar(frame.colour = "black",
                       ticks.colour = "black",
                       barheight = 35)


  hgc <- guide_colorbar(frame.colour = "black",
                       ticks.colour = "black",
                       barwidth = 35)

  plot_trans <- "identity"

  # filter results ----------------------------------------------------------


  bad_sims <- processed_grid %>%
    select(-fishery_effect,-density_ratio,-mpa_size) %>%
    unnest(cols = c(mpa_effect)) %>%
    group_by(experiment) %>%
    mutate(year = 1:length(year)) %>%
    mutate(years_protected = year - year_mpa + 1) %>%
    filter(years_protected <= 0) %>%
    mutate(b0 = `no-mpa`[year == min(year)]) %>%
    mutate(depletion = pmax(0,1 - `no-mpa` / b0)) %>%
    mutate(pop_effect = pmin(1, (`with-mpa` - `no-mpa`) / b0)) %>%
    summarise(bad = any(depletion > 0.95 |
                          abs(pop_effect) > 1e-3)) %>%
    filter(bad == TRUE)

  processed_grid <- processed_grid %>%
    filter(!experiment %in% unique(bad_sims$experiment))

  write_rds(processed_grid, file.path(run_dir,"filtered_processed_grid.rds"))


# create summary of results -----------------------------------------------

  simmed_fish_life <-
    processed_grid$scientific_name %>% str_split(' ', simplify = T) %>%
    as_data_frame() %>%
    select(1:2) %>%
    set_names(c('genus', 'species')) %>%
    unique()

  simmed_fish_life <- simmed_fish_life %>%
    mutate(life_traits = future_map2(genus, species, safely(get_fish_life)))

  simmed_fish_life <- simmed_fish_life %>%
    mutate(fish_life_worked = map(life_traits, 'error') %>% map_lgl(is.null)) %>%
    filter(fish_life_worked) %>%
    mutate(life_traits = map(life_traits, 'result')) %>%
    unnest(cols = life_traits) %>%
    mutate(taxa = paste(genus,species)) %>%
    set_names(tolower)


  outcomes <- processed_grid %>%
    left_join(simmed_fish_life %>% select(taxa, m), by = c("scientific_name" = "taxa")) %>%
    select(-fishery_effect, -density_ratio) %>%
    unnest() %>%
    group_by(experiment) %>%
    mutate(year = 1:length(year)) %>%
    ungroup() %>%
    mutate(years_protected = year - year_mpa + 1) %>%
    mutate(mpa_effect = pmax(-.5,pmin(mpa_effect,1))) %>%
    group_by(experiment) %>%
    mutate(b0 = `no-mpa`[year == min(year)]) %>%
    mutate(depletion = 1 - `no-mpa`/b0) %>%
    mutate(final_depletion = depletion[year == max(year)],
           f = f_v_m * m) %>%
    mutate(u = 1 - exp(-f)) %>%
    mutate(final_u = u[year == max(year)]) %>%
    ungroup() %>%
    mutate(pop_effect = pmin(1,(`with-mpa` - `no-mpa`) / b0))

  write_rds(outcomes, file.path(run_dir,"outcomes.rds"))

  density_ratios <- processed_grid %>%
    select(-fishery_effect) %>%
    unnest() %>%
    group_by(experiment) %>%
    mutate(year = 1:length(year)) %>%
    ungroup() %>%
    mutate(years_protected = year - year_mpa + 1) %>%
    left_join(outcomes %>% select(experiment, year, depletion), by = c("experiment", "year"))

  write_rds(density_ratios, file.path(run_dir,"density_ratios.rds"))


  # make examples -----------------------------------------------------------

  fish <- create_fish(r0 = 100)

  fleet <- create_fleet(fleet_model = "constant-effort", fish = fish,
                        initial_effort = 1, q  = 1)


  sizes <- tibble(mpa_size = seq(0,0.75, by = 0.25))

  year_mpa <- 20
  sizes <- sizes %>%
    mutate(foo = map(mpa_size, ~ run_mpa_experiment(fish = fish,
                                                    fleet = fleet,
                                                    year_mpa = year_mpa,
                                                    mpa_size = .x,
                                                    sim_years = 50,
                                                    num_patches = 25,
                                                    burn_years = 1,
                                                    sprinkler = FALSE,
                                                    mpa_habfactor = 100,
                                                    enviro = NA,
                                                    enviro_strength = NA,
                                                    rec_driver = 'stochastic',
                                                    simseed = 42)))

  compare <- sizes %>%
    mutate(delta = map(foo,"outcomes")) %>%
    select(-foo) %>%
    unnest() %>%
    mutate(years_protected = year - year_mpa) %>%
    filter(years_protected <=short_frame)

  ribbon <- compare %>%
    select(years_protected, mpa_size, experiment, biomass) %>%
    spread(experiment, biomass) %>%
    group_by(mpa_size) %>%
    mutate(delta = mean(`with-mpa` - `no-mpa`))

  labelfoo <- function(x){
    paste0("MPA:", percent(as.numeric(x)))

  }


  effect_example_plot <- compare %>%
    ggplot() +
    geom_ribbon(
      data = ribbon,
      aes(
        years_protected,
        ymin = `no-mpa`,
        ymax = `with-mpa`,
        fill = delta
      ),
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(years_protected , biomass, color = experiment), size = 1.5) +
    geom_vline(aes(xintercept = 0), linetype = 2, color = "red") +
    facet_wrap(~ mpa_size, labeller = labeller(mpa_size = labelfoo)) +
    labs(x = "Years with MPA", y = "Biomass") +
    scale_color_aaas(labels = c("Without MPAs", "With MPAs"), name  = '') +
    scale_fill_gradient(low = "white", high = "steelblue",
                        guide = gc)



  # calculate response ratios --------------------------------------------------------

  base_run <- model_runs %>%
    filter(var_names == "pisco_a", data_to_use == "all", center_scale == TRUE)

  mpa_run <- model_runs %>%
    filter(data_source == "pisco",
           var_names == "pisco_a",
           data_to_use == "mpa_only")

  
  fished_run <- model_runs %>%
    filter(data_source == "pisco",
           var_names == "pisco_a",
           data_to_use == "fished_only")
  
  kfm_run <- model_runs %>%
    filter(data_source == "kfm")

  site_data <- read_csv(here::here("data",'Final_Site_Table_UCSB.csv')) %>%
    magrittr::set_colnames(., tolower(colnames(.))) %>%
    select(
      site,
      side,
      mpagroup,
      mpa_status,
      reserve,
      region,
      year_mpa,
      mpaareanm2,
      lat_wgs84,
      lon_wgs84
    ) %>%
    rename(lat_wgs84 = lon_wgs84,
           lon_wgs84 = lat_wgs84) %>%
    unique() %>%
    mutate(eventual_mpa = (year_mpa > 0))

  top_species <- pisco_abundance_data$classcode %>% unique()

  cip_data <- pisco_data %>%
    left_join(site_data, by = c("site","side")) %>%
    filter(is.na(eventual_mpa) == F) %>%
    filter(region %in% c("ANA", "SCI","SRI",'SMI'),
           classcode %in% top_species)


  processed_grid$adult_movement <- (2 * processed_grid$adult_movement) / num_patches


  ## empirical response ratios
  
  #All biomass estimates were converted to metric tons per hectare (t ha−1) to facilitate comparisons with other studies in California and more globally.
  # For each fish species, we summed biomass over the different levels in the water column and calculated 
  # the mean biomass per site per year (averaging across all transects at the site). 
  # These estimates represented the lowest level of replication for all analyses.
  
  
  # pisco_data %>% 
  #   group_by(year, site_side) %>% 
  #   summarise(nl = n_distinct(level),
  #             nt = n_distinct(transect))
  
  consistent_sites <- pisco_abundance_data %>% 
    # filter(classcode %in% unique(pisco_abundance_data$classcode)) %>% 
    select(site_side, year) %>% 
    unique() %>% 
    group_by(site_side) %>% 
    mutate(has_all = all((2003:2012) %in% year)) %>% # results change dramatically if you filter until present time... 
    ungroup() %>% 
    filter(has_all) %>% 
    left_join(site_data %>% unite("site_side", site, side, sep  = "-"), by = c("site_side")) %>% 
    filter(year_mpa == 0 | year_mpa > 2000)
  
  biomass_density  <- pisco_abundance_data %>%
    filter(year >= 2003,
           site_side %in% unique(consistent_sites$site_side)) %>% 
    group_by(year, site,side, region, zone, transect,eventual_mpa, classcode, targeted) %>%
    summarise(total_classcode_density = sum(exp(log_density))) %>%  # sum density across all levels of a transect
    group_by(year, site,side, region,eventual_mpa, classcode, targeted) %>%
    summarise(md = mean(total_classcode_density)) %>% # calculate mean density per year site, side, species, averaging over zone, transect
    group_by(year, site,side,region, eventual_mpa, targeted) %>% 
    summarise(total_biomass_density = (sum(md) / 1e6) * 10000, # calculate total and mean biomass densities across all species per year site side
              mean_biomass_density = (mean(md) / 1e6) * 10000) %>% 
    ungroup() %>% 
    mutate(eventual_mpa = as.numeric(eventual_mpa)) %>% 
    mutate(mpa_location = case_when(eventual_mpa == 1 ~ "IN", TRUE ~ "OUT")) %>% 
    mutate(fyear = factor(year)) 
  
  
  # look at response ratios
  
  bd_trend_plot <- biomass_density %>% 
    ggplot(aes(year, total_biomass_density, color = mpa_location, fill = mpa_location)) + 
    geom_smooth() + 
    facet_wrap(~targeted) + 
    labs(x = '', y = "Mean Total Biomass Density (MT/Hectare)")
  
  targ_rr_fit = stan_glmer(
    total_biomass_density ~ (fyear - 1|eventual_mpa) + (1 |
                                                  region),
    data = biomass_density %>% filter(targeted == 1),
    family = Gamma(link = "log"),
    cores = 4,
    iter = 2000
  )
  
  targ_rr_coefs <- tidybayes::spread_draws(targ_rr_fit, `(Intercept)`,b[year,mpa]) %>% 
    ungroup() %>% 
    filter(!str_detect(mpa,"region")) %>% 
    mutate(year = as.integer(str_remove_all(year,"fyear")),
           mpa = as.integer(str_remove_all(mpa,"eventual_mpa:"))) %>%
    mutate(mean_density = exp(b + `(Intercept)`)) %>% 
    group_by(mpa,year) %>% 
    mutate(prank = percent_rank(mean_density)) %>% 
    filter(prank > 0.05, prank < 0.95) %>% 
    ungroup()
  
  raw_targ_rr_plot <-  targ_rr_coefs %>%
    ggplot(aes(
      mean_density,
      year,
      fill = mpa == 1,
      group = interaction(year, mpa)
    )) +
    ggridges::geom_density_ridges(alpha = 0.75, color = "transparent")
    
  targ_rr_plot <-   targ_rr_coefs %>%
    select(-b, -prank) %>%
    pivot_wider(names_from = mpa, values_from = mean_density) %>%
    mutate(response_ratio = `1` / `0`) %>%
    # group_by(year) %>%
    # mutate(prank = percent_rank(response_ratio)) %>%
    # ungroup() %>%
    # filter(prank > 0.05, prank < 0.95) %>%
    ggplot(aes(response_ratio, year, group = year)) +
    geom_vline(aes(xintercept = 1), color = "red", linetype = 2) + 
    ggridges::geom_density_ridges(alpha = 0.75, color = "transparent") + 
    scale_x_continuous(name = "Response Ratio") + 
    scale_y_continuous(name = element_blank()) + 
    labs(title = "Targeted")
  
  # repeat for non-targeted
  
  
  nontarg_rr_fit = stan_glmer(
    total_biomass_density ~ (fyear - 1|eventual_mpa) + (1 |
                                                          region),
    data = biomass_density %>% filter(targeted == 0),
    family = Gamma(link = "log"),
    cores = 4,
    iter = 2000
  )
  
  nontarg_rr_coefs <- tidybayes::spread_draws(nontarg_rr_fit, `(Intercept)`,b[year,mpa]) %>% 
    ungroup() %>% 
    filter(!str_detect(mpa,"region")) %>% 
    mutate(year = as.integer(str_remove_all(year,"fyear")),
           mpa = as.integer(str_remove_all(mpa,"eventual_mpa:"))) %>%
    mutate(mean_density = exp(b + `(Intercept)`)) %>% 
    group_by(mpa,year) %>% 
    mutate(prank = percent_rank(mean_density)) %>% 
    filter(prank > 0.05, prank < 0.95) %>% 
    ungroup()
  
  raw_nontarg_rr_plot <-  nontarg_rr_coefs %>%
    ggplot(aes(
      mean_density,
      year,
      fill = mpa == 1,
      group = interaction(year, mpa)
    )) +
    ggridges::geom_density_ridges(alpha = 0.75, color = "transparent") 
    
  nontarg_rr_plot <-   nontarg_rr_coefs %>%
    select(-b, -prank) %>%
    pivot_wider(names_from = mpa, values_from = mean_density) %>%
    mutate(response_ratio = `1` / `0`) %>%
    # group_by(year) %>%
    # mutate(prank = percent_rank(response_ratio)) %>%
    # ungroup() %>%
    # filter(prank > 0.05, prank < 0.95) %>%
    ggplot(aes(response_ratio, year, group = year)) +
    geom_vline(aes(xintercept = 1), color = "red", linetype = 2) + 
    ggridges::geom_density_ridges(alpha = 0.75, color = "transparent") + 
    scale_x_continuous(name = "Response Ratio") + 
    scale_y_continuous(name = element_blank()) + 
    labs("Non-Targeted")
  
  
  nontarg_rr <-   nontarg_rr_coefs %>%
    select(-b, -prank) %>%
    pivot_wider(names_from = mpa, values_from = mean_density) %>%
    mutate(nontarg_response_ratio = `1` / `0`) %>% 
    select(.draw,contains("_ratio"), year)
  
  
  targ_rr <-   targ_rr_coefs %>%
    select(-b, -prank) %>%
    pivot_wider(names_from = mpa, values_from = mean_density) %>%
    mutate(targ_response_ratio = `1` / `0`) %>% 
    select(.draw,contains("_ratio"), year)
  
  targ_v_nontarg_rr <- nontarg_rr %>% 
    left_join(targ_rr, by = c(".draw","year")) %>% 
    mutate(tarv_v_nontarg = targ_response_ratio / nontarg_response_ratio)
  
  targ_v_nontarg_plot <-   targ_v_nontarg_rr %>%
    ggplot(aes(tarv_v_nontarg, year, group = year)) +
    geom_vline(aes(xintercept = 1), color = "red", linetype = 2) + 
    ggridges::geom_density_ridges(alpha = 0.75, color = "transparent") + 
    scale_x_continuous(name = "Targeted to Non-Targeted Response Ratio") + 
    scale_y_continuous(name = element_blank())

    # process simulation outcomes

  eqo <- outcomes %>%
    filter(year == max(year))

  fishery_outcomes <- processed_grid %>%
    select(-mpa_effect, -density_ratio) %>%
    unnest() %>%
    group_by(experiment) %>%
    mutate(year = 1:length(year)) %>%
    ungroup() %>%
    mutate(years_protected = year - year_mpa + 1) %>%
    mutate(fishery_effect = pmin(mpa_effect,1)) %>%
    left_join(outcomes %>% select(experiment, year, depletion), by = c("experiment", "year")) %>%
    filter(fleet_model != "constant-catch") %>%
    mutate(msy_effect = pmin(1,(`with-mpa` - `no-mpa`) / msy))


  fishery_eqo <- fishery_outcomes %>%
    filter(year == max(year))


  facet_labels <- c(
    mpa_size = "Range in MPA",
    depletion = "Pre-MPA Depletion"
  )




  ## ------------------------------------------------------------------------
  short_term <- outcomes %>%
    filter(years_protected == short_frame, years_protected >= 0)


  fishery_short_term <- fishery_outcomes %>%
    filter(years_protected <= short_frame, years_protected >= 0)

  median_effect <- median(short_term$mpa_effect[short_term$year == max(short_term$year)])

  median_pop_effect <- median(short_term$pop_effect[short_term$year == max(short_term$year)])

  median_fish_effect <- median(fishery_short_term$mpa_effect[short_term$year == max(short_term$year)], na.rm = TRUE)

  median_fish_effect <- median(fishery_short_term$msy_effect[short_term$year == max(short_term$year)], na.rm = TRUE)


  pos_runs <- round(mean(short_term$pop_effect[short_term$year == max(short_term$year)] > 0),2)

  pos_fish_runs <- round(mean(fishery_short_term$msy_effect[fishery_short_term$year == max(fishery_short_term$year)] > 0, na.rm = TRUE),2)

  pos_fish_effect <- round(mean(short_term$f[short_term$year == max(short_term$year)] > 0),2)

  small_effect <- round( 100* median(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

  min_small_effect <- round( 100* min(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

  max_small_effect <- round( 100* max(short_term$mpa_effect[short_term$mpa_size <= 0.25 & short_term$depletion <= 0.5]))

  medium_effect <- round(100*median(short_term$mpa_effect[short_term$depletion > 0.5 & short_term$depletion <= 0.75]))

  big_effect <- round(100*median(short_term$mpa_effect[short_term$depletion > 0.5 & short_term$depletion > 0.75]))


  ## ----conservation-effect, fig.cap = "Median (A) and range (B) regional MPA conservation effect (expressed as percent changes in biomass with MPAs relative to what would have occured without MPAs) after 15 years of protection. For (A), X-axes indicate the pre-MPA depletion of the fishery, where depletion is the percentage of unfished biomass that has been removed from the population, and Y-axes is the percent of the population's range encompasssed inside an MPA. For B), y-axes show the regional conservation effect.", include = FALSE----


  pop_depletion_plot <- outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(depletion, mpa_effect), binwidth = c(.05, .05),
               show.legend = FALSE) +
    scale_fill_viridis(
      option = "A",
      trans = plot_trans,
      guide = gc,
      name = "Median Effect"
    ) +
    scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
    scale_y_continuous(labels = percent, name = "")

  fleet_importance_plot <- outcomes %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(mpa_effect, fill = fleet_model)) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    geom_histogram(show.legend = FALSE, color = "black") +
    facet_wrap(~fleet_model) +
    scale_fill_npg() +
    labs(y = "# of Sims") +
    scale_x_percent(name = "MPA Conservation Effect")


  # fleet_importance_plot <- outcomes %>%
  #   group_by(experiment) %>%
  #   filter(years_protected == max(years_protected)) %>%
  #   ggplot(aes(mpa_effect, fill = fleet_model)) +
  #   geom_vline(aes(xintercept = 0), linetype = 2) +
  #   geom_histogram(show.legend = FALSE, color = "black", aes(alpha = mpa_effect)) +
  #   facet_wrap(~fleet_model) +
  #   scale_fill_npg() +
  #   labs(y = "# of Sims") +
  #   scale_x_percent(name = "MPA Conservation Effect")
  #


  cc_outcomes <- outcomes %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected), fleet_model == "constant-catch") %>%
    mutate(conservation_loss = as.factor(mpa_effect < 0)) %>%
    select(-mpa_effect,-pop_effect)


  # negative_tree = train(
  #   conservation_loss ~  adult_movement + larval_movement + size_limit + f_v_m + mpa_size + factor(density_dependence_form) + density_movement_modifier,
  #   method = "rpart",
  #   data =  cc_outcomes
  # )


  # class(rpart.plot::rpart.plot(negative_tree$finalModel))

  cc_importance_plot <- outcomes %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected), fleet_model == "constant-catch", mpa_effect < 0) %>%
    ggplot(aes(mpa_size, adult_movement)) +
    geom_point(alpha = 0.5, size = 4) +
    geom_density2d(size = 2, color = "tomato",alpha = 0.75) +
    scale_y_percent(name = "Adult Movement") +
    scale_x_percent(name = "MPA Size")


  pop_size_plot <- outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(mpa_size, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      trans = plot_trans,
      guide = hgc,
      name = "# of Sims"
    ) +
    scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Conservation Effect") +
    theme(legend.position = "top")


  pop_facet_effect_plot <- outcomes %>%
    filter(years_protected == short_frame) %>%
    select(mpa_effect, mpa_size, depletion) %>%
    gather(measure, value, -mpa_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      trans = plot_trans,
      guide = hgc,
      name = "# of Sims"
    ) +
    facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Conservation Effect")

  pop_log10_facet_effect_plot <- outcomes %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ungroup() %>%
    select(mpa_effect, mpa_size, depletion) %>%
    gather(measure, value, -mpa_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, mpa_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      trans = "log10",
      guide = gc,
      name = "# of Sims"
    ) +
    facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Conservation Effect")


  pop_f_facet_effect_plot <- outcomes %>%
    filter(years_protected == short_frame) %>%
    select(mpa_effect, mpa_size, f_v_m) %>%
    gather(measure, value,-mpa_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, mpa_effect),
               binwidth = c(.05, .05),
               show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      trans = plot_trans,
      guide = gc,
      name = "# of Sims"
    ) +
    facet_wrap( ~ measure,
                strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Conservation Effect")



  pop_depletion_and_size_plot <- outcomes %>%
    filter(years_protected >= 0) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(mpa_effect)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_tile() +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = hgc,
                       name = "Median Effect") +
    labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent) +
    theme(legend.position = "top")

  eq_pop_depletion_and_size_plot <- outcomes %>%
    filter(years_protected >= 0) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == max(years_protected)], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == max(years_protected)) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(mpa_effect)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_tile() +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = hgc,
                       name = "Median Pop. Effect") +
    labs(x = "Equilibrium Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent)


  eq_pop_f_and_size_plot <- outcomes %>%
    filter(years_protected >= 0) %>%
    group_by(experiment) %>%
    mutate(f = plyr::round_any(f_v_m[years_protected == 0], .1),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == max(years_protected)) %>%
    group_by(f, mpa_size) %>%
    summarise(median_mpa_effect = median(mpa_effect)) %>%
    ggplot(aes(f, mpa_size, fill = median_mpa_effect)) +
    geom_tile() +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = gc,
                       name = "Median Effect") +
    labs(x = "EQ F/M", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous()



  pop_depletion_and_size_plot <- outcomes %>%
    filter(years_protected >= 0) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(mpa_effect)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_tile() +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_viridis(labels = percent,
                       guide = hgc,
                       name = "Median Effect") +
    labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent) +
    theme(legend.position = "top")

  pop_combo_plot <-
    (pop_depletion_and_size_plot + labs(title = "A")) + ((pop_size_plot + labs(title = "B")) / pop_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), units = "lines"),
      axis.text.x = element_text(size = 8),
      legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.y = element_text(size = 10)
    )


  ## ----fishery-effects,fig.cap = "Median (A) and range (B) MPA fishery effects, expressed as the difference in catch with and without MPAs  as a proportion of MSY, after 15 years of protection. For (A), X-axes indicate the pre-MPA depletion of the fishery, where depletion is the percentage of unfished biomass that has been removed from the population, and Y-axes is the percent of the population's range encompasssed inside an MPA. For B), y-axes show the regional conservation effect. Constant-catch scenarios are not included in this plot since by definition catches are equal with or without MPAs", include = FALSE----

  msy_depletion_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(depletion, msy_effect), binwidth = c(.05, .05),
               show.legend = FALSE) +
    scale_fill_viridis(
      option = "A",
      guide = hgc,
      name = "Median Effect"
    ) +
    scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
    scale_y_continuous(labels = percent, name = "Fishery Effect")


  msy_size_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(mpa_size, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = gc,
      name = "# of Sims"
    ) +
    scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Fishery Effect") +
    theme(legend.position = "right")

  msy_facet_effect_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    select(msy_effect, mpa_size, depletion) %>%
    gather(measure, value, -msy_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, msy_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = hgc,
      name = "# of Sims"
    ) +
    facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Fishery Effect")


  msy_depletion_and_size_plot <- fishery_outcomes %>%
    filter(years_protected >= 0,!is.na(msy)) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(msy_effect, na.rm = T)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_raster(interpolate = TRUE) +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_gradient2(midpoint = 0,
                         low = "tomato",
                         high = "steelblue",
                         mid = "white",
                         labels = percent,
                         guide = hgc,
                         name = "Median Effect") +
    labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent, limits = c(0,NA))

  msy_fishery_combo_plot <-
    (msy_depletion_and_size_plot + labs(title = "A")) + ((msy_size_plot + labs(title = "B")) / msy_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
      plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.x = element_text(size = 8),
      legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.y = element_text(size = 10)
    )



  catch_depletion_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(depletion, fishery_effect), binwidth = c(.05, .05),
               show.legend = FALSE) +
    scale_fill_viridis(
      option = "A",
      guide = hgc,
      name = "Median Effect"
    ) +
    scale_x_continuous(labels = percent, name = "Pre-MPA Depletion") +
    scale_y_continuous(labels = percent, name = "Fishery Effect")


  catch_size_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    ggplot() +
    geom_bin2d(aes(mpa_size, fishery_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = hgc,
      name = "# of Sims"
    ) +
    scale_x_continuous(labels = percent, name = "Range in MPA") +
    scale_y_continuous(labels = percent, name = "Fishery Effect") +
    theme(legend.position = "right")

  catch_facet_effect_plot <- fishery_outcomes %>%
    filter(years_protected == short_frame) %>%
    select(fishery_effect, mpa_size, depletion) %>%
    gather(measure, value, -fishery_effect) %>%
    ggplot() +
    geom_bin2d(aes(value, fishery_effect), binwidth = c(.05, .05), show.legend = TRUE) +
    scale_fill_viridis(
      option = "A",
      guide = hgc,
      name = "# of Sims"
    ) +
    facet_wrap(~measure, labeller = labeller(measure = facet_labels), strip.position = "bottom") +
    scale_x_percent(name = "") +
    scale_y_percent(name = "Fishery Effect")


  catch_depletion_and_size_plot <- fishery_outcomes %>%
    filter(years_protected >= 0,!is.na(msy)) %>%
    group_by(experiment) %>%
    mutate(depletion = plyr::round_any(depletion[years_protected == 0], .05),
           mpa_size = plyr::round_any(mpa_size, .05)) %>%
    filter(years_protected == short_frame) %>%
    group_by(depletion, mpa_size) %>%
    summarise(median_mpa_effect = median(fishery_effect, na.rm = T)) %>%
    ggplot(aes(depletion, mpa_size, fill = median_mpa_effect)) +
    geom_raster(interpolate = FALSE) +
    # geom_contour(aes(z = median_mpa_effect)) +
    scale_fill_gradient2(midpoint = 0,
                         low = "tomato",
                         high = "steelblue",
                         mid = "white",
                         labels = percent,
                         guide = hgc,
                         name = "Median Effect") +
    labs(x = "Pre-MPA Depletion", y = "Range in MPA") +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent, limits = c(0,NA))

  catch_fishery_combo_plot <-
    (catch_depletion_and_size_plot + labs(title = "A")) + ((catch_size_plot + labs(title = "B")) / catch_depletion_plot)  + plot_layout(widths = c(1.5, 1)) & theme(
      plot.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.x = element_text(size = 8),
      legend.box.margin = unit(c(0, 0, 0, 0), units = "lines"),
      axis.text.y = element_text(size = 10)
    )




  ## ----density-ratio-plot--------------------------------------------------


  unbiased_dr <- density_ratios %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(1 - depletion, pmin(10,true_density_ratio - 1))) +
    geom_point(size = 2, alpha = 0.75)  +
    labs(x = "Percent of Unfished Biomass", y = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density of identical scenario without MPAs",
         title = "'Unbiased' Density Ratio") +
    scale_x_percent() +
    scale_y_percent()

  biased_dr <- density_ratios %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(1 - depletion, pmin(10,biased_density_ratio - 1))) +
    geom_point() +
    labs(x = "Percent of Unfished Biomass", y = "Percent Difference in Inside/Outside MPA Densities", caption = "Density ratios calculated as density inside MPAs relative to density outside MPAs, weighted by distance from MPA border",
         title = "'Biased' Density Ratio") +
    scale_x_percent() +
    scale_y_percent()




  ## ----dr-vs-mpa-effect-plot-----------------------------------------------

  unbiased_dr_plot <- density_ratios %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(y = pmin(4,mpa_effect), x = pmin(4,true_density_ratio - 1))) +
    geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
    labs(y = "True % Biomass Difference", x = "% Difference in Inside/Outside Biomass") +
    scale_x_percent() +
    scale_y_percent() +
    geom_smooth(method = "lm") +
    geom_abline(aes(slope = 1, intercept = 0)) +
    scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)


  baci_plot <- density_ratios %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(y = mpa_density - nompa_density, x = baci)) +
    geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
    labs(y = "True Change in Density", x = "BACI Estimate of Density Change") +
    geom_smooth(method = "lm") +
    geom_abline(aes(slope = 1, intercept = 0)) +
    scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)

  biased_dr_plot <- density_ratios %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ggplot(aes(x =  pmin(4,biased_density_ratio - 1),y = pmin(4, mpa_effect))) +
    geom_point(aes(color = adult_movement),size = 2, alpha = 0.75)  +
    labs(y = "True % Biomass Difference", x = "% Difference in Inside/Outside Biomass") +
    scale_x_percent() +
    scale_y_percent() +
    geom_smooth(method = "lm") +
    geom_abline(aes(slope = 1, intercept = 0)) +
    scale_color_viridis(guide = gc, name = "Adult Movement", labels = percent)




  ## ------------------------------------------------------------------------
  fish_v_fishing <- processed_grid %>%
    select(-density_ratio) %>%
    unnest() %>%
    group_by(experiment) %>%
    mutate(pop_effect = pmin(1,(`with-mpa` - `no-mpa`) / `no-mpa`[year == min(year)])) %>%
    rename(fishery_effect = mpa_effect1) %>%
    group_by(experiment) %>%
    mutate(year = 1:length(year)) %>%
    ungroup() %>%
    mutate(years_protected = year - year_mpa + 1) %>%
    left_join(outcomes %>% select(experiment, year, depletion),
              by = c("experiment", "year")) #%>%
  # left_join(fishery_outcomes %>% select(experiment, year, msy_effect))

  fish_v_fishing %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected)) %>%
    ungroup() %>%
    ggplot(aes(pmin(2,fishery_effect), pmin(2,mpa_effect))) +
    geom_bin2d() +
    scale_fill_viridis(guide = gc)


  labelfoo <- function(x){

    paste0("%MPA:",x)

  }

  fish_v_fishing_plot <- fish_v_fishing %>%
    group_by(experiment) %>%
    filter(years_protected == max(years_protected),
           fleet_model != "constant-catch") %>%
    ungroup() %>%
    # mutate(mpa_bins = cut_width(100*mpa_size, width = 25)) %>%
    mutate(mpa_bins = cut(100*mpa_size, breaks = seq(0,100, by = 25))) %>%
    ggplot(aes(x =pmin(1, pop_effect), y = pmin(1, fishery_effect))) +
    geom_hline(aes(yintercept = 0), size = 1) +
    geom_vline(aes(xintercept = 0), size = 1) +
    # geom_smooth(method = "lm") +
    geom_point(alpha = 0.25, aes(fill = mpa_bins), shape = 21, show.legend = FALSE, size = 2) +
    # geom_hex(bins = 10) +
    # geom_abline(aes(slope = 1, intercept = 0), color = "red", size = 1) +
    geom_abline(aes(slope = -1, intercept = 0), color = "steelblue", size = 1) +
    geom_smooth(aes(color = mpa_bins),method = "lm", linetype = 1, se = TRUE) +

    scale_x_percent(name = "% of Unfished Biomass Recovered") +
    scale_y_percent(name = "% Change in Catch") +
    # scale_(labels = percent, guide = gc) +
    scale_color_brewer(palette = "Reds",labels = percent, guide = gc) +
    scale_fill_brewer(palette = "Reds",labels = percent, guide = gc) +
    facet_wrap(~mpa_bins,
               labeller = labeller(mpa_bins = labelfoo))




  ## ----ci-map, fig.cap = "Map of study region and sampling locations. Shaded polygons indicate location of MPAs. Points represent sampling locations, and color indicates the number of observations recorded at a given point", include = FALSE----

  ci_map <- pisco_data %>%
    left_join(site_data, by = c("site","side")) %>%
    filter(is.na(eventual_mpa) == F) %>%
    filter(region %in% c("ANA", "SCI","SRI",'SMI'),
           classcode %in% top_species) %>%
    group_by(lat_wgs84, lon_wgs84, site, region, eventual_mpa) %>%
    count()



  ci_map <-  ci_map %>%
    dplyr::mutate(geometry = purrr::map2(lon_wgs84, lat_wgs84, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
    ungroup() %>%
    mutate(geometry = sf::st_sfc(geometry, crs = 4326)) %>%
    sf::st_sf()


  bbox <- sf::st_bbox(ci_map)

  ci_mpa_plot <-  ggmap(channel_islands) +
    geom_sf(data = ca_mpas, inherit.aes = FALSE, alpha = 0.5) +
    coord_sf(xlim = c(bbox['xmin'] - .3, bbox['xmax'] + .06),
             ylim = c(bbox['ymin'] - .1, bbox['ymax'] + .4)) +
    labs(x = "Longitude", y = "Latitude")


  pisco_ci_map_plot <-  ggmap(channel_islands) +
    geom_sf(data = ca_mpas, inherit.aes = FALSE, alpha = 0.5) +
    geom_sf(data = ci_map,inherit.aes = FALSE,
            aes(color = n),
            alpha = 0.75,
            size = 1) +
    coord_sf(xlim = c(bbox['xmin'] - .3, bbox['xmax'] + .06),
             ylim = c(bbox['ymin'] - .1, bbox['ymax'] + .4)) +
    scale_color_gradient(low = "white", high = "orangered",
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth = 1,
                                                barheight = unit(13, "lines")),
                         name = "Samples")  +
    labs(x = "Longitude", y = "Latitude")

  # a = outcomes %>%
  #   group_by(experiment) %>%
  #   filter(years_protected == max(years_protected)) %>%
  #   ungroup() %>%
  #   select(
  #     -experiment,
  #     -msy,
  #     -b_msy,
  #     -year,
  #     -mpa_size1,-`no-mpa`,
  #     -`with-mpa`,
  #     -years_protected,
  #     -depletion,
  #     -scientific_name,
  #     -pop_effect,
  #     -b0    ) %>%
  #   na.omit()
  #
  #
  # d <- recipes::recipe(mpa_effect ~ ., data = a) %>%
  #   step_dummy(all_nominal())
  #
  # d <- prep(d, data = a) %>%
  #   juice() %>%
  #   mutate(rando = rnorm(nrow(.)))
  #
  #
  # b <- caret::train(mpa_effect ~ .,
  #                   data = d,
  #                   importance = "impurity_corrected",
  #                   method = "ranger")
  #
  #
  # imp <- b$finalModel$variable.importance
  #
  # i <- imp %>%
  #   broom::tidy() %>%
  #   mutate(names = fct_reorder(names, x))
  #
  #
  # importance_plot <- i %>%
  #   filter(x >= x[names == "rando"]) %>%
  #   top_n(15,x) %>%
  #   ggplot(aes(names,x)) +
  #   geom_col(fill = "steelblue", color = "black") +
  #   coord_flip() +
  #   labs(x = "", y = "Importance to MPA Conservation Effect")

  # ggsave(
  #   filename =  file.path(fig_dir, "channel-islands.pdf"),
  #   plot  = pisco_ci_map_plot,
  #   width = fig_width,
  #   height = fig_height,
  #   useDingbats = TRUE
  # )

  # pdf(
  #   file =  file.path(fig_dir, "channel-islands.pdf")
  #   ,
  #   width = 8,
  #   height = 5,
  #   useDingbats = TRUE
  # )
  # print(pisco_ci_map_plot)
  # dev.off()



  ## ------------------------------------------------------------------------
  ci_effects <- short_term %>%
    filter(years_protected == max(years_protected),
           mpa_size  <= 0.25, mpa_size >= 0.15)

  med_ci_effect <- median(ci_effects$mpa_effect)

  # rethinking::HPDI(ci_effects$mpa_effect, prob = 0.5)

  quart_ci_effect <- quantile(ci_effects$mpa_effect)



  # plot the actual data ---------------------------------------------------------------

  used_data <- pisco_abundance_data
  # used_data <- base_run$data[[1]] %>%
  #   left_join(life_history_data %>% select(classcode, commonname), by = "classcode")
  # 

  species_plot <- used_data %>%
    group_by(commonname) %>%
    summarise(obs = sum(any_seen),
              targeted = unique(targeted)) %>%
    mutate(commonname = fct_reorder(commonname, obs)) %>%
    ggplot(aes(commonname,obs, fill = targeted == 1)) +
    geom_col(color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "") +
    labs(x= "", y = "Observations")



  ## ----raw-trend, fig.cap="Centered and scaled mean annual density of included species (faded lines) and smoothed means of targeted and non-targeted groups, and mean (darker lines) and 95% confidence interval of the mean (ribbon) over time", include = FALSE----
  raw_did_plot <- used_data %>%
    select(year, targeted,classcode, log_density, any_seen, region) %>%
    mutate(density = exp(log_density) * any_seen) %>%
    # group_by(region, classcode) %>%
    # mutate(density = scale(density)) %>%
    group_by(classcode,targeted, year) %>%
    summarise(md = mean(density)) %>%
    group_by(classcode,targeted) %>%
    mutate(smd = (md - mean(md)) / sd(md)) %>%
    ggplot(aes(year, smd, color = targeted == 1, fill = targeted == 1)) +
    geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
    geom_line(aes(group = interaction(targeted, classcode)), alpha = 0.5) +
    geom_smooth(size = 2) +
    labs(y = "Centered and Scaled Mean Density", x = "Year") +
    scale_color_manual(values = c("steelblue", "tomato") ,name = "Targeted?") +
    scale_fill_manual(values = c("steelblue", "tomato") ,name = "Targeted?") +
    labs(y = "Centered and Scaled Mean Density", x = "Year")


  raw_did_plot <- raw_did_plot +
    scale_color_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "") +
    scale_fill_manual(values = c("steelblue", "tomato"), labels = c("Non-Targeted", "Targeted"), name = "")


  # ggsave(
  #   filename =  file.path(fig_dir, "pop-trends.pdf"),
  #   plot  = raw_did_plot,
  #   width = fig_width,
  #   height = fig_height,
  #   useDingbats = TRUE
  # )

  # pdf(
  #   file =  file.path(fig_dir, "pop-trends.pdf"),
  #   width = 8,
  #   height = 5,
  #   useDingbats = TRUE
  # )
  # print(raw_did_plot)
  # dev.off()




  ## ----did-plot, fig.cap = "Estimated divergence in biomass densities of targeted and non-targeted fishes throughout the Channel Islands (i.e. integrated across inside and outside of MPAs). MPAs are implemented in 2003 (red dashed line). Estimates are from a regression on log(abundance index), so estimated effects roughly correspond to percentage changes", include = FALSE----

  # did_plot <- did_plot +
  #   scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent, limits = c(-.75,1))

  # doh_did_plot <- doh_did_plot +
  #   scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent, limits = c(-.75,1))

  # EDM ---------------------------------------------------------------------



  cip_data <- pisco_data %>%
    left_join(site_data, by = c("site","side")) %>%
    filter(is.na(eventual_mpa) == F) %>%
    filter(region %in% c("ANA", "SCI","SRI",'SMI'),
           classcode %in% top_species)

  dat <- cip_data %>%
    group_by(broadtrophic, region, site, side, transect, year, eventual_mpa) %>%
    summarise(density = sum(density_g_m2)) %>%
    group_by(broadtrophic, region, site, year) %>%
    summarise(mean_density = mean(density, na.rm = T)) %>%
    group_by(broadtrophic,region, year) %>%
    summarise(regional_density = mean(mean_density, na.rm = T)) %>%
    ungroup() %>%
    spread(broadtrophic, regional_density) %>%
    group_by(region) %>%
    arrange(year) %>%
    mutate(lag4_carnivore = lag(carnivore,4),
           lag4_piscovore = lag(piscivore,4)) %>%
    na.omit() %>%
    ungroup() %>%
    filter(year > 1999)

  dat <- dat %>%
    select(-lag4_carnivore,-lag4_piscovore)%>%
    group_by(region)%>%
    mutate_at(vars(carnivore:planktivore),funs(norm=(.-mean(.))/sd(.)))->dat

  group_trends_plot <- dat %>%
    select(region,year,contains("norm"))%>%
    gather("group","density",-region,-year)%>%
    ggplot(aes(year,density,col=group))+
    geom_line()+
    scale_color_locuszoom(name="Trophic Troup",labels=c("carnivore","herbivore","piscivore","planktivore"))+
    labs(x="year",y="Centered and Scaled Density")+
    facet_wrap(~region)

  ## split data by region to analye timeseries from each island separately
  ana.dat <- dat %>% filter(region=="ANA")
  sci.dat <- dat %>% filter(region=="SCI")
  smi.dat <- dat %>% filter(region=="SMI")
  sri.dat <- dat %>% filter(region=="SRI")

  datnest <- dat %>% group_by(region) %>% nest()

  datnorm <- dat %>% select(region,year,contains('norm')) %>% ungroup()

  # have to record the segments corresponding to each "replicate" so simplex algorithm does not try to make predictions crossing time barriers
  segs <- datnorm %>% mutate(ind=row_number()) %>% group_by(region) %>% summarise(first=first(ind),last=last(ind)) %>%
    select(-region)

  var_names <- c("carnivore_norm","herbivore_norm","piscivore_norm","planktivore_norm")


  regions.combined.simp.list <- map(var_names,function(x){
    temp <- datnorm %>% ungroup() %>% select(matches(x)) %>% as.data.frame()
    out <- simplex(as.numeric(temp[,1]),E=1:10,lib=as.matrix(segs),silent=T) %>%
      mutate(trophic=x)
    out
  })


  embed_plot <- bind_rows(regions.combined.simp.list) %>%
    ggplot(aes(E,rho,color=trophic))+
    geom_line(size=2)+
    facet_wrap(~trophic,nrow=2,scales="free_y")+
    geom_hline(yintercept = 0,color="black")+
    labs(x="Embedding Dimension (E)",y=expression(paste("Skill, ",rho)))+
    scale_x_continuous(breaks=seq(0,12,by=2))+
    scale_color_locuszoom()+
    guides(color=F)

  tempE <- 8

  temp <- ccm(datnorm,lib=as.matrix(segs),pred=as.matrix(segs),E=tempE,lib_column= 'carnivore_norm',target_column = 'herbivore_norm',lib_sizes = c(10,25,50,75),num_samples=100,replace=T,silent=T,RNGseed = 41389)

  c_xmap_h <- temp %>%
    group_by(lib_size)%>%
    summarise(rhomean=mean(rho,na.rm=T),upper=quantile(rho, 0.975),lower=quantile(rho, 0.025))%>%
    ungroup()%>%
    ggplot(aes(lib_size,rhomean))+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")+
    geom_line(color="darkorchid3")+
    geom_hline(aes(yintercept = 0), linetype = 2) +
    labs(x="",y=expression(paste(rho, " (predictive skill)")),title="Carnivores on Herbivores")

  # cross map herbivore to carnivore
  # inspect the output of simplex from the previous step and use the best embedding dimension (highest rho) for the carnivore time series

  tempE <- 6

  temp <- ccm(datnorm,lib=as.matrix(segs),pred=as.matrix(segs),E=tempE,lib_column= 'herbivore_norm',target_column = 'carnivore_norm',lib_sizes = c(10,25,50,75),num_samples=100,replace=T,silent=T,RNGseed = 41389)

  h_xmap_c <- temp %>%
    group_by(lib_size)%>%
    summarise(rhomean=mean(rho,na.rm=T),upper=quantile(rho, 0.975),lower=quantile(rho, 0.025))%>%
    ungroup()%>%
    ggplot(aes(lib_size,rhomean))+
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")+
    geom_line(color="darkorchid3")+
    geom_hline(aes(yintercept = 0), linetype = 2) +
    labs(x="Library Size",y=expression(paste(rho, " (predictive skill)")),title="Herbivores on Carnivores")

  cross_map_plot <- h_xmap_c + c_xmap_h



  # environment effects and MPA only-------------------------------------------------------------------------

  facet_labels <- c(
    `TRUE` = "MPA",
    `FALSE` = "Non-MPA"
  )


  species_mpa_trend_plot <- cip_data %>%
    group_by(year, classcode, eventual_mpa, targeted) %>%
    summarise(mean_density = mean(density_g_m2, na.rm = T)) %>%
    ungroup() %>%
    mutate(targeted = targeted == 1) %>%
    group_by(classcode, eventual_mpa) %>%
    mutate(cs_density = (mean_density - mean(mean_density)) / sd(mean_density)) %>%
    ungroup() %>%
    rename(`Eventual MPA?` = eventual_mpa) %>%
    ggplot(aes(year, cs_density, color = targeted, fill = targeted)) +
    geom_vline(aes(xintercept = 2003), linetype = 2, color = "red") +
    geom_line(aes(group = interaction(classcode, targeted)),alpha = 0.5) +
    geom_smooth() +
    scale_color_manual(values = c("steelblue", "tomato") ,name = "", labels = c("Non-Targeted","Targeted")) +
    scale_fill_manual(values = c("steelblue", "tomato") ,name = "",labels = c("Non-Targeted","Targeted")) +

    facet_wrap(~`Eventual MPA?`, labeller = labeller(`Eventual MPA?` = facet_labels)) +
    labs(y = "Centered and Scaled Mean Density", x = "Year")

  species_mpa_trend_plot



  # mpa_did_plot <- mpa_run$did_plot[[1]] +
  #   labs(x = "Year", y = "Estimate of Regional Effect",
  #        caption = "")

  temp_trends_plot <- cip_data %>%
    mutate(year_month = year + (month / 12 - .1)) %>%
    filter(region %in% c("ANA", "SCI","SRI","SMI")) %>%
    group_by(year,region) %>%
    summarise(mean_temperature = mean(mean_temp)) %>%
    ggplot(aes(year, mean_temperature, color = region)) +
    geom_line(size = 1) +
    scale_color_locuszoom(name = "") +
    scale_y_continuous(name =bquote("Mean Temperature"~(degree~ "C"))) +
    labs(x = "Year")



  # kfm example -------------------------------------------------------------

  # kfm_did_plot <- kfm_run$did_plot[[1]] +
  #   scale_y_continuous(name = "~Divergence from Non-Targeted", labels = percent) +
  #   labs(caption = "")


  # regional catches --------------------------------------------------------

  cdfw_catches <-
    read_csv(file = here::here("data", 'cdfw-catches.csv')) %>%
    group_by(sci_name, year) %>%
    summarise(catch = sum(pounds_caught, na.rm = T)) %>%
    left_join(
      life_history_data %>% select(taxa, classcode) %>% mutate(taxa = tolower(taxa)),
      by = c("sci_name" = "taxa")
    ) %>%
    filter(!is.na(classcode)) %>%
    left_join(life_history_data %>% select(classcode, commonname), by = "classcode") %>%
    filter(classcode %in% abundance_data$data[[1]]$classcode)

  catch_trend_plot <- cdfw_catches %>%
    group_by(classcode) %>%
    mutate(catch = catch * 0.000453592) %>%
    mutate(scatch = scale(catch),
           ncatch = length(catch),
           mcatch = mean(catch)) %>%
    ungroup() %>%
    filter(ncatch > 10) %>%
    ggplot(aes(year, scatch)) +
    geom_vline(aes(xintercept = 2003), color = 'red', linetype = 2) +
    geom_line(show.legend = TRUE, aes(size = mcatch)) +
    geom_point(show.legend = TRUE, aes(size = mcatch)) +
    facet_wrap(~commonname) +
    theme_minimal() +
    labs(y = "Centered and Scaled Catch", x = "Year") +
    scale_size_continuous(name = "Mean Catch (tons)")




  # simulation testing ------------------------------------------------------

  pisco_sim_plot <- pisco_performance$mixed_effect_did +
    labs(title = "", x= "Year") +
    scale_y_percent(name = "MPA Effect")


  simple_sim_plot <- simple_performance$mixed_effect_did +
    labs(title = "", x= "Year") +
    scale_y_percent(name = "MPA Effect")



  # case studies ------------------------------------------------------------

  fish <- create_fish(r0 = 100, adult_movement = 2, larval_movement = 10)

  fleet <- create_fleet(fish = fish,
                        initial_effort = 10,
                        q = 1,
                        target_catch = 25,
                        fleet_model = "constant-effort")

  experiment <- spasm::mpa_counterfactual(
    fish = fish,
    fleet = fleet,
    year_mpa = 25,
    mpa_size = .5,
    sim_years = 50,
    burn_years = 1,
    num_patches = 50,
    random_mpas = FALSE,
    min_size = 0.1
  )

  raw <- experiment$raw_outcomes %>%
    select(year,
           patch,
           biomass,
           biomass_caught,
           profits,
           effort,
           experiment,
           mpa) %>%
    rename(catch = biomass_caught) %>%
    gather(metric, value, -year,-patch,-experiment,-mpa) %>%
    group_by(year, patch, metric, experiment) %>%
    summarise(value = sum(value),
              mpa = unique(mpa)) %>%
    group_by(year, patch) %>%
    mutate(mpa = unique(mpa[experiment == "with-mpa"])) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(value = value / max(value)) %>%
    ungroup() %>%
    spread(experiment, value) %>%
    ungroup() %>%
    mutate(delta = `with-mpa` / `no-mpa`,
           ref = 1)

  doh <- raw %>%
    filter(metric %in% c("biomass","effort")) %>%
    ggplot() +
    geom_col(
      aes(patch, `with-mpa`, fill = mpa),
      width = 1
    ) +
    geom_line(aes(patch, `no-mpa`), size = 2,
              color = "black") +
    facet_wrap(~ metric) +
    coord_polar() +
    gganimate::transition_time(year) +
    gganimate::ease_aes('linear') +
    labs(title = 'Year: {frame_time}', x = "'",  y = "") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    scale_fill_npg(labels = c("Witout MPAs","With MPAs"), name = "")


  file.path(fig_dir,"simple_cs.gif")

  gganimate::anim_save(animation = doh,filename = file.path(fig_dir,"simple_cs.gif")
  )

  raw_sum <- raw %>%
    group_by(year, metric) %>%
    summarise(no_mpa = sum(`no-mpa`),
              with_mpa = sum(`with-mpa`)) %>%
    ungroup() %>%
    gather(experiment, value, no_mpa:with_mpa) %>%
    filter(year > 20)

  ribbon <- raw_sum %>%
    group_by(year, metric) %>%
    summarise(ymin = min(value),
              ymax = max(value),
              delta = value[experiment == "with_mpa"] - value[experiment == "no_mpa"]) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(max_delta = last(delta)) %>%
    ungroup()

  simple_cs_plot <- raw_sum %>%
    ungroup() %>%
    ggplot() +
    geom_line(aes(year, value, color = experiment), size = 1.5) +
    geom_ribbon(data = ribbon, aes(
      year,
      ymin = ymin,
      ymax = ymax,
      fill = max_delta
    ),
    alpha = 0.5,
    show.legend = FALSE) +
    facet_wrap( ~ metric, scales = "free_y",
                labeller = labeller(metric = tools::toTitleCase)) +
    scale_fill_gradient(low = "tomato", high = "steelblue") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(1,"lines")) +
    labs(x = "Year") +
    scale_color_npg(labels = c("Without MPAs","With MPAs"), name = '')


  # trickier case study

  fish <-
    create_fish(r0 = 100,
                adult_movement = 2,
                larval_movement = 10,
                density_dependence_form = 1)

  fleet <- create_fleet(fish = fish,
                        initial_effort = 10,
                        q = 1,
                        target_catch = 25,
                        fleet_model = "open-access",
                        max_cr_ratio = 0.5,
                        effort_allocation = "profit-gravity"
  )

  experiment <- spasm::mpa_counterfactual(
    fish = fish,
    fleet = fleet,
    year_mpa = 25,
    mpa_size = .5,
    sim_years = 50,
    burn_years = 1,
    num_patches = 50,
    random_mpas = TRUE,
    min_size = 0.1,
    sprinkler = FALSE,
    mpa_habfactor = 1
  )

  raw <- experiment$raw_outcomes %>%
    select(year,
           patch,
           biomass,
           biomass_caught,
           profits,
           effort,
           experiment,
           mpa) %>%
    rename(catch = biomass_caught) %>%
    gather(metric, value, -year,-patch,-experiment,-mpa) %>%
    group_by(year, patch, metric, experiment) %>%
    summarise(value = sum(value),
              mpa = unique(mpa)) %>%
    group_by(year, patch) %>%
    mutate(mpa = unique(mpa[experiment == "with-mpa"])) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(value = value / max(value)) %>%
    ungroup() %>%
    spread(experiment, value) %>%
    ungroup() %>%
    mutate(delta = `with-mpa` / `no-mpa`,
           ref = 1)

  doh <- raw %>%
    filter(metric %in% c("biomass","effort")) %>%
    ggplot() +
    geom_col(
      aes(patch, `with-mpa`, fill = mpa),
      width = 1
    ) +
    geom_line(aes(patch, `no-mpa`), size = 2,
              color = "black") +
    facet_wrap(~ metric) +
    coord_polar() +
    gganimate::transition_time(year) +
    gganimate::ease_aes('linear') +
    labs(title = 'Year: {frame_time}', x = "'",  y = "") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    scale_fill_npg(labels = c("Fished","MPA"), name = "")


  gganimate::anim_save(animation = doh,filename = file.path(fig_dir,"complex_cs.gif")
  )

  raw_sum <- raw %>%
    group_by(year, metric) %>%
    summarise(no_mpa = sum(`no-mpa`),
              with_mpa = sum(`with-mpa`)) %>%
    ungroup() %>%
    gather(experiment, value, no_mpa:with_mpa) %>%
    filter(year > 20)

  ribbon <- raw_sum %>%
    group_by(year, metric) %>%
    summarise(ymin = min(value),
              ymax = max(value),
              delta = value[experiment == "with_mpa"] - value[experiment == "no_mpa"]) %>%
    ungroup() %>%
    group_by(metric) %>%
    mutate(max_delta = last(delta)) %>%
    ungroup()

  complex_cs_plot <- raw_sum %>%
    ungroup() %>%
    ggplot() +
    geom_line(aes(year, value, color = experiment), size = 1.5) +
    geom_ribbon(data = ribbon, aes(
      year,
      ymin = ymin,
      ymax = ymax,
      fill = max_delta
    ),
    alpha = 0.5,
    show.legend = FALSE) +
    facet_wrap( ~ metric, scales = "free_y",
                labeller = labeller(metric = tools::toTitleCase)) +
    scale_fill_gradient(low = "tomato", high = "steelblue") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(1,"lines")) +
    labs(x = "Year") +
    scale_color_npg(labels = c("Without MPAs","With MPAs"), name = '')



  # save plots --------------------------------------------------------------


  rm(diagnostic_plots)

  plots <- ls()[str_detect(ls(),"_plot")]



  savefoo <- function(fig,
                      device = "pdf",
                      fig_width = 6,
                      fig_height = 5) {
    ggsave(
      filename =  file.path(fig_dir, paste(fig, device, sep = '.')),
      plot  = get(fig),
      width = fig_width,
      height = fig_height
    )

  }

  walk(plots, safely(savefoo), device = device, fig_height = fig_height,
       fig_width = fig_width)

}


# knit paper --------------------------------------------------------------
if (knit_paper == TRUE){

  rmarkdown::render(here::here("documents","ovando-regional-effects-of-mpas.Rmd"), params = list(run_name = run_name))

  rmarkdown::render(here::here("documents","ovando-regional-effects-of-mpas-si.Rmd"), params = list(run_name = run_name))

}

