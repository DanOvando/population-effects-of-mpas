make_sim_c_worthy <- function(subdata,
                          non_nested_vars,
                          center_and_scale = T,
                          fixed_did = T,
                          fixed_regions = T,
                          include_intercept = T,
                          numeric_species_key){


  x_non_nested <- subdata %>%
    select(non_nested_vars) %>%{
      if (include_intercept == T){
        mutate(.,intercept = 1)
      } else{
        .
      }
    } %>%
    spread_factors(drop_one = T) %>%{
      if (center_and_scale == T){
        purrrlyr::dmap(.,center_scale)

      } else
      {
        .
      }
    } %>%
    # mutate(cumulative_n_obs_2 = cumulative_n_obs ^ 2) %>%
    mutate(targeted = subdata$targeted)


 # x_enso <- subdata %>%
 #    select(classcode, enso) %>%
 #    mutate(index = 1:nrow(.)) %>%
 #    spread(classcode, enso, fill = 0) %>%
 #    arrange(index) %>%
 #    select(-index) %>%
 #    set_names(paste0("enso_",colnames(.)))

 # x_pdo <- subdata %>%
 #   select(classcode, pdo) %>%
 #   mutate(index = 1:nrow(.)) %>%
 #   spread(classcode, pdo, fill = 0) %>%
 #   arrange(index) %>%
 #   select(-index) %>%
 #   set_names(paste0("pdo_",colnames(.)))

 # x_non_nested <- x_non_nested %>%
 #   bind_cols(x_enso,
 #             x_pdo)

 x_did <- subdata %>%
   select(targeted, factor_year) %>%
   mutate(index = 1:nrow(.)) %>%
   spread(factor_year, targeted, fill = 0) %>%
   arrange(index) %>% {
     if (fixed_did == T){
       select(.,-(1:2))
     } else{
       select(., -1)
     }
   }


 x_year_species <- subdata %>%
   mutate(year_classcode = paste(classcode, year, sep = "-")) %>%
   select(year_classcode) %>%
   spread_factors(drop_one = F)

 species_index <- subdata$numeric_classcode

 year_species_index <-
   data_frame(classcode = colnames(x_year_species)) %>%
   mutate(classcode = str_extract_all(classcode, "(?<=-).*(?=-)", "")) %>%
   left_join(numeric_species_key, by = "classcode") %>%
   {
     .$numeric_classcode
   }

 x_region_cluster <- subdata %>%
   select(region, geographic_cluster) %>%
   mutate(region_cluster = paste(geographic_cluster, region, sep = "-")) %>%
   select(region_cluster) %>%
   spread_factors(drop_one = fixed_regions)

 region_cluster_index <-
   data_frame(region = colnames(x_region_cluster)) %>%
   mutate(region = str_replace_all(region, "(\\D)", "") %>% as.factor() %>% as.numeric()) %>%
   {
     .$region
   }



return(
  list(x_non_nested = x_non_nested,
       x_did = x_did,
       x_year_species = x_year_species,
       species_index = species_index,
       year_species_index = year_species_index,
       x_region_cluster = x_region_cluster,
       region_cluster_index = region_cluster_index
       )
)


}