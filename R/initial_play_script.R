# initial script to look at and get a feel for the data

# packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sf)
library(tmap)
library(concaveman)
library(mobr)

# read data in
dat <- read_csv("Data/UMBC_Swan_20211029_w_City.csv")


# total number of species
length(unique(dat$scientific_name))

# does scientific name and common name match
length(unique(dat$scientific_name))==length(unique(dat$common_name))


# number of landuses
unique(dat$fia_lu_abbr)

# cities
table(dat$evalid)

# lets visualize the spatial distribution of
# the experimental treatment within each city
# first convert to sf object
# and name everything that isn't forest to "other"
dat_sf <- dat %>%
  mutate(treatment=ifelse(fia_lu_abbr=="Forest", fia_lu_abbr, "Other")) %>%
  st_as_sf(coords=c("lon", "lat"), crs=3426)

tm_shape(dat_sf)+
  tm_dots(size=1, col="treatment")+
  tm_facets(by="evalid")

# see how many sites per each city
dat %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  left_join(., dat %>%
              unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
              dplyr::select(site_coords) %>%
              distinct() %>%
              mutate(site_id=1:nrow(.)), by="site_coords") %>%
  mutate(treatment=ifelse(fia_lu_abbr=="Forest", fia_lu_abbr, "Other")) %>%
  group_by(evalid, treatment) %>%
  summarize(N=length(unique(site_id)))

# let's see how different the cities are in size
city_size_function <- function(city_name){
  
  area <- dat_sf %>%
    dplyr::filter(evalid == city_name) %>%
    concaveman(., concavity=1) %>%
    st_area() %>%
    as.numeric()
  
  out <- data.frame(`area_us_survey_foot^2`=area,
                    city=city)
  
  return(out)
  
}

city_area <- bind_rows(lapply(unique(dat$evalid), city_size_function))


# data for analysis
analysis_dat <- dat %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  left_join(., dat %>%
              unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
              dplyr::select(site_coords) %>%
              distinct() %>%
              mutate(site_id=1:nrow(.)), by="site_coords") %>%
  mutate(site_id=as.character(site_id)) %>%
  mutate(treatment=ifelse(fia_lu_abbr=="Forest", fia_lu_abbr, "Other")) %>%
  dplyr::filter(!evalid %in% c("SanDiego2017Curr", "WashingtonDC2018Curr"))

unique(analysis_dat$evalid)

# a potential method to get species x site matrix and calculate mob stuff
# A resampling approach
# it takes each city separately (evalid)
# and then depending on the total number of observations of the 'treatment'
# this is used for the resampling
# Forest is always < OTHER
# so we can pick a number of sites to randomly sample
resampling_beta_method <- function(city_name, number_of_points){
  
  message(paste0("City identifier: ", city_name))
  
  dat_tmp <- analysis_dat %>%
    dplyr::filter(evalid==city_name)
  
  # sample a number of points from each forest and "other"
  aggregate_samples_function <- function(number_of_points, draw_number){
    
    message(paste0("Compiling draw number ", draw_number))
    
    # sample sites
    site_sample <- dat_tmp %>%
      dplyr::select(site_id, treatment) %>%
      distinct() %>%
      group_by(treatment) %>%
      sample_n(number_of_points)
    
    nrow(site_sample)==length(unique(site_sample$site_id))
    
    # get coords/shapefile of those random sites
    sample_coords <- dat_tmp %>%
      dplyr::filter(site_id %in% local(site_sample$site_id)) %>%
      dplyr::select(lat, lon, site_id, treatment) %>%
      distinct()
    
    area_all_sites <- sample_coords %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    area_forest <- sample_coords %>%
      dplyr::filter(treatment=="Forest") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    area_other <- sample_coords %>%
      dplyr::filter(treatment=="Other") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    centroid_forest <- sample_coords %>%
      dplyr::filter(treatment=="Forest") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_centroid() %>%
      st_coordinates() %>%
      data.frame()
    
    centroid_other <- sample_coords %>%
      dplyr::filter(treatment=="Other") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_centroid() %>%
      st_coordinates() %>%
      data.frame()
      
    # get the area and coords of the forest and 'other' treatments
    # into a single dataframe to read in below
    attributes_df <- data.frame(treatment=c("Forest", "Other"),
                                area_m2=c(area_forest, area_other),
                                centroid_lat=c(centroid_forest$Y, centroid_other$Y),
                                centroid_lon=c(centroid_forest$X, centroid_other$X))
    
    # then get the data from those randomly sampled sites
    # and turn it into a matrix for mobr format
    dat_sample <- dat_tmp %>%
      dplyr::filter(site_id %in% local(site_sample$site_id)) %>%
      group_by(treatment, scientific_name) %>%
      summarize(abund=sum(count_tree)) %>%
      ungroup() %>%
      mutate(sample_number=draw_number) %>%
      mutate(number_of_points=number_of_points) %>%
      left_join(., attributes_df) %>%
      mutate(equal_number_of_points=nrow(site_sample)==length(unique(site_sample$site_id)))
    
    return(dat_sample)
    
  }
  
  # now apply the above function 100 times
  # for now just for '5' sites aggregated at a time
  sampled_community_aggregated <- bind_rows(lapply(c(1:50), function(x){aggregate_samples_function(number_of_points, x)}))
  
  sampled_community_aggregated <- sampled_community_aggregated %>%
    ungroup() %>%
    left_join(., sampled_community_aggregated %>%
                dplyr::select(treatment, sample_number) %>%
                distinct() %>%
                mutate(agg_site_id=1:nrow(.)), by=c("treatment", "sample_number"))
  
  # now get sampled_community_aggregated ready for analysis in mob
  # get the community matrix
  env <- sampled_community_aggregated %>%
    dplyr::select(agg_site_id, centroid_lat, centroid_lon, treatment, area_m2) %>%
    distinct()
  
  comm <- sampled_community_aggregated %>%
    dplyr::select(agg_site_id, scientific_name, abund) %>%
    pivot_wider(names_from=scientific_name,
              values_from=abund,
              values_fill=0)
  
  row.names(comm) <- comm$agg_site_id
  comm <- comm[ , -1]
  comm[1:5, 1:5]
    
  tree_mob <- make_mob_in(comm, env, coord_names = c('centroid_lon', 'centroid_lat'),
                          latlong = TRUE)
  tree_mob
  
  ## multi-metric MoB analysis 
  stats <- get_mob_stats(tree_mob, group_var = 'treatment', n_perm = 1)
  
  final_summary_df <- stats$samples_stats %>%
    mutate(level="sample") %>%
    bind_rows(stats$groups_stats %>%
                mutate(level="groups")) %>%
    mutate(number_of_sites_aggregated=number_of_points) %>%
    mutate(city=city_name)
  
  return(final_summary_df)
    
  
}

number_points_5 <- bind_rows(lapply(unique(analysis_dat$evalid), function(x){resampling_beta_method(x, 5)}))


# first attempt at a plot
number_points_5 %>%
  dplyr::filter(level=="sample") %>%
  dplyr::filter(index %in% c("beta_S", "beta_S_n", "beta_S_PIE")) %>%
  ggplot(., aes(x=group, y=value, fill=group))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  coord_flip()+
  facet_wrap(index~city, scales="free")+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Beta diversity")+
  xlab("")

ggsave("temp_fig.png", width=8.5, height=6.6, units="in")




# QUESTIONS
#1) These coords for Austin are classified as both forest AND "other"
# 30.23515, -97.89663
b <- dat %>%
  dplyr::filter(evalid=="Austin2017Curr") %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  dplyr::filter(site_coords=="30.235153_-97.896629")

# 2) What does count_tree==0 mean? Only when saplings are recorded?

