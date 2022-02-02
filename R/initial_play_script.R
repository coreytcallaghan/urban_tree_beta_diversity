# initial script to look at and get a feel for the data

# packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sf)
library(tmap)
library(concaveman)

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
city_size_function <- function(city){
  
  area <- dat_sf %>%
    dplyr::filter(evalid == city) %>%
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
resampling_beta_method <- function(city_name){
  
  dat_tmp <- analysis_dat %>%
    dplyr::filter(evalid==city_name)
  
  # sample a number of points from each forest and "other"
  sample_function <- function(number_of_points){
    
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
      
    # then get the data from those randomly sampled sites
    # and turn it into a matrix for mobr format
    dat_sample <- dat_tmp %>%
      dplyr::filter(site_id %in% local(site_sample$site_id)) %>%
      group_by(treatment, scientific_name) %>%
      summarize(abund=sum(count_tree)) %>%
      ungroup() %>%
      pivot_wider(names_from=scientific_name,
                  values_from=abund,
                  values_fill=0)
    
  }
  
  
  
}





# QUESTIONS
#1) These coords for Austin are classified as both forest AND "other"
# 30.23515, -97.89663
b <- dat %>%
  dplyr::filter(evalid=="Austin2017Curr") %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  dplyr::filter(site_coords=="30.235153_-97.896629")



