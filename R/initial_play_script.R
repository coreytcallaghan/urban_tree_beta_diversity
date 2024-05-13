# initial script to look at and get a feel for the data

# packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sf)
library(tmap)
library(tmaptools)
library(concaveman)
library(mobr)
library(lmerTest)
library(lme4)
library(tidyr)
library(tibble)

#turn off the s2 processing via 
sf::sf_use_s2(FALSE)

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
                    city=city_name)
  
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
  dplyr::filter(!evalid %in% c("SanDiego2017Curr", "WashingtonDC2018Curr")) %>%
  group_by(site_coords) %>%
  mutate(number_landuse=length(unique(fia_lu_abbr))) %>%
  dplyr::filter(number_landuse==1)

unique(analysis_dat$evalid)


# summarize the data a little bit
analysis_dat %>%
  group_by(evalid) %>%
  summarize(number_sites=length(unique(site_id)),
            number_species=length(unique(scientific_name)))

analysis_dat %>%
  group_by(evalid, treatment) %>%
  summarize(number_sites=length(unique(site_id)),
            number_species=length(unique(scientific_name)))

# number of individuals per site_id
individuals_per_site <- analysis_dat %>%
  group_by(evalid, site_id, treatment) %>%
  summarize(number_individuals_sapling=sum(count_sapling),
            number_individuals_tree=sum(count_tree)) %>%
  mutate(total_individuals=number_individuals_sapling+number_individuals_tree)

mean(individuals_per_site$total_individuals)
median(individuals_per_site$total_individuals)
sd(individuals_per_site$total_individuals)

# make a map of the study sites?
analysis_dat_sf <- analysis_dat %>%
  dplyr::select(evalid, site_id, treatment, lon, lat) %>%
  distinct() %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326)

tm_shape(analysis_dat_sf)+
  tm_dots(size=1, col="treatment")+
  tm_facets(by="evalid")

port <- analysis_dat_sf %>%
  dplyr::filter(evalid=="PortlandOR2018Curr")

port_osm <- read_osm(port, ext=1.1)

portland_map <- tm_shape(port_osm)+
  tm_rgb()+
  tm_shape(port)+
  tm_dots(size=0.5, col="treatment", palette=c(Other='cyan', Forest='green'),
          legend.show=FALSE)

portland_map

austin <- analysis_dat_sf %>%
  dplyr::filter(evalid=="Austin2017Curr")

austin_osm <- read_osm(austin, ext=1.1)

austin_map <- tm_shape(austin_osm)+
  tm_rgb()+
  tm_shape(austin)+
  tm_dots(size=0.5, col="treatment", palette=c(Other='cyan', Forest='green'),
          legend.show=FALSE)

austin_map

houston <- analysis_dat_sf %>%
  dplyr::filter(evalid=="Houston2017Curr")

houston_osm <- read_osm(houston, ext=1.1)

houston_map <- tm_shape(houston_osm)+
  tm_rgb()+
  tm_shape(houston)+
  tm_dots(size=0.5, col="treatment", palette=c(Other='cyan', Forest='green'),
          legend.show=FALSE)

houston_map

sanan <- analysis_dat_sf %>%
  dplyr::filter(evalid=="SanAntonio2018Curr")

sanan_osm <- read_osm(sanan, ext=1.1)

sanan_map <- tm_shape(sanan_osm)+
  tm_rgb()+
  tm_shape(sanan)+
  tm_dots(size=0.5, col="treatment", palette=c(Other='cyan', Forest='green'),
          legend.show=FALSE)

sanan_map

map_all <- tmap_arrange(sanan_map, houston_map, austin_map, portland_map, ncol=2, nrow=2)

tmap_save(map_all, "Figures/map_of_cities.png", width=7.5, height=6.9, units="in")

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
      as.numeric() %>%
      unique()
    
    area_forest <- sample_coords %>%
      dplyr::filter(treatment=="Forest") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric() %>%
      unique()
    
    area_other <- sample_coords %>%
      dplyr::filter(treatment=="Other") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric() %>%
      unique()
    
    centroid_forest <- sample_coords %>%
      dplyr::filter(treatment=="Forest") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_centroid() %>%
      st_coordinates() %>%
      data.frame() %>%
      distinct()
    
    centroid_other <- sample_coords %>%
      dplyr::filter(treatment=="Other") %>%
      st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_centroid() %>%
      st_coordinates() %>%
      data.frame() %>%
      distinct()
      
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
      summarize(abund=sum(count_tree),
                biomass=sum(biomass_tree),
                carbon=sum(carbon_tree)) %>%
      ungroup() %>%
      mutate(sample_number=draw_number) %>%
      mutate(number_of_points=number_of_points) %>%
      left_join(., attributes_df) %>%
      mutate(equal_number_of_points=nrow(site_sample)==length(unique(site_sample$site_id)))
    
    return(dat_sample)
    
  }
  
  # now apply the above function 100 times
  # for now just for '5' sites aggregated at a time
  sampled_community_aggregated <- bind_rows(lapply(c(1:1000), function(x){aggregate_samples_function(number_of_points, x)}))
  
  sampled_community_aggregated <- sampled_community_aggregated %>%
    ungroup() %>%
    left_join(., sampled_community_aggregated %>%
                dplyr::select(treatment, sample_number) %>%
                distinct() %>%
                mutate(agg_site_id=1:nrow(.)), by=c("treatment", "sample_number"))
  
  # now get sampled_community_aggregated ready for analysis in mob
  # get the community matrix
  env <- sampled_community_aggregated %>%
    group_by(agg_site_id, treatment) %>%
    summarize(biomass=sum(biomass),
              carbon=sum(carbon)) %>%
    left_join(., sampled_community_aggregated %>%
                dplyr::select(agg_site_id, centroid_lat, centroid_lon, treatment, area_m2) %>%
                distinct())
  
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
    mutate(scale="alpha") %>%
    mutate(agg_site_id=rep(env %>% arrange(treatment) %>% .$agg_site_id, 7)) %>%
    bind_rows(stats$groups_stats %>%
                mutate(scale="gamma")) %>%
    mutate(number_of_sites_aggregated=number_of_points) %>%
    mutate(city=city_name) %>%
    mutate(scale=ifelse(grepl("beta", index)==TRUE, "beta", scale)) %>%
    left_join(., tree_mob$env)
  
  return(final_summary_df)
    
}

number_points_5 <- bind_rows(lapply(unique(analysis_dat$evalid), function(x){resampling_beta_method(x, 5)}))
saveRDS(number_points_5, "intermediate_results/number_points_5_analysis.RDS")

number_points_10 <- bind_rows(lapply(unique(analysis_dat$evalid), function(x){resampling_beta_method(x, 10)}))
saveRDS(number_points_10, "intermediate_results/number_points_10_analysis.RDS")

number_points_15 <- bind_rows(lapply(unique(analysis_dat$evalid), function(x){resampling_beta_method(x, 15)}))
saveRDS(number_points_15, "intermediate_results/number_points_15_analysis.RDS")

number_points_20 <- bind_rows(lapply(unique(analysis_dat$evalid), function(x){resampling_beta_method(x, 20)}))
saveRDS(number_points_20, "intermediate_results/number_points_20_analysis.RDS")


number_points_5 <- readRDS("intermediate_results/number_points_5_analysis.RDS")
number_points_10 <- readRDS("intermediate_results/number_points_10_analysis.RDS")
number_points_15 <- readRDS("intermediate_results/number_points_15_analysis.RDS")
number_points_20 <- readRDS("intermediate_results/number_points_20_analysis.RDS")

# first attempt at a plot
number_points_5 %>%
  dplyr::filter(level=="sample") %>%
  #dplyr::filter(index %in% c("N", "beta_S", "beta_S_n", "beta_S_PIE")) %>%
  ggplot(., aes(x=group, y=value, fill=group))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  coord_flip()+
  facet_wrap(index~city, scales="free", ncol=4)+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Beta diversity")+
  xlab("")

ggsave("Figures/all_results.png", width=8.5, height=6.6, units="in")

# a figure just for N
number_points_5 %>%
  dplyr::filter(level=="sample") %>%
  dplyr::filter(index=="N") %>%
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

# run a model
# separately for each index
beta_s_n <- number_points_5 %>%
  dplyr::filter(index=="beta_S_n")

mod_beta_s_n <- lmer(value ~ group + (1|city), data=beta_s_n)
summary(mod_beta_s_n)
confint(mod_beta_s_n)

beta_s <- number_points_5 %>%
  dplyr::filter(index=="beta_S")

mod_beta_s <- lmer(value ~ group + (1|city), data=beta_s)
summary(mod_beta_s)
confint(mod_beta_s)

beta_s_pie <- number_points_5 %>%
  dplyr::filter(index=="beta_S_PIE")

mod_beta_s_pie <- lmer(value ~ group + (1|city), data=beta_s_pie)
summary(mod_beta_s_pie)
confint(mod_beta_s_pie)

# make a figure of the model results
summary(mod_beta_s_n)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column(var="term") %>%
  mutate(lwr_95=confint(mod_beta_s_n)[4]) %>%
  mutate(upr_95=confint(mod_beta_s_n)[4, 2]) %>%
  mutate(model="Beta_S_n") %>%
  bind_rows(summary(mod_beta_s)$coefficients %>%
              as.data.frame() %>%
              rownames_to_column(var="term") %>%
              mutate(lwr_95=confint(mod_beta_s)[4]) %>%
              mutate(upr_95=confint(mod_beta_s)[4, 2]) %>%
              mutate(model="Beta_S")) %>%
  bind_rows(summary(mod_beta_s_pie)$coefficients %>%
              as.data.frame() %>%
              rownames_to_column(var="term") %>%
              mutate(lwr_95=confint(mod_beta_s_pie)[4]) %>%
              mutate(upr_95=confint(mod_beta_s_pie)[4, 2]) %>%
              mutate(model="Beta_S_PIE")) %>%
  dplyr::filter(term=="groupOther") %>%
  mutate(term2="Other") %>%
  ggplot(., aes(x=model, y=Estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lwr_95, ymax=upr_95), width=0.5)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  coord_flip()+
  xlab("")+
  ylab("The beta diversity difference between other and forest")+
  geom_hline(yintercept=0, color="red", linetype="dashed")

ggsave("Figures/model_results.png", width=5.6, height=4.8, units="in")









# QUESTIONS
#1) These coords for Austin are classified as both forest AND "other"
# 30.23515, -97.89663
b <- dat %>%
  dplyr::filter(evalid=="Austin2017Curr") %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  dplyr::filter(site_coords=="30.235153_-97.896629")

# 2) What does count_tree==0 mean? Only when saplings are recorded?

# 3) are all the fia plots the same size?
