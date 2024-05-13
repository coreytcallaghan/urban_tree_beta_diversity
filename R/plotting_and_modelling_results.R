# script for plotting and modelling results

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

# how many sampling sites had <5 species?
analysis_dat %>%
  group_by(site_coords) %>%
  summarize(number_species=length(unique(scientific_name))) %>%
  dplyr::filter(number_species<3) %>%
  nrow(.)

# how many total sampling sites
length(unique(analysis_dat$site_coords))


# make a summary figure of the total overall data
analysis_dat %>%
  mutate(treatment=gsub("Other", "Built", treatment)) %>%
  group_by(scientific_name, evalid, treatment) %>%
  summarize(N=n()) %>%
  mutate(N=1) %>%
  ungroup() %>%
  pivot_wider(names_from=treatment, values_from=N, values_fill=0) %>%
  ungroup() %>%
  mutate(Both=rowSums(across(where(is.numeric)))) %>%
  mutate(Both=ifelse(Both==2, 1, 0)) %>%
  pivot_longer(!c("scientific_name", "evalid"), names_to="treatment", values_to="N") %>%
  group_by(evalid, treatment) %>%
  summarize(species_richness=sum(N)) %>%
  ggplot(., aes(fill=treatment, x=evalid, y=species_richness))+
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Species richness")+
  scale_fill_brewer(palette="Set1")

# get unique species in forest land use
# make a summary figure of the total overall data
analysis_dat %>%
  group_by(scientific_name, evalid, treatment) %>%
  summarize(N = n(), .groups = 'drop') %>%
  pivot_wider(names_from = treatment, values_from = N, values_fill = list(N = 0)) %>%
  mutate(Unique_Forest = ifelse(Forest > 0 & Other == 0, 1, 0),
         Unique_Other = ifelse(Other > 0 & Forest == 0, 1, 0)) %>%
  filter(Unique_Forest == 1 | Unique_Other == 1) %>%
  group_by(evalid) %>%
  summarize(Forest_Species = sum(Unique_Forest),
            Other_Species = sum(Unique_Other), .groups = 'drop') %>%
  pivot_longer(cols = c(Forest_Species, Other_Species), names_to = "treatment", values_to = "species_richness") %>%
  left_join(., analysis_dat %>%
              group_by(evalid) %>%
              summarize(N=length(unique(scientific_name)))) %>%
  mutate(percent_unique=(species_richness/N)*100)


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

# SUMMARIZE MODEL RESULTS
number_points_5 <- readRDS("intermediate_results/number_points_5_analysis.RDS")
number_points_10 <- readRDS("intermediate_results/number_points_10_analysis.RDS")
number_points_15 <- readRDS("intermediate_results/number_points_15_analysis.RDS")
number_points_20 <- readRDS("intermediate_results/number_points_20_analysis.RDS")

# first attempt at a plot
number_points_5 %>%
  dplyr::filter(scale=="alpha") %>%
  #dplyr::filter(index %in% c("N", "beta_S", "beta_S_n", "beta_S_PIE")) %>%
  ggplot(., aes(x=group, y=value, fill=group))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  coord_flip()+
  facet_wrap(index~city, scales="free", ncol=4)+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Alpha diversity")+
  xlab("")

ggsave("Figures/alpha_diversity.png", width=8.5, height=6.6, units="in")

# a figure just for N
number_points_5 %>%
  dplyr::filter(scale=="beta") %>%
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

ggsave("Figures/beta_diversity.png", width=8.5, height=6.6, units="in")

number_points_5 %>%
  dplyr::filter(scale=="gamma") %>%
  #dplyr::filter(index %in% c("N", "beta_S", "beta_S_n", "beta_S_PIE")) %>%
  ggplot(., aes(x=group, y=value, fill=group))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  coord_flip()+
  facet_wrap(index~city, scales="free", ncol=4)+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Gamma diversity")+
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


# differences between biomass and carbon for each of the cities
number_points_5 %>%
  dplyr::select(8:12) %>%
  dplyr::filter(complete.cases(.)) %>%
  mutate(city=c(rep("Austin", 2000), rep("Houston", 2000), rep("Portland", 2000), rep("San Antonio", 2000))) %>%
  dplyr::select(treatment, city, biomass, carbon) %>%
  pivot_longer(!c("treatment", "city"), names_to="attribute", values_to="value") %>%
  ggplot(., aes(x=treatment, y=value, fill=treatment))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  coord_flip()+
  facet_wrap(attribute~city, scales="free", ncol=4)+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Attribute")+
  xlab("")
