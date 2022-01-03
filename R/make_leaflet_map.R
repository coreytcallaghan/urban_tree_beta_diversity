# script to make a leaflet map for initial exploration of where
# the data are


# packages
library(leaflet)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sf)
library(htmlwidgets)

# read data in
dat <- read_csv("Data/UMBC_Swan_20211029_w_City.csv")

dat2 <- dat %>%
  unite(lon_lat, lon, lat, sep="_", remove=FALSE) %>%
  group_by(lon_lat) %>%
  summarize(species_richness=length(unique(scientific_name))) %>%
  mutate(point_id=1:nrow(.)) %>%
  left_join(., dat %>%
              unite(lon_lat, lon, lat, sep="_", remove=FALSE) %>%
              dplyr::select(lon_lat, lon, lat) %>%
              distinct()) %>%
  dplyr::select(-lon_lat) %>%
  dplyr::select(point_id, lon, lat, species_richness)

ggplot(dat2, aes(x=species_richness))+
  geom_histogram(fill="gray80", color="black")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Species richness")+
  ylab("Number of points")

# points sf
points_sf <- dat2 %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326)

plot(points_sf)

map <- leaflet(points_sf) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(clusterOptions = markerClusterOptions())

map

# but make it a bit better now.

# first cut the continuous variable into bins
# these bins are now factors
pal <- colorFactor(
  palette = colorRampPalette(rainbow(10))(length(points_sf$species_richness)), 
  domain = points_sf$species_richness)


map <- leaflet(points_sf) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(clusterOptions = markerClusterOptions(),
                   color=~pal(species_richness)) %>% 
  addLegend("bottomright", pal = pal, values = ~species_richness,
            title = "Species richness")

map


saveWidget(map, file="map_of_points.html")
