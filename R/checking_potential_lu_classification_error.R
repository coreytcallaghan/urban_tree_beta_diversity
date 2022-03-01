# packages
library(dplyr)
library(tidyr)
library(readr)

# read data in
dat <- read_csv("Data/UMBC_Swan_20211029_w_City.csv")

dat2 <- dat %>%
  unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
  left_join(., dat %>%
              unite(site_coords, lat, lon, remove=FALSE, sep="_") %>%
              dplyr::select(site_coords) %>%
              distinct() %>%
              mutate(site_id=1:nrow(.)), by="site_coords")

length(unique(dat2$site_id))

test <- dat2 %>%
  group_by(evalid, site_coords, site_id) %>%
  summarize(number_landuse=length(unique(fia_lu_abbr)))

data_to_check <- dat2 %>%
  dplyr::filter(site_id %in% local(test %>%
                  dplyr::filter(number_landuse==2) %>%
                  .$site_id))
