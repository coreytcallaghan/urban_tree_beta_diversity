# initial script to look at and get a feel for the data

# packages
library(dplyr)
library(ggplo2)
library(tidyr)
library(readr)

# read data in
dat <- read_csv("Data/UMBC_Swan_20211029_w_City.csv")


# total number of species
length(unique(dat$scientific_name))


# number of landuses
unique(dat$fia_lu_abbr)

# cities
table(dat$evalid)
