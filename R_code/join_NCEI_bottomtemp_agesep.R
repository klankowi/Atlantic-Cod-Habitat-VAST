# Clear workspace
rm(list=ls())

# Load libraries
library(tidyverse)
library(sf)
library(raster)
library(here)
library(sp)
library(stars)

# Load survey location data
survs <- read.csv(here('Data/VAST_input/cod_agesep_VASTdata.csv'))

# Add month
survs$MONTH <- month(as.POSIXct(survs$DATE,
                                format="%m/%d/%Y"),
                     abbr=TRUE, label=TRUE)
table(survs$MONTH)

# Add decade
survs$DECADE[survs$YEAR >=1974 & survs$YEAR <=1984] <- 'Y1975-1984'
survs$DECADE[survs$YEAR >=1985 & survs$YEAR <=1994] <- 'Y1985-1994'
survs$DECADE[survs$YEAR >=1995 & survs$YEAR <=2004] <- 'Y1995-2004'
survs$DECADE[survs$YEAR >=2005 & survs$YEAR <=2022] <- 'Y2005-2017'

# Add decade-month, order
survs$DECMON <- paste0(survs$DECADE, "_", survs$MONTH)
head(survs)

# Split survs
survslist <- split(survs, f=survs$DECMON)

# Load raster data
files.list <- list.files(here('Data/Density_Covariates/Bottom_temp/rasters'))
tiflist <- vector(mode='list', length=length(files.list))
for(i in 1:length(files.list)){
  tiflist[[i]] <- read_stars(paste0(here(
    'Data/Density_Covariates/Bottom_temp/rasters'),
                                    '/', files.list[i]))
}
names(tiflist) <- files.list

# Loop, extracting data
for(i in 1:length(survslist)){
  bottem.stars <- st_extract(tiflist[[i]],
                               st_as_sf(survslist[[i]], coords=c('LON', 'LAT'),
                                        crs='EPSG:4326'))
  colnames(bottem.stars) <- c('bottem', 'geometry')
  survslist[[i]]$BOTTEMP <- bottem.stars$bottem
}

survs <- do.call(rbind, survslist)  
rownames(survs) <- NULL  
head(survs)  
summary(survs$BOTTEMP)
