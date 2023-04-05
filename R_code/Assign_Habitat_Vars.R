# Code to append temporally-stationary (more or less) habitat variables to
# survey stations. Data include bathymetry, bottom sediment probabilities,
# and rugosity condition (rough/smooth). 
# Katie Lankowicz
# 15 Dec 2022

# Clear workspace
rm(list=ls())

#### Load libraries ####
library(tidyverse, quietly=T,verbose=F)
library(sf, quietly=T,verbose=F)
library(raster, quietly=T, verbose=F)
library(rgeos, quietly=T, verbose=F)
library(here, quietly=T, verbose=F)
library(viridis, quietly=T, verbose=F)
library(RColorBrewer, quietly=T, verbose=F)
library(spatialEco, quietly=T, verbose=F)
library(marmap, quietly=T, verbose=F)
library(beepr, quietly=T, verbose=F)

# Load data
# Load rugosity shapefile
load(here("data/RData_Storage/rugosity_wgs.RData"))

# Load sediment shapefile
sed <- st_read(here("data/GIS/Sediment_Krig_1K_Polygons.shp"),
               quiet=TRUE)

# Convert to sf
sed_sf   <- st_as_sf(sed)

# Fix slight self-intersection issue 
rugos_sf <- st_make_valid(rugos_sf_wgs)
sed_sf <- st_make_valid(sed_sf)

# Transform to unprojected lat-lon
st_crs(rugos_sf); st_crs(sed_sf)
sed_sf <- st_transform(sed_sf, crs = st_crs(rugos_sf))

# Add survey stations
survs <- read.csv(here("data/Dataframes/Bio_Data_Agesep.csv"))

# Remove NA values (cannot be sf object), Convert to sf
survs.n0 <- survs[is.na(survs$LON)==FALSE,]
survs.n0 <- survs.n0[is.na(survs.n0$LAT)==FALSE,]
survs_sf <- st_as_sf(survs.n0, coords = c('LON', 'LAT'))

# Set CRS
st_crs(survs_sf) <- st_crs(rugos_sf)

# Add stock areas
strat_sf <- st_read(here("data/GIS/codstox.shp"),
                 quiet=T)
strat_sf <- st_cast(strat_sf, "POLYGON")
strat_sf <- st_make_valid(strat_sf)
ggplot() +
  geom_sf(data=strat_sf, aes(fill=STOCK))

# Remove points on land
survs_sf <- erase.point(survs_sf, strat_sf, inside=F)
survs_sf <- st_as_sf(survs_sf)

# Pull NOAA bathymetry data
Bathy <- getNOAA.bathy(lon1 = -80, lon2 = -60,
                       lat1 = 30, lat2 = 50, resolution = 1)

# Convert data to raster
Bathy_Raster <- marmap::as.raster(Bathy)

# Remove intermediates
rm(rugos, sed, survs, survs.n0, Bathy)

# Join with characteristics of polygon datasets
survs_sf <- st_join(survs_sf, left=TRUE, sed_sf[,1:10])
beep(1)
names(rugos_sf) <- c('rugos', 'geometry')
survs_sf <- st_join(survs_sf, left=TRUE, rugos_sf[,1])
beep(1)

# Join with characteristics of bathy raster
survs.spdf <- as(survs_sf, "Spatial")
survs.spdf$value <- raster::extract(Bathy_Raster,
                                    survs.spdf)
survs.df <- as.data.frame(survs.spdf)
survs.df <- dplyr::select(survs.df, HAUL_ID, value)
colnames(survs.df) <- c('HAUL_ID', 'BATHY.DEPTH')
survs_sf <- merge(survs_sf, survs.df, by="HAUL_ID")
survs_sf <- survs_sf[with(survs_sf,
                          order(DATE, HAUL_ID)),]
rownames(survs_sf) <- NULL
head(survs_sf)

# Great. Save this.
save(survs_sf, 
     file = here("data/RData_Storage/surveys_habitat_agesep2.RData"))
beep(2)
