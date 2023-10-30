## Join all data together
rm(list=ls())

library(here)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(MultiscaleDTM)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Station data (lat-lon, reported env, cod caught per trawl) ####
stations <- read.csv(here('Data/Survey_Data/Survey_Data.csv'))

# Cleaning
# Remove unnecessary columns
stations <- dplyr::select(stations, -STOCK, -AREA, -STRATUM,
                          -SEASON, -BOTTOM.TEMP, -SALINITY,
                          -BOTTOM.TYPE, -TRUE_SEASON)

# Force datestring to posixct
stations$DATE <- as.POSIXct(stations$DATE, format="%m/%d/%Y")
stations[is.na(stations$DATE),]

# Remove stations prior to 1982
stations <- stations[stations$YEAR >=1982,]

# Remove Sentinel survey
stations <- stations[stations$SURVEY != "Sentinel",]

# Remove surveys that took place Jan-Feb 1982 (belongs to 1981 fall season)
stations$yrmo <- paste0(stations$YEAR, "-", 
                        str_pad(lubridate::month(stations$DATE),
                                side="left", width=2, pad="0"))

stations <- stations[stations$yrmo %notin% c('1982-01',
                                             '1982-02',
                                             '2022-03', 
                                             '2022-04', 
                                             '2022-05'),]

# Add categorical season
stations$SEASON <- NA
stations$SEASON[lubridate::month(stations$DATE) %in% c(3:8)] <- 'A.SPRING'
stations$SEASON[lubridate::month(stations$DATE) %in% c(1,2,9,10,11,12)] <- 'B.FALL'

# Keep fall season temporally continuous (Jan and Feb go with previous year)
stations$YEAR.SEASON <- stations$YEAR
for(i in 1:nrow(stations)){
  if(lubridate::month(stations$DATE[i]) %in% c(1,2)){
    stations$YEAR.SEASON[i] <- stations$YEAR.SEASON[i]-1
  }
}

# Add timestep
stations$TIME <- as.numeric(as.factor(paste0(stations$YEAR.SEASON, 
                                             "-",
                                             stations$SEASON)))

# Remove stations without valid spatial information
stations <- stations[!is.na(stations$LAT),]
stations <- stations[!is.na(stations$LON),]

# Force correct depth units (RIDEM reports in feet, DFO in fathoms)
# Everything should be in meters
stations$DEPTH[stations$SURVEY == 'RIDEM Trawl'] <- 
  stations$DEPTH[stations$SURVEY == 'RIDEM Trawl']*0.3048

stations$DEPTH[stations$SURVEY == 'DFO Trawl'] <- 
  stations$DEPTH[stations$SURVEY == 'DFO Trawl']*1.8288

# Remove points on land
stations_sf <- st_as_sf(stations, coords=c('LON', 'LAT'))
st_crs(stations_sf) <- 'EPSG:4326'
codstox <- st_read(here('Data/GIS/codstox.shp'))
codstox <- st_transform(codstox, st_crs(stations_sf))
codstox <- st_make_valid(codstox)
stations_sf <- st_intersection(stations_sf, codstox)
stations <- sfheaders::sf_to_df(stations_sf, fill=TRUE)
stations <- dplyr::select(stations, -sfg_id, -point_id)
stations <- stations %>% 
  rename('LON' = 'x',
         'LAT' = 'y')
  
# Order
stations <- stations[with(stations, order(DATE)),]
rownames(stations) <- NULL

# View
head(stations)

# Remove unnecessary columns
stations <- dplyr::select(stations, -yrmo, -SEASON, -YEAR, 
                          -YEAR.SEASON, -OBJECTID, -STOCK,
                          -Shape_Leng, -Shape_Area)

#### Pull depth from GEBCO bathymetry ####
# Pull GEBCO bathymetry data
Bathy_Raster <- rast(here('Data/Bathymetry/gebco_2023.tif'))

# Join with characteristics of bathy raster
survs_sf <- st_as_sf(stations, coords=c('LON', 'LAT'))
st_crs(survs_sf) <- "EPSG:4326"
bathy <- extract(Bathy_Raster, survs_sf, bind=T)
bathy <- as.data.frame(bathy)
bathy <- dplyr::select(bathy, -DATE)
stations <- left_join(stations, bathy,
                      by=c('INDEX_NAME', 'SURVEY', 'HAUL_ID',
                           'COD_N', 'COD_KG', 'DEPTH', 'SURFACE.TEMP',
                           'TIME'))
stations$gebco_2023 <- stations$gebco_2023 * -1

# Find stations with gebco depth on land
s.0 <- stations[stations$gebco_2023 <= 0,]
table(s.0$SURVEY)
# All inshore stations, could be in pixels that include land
# Replace with reported depth value
stations$BATHY.DEPTH <- stations$gebco_2023
stations$BATHY.DEPTH[stations$BATHY.DEPTH <=0] <- 
  stations$DEPTH[stations$BATHY.DEPTH <=0]

# Rename and remove columns
stations <- dplyr::select(stations, -DEPTH, -gebco_2023)

# Remove intermediates
rm(bathy, Bathy_Raster, codstox, s.0, stations_sf, survs_sf, i)

#### Add rugosity data ####
stations_sf <- st_as_sf(stations, coords=c('LON', 'LAT'))
st_crs(stations_sf) <- 'EPSG:4326'

# Pull rugosity data
load(here('Data/Density_Covariates/Rugosity/rast_rugosity.RData'))
Rugosity <- rast(masked.raster)
rm(masked.raster)
rugos <- extract(Rugosity, stations_sf, bind=T)
rugos <- as.data.frame(rugos)
rugos <- dplyr::select(rugos, -DATE)
stations_sf <- left_join(stations_sf, rugos,
                      by=c('INDEX_NAME', 'SURVEY', 'HAUL_ID',
                           'COD_N', 'COD_KG', 'SURFACE.TEMP',
                           'TIME', 'BATHY.DEPTH'))
stations_sf <- stations_sf %>% 
  rename('rugos' = 'band1')
table(stations_sf$SURVEY[is.na(stations_sf$rugos)])
# Problem: rugosity values not recorded for 796 tows (728 from RIDEM in NB)
# Solution: remove? no evidence that this is important habitat anyway.
rm(rugos, Rugosity, stations)

#### Add sediment data ####
sed <- st_read(here('Data/GIS/Sediment_Krig_1K_Polygons.shp'))

# Convert to sf
sed_sf   <- st_as_sf(sed)

# Fix self-intersection
sed_sf <- st_make_valid(sed_sf)

# Remove variance
sed_sf <- dplyr::select(sed_sf, -cobble_V, -gravel_V, -mud_V,
                        -rock_V, sand_V,)

# Transform
stations_sf <- st_transform(stations_sf, st_crs(sed_sf))

# Join
stations_sf <- st_join(stations_sf, left=TRUE, sed_sf[,1:5])
stations_sf <- st_transform(stations_sf, 'EPSG:4326')

# Remove intermediates
rm(sed, sed_sf)

#### Add climate index data ####
# Load indices
nao <- read.csv(here('Data/Climate_Indices/daily_NAO.csv'))
amo <- read.csv(here('Data/Climate_Indices/monthly_AMO.csv'))

# Adjust NAO structure
nao$DATE <- as.POSIXct(paste0(nao$year, '-', nao$month, '-', 
                              nao$day),
                       format='%Y-%m-%d')
nao <- nao %>% 
  dplyr::select(-year, -month, -day) %>% 
  rename('nao' = 'nao_index_cdas')

# Merge NAO
stations_sf <- left_join(stations_sf, nao, by=c('DATE'))

# Adjust stations month and year structure
stations_sf$Month <- lubridate::month(stations_sf$DATE,
                                      label=TRUE, abbr = TRUE)
stations_sf$Year <- lubridate::year(stations_sf$DATE)

# Adjust AMO struture
amo <- amo %>% 
  rename('amo' = 'Value')

# Merge AMO
stations_sf <- left_join(stations_sf, amo, by=c('Year', 'Month'))

stations_sf <- dplyr::select(stations_sf, -Year, -Month)

# Remove intermediates
rm(amo, nao)

#### Add area swept ####
# Load area swept df
as <- read.csv(here('Background_Info/area_swept_perhaul.csv'))
as <- dplyr::select(as, HAUL_ID, AREA_SWEPT)

# Merge
stations_sf <- left_join(stations_sf, as, by=c("HAUL_ID"))

# Remove intermediates
rm(as)

#### Add bottom temperature
# Sf to df
stations <- sfheaders::sf_to_df(stations_sf, fill=TRUE)

stations <- stations %>% 
  rename('LON' = 'x',
         'LAT' = 'y') %>% 
  dplyr::select(-point_id, -sfg_id)

# Format date
stations$DATE <- as.Date(stations$DATE,
                         format='%Y-%m-%d')
stations$YEAR <- lubridate::year(stations$DATE)
stations$matching <- paste0('X',
                            stations$YEAR, 
                            '.',
                            str_pad(lubridate::month(stations$DATE), 2, 'left', '0'),
                            '.',
                            str_pad(lubridate::day(stations$DATE), 2, 'left', '0'))
stations$h_bt <- NA

# Set years
years <- 1982:2020

# Loop through years
library(raster)
for(i in years){
  brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                             'hubert', 'extended_grd_krig'),
                        '/', i, '.grd'))
  
  dates <- seq(as.Date(paste0(i, "-01-01")),
               as.Date(paste0(i, "-12-31")),by="1 day")
  
  if(length(dates) < nlayers(brick)){
    dates <- c(dates, NA)
  }
  
  names(brick) <- dates
  
  stations.yr <- subset(stations, stations$YEAR == paste0(i))
  stations.yr <- dplyr::select(stations.yr,
                               HAUL_ID, YEAR, matching, LON, LAT)
  
  
  for(q in 1:nrow(stations.yr)){
    print(round(q/ nrow(stations.yr) * 100),1)
    stations.yr$h_bt[q] <- extract(brick[[stations.yr[q,"matching"]]], 
                                   cbind(stations.yr[q,"LON"], 
                                         stations.yr[q,"LAT"]))
  }
  #head(stations.yr)
  
  #stations.yr <- dplyr::select(stations.yr, HAUL_ID, h_bt)
  
  stations$h_bt[stations$HAUL_ID %in% stations.yr$HAUL_ID] <- 
    stations.yr$h_bt
  
  rm(stations.yr, dates, q, brick)
  
}
head(stations)
ns <- stations[is.na(stations$h_bt),]
table(ns$YEAR)

# Copy 2020 to 2022
brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                           'hubert', 'extended_grd_krig'),
                      '/', 2020, '.grd'))

dates <- seq(as.Date(paste0(2022, "-01-01")),
             as.Date(paste0(2022, "-12-31")),by="1 day")

if(length(dates) < nlayers(brick)){
  dates <- c(dates, NA)
}

names(brick) <- dates

stations.yr <- subset(stations, stations$YEAR == paste0(2022))
stations.yr <- dplyr::select(stations.yr,
                             HAUL_ID, YEAR, matching, LON, LAT)


for(q in 1:nrow(stations.yr)){
  print(round(q/ nrow(stations.yr) * 100),1)
  stations.yr$h_bt[q] <- extract(brick[[stations.yr[q,"matching"]]], 
                                 cbind(stations.yr[q,"LON"], 
                                       stations.yr[q,"LAT"]))
}

stations$h_bt[stations$HAUL_ID %in% stations.yr$HAUL_ID] <- 
  stations.yr$h_bt

rm(stations.yr, dates, q, brick)
head(stations)
ns <- stations[is.na(stations$h_bt),]
table(ns$YEAR)

# Copy 2020 to 2021
brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                           'hubert', 'extended_grd_krig'),
                      '/', 2020, '.grd'))

dates <- seq(as.Date(paste0(2021, "-01-01")),
             as.Date(paste0(2021, "-12-31")),by="1 day")

if(length(dates) < nlayers(brick)){
  dates <- c(dates, NA)
}

names(brick) <- dates

stations.yr <- subset(stations, stations$YEAR == paste0(2021))
stations.yr <- dplyr::select(stations.yr,
                             HAUL_ID, YEAR, matching, LON, LAT)

for(q in 1:nrow(stations.yr)){
  print(round(q/ nrow(stations.yr) * 100),1)
  stations.yr$h_bt[q] <- extract(brick[[stations.yr[q,"matching"]]], 
                                 cbind(stations.yr[q,"LON"], 
                                       stations.yr[q,"LAT"]))
}

stations$h_bt[stations$HAUL_ID %in% stations.yr$HAUL_ID] <- 
  stations.yr$h_bt

rm(stations.yr, dates, q, brick)
head(stations)
ns <- stations[is.na(stations$h_bt),]
table(ns$YEAR)

# Some msising values in 1996 and 2018. Just too close to shore.
# Will have to remove.

# Remove intermediates
rm(ns, stations_sf, i, years)

# Save
write.csv(stations,
          row.names = F,
          here('Data/VAST_input/cod_rawage_VASTdata.csv'))
