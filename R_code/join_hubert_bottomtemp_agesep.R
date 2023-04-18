# Clear workspace
rm(list=ls())

library(raster)
library(tidyverse)
library(here)

# Load surveys
stations <- read.csv(here('Data/VAST_input/cod_agesep_VASTdata.csv'))
stations <- subset(stations, AGEGROUP == 'Age0-2')
stations$DATE <- as.Date(stations$DATE,
                            format='%m/%d/%Y')
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
for(i in years){
  brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                             'hubert', 'extended_grd'),
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


# Copy 2020 to 2021
brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                           'hubert', 'extended_grd'),
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
#head(stations.yr)

#stations.yr <- dplyr::select(stations.yr, HAUL_ID, h_bt)

stations$h_bt[stations$HAUL_ID %in% stations.yr$HAUL_ID] <- 
  stations.yr$h_bt

rm(stations.yr, dates, q, brick)
head(stations)
ns <- stations[is.na(stations$h_bt),]
table(ns$YEAR)

# Rebind
fullst <- read.csv(here('Data/VAST_input/cod_agesep_VASTdata.csv'))
hbt <- dplyr::select(stations, HAUL_ID, h_bt)

test <- merge(fullst, hbt, by=c('HAUL_ID'))
nas <- test[is.na(test$h_bt),]
table(nas$YEAR)

# Save
write.csv(test,
          row.names = F,
          here('Data/VAST_input/cod_agesep_VASTdata.csv'))
