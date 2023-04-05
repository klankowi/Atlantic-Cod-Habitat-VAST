# Adapted from code by Sarah Gaichas
# Joins existing field station SST data and OISST satellite SST data, tests
# relationship between the two, visualizes OISST inputs at stations missing
# empirical data.

# Clear workspace
rm(list=ls())

# Load libraries
library(tidyverse)
library(here)
library(DT)
library(FishStatsUtils)
library(sf)
library(raster)
library(terra)
library(nngeo)
library(data.table)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))


#list of SST dataframes
SSTdfs <- list.files(here("data-raw/gridded/sst_data/"), pattern = "*.rds")

# Create empty tibble to fill
stn_OISST <- tibble()

# Load station data
load(here("data/RData_Storage/surveys_habitat.RData"))
stations <- sfheaders::sf_to_df(survs_sf, fill=TRUE)
stations <- dplyr::select(stations, -c(sfg_id, point_id))
head(stations)  
stations <- rename(stations, LON=x)
stations <- rename(stations, LAT=y)
rm(survs_sf)

# Create datestring to merge SST and OISST
stations$DATE <- as.POSIXct(stations$DATE, 
                            format = "%Y-%m-%d %H:%M:%S")
stations$month <- substr(stations$DATE, start=6, stop=7)
stations$day <- substr(stations$DATE, start=9, stop=10)
stations$yrmody <- paste0(stations$YEAR, stations$month, stations$day)
stations <- st_as_sf(stations, coords=c("LON", "LAT"),
                     na.fail=T)
stations <- subset(stations, YEAR > 1981)

# Initialize progress bar
pb <- txtProgressBar(min=0, max=length(SSTdfs), initial=0, char="=", style=3)

# Merge in year-based loops
for(df in SSTdfs){
  # Call the number of the current .rds file in the vector of .rds files
  t <- match(df, SSTdfs)
  
  # Update progress bar
  setTxtProgressBar(pb, t)
  getTxtProgressBar(pb)
  
  # Call annual .rds OISST file
  sstdf <- readRDS(paste0(here("data-raw/gridded/sst_data/", df)))
  
  # Create string of sampled years (removes cod data prior to 1982)
  stationsyr <- stations %>%
    filter(YEAR == unique(sstdf$year))
  
  # Filter out OISST days not sampled for cod, convert to sf object
  sstdf_survdays <- sstdf %>%
    dplyr::mutate(yrmody = as.numeric(paste0(year, month, day)) )%>%
    dplyr::filter(yrmody %in% unique(stationsyr$yrmody)) %>%
    dplyr::mutate(year = as.numeric(year),
                  month = as.numeric(month),
                  day = as.numeric(day),
                  declon = Lon,
                  declat = Lat) %>%
    dplyr::select(-Lon, -Lat) %>%
    sf::st_as_sf(coords=c("declon","declat"), crs=4326, remove=FALSE)  
  
  # Join by nearest neighbor and date
  yrdietOISST <- suppressMessages(do.call('rbind', 
                         lapply(split(stationsyr, 
                                      1:nrow(stationsyr)), 
                  function(x) {
                        st_join(x, 
                                sstdf_survdays[sstdf_survdays$yrmody == unique(x$yrmody),],
                                join = st_nn, k = 1, progress = TRUE)
                               }
  )))
  
  # Bind one year of the loop to the initialized tibble
  stn_OISST <- rbind(stn_OISST, yrdietOISST)
  
  # Close segment of progress bar every loop
  close(pb)
  
}

# Save output
saveRDS(stn_OISST, here("data/RData_Storage/stn_OISST.rds"))

# Read in station data and station-OISST
stn_OISST_merge <- stn_OISST %>%
  dplyr::select(-month.y,
                -day.y,
                -yrmody.y) %>% 
  dplyr::rename(month = month.x,
                day = day.x,
                yrmody = yrmody.x,
                oisst = sst) %>%
  dplyr::select(HAUL_ID, oisst, declon, declat) %>%
  sf::st_drop_geometry()

# Merge
agg_stn_all_OISST <- left_join(stations, stn_OISST_merge)

# Save output
saveRDS(agg_stn_all_OISST, here("data/RData_Storage/agg_stn_all_OISST.rds"))
#agg_stn_all_OISST <- readRDS(here("data/agg_stn_all_OISST.rds"))

# Create dataset of comparisons
comparesst <- agg_stn_all_OISST %>%
  dplyr::filter(YEAR>1981)%>%
  dplyr::select(SURFACE.TEMP, oisst, declon, declat, SURVEY) %>%
  na.omit()

# Plot
ggplot2::ggplot(comparesst, aes(x=SURFACE.TEMP, y=oisst, color=SURVEY, fill=SURVEY)) +
  geom_point(alpha=0.8, pch=16)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Field surface temperature (deg C)",
       y = "OISST (deg C)") +
  theme_bw() 

# Map annual shifts
mapsst <- agg_stn_all_OISST %>%
  dplyr::filter(YEAR>1981) %>%
  #dplyr::filter(SEASON == 'FALL' | SEASON=='SPRING') %>% 
  dplyr::mutate(sstdiff = SURFACE.TEMP-oisst) %>%
  dplyr::select(HAUL_ID, YEAR, SEASON, TRUE_SEASON, 
                declon, declat, SURFACE.TEMP, oisst, sstdiff) 

st_crs(mapsst) <- "EPSG:4326"

yrmap <- function(mapyr){
  ggplot2::ggplot(mapsst%>%filter(YEAR==mapyr)) +
    geom_sf(data = ecodata::coast) +
    coord_sf(xlim = c(-77, -65), ylim = c(35, 45)) + 
    geom_point(aes(x=declon, y=declat, colour=sstdiff)) +
    scale_color_gradient2(low = "blue",
                          mid = "yellow",
                          high = "red",
                          midpoint = 0,
                          limits=c(-25, 15),
                          na.value = "gray") +
    theme_bw() +
    facet_wrap(~TRUE_SEASON) +
    ggtitle(paste("SST difference survey-OISST:", mapyr, sep = " "))
}

for(mapyr in 2019:2022){
  
  #cat("  \n####",  as.character(mapyr),"  \n")
  print(yrmap(mapyr)) 
  #cat("  \n")   
  
}
