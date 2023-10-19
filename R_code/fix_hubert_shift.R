rm(list=ls())

library(here)
library(ncdf4)
library(raster)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(viridis)
library(ggnewscale)
# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=11),
                axis.text.y=element_text(size=11),
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

goodyear <- 1982:1992
badyear <- 1994:2021
codpoly <- st_read(here('Data', 'GIS', 'cod_region_wgs.shp'))
crop_area <- extent(codpoly)


good.brick <- brick(paste0(here('Data', 'Density_Covariates',
                           'Bottom_temp', 'Old', 'hubert', 'annual_nc'),
                      '/', goodyear[1], '.grd'))
good.brick <- extend(good.brick, crop_area)

if(lubridate::leap_year(goodyear[1])==FALSE){
  good.brick <- dropLayer(good.brick, 366)
}

good.layer <- good.brick[[1]]
good.layer[!is.na(good.layer)] <- 1

bad.brick <- brick(paste0(here('Data', 'Density_Covariates',
                                'Bottom_temp', 'Old', 'hubert', 'annual_nc'),
                           '/', badyear[1], '.grd'))
bad.brick <- extend(bad.brick, crop_area)

if(lubridate::leap_year(badyear[1])==FALSE){
  bad.brick <- dropLayer(bad.brick, 366)
}

bad.layer <- bad.brick[[1]]
bad.layer[!is.na(bad.layer)] <- 1

good.rast <- terra::rast(good.layer)
bad.rast <- terra::rast(bad.layer)

ggplot() +
  geom_spatraster(data=good.rast, alpha=0.7)+
  scale_fill_viridis(na.value=NA, option='inferno', 'good')+
  new_scale_fill()+
  geom_spatraster(data=bad.rast, alpha=0.7)+
  scale_fill_viridis(na.value=NA, option='viridis', 'bad')

bad.matrix <- raster::as.matrix(bad.layer)
good.matrix <- raster::as.matrix(good.layer)

fix.matrix <- bad.matrix
fix.matrix[,1] <- bad.matrix[,ncol(bad.matrix)]
fix.matrix[,2:(ncol(bad.matrix))] <- bad.matrix[,1:ncol(bad.matrix)-1]

fix.layer <- raster::raster(fix.matrix)
extent(fix.layer) <- extent(good.layer)
crs(fix.layer) <- crs(good.layer)

fix.matrix2 <- fix.matrix
#fix.matrix2[8,115:123] <- fix.matrix[38,115:123]
fix.matrix2[9:46,115:123] <- fix.matrix[8:45,115:123]
fix.matrix2[8,115] <- NA

fix.layer2 <- raster::raster(fix.matrix2)
extent(fix.layer2) <- extent(good.layer)
crs(fix.layer2) <- crs(good.layer)


fix.rast <- terra::rast(fix.layer)
fix.rast2 <- terra::rast(fix.layer2)
# 8:38, 115:123


ggplot() +
  geom_spatraster(data=good.rast, alpha=0.9)+
  scale_fill_viridis(na.value=NA, option='viridis', 'good')+
  new_scale_fill()+
  #geom_spatraster(data=fix.rast[,110:123], alpha=0.7)+
  #scale_fill_viridis(na.value=NA, option='viridis', 'fix')
  #new_scale_fill()+
  geom_spatraster(data=fix.rast2, alpha=0.4)+
  scale_fill_viridis(na.value=NA, option="inferno", "fix2")
