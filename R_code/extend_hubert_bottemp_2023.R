rm(list=ls())

library(here)
library(ncdf4)
library(raster)
library(tidyverse)
library(sf)

fillRaster <- function(r){
  
  xyV = as.data.frame(m,xy=TRUE)
  colnames(xyV) <- c('x', 'y', 'layer')
  sp::coordinates(xyV)=~x+y
  
  miss = xyV$layer==999
  miss[is.na(miss)] <- FALSE
  
  m = automap::autoKrige(
    layer~1,
    input_data = xyV[xyV$layer!=999 & !is.na(xyV$layer),],
    new_data=xyV[miss,])
  
  rfill = raster(r)
  rfill[] = r[]
  rfill[miss] = m$krige_output$var1.pred
  
  return(rfill)
}

years <- 1982:2020
codpoly <- st_read(here('Data', 'GIS', 'cod_region_wgs.shp'))
crop_area <- extent(codpoly)

options(warn=-1)

for(i in years){
  message(i)
  brick <- brick(paste0(here('Data', 'Density_Covariates',
                             'Bottom_temp', 'Old', 'hubert', 'annual_nc'),
                        '/', i, '.grd'))
  brick <- extend(brick, crop_area)
  
  
  for(j in 1:1){
    #brick[[j]] <- extend(brick[[j]], crop_area, snap='out')
    print(j)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    #plot(crop_area, add=T)
    
    r <- brick[[j]]
    r <- extend(r, extent(codpoly), snap="out")
    
    #plot(r)
    #plot(codpoly, add=T, col=NA)
    #plot(crop_area, add=T)
    
    r@data@values[is.na(r@data@values)] <- 999
    m <- mask(r, codpoly)
    
    rfill <- fillRaster(m)
    rfill@data@values[rfill@data@values==999] <- NA
    
    blank.brick <- rfill
    
    rm(r, m, rfill)
    
  }
  
  for(j in 2:nlayers(brick)){
    #brick[[j]] <- extend(brick[[j]], crop_area, snap='out')
    print(j)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    #plot(crop_area, add=T)
    
    r <- brick[[j]]
    r <- extend(r, extent(codpoly), snap="out")
    
    #plot(r)
    #plot(codpoly, add=T, col=NA)
    #plot(crop_area, add=T)
    
    r@data@values[is.na(r@data@values)] <- 999
    m <- mask(r, codpoly)
    
    rfill <- fillRaster(m)
    rfill@data@values[rfill@data@values==999] <- NA
    
    blank.brick <- addLayer(blank.brick, rfill)
    
  }
  
  raster::writeRaster(blank.brick, 
                      filename=paste0(here('Data', 'Density_Covariates',
                                           'Bottom_temp', 'hubert',
                                           'extended_grd_krig'),
                                      '/', i, '.grd'), 
                      bandorder='BIL', 
                      overwrite=TRUE)
  
}

options(warn=0)
