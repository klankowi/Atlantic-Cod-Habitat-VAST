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
goodyears <- 1982:1992
badyears <- 1994:2020
transyears <- 1993
codpoly <- st_read(here('Data', 'GIS', 'cod_region_wgs.shp'))
# Remove islands
codpoly2 <- smoothr::fill_holes(codpoly, threshold=30000000000)
# Add buffer (make sure we get all the way inshore)
codpoly2 <- st_buffer(codpoly2, 5000)
crop_area <- extent(codpoly2)
codpoly <- codpoly2
rm(codpoly2)

options(warn=-1)

for(i in years){
  message(i)
  
  brick <- brick(paste0(here('Data', 'Density_Covariates',
                             'Bottom_temp', 'Old', 'hubert', 'annual_nc'),
                        '/', i, '.grd'))
  brick <- extend(brick, crop_area)
  
  if(lubridate::leap_year(i)==FALSE){
    brick <- dropLayer(brick, 366)
  }
  
  if(i %in% goodyears){
    for(j in 1:1){
      print(j)
      
      r <- brick[[j]]
      r <- extend(r, extent(codpoly), snap="out")
      
      r@data@values[is.na(r@data@values)] <- 999
      m <- mask(r, codpoly)
      
      rfill <- fillRaster(m)
      rfill@data@values[rfill@data@values==999] <- NA
      
      blank.brick <- rfill
      
      rm(r, m, rfill)
      
    }
    
    for(j in 2:nlayers(brick)){
      print(j)
      
      r <- brick[[j]]
      r <- extend(r, extent(codpoly), snap="out")
      
      r@data@values[is.na(r@data@values)] <- 999
      m <- mask(r, codpoly)
      
      rfill <- fillRaster(m)
      rfill@data@values[rfill@data@values==999] <- NA
      
      blank.brick <- addLayer(blank.brick, rfill)
      
    }
  }
  
  if(i %in% badyears){
    for(j in 1:1){
      print(j)
      
      r <- brick[[j]]
      bad.matrix <- raster::as.matrix(r)
      fix.matrix <- bad.matrix
      fix.matrix[,1] <- bad.matrix[,ncol(bad.matrix)]
      fix.matrix[,2:(ncol(bad.matrix))] <- bad.matrix[,1:ncol(bad.matrix)-1]
      fix.matrix2 <- fix.matrix
      #fix.matrix2[8,115:123] <- fix.matrix[38,115:123]
      fix.matrix2[9:46,115:123] <- fix.matrix[8:45,115:123]
      fix.matrix2[8,115] <- NA
      fix.layer2 <- raster::raster(fix.matrix2)
      extent(fix.layer2) <- extent(r)
      crs(fix.layer2) <- crs(r)
      r <- fix.layer2
      rm(bad.matrix, fix.layer2, fix.matrix, fix.matrix2)
      
      r <- extend(r, extent(codpoly), snap="out")
      
      r@data@values[is.na(r@data@values)] <- 999
      m <- mask(r, codpoly)
      
      rfill <- fillRaster(m)
      rfill@data@values[rfill@data@values==999] <- NA
      
      blank.brick <- rfill
      
      rm(r, m, rfill)
      
    }
    
    for(j in 2:nlayers(brick)){
      print(j)
      
      r <- brick[[j]]
      bad.matrix <- raster::as.matrix(r)
      fix.matrix <- bad.matrix
      fix.matrix[,1] <- bad.matrix[,ncol(bad.matrix)]
      fix.matrix[,2:(ncol(bad.matrix))] <- bad.matrix[,1:ncol(bad.matrix)-1]
      fix.matrix2 <- fix.matrix
      #fix.matrix2[8,115:123] <- fix.matrix[38,115:123]
      fix.matrix2[9:46,115:123] <- fix.matrix[8:45,115:123]
      fix.matrix2[8,115] <- NA
      fix.layer2 <- raster::raster(fix.matrix2)
      extent(fix.layer2) <- extent(r)
      crs(fix.layer2) <- crs(r)
      r <- fix.layer2
      rm(bad.matrix, fix.layer2, fix.matrix, fix.matrix2)
      
      r <- extend(r, extent(codpoly), snap="out")
      
      r@data@values[is.na(r@data@values)] <- 999
      m <- mask(r, codpoly)
      
      rfill <- fillRaster(m)
      rfill@data@values[rfill@data@values==999] <- NA
      
      blank.brick <- addLayer(blank.brick, rfill)
      
    }
  }
    
    if(i %in% transyears){
      for(j in 1:1){
        print(j)
        
        r <- brick[[j]]
        r <- extend(r, extent(codpoly), snap="out")
        
        r@data@values[is.na(r@data@values)] <- 999
        m <- mask(r, codpoly)
        
        rfill <- fillRaster(m)
        rfill@data@values[rfill@data@values==999] <- NA
        
        blank.brick <- rfill
        
        rm(r, m, rfill)
        
      }
      
      for(j in 2:50){
        print(j)
        
        r <- brick[[j]]
        r <- extend(r, extent(codpoly), snap="out")
        
        r@data@values[is.na(r@data@values)] <- 999
        m <- mask(r, codpoly)
        
        rfill <- fillRaster(m)
        rfill@data@values[rfill@data@values==999] <- NA
        
        blank.brick <- addLayer(blank.brick, rfill)
        
        rm(r, m, rfill)
        
      }
      
      for(j in 51:51){
        blank.brick <- addLayer(blank.brick, blank.brick[[50]])
      }
      
      for(j in 52:nlayers(brick)){
        print(j)
        
        r <- brick[[j]]
        bad.matrix <- raster::as.matrix(r)
        fix.matrix <- bad.matrix
        fix.matrix[,1] <- bad.matrix[,ncol(bad.matrix)]
        fix.matrix[,2:(ncol(bad.matrix))] <- bad.matrix[,1:ncol(bad.matrix)-1]
        fix.matrix2 <- fix.matrix
        #fix.matrix2[8,115:123] <- fix.matrix[38,115:123]
        fix.matrix2[9:46,115:123] <- fix.matrix[8:45,115:123]
        fix.matrix2[8,115] <- NA
        fix.layer2 <- raster::raster(fix.matrix2)
        extent(fix.layer2) <- extent(r)
        crs(fix.layer2) <- crs(r)
        r <- fix.layer2
        rm(bad.matrix, fix.layer2, fix.matrix, fix.matrix2)
        
        r <- extend(r, extent(codpoly), snap="out")
        
        r@data@values[is.na(r@data@values)] <- 999
        m <- mask(r, codpoly)
        
        rfill <- fillRaster(m)
        rfill@data@values[rfill@data@values==999] <- NA
        
        blank.brick <- addLayer(blank.brick, rfill)
        
      }
    
    
  }

  raster::writeRaster(blank.brick, 
                      filename=paste0(here('Data', 'Density_Covariates',
                                           'Bottom_temp', 'hubert',
                                           'extended_grd_krig'),
                                      '/', i, '.grd'), 
                      bandorder='BIL', 
                      overwrite=TRUE)
  
  rm(brick, blank.brick, m, r, rfill)
  
}

options(warn=0)
