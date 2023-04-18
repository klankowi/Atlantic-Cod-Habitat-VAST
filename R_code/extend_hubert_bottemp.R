rm(list=ls())

library(here)
library(ncdf4)
library(raster)
library(tidyverse)

years <- 1982:2020
codpoly <- st_read(here('Data', 'GIS', 'cod_region_wgs.shp'))
crop_area <- extent(codpoly)

for(i in years){
  message(i)
  brick <- brick(paste0(here('Data', 'Density_Covariates',
                             'Bottom_temp', 'hubert', 'annual_nc'),
                        '/', i, '.grd'))
  brick <- extend(brick, crop_area)
  
  for(j in 1:nlayers(brick)){
    print(j)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
    
    brick[[j]] <- focal(brick[[j]],
                        w=matrix(1,nrow=3, ncol=3), fun=mean, 
                        NAonly=TRUE, na.rm=TRUE)
    #plot(brick[[j]])
    #plot(codpoly, add=T, col=NA)
  }
  
  raster::writeRaster(brick, 
                      filename=paste0(here('Data', 'Density_Covariates',
                                           'Bottom_temp', 'hubert',
                                           'extended_grd'),
                                      '/', i, '.grd'), 
                      bandorder='BIL', 
                      overwrite=TRUE)
  
}


