# Pull, clean, and compile bottom temperature data from NOAA NCEI temperature
# data, statistical mean on a 1/10deg grid. Do so for 1982-2022.

# Clear workspace
rm(list=ls())

# Load libraries
library(ncdf4)
library(raster)
library(tidyverse)
library(here)
library(sp)
library(rgdal)
library(chron)
library(viridis)
library(intamap)
library(lubridate)
library(dplyr)
library(sf)

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
                plot.caption=element_text(hjust=0, face='italic', size=12),
                strip.text=element_text(size=12)))

# Load coast shapefile
coast <- ecodata::coast
coast <- st_transform(coast, crs='EPSG:4326')

# Set lat-lon bounds
minlat <- 35; minlon <- -78
maxlat <- 46; maxlon <- -65

# Pull index of folder names (decades)
folds <- list.dirs(here('Data/Density_Covariates/Bottom_temp/nc_files'))
folds <- folds[2:5]

# Name decades
decade_names <- c('Y1975-1984', 'Y1985-1994', 'Y1995-2004', 'Y2005-2017')

# Name months
month_names <- month.abb

# Load VAST polygon shapefile
vpoly <- st_read(here('Data/GIS/NWAtlantic.shp'))
vpoly <- st_transform(vpoly, crs='EPSG:4326')

# Create blank to collect data over decades
overband <- data.frame(
  lon=NA,
  lat=NA,
  tempc=NA,
  depth=NA,
  decade=NA,
  month=NA
)

for(i in 1:
    length(folds)
    #1 # for testing
    ){
  # Pull index of file names (months)
  mons <- list.files(folds[i])
  # Print name of decade (to visualize progress)
  message(decade_names[i])
  
  # Create blank to collect data over months
  innerband <- data.frame(
    lon=NA,
    lat=NA,
    tempc=NA,
    depth=NA,
    decade=NA,
    month=NA
  )
  
  for(j in 1:
      length(mons)
      #3 # for testing
      ){
    # Print name of month (to visualize progress)
    print(month_names[j])
    # Open NC file
    nc_data <- nc_open(paste0(folds[i], '/', mons[j]))
    names(nc_data$var)
    
    # Load statistical mean temperature variable
    t_mn <- ncvar_get(nc_data, 't_mn')
    
    # Get lat-lon
    lat <- ncvar_get(nc_data, 'lat')
    lat <- as.vector(lat)
    lon <- ncvar_get(nc_data, 'lon')
    lon <- as.vector(lon)
    
    # Get depths
    depths <- ncvar_get(nc_data, 'depth')
    depths <- as.vector(depths)

    # Close netcdf file
    nc_close(nc_data)
    
    # Check dims of t_mn
    dim(t_mn)
    
    # Create blank
    totbands <- data.frame(
      lon=NA,
      lat=NA,
      tempc=NA,
      depth=NA
    )
    
    # Loop through depths 
    for(k in 1:dim(t_mn)[3]){
      # Pull array slice of depthband i
      depband <- t_mn[,,k]
      
      # Assign row and column names as lat-lon locations
      rownames(depband) <- lon
      colnames(depband) <- lat
      
      # Convert to df
      depband <- as.data.frame(as.table(depband))
      colnames(depband) <- c('lon', 'lat', 'tempc')
      
      # Set depth
      depband$depth <- depths[k]
      
      # Remove NA values (on land, not collected, etc.)
      depband <- depband[complete.cases(depband),]
      
      # Change structure
      depband$lat <- as.numeric(as.character(depband$lat))
      depband$lon <- as.numeric(as.character(depband$lon))
      
      # Remove locations outside of lat-lon interest bounds
      depband <- depband %>% 
        filter(lat >= minlat & lat <= maxlat) %>% 
        filter(lon >= minlon & lon <= maxlon)
      
      # Save
      totbands <- rbind(totbands, depband)
      
    }
    
    # Remove that first row
    totbands <- totbands[-1,]
    
    # For each location, save deepest temperature estimation
    totbands$latlon <- paste0(totbands$lat, totbands$lon)
    maxdep <- totbands %>% 
      group_by(latlon) %>% 
      slice_max(depth)
    
    # Clean, re-order
    maxdep <- as.data.frame(maxdep)
    maxdep <- maxdep[with(maxdep, order(lat, lon, depth)),]
    rownames(maxdep) <- NULL
    
    # Clip to VAST polygon
    maxsf <- st_as_sf(maxdep, coords=c('lon', 'lat'))
    st_crs(maxsf) <- 'EPSG:4326'
    aoi <- st_intersection(maxsf, vpoly)
    
    mostin.save <- sfheaders::sf_to_df(aoi, fill=T)
    mostin.save <- dplyr::select(mostin.save,
                                 x, y, depth, tempc)
    
    mostin.save$decade <- decade_names[i]
    mostin.save$month <- j
    colnames(mostin.save) <- c('lon', 'lat', 'depth', 'tempc', 'decade', 
                               'month')
    
    # Remove intermediates
    rm(aoi, depband, maxdep, maxsf, nc_data, totbands, 
       lat, lon, depths, t_mn)
    
    # Bind
    innerband <- rbind(innerband, mostin.save)
    
    rm(mostin.save)
  }
  
  innerband <- innerband[-1,]
  # Bind
  overband <- rbind(overband, innerband)
  
  rm(innerband)
  
}

overband <- overband[-1,]
overband$decade <- factor(overband$decade)
summary(overband)
head(overband)
table(overband$decade)
table(overband$month)
table(overband$depth)

overband.sf <- st_as_sf(overband, coords=c('lon', 'lat'))
st_crs(overband.sf) <- 'EPSG:4326'
for(i in 1:nrow(overband.sf)){
  overband.sf$moab[i] <- month.abb[overband.sf$month[i]]
}
overband.sf$moab <- factor(overband.sf$moab,
                           levels=c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))


declist <- split(overband.sf, f=overband.sf$decade)

for(i in 1:length(declist)){
  # Plot
  ggplot() +
    geom_sf(data=coast) +
    geom_sf(data=declist[[i]], aes(col=depth)) +
    scale_color_viridis_c(option='mako',
                          direction = -1,
                          na.value=NA) +
    coord_sf(xlim=c(-78, -65),
             ylim=c(35, 46)) +
    facet_wrap(vars(moab)) +
    ggtitle(paste0('Maximum depth estimated by month, ',
                   decade_names[i]))
  
  ggplot() +
    geom_sf(data=coast) +
    geom_sf(data=declist[[i]], aes(col=tempc)) +
    scale_color_viridis_c(option='inferno',
                          direction = -1,
                          na.value=NA,
                          begin=0.25, end=1) +
    coord_sf(xlim=c(-78, -65),
             ylim=c(35, 46)) +
    facet_wrap(vars(moab)) +
    ggtitle(paste0('Bottom temperature at maximum depth by month, ',
                   decade_names[i]))
}
  
for(i in 1:length(declist)){
  p <- ggplot() +
    geom_point(data=declist[[i]],
               aes(x=depth, y=tempc)) +
    geom_smooth(data=declist[[i]],
                aes(x=depth, y=tempc),
                method='loess',
                span=0.9) +
    xlim(0,300) +
    ylim(-2, 32) +
    facet_wrap(vars(moab)) +
    xlab('Depth (m)') +
    ylab('Temperature (C)') +
    ggtitle(paste0('Monthly Depth-Temperature relationships, ',
            decade_names[i]))
  plot(p)
}

# Set number of rows and columns
n.row <- length(table(overband$lat)) 
n.col <- length(table(overband$lon))

# Create matrix of desired dimensions, rasterize
r <- raster(matrix(0, n.row, n.col))

# Set extent to match points retained from SCHISM manipulation
# Make sure extent is exactly 2500m by 500m so cell resolution is 10m
extent(r) <- c(st_bbox(vpoly)[1], st_bbox(vpoly)[3],
               st_bbox(vpoly)[2], st_bbox(vpoly)[4])

# Set projection
# NAD83 HARN StatePlane Maryland FIPS 1900 (m)
r@crs <- CRS("EPSG:4326")

# Cut to vpoly
vpoly <- readOGR(here('Data/GIS/cod_region_rmislands_buf3.shp'))
vpoly <- spTransform(vpoly, CRS("+init=epsg:4326"))
r <- crop(r, vpoly)
r <- mask(r, vpoly)
plot(r)
# Convert raster to SpatialPixels
rpix <- as(r, "SpatialPixels")

declist <- split(overband, f=overband$decade)
for(i in 1:length(declist)){
  declist[[i]] <- split(declist[[i]], f=declist[[i]]$month)
}

# Final cutting polygon
#closepoly <- readOGR(here('Data/GIS/cod_region_rmislands.shp'))
#closepoly <- spTransform(closepoly, CRS("+init=epsg:4326"))

# Set k, typically 5 or 10
k <- 10

# Set color scheme
cols <- viridis(256)

for(i in 2:
    length(declist)
    #1 # testing 
    ){
  message(decade_names[i])
  for(j in 1:
      length(declist[[i]])
      # 3 #testing
      ){
    print(month_names[j])
    # Subset to decade-month combo
    test <- declist[[i]][[j]]
    # Points to spatial points data frame
    spdat <- sp::SpatialPointsDataFrame(
      # Coordinates set to lon and lat (will be same for all timesteps)
      coords=test[,c("lon", "lat")],
      # Data set to variable at timestep i
      data  =test[,c('depth', 'tempc')],
      # CRS set to NAD83 HARN StatePlane Maryland FIPS 1900 (m)
      proj4string = CRS("+init=epsg:4326")
    )
    # Rename dep to value, artifact of intamap package
    names(spdat@data) <- c('idx', 'value')
    
    # Suppress warning associated with CRS comments
    options(warn=-1)
    
    # Set up idwObject
    idwObject <-  createIntamapObject(
      # Set observations to SpatialPointsDataFrame
      observations        = spdat,
      # Set formula to universal kriging
      formulaString       = as.formula(value~1),
      # Set prediction locations to blank raster of model space
      predictionLocations = rpix,
      # Set method to inverse distance weighting
      class               = "idw"
    )
    
    # k-fold brute-force Cross-validation selection of IDP
    # Estimate best IDP for IDW interpolation
    idwObject <-  estimateParameters(
      idwObject,
      # Set range of possible inverse distance power values
      # min=0, max=5, step=0.05. All steps will be evaluated.
      idpRange = seq(0.00,5.00,0.05),
      # Call k-fold value
      nfold=k
    )
    
    # Interpolate using IDP chosen above
    fitmax <- gstat::gstat(
      # Again, universal kriging
      formula=value ~ 1, 
      # Call spatial points dataframe with 
      data=spdat, 
      # Set number of neighbors to consider. Bigger=smoother               
      nmax=4, 
      # Set IDW estimated above
      set=list(idp=idwObject$inverseDistancePower)
    )
    
    # Rasterize interpolation based on previous steps
    maxint <- raster::interpolate(r, model=fitmax, ext=r)
    maxint <- crop(maxint, vpoly)
    maxint <- mask(maxint, vpoly)
    
    # Save minimum and maximum values for plotting (maximizes color stretch)
    mindep <- floor(min(maxint@data@values))
    maxdep <- ceiling(max(maxint@data@values))
    
    # Plot raster
    raster::plot(maxint, asp=1, zlim=c(mindep,maxdep), col=cols,
                 legend=T, ylab='Latitude (m)', 
                 xlab='Longitude (m)',
                 legend.args=list(text='Bottom temp (C)',side=4, 
                                  font=2, line=2.5, cex=0.8),
                 xlim=c(maxint@extent[1], maxint@extent[2]))
    
    # Save raster
    writeRaster(maxint,
                filename=paste0(here('Data/Density_Covariates/Bottom_temp/rasters/'),
                                '/', decade_names[i], '_',
                                month_names[j], '.tif'))
    
    rm(fitmax, idwObject, maxint, spdat, test, maxdep, mindep)
    
    
  }
}
