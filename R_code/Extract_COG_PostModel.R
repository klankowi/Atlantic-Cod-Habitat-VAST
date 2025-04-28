rm(list=ls())          

# Load packages
library(here)
library(tidyverse)
library(VAST)
library(sf)
library(centr)

# Install unitless
#install_unit(symbol='unitless', def='unitless', name='unitless')

# Load VAST fit data
load(here('VAST_runs/large/Overall_BC/ALL_Catchability2/Overall_BC_largecod_allstrat_natsplin_fsON_ALL_Catchability_SVC.RData'))
rm(list=setdiff(ls(), c('fit')))

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Extract Data
Y_gt = fit$Report$D_gct[,1,]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
panel_labels = fit$year_labels
file_name = "density"
working_dir = here('VAST_runs/large/Overall_BC/ALL_Catchability2')
setwd(working_dir)
fun = mean

# Call data
zlim = range(Y_gt, na.rm = TRUE)
MapSizeRatio = map_list$MapSizeRatio
Y_gt = Y_gt[map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                    0), "x2i"], , drop = FALSE]
n_cells = nrow(Y_gt)

# Call location list
loc_g = map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                0), c("Lon", "Lat")]
# Call projections
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS("EPSG:26919")

# Make list to append data
Big_Data <- vector("list", length=length(panel_labels))

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
region_shape<- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
# First, we need our region shapefile
region_shape$OBJECTID <- 'ALL'
colnames(region_shape) <- c("Region", 'geometry')
# Second, get our index area shapefile
# We could just use this same shapefile in the "index_shapes" argument, but to 
# show off the new functions we wrote, we will also want to have a sf multipolygon
# shapefiles with areas defined within this general region
index_areas<- c("WGOM", "EGOM", "GBK", "SNE")
for(i in seq_along(index_areas)){
  index_area_temp<- st_read(paste0(here("Data/GIS"), '/', 
                                   index_areas[i], "_UTM.shp"), quiet=T)
  index_area_temp <- st_transform(index_area_temp, "EPSG:4326")
  index_area_temp <- st_make_valid(index_area_temp)
  index_area_temp <- dplyr::select(index_area_temp, STOCK, geometry)
  colnames(index_area_temp) <- c("Region", 'geometry')
  
  if(i == 1){
    index_area_shapes<- index_area_temp
  } else {
    index_area_shapes<- bind_rows(index_area_shapes, index_area_temp)
  }
}
index_area_shapes <- bind_rows(region_shape, index_area_shapes)
index_area_shapes <- st_make_valid(index_area_shapes)
index_area_shapes <- st_transform(index_area_shapes, 
                                  "EPSG:26919")

for (tI in 1: ncol(Y_gt)) {
  print(tI)
  Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                           data = data.frame(y = Y_gt[, tI]), 
                                           proj4string = CRS_orig)
  Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
  Points_proj = sp::spTransform(Points_orig, CRS_proj)
  Zlim = zlim
  xlim = Points_proj@bbox[1, ]
  ylim = Points_proj@bbox[2, ]
  
  cell.size = mean(diff(Points_proj@bbox[1, ]), diff(Points_proj@bbox[2, 
  ]))/floor(sqrt(n_cells))
  Points_sf = sf::st_as_sf(Points_proj)
  grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
  grid_i = sf::st_intersects(Points_sf, grid)
  grid = sf::st_sf(grid, y = tapply(Points_sf$y, INDEX = factor(as.numeric(grid_i), 
                                                                levels = 1:length(grid)), FUN = mean, na.rm = TRUE))
  
  grid <- st_intersection(grid, index_area_shapes)
  
  grid <- split(grid, f=grid$Region)
  for(q in 1:length(grid)){
    grid[[q]] <- grid[[q]] %>% 
      filter(!is.na(y)) %>% 
      mean_center(group = 'Region', weight='y') %>% 
      sfheaders::sf_to_df(fill=T) %>% 
      mutate(easting = as.numeric(x),
             northing = as.numeric(y)) %>% 
      dplyr::select(-sfg_id, -point_id, -x, -y)
  }
  grid <- do.call(rbind, grid)
  grid$TS <- panel_labels[tI]
  Big_Data[[tI]] <- grid
  
  rm(grid, grid_i, Points_sf, cell.size, ylim, xlim, Zlim, Points_proj,
     Points_LongLat, Points_orig)
}

allData <- do.call("rbind", Big_Data)

allData <- allData %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

write.csv(allData,
          'COG_New_Process.csv')
