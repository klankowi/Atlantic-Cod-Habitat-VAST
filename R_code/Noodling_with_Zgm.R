# Create Z_gm matrix dataset to calculate center of gravity shifts specific to
# each stock area.

# Would typically run the model for the whole spatial domain first.

# Prepare workspace
rm(list=ls())

# Load libraries
library(VAST)
library(tidyverse)
library(sf)
library(here)
library(ggpubr)
library(ggnewscale)
library(patchwork)
library(splines)  # Used to include basis-splines
library(INLAspacetime) # used for mesh-building

# Load functions
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
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

# Load VAST fit data
load(here('VAST_runs/add_climate_aja3/add_climate_aja3.Rdata'))
rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs", "cat.labs")))
fitold <- fit
rm(fit)

strata_names <- c('EGOM', 'GBK', 'SNE', 'WGOM')

source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape<- st_read(here::here("", "Data/GIS/cod_region_utm.shp"))
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
region_shape$OBJECTID <- 'ALL'
colnames(region_shape) <- c("Region", 'geometry')

# Second, get our index area shapefile
# We could just use this same shapefile in the "index_shapes" argument, but to 
# show off the new functions we wrote, we will also want to have a sf multipolygon
# shapefiles with areas defined within this general region
index_areas<- c("EGOM", "GBK", "SNE", "WGOM")

for(i in seq_along(index_areas)){
  index_area_temp<- st_read(paste0(here::here("", "Data/GIS"), '/', 
                                   index_areas[i], "_UTM.shp"))
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

index_area_shapes <- bind_rows(index_area_shapes, region_shape)

index_area_shapes <- st_make_valid(index_area_shapes)

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("EGOM", "GBK", "SNE", "WGOM", "ALL"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         #strata.limits = strata_use, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c("fitold", "%notin%", "year.labs", 'vast_extrap_grid',
                        'strata_names', "cat.labs")))

# Load 2000 cell, strata-specific Z_gm data
cell.strata <- read.csv(here('Data/VAST_input/Z_gm_2000cells.csv'))

for(i in 1:length(strata_names)){
  z_gm_use <- subset(cell.strata, STRATA == strata_names[i])
  z_gm_use <- as.matrix(cbind(z_gm_use$E_km, z_gm_use$N_km))
  colnames(z_gm_use) <- colnames(fitold$data_list$Z_gm)
  
  working_dir <- here::here(paste0("VAST_runs/add_climate_aja3/", 
                                   '/',
                                   strata_names[i]))
  # Create working directory if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
  }
  
  setwd(working_dir)
  
  fit = fit_model( 
    
    # Call settings
    settings = fitold$settings, 
    
    # Call survey data info
    Lat_i = fitold$data_frame$Lat_i,
    Lon_i = fitold$data_frame$Lon_i,
    t_i = fitold$data_frame$t_i,
    b_i = fitold$data_frame$b_i,
    a_i = fitold$data_frame$a_i,
    v_i = fitold$data_frame$v_i,
    c_iz = fitold$data_frame$c_iz,
    
    # Call covariate info
    X1_formula =  fitold$X1_formula,
    X2_formula =  fitold$X2_formula,
    covariate_data = fitold$covariate_data,
    
    # Call spatial 
    input_grid=vast_extrap_grid,
    "Extrapolation_List" = fitold$extrapolation_list,
    Z_gm = z_gm_use,
    
    # Set naming conventions
    category_names = cat.labs,
    year_labels = year.labs,
    
    # Tell model to run
    build_model = TRUE,
    run_model = FALSE
    
  )
  
  
  beep(8)
    
  )
  
  
  
  
}
