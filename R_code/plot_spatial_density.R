### A quick demonstration of how to extract map quantities and
### plot them externally. Cole Monnahan | May 2021
rm(list=ls())

# Load libraries
library(VAST)                           # 3.8.0
library(ggplot2)                        # 2.10.0
library(dplyr)
library(tidyr)
library(here)
library(sf)
library(ggpubr)
library(ggpattern)
library(ggnewscale)

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.title = element_text(size=10),
                legend.text = element_text(size=10),
                legend.background = element_blank(),
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=16, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load VAST fit data
load(here("VAST_runs/tuna10_both/tuna10_both_allfirst.Rdata"))
rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs")))

# Remake map list locally for recreating plots
mdl <- make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

# Load outline of mapped area
area.outline <- st_read(here('Data/GIS/Combined_Bluefin2.shp'))
area.outline <- st_transform(area.outline, "EPSG:4326")

## Get the model estimate of density for each category and year;
# link it spatially to a lat/lon extrapolation point.
D_gt <- fit$Report$D_gct[,,] # drop the category
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=year.labs)
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')

# Adjust data to log abundance, strip units
D$D <- strip_units(D$D)
D$logD <- log(D$D)

# Set CRS 
projargs <- fit$extrapolation_list$projargs
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS(projargs)

# Load spatial information
coast <- ecodata::coast
coast <- st_transform(coast, "EPSG:4326")

us <- st_read(here('Data/GIS/NWAtlantic.shp'))
us <- st_transform(us, st_crs(coast))
us$STOCK <- 'US'
us <- dplyr::select(us, -FID)

can <- st_read(here('Data/GIS/CanadaEEZ.shp'))
can <- st_transform(can, st_crs(coast))
can$STOCK <- 'Canada'
can <- dplyr::select(can, -OBJECTID, -Id, -Shape_Leng, -Shape_Area)

stocks <- rbind(us, can)
stocks <- st_transform(stocks, st_crs(coast))
stocks <- st_make_valid(stocks)

new_bb <- st_bbox(stocks)

# Outer loop: Years
for(i in 1:length(year.labs)){
  Cat.sub <- D[D$Year == paste0(year.labs[i]),]
  
  # Set center of cells
  loc_g <- cbind(Cat.sub$Lon, 
                 Cat.sub$Lat)
    
  n_cells <- dim(loc_g)[1]
    
  Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                             data = Cat.sub, 
                                             proj4string = CRS_orig)
    
  Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
    
  Points_proj = sp::spTransform(Points_orig, CRS_proj)
    
  cell.size = mean(diff(Points_proj@bbox[1, ]), 
                   diff(Points_proj@bbox[2,]))/floor(sqrt(n_cells))
  
  Points_sf = sf::st_as_sf(Points_proj)
    
  grid_fall = sf::st_make_grid(Points_sf, cellsize = cell.size)
  grid_i_fall = sf::st_intersects(Points_sf, grid_fall)
  grid_fall = sf::st_sf(grid_fall, Abundance = tapply(Points_sf$logD, 
                                    INDEX = factor(as.numeric(grid_i_fall),
                                                   levels = 1:length(grid_fall)), 
                                    FUN = mean, na.rm = TRUE))
  grid_fall <- st_transform(grid_fall, "EPSG:4326")

  # Plot fall
  fall <- ggplot()+
      geom_sf(data=grid_fall, aes(fill=Abundance, col=Abundance)) +

      geom_sf(data=area.outline, fill=NA, col='black', pch=19, cex=0.5)+
      geom_sf(data=coast, fill='gray')+
      
      coord_sf(xlim=c(new_bb[1], new_bb[3]),
               ylim=c(new_bb[2], new_bb[4]),
               crs="EPSG:4326")+
      
      theme(legend.position = c(0.85,0.22),
            legend.background = element_rect(fill='white', linewidth = 0)) +
    
      guides(fill = guide_legend(ncol = 2),
             color= guide_legend(ncol=2)) +
    
      theme(legend.key.size = unit(0.2, 'in'))
      
  fall
    
    # Save
    ggsave(fall, 
           filename = 
             paste0(here(),
                    '/Plot_output/',
                    names(D.list)[i], "_total_spatialdensity.png"),
           device="png",
           width = 11, height = 8.5, units='in'
           )
  
}
