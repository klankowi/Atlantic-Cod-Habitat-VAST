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
load(here("VAST_runs/herring_graham/herring1.Rdata"))
rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs")))

# Remake map list locally for recreating plots
mdl <- make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

# Load outline of mapped area
area.outline <- st_read(here('Data/GIS/Prey_Management/Herring_Management_Areas.shp'),
                        quiet=T)
area.outline <- st_transform(area.outline, "EPSG:4326")

## Get the model estimate of density for each category and year;
# link it spatially to a lat/lon extrapolation point.
D_gt <- as.data.frame(strip_units(fit$Report$D_gct[,1,])) # drop the category
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=year.labs)
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')

# Adjust data to log abundance, strip units
D$D <- log(strip_units(D$D))

# Force into seasons
D <- D %>% 
  mutate(OldLab = Year) %>% 
  separate(Year, into=c('Year', 'Season'))

D$Year <- as.numeric(D$Year)
D$Season[D$Season == 'Jan'] <- 'Winter-Spring Feeding'
D$Season[D$Season == 'Jul'] <- 'Summer Feeding-Spawning'
D$Season <- factor(D$Season,
                   levels=c(
                     'Winter-Spring Feeding',
                     'Summer Feeding-Spawning'
                   ))

# Set CRS 
projargs <- fit$extrapolation_list$projargs
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS(projargs)

# Load spatial information
coast <- ecodata::coast
coast <- st_transform(coast, "EPSG:4326")

stocks <- st_make_valid(area.outline)

new_bb <- st_bbox(stocks)

Dlims <- data.frame(
  Year=NA,
  Season=NA,
  Min=NA,
  Max=NA
)

# Scale colors
for(i in 1:length(year.labs)){
  Cat.sub <- D[D$OldLab == paste0(year.labs[i]),]
  Inlims <- data.frame(
    Year=NA,
    Season=NA,
    Min=NA,
    Max=NA
  )
  Inlims$Year <- Cat.sub$Year[1]
  Inlims$Season <- Cat.sub$Season[1]
  Inlims$Min <- min(Cat.sub$D, na.rm=T)  
  Inlims$Max <- max(Cat.sub$D, na.rm=T)
  
  Dlims <- rbind(Dlims, Inlims)
  
  rm(Cat.sub, Inlims)
}
Dlims <- Dlims[!is.na(Dlims$Year),]
summary(Dlims)

summary(Dlims[Dlims$Season == 'Winter-Spring Feeding',])
summary(Dlims[Dlims$Season == 'Summer Feeding-Spawning',])

# Outer loop: Years
for(i in 1:length(year.labs)){
  Cat.sub <- D[D$OldLab == paste0(year.labs[i]),]
  
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
  grid_fall = sf::st_sf(grid_fall, Abundance = tapply(Points_sf$D, 
                                    INDEX = factor(as.numeric(grid_i_fall),
                                                   levels = 1:length(grid_fall)), 
                                    FUN = mean, na.rm = TRUE))
  grid_fall <- st_transform(grid_fall, "EPSG:4326")

  # Plot fall
  fall <- ggplot()+
      geom_sf(data=grid_fall, aes(fill=Abundance, col=Abundance)) +
    
    scale_color_viridis_c(limits=c(min(Dlims$Min), max(Dlims$Max)),
                          breaks=c(min(Dlims$Min), max(Dlims$Max)),
                          labels=c('Low',
                                   'High'),
                          
                          na.value = 'transparent',
                          option=('viridis'),
                          direction = 1) +
    
    scale_fill_viridis_c(limits=c(min(Dlims$Min), max(Dlims$Max)),
                         breaks=c(min(Dlims$Min), max(Dlims$Max)),
                         labels=c('Low', 
                                  'High'),
                         
                         na.value = 'transparent',
                         option=('viridis'),
                         direction = 1) +

      geom_sf(data=area.outline, fill=NA, col='black', pch=19, cex=0.5)+
      geom_sf(data=coast, fill='gray')+
      
      coord_sf(xlim=c(new_bb[1], new_bb[3]),
               ylim=c(new_bb[2], new_bb[4]),
               crs="EPSG:4326")+
    
      labs(color='Log(Abund)', fill='Log(Abund)') +
      
      theme(legend.position = 'right',
            legend.background = element_rect(fill='white', linewidth = 0)) +
    
      theme(legend.key.size = unit(0.2, 'in')) +
      
      ggtitle(paste0(Cat.sub$Season[1], ' ', Cat.sub$Year[1]))
      
  fall
    
    # Save
    ggsave(fall, 
           filename = 
             paste0(here(),
                    '/VAST_runs/herring_graham/annual_spatdens/',
                    Cat.sub$Season[1], '_',
                    Cat.sub$Year[1],
                    ".png"),
           device="png",
           width = 11, height = 8.5, units='in'
           )
  
}
