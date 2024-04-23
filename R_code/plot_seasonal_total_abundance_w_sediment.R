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
load(here("VAST_runs/add_climate_aja4/ALL/add_climate_aja4_all.RData"))

# Load sediment grids
load(here('Data/Density_Covariates/Sediment/sediment_grids.RData'))
sedgrid <- grid
rm(c, cent.vals, closed.areas, coast, grid, grid_noNA, 
   hospital_names, in.hard, in.mix, in.soft, surveys, temp, i, substrates,
   wss, hospital_labeller, data.cats, data.f, data.list, pb, seds, test,
   use.cell, use.seds, cat.labs)
hard.agg <- st_union(hard)
hard.agg <- st_transform(hard.agg, "EPSG:4326")
mix.agg <- st_union(mix)
mix.agg <- st_transform(mix.agg, "EPSG:4326")

# Load WEA
wea <- st_read(here('Data/GIS/Final_WEA_Poly.shp'))
wea <- st_transform(wea, st_crs(sedgrid))

# Pull vector of years
years <- fit$year_labels

# Remake map list locally for recreating plots
mdl <- make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)

# Load outline of mapped area
area.outline <- st_read(here('Data/GIS/cod_region_wgs.shp'))
area.outline <- st_transform(area.outline, "EPSG:4326")

# Crop sediment types to area outline
hard.agg <- st_intersection(hard.agg, area.outline)
mix.agg <- st_intersection(mix.agg, area.outline)
test <- c(hard.agg, mix.agg)
test <- st_as_sf(test)
test$Substrate <- c('Hard', 'Mix')
test$Substrate <- as.factor(test$Substrate)

# Call names
category_names <- fit$category_names
size_dist <- c('< 40', '40-70', '> 70')

## Get the model estimate of density for each category and year;
# link it spatially to a lat/lon extrapolation point.

# Category 1
D_gt.1 <- fit$Report$D_gct[,1,] # drop the category
dimnames(D_gt.1) <- list(cell=1:nrow(D_gt.1), year=years)
D_gt.1 <- D_gt.1 %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D.1 <- merge(D_gt.1, mdl$PlotDF, by.x='cell', by.y='x2i')
#D.1 <- separate(D.1, Year, into = c("Year", "Season"), sep = " (?=[^ ]+$)")
D.1$Cat <- 1

# Category 2
D_gt.2 <- fit$Report$D_gct[,2,] # drop the category
dimnames(D_gt.2) <- list(cell=1:nrow(D_gt.2), year=years)
D_gt.2 <- D_gt.2 %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D.2 <- merge(D_gt.2, mdl$PlotDF, by.x='cell', by.y='x2i')
#D.2 <- separate(D.2, Year, into = c("Year", "Season"), sep = " (?=[^ ]+$)")
D.2$Cat <- 2

# Category 3
D_gt.3 <- fit$Report$D_gct[,3,] # drop the category
dimnames(D_gt.3) <- list(cell=1:nrow(D_gt.3), year=years)
D_gt.3 <- D_gt.3 %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D.3 <- merge(D_gt.3, mdl$PlotDF, by.x='cell', by.y='x2i')
#D.3 <- separate(D.3, Year, into = c("Year", "Season"), sep = " (?=[^ ]+$)")
D.3$Cat <- 3

# Rebind to new shape
D <- rbind(D.1, D.2, D.3)

# Adjust data to log abundance, strip units
D$D <- strip_units(D$D)
D$logD <- log(D$D)

#Rebind to list
D.list <- split(D, f=D$Cat)
names(D.list) <- c("Small", "Medium", "Large")

# Set CRS 
projargs <- fit$extrapolation_list$projargs
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS(projargs)

# Outer loop: Categories
for(i in c(1:3)){
  Cat.sub <- D.list[[i]]
  
  # Sum for total abundance
  Cat.sub$loc <- paste0(Cat.sub$Lat, Cat.sub$Lon)
  Cat.sub <- Cat.sub %>% 
    separate(Year, c("Year", "Season"), " ")
  Cat.sub$loc <- paste0(Cat.sub$loc, "-", Cat.sub$Season)
  
  Cat.sub <- Cat.sub %>% 
    filter(Year >=2012)
  
  Sum <- Cat.sub %>% 
    group_by(loc) %>% 
    summarise(Sum=sum(D, na.rm=T))
  Sum <- Sum %>% as.data.frame()
  
  # Merge
  tot <- merge(Cat.sub, Sum, by=c('loc'))
  tot <- dplyr::select(tot, Lat, Lon, Season, Sum)
  tot <- unique(tot)
    
    loc_g <- cbind(tot$Lon, 
                   tot$Lat)
    
    n_cells <- dim(loc_g)[1]
    
    Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                             data = tot, 
                                             proj4string = CRS_orig)
    
    Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
    
    Points_proj = sp::spTransform(Points_orig, CRS_proj)
    
    cell.size = mean(diff(Points_proj@bbox[1, ]), 
                     diff(Points_proj@bbox[2,]))/floor(sqrt(n_cells))
    Points_sf = sf::st_as_sf(Points_proj)
    
    Points_fall = subset(Points_sf, Season == "Fall")
    Points_spring = subset(Points_sf, Season == 'Spring')
    
    grid_fall = sf::st_make_grid(Points_fall, cellsize = cell.size)
    grid_i_fall = sf::st_intersects(Points_fall, grid_fall)
    grid_fall = sf::st_sf(grid_fall, Abundance = tapply(Points_fall$Sum, 
                                      INDEX = factor(as.numeric(grid_i_fall),
                                                     levels = 1:length(grid_fall)), 
                                      FUN = mean, na.rm = TRUE))
    grid_fall <- st_transform(grid_fall, "EPSG:4326")
    coast <- ecodata::coast
    coast <- st_transform(coast, "EPSG:4326")

    max.D_fall <- max(grid_fall$Abundance, na.rm=T)
    max.D_fall <- ceiling(max.D_fall)
    
    grid_spring = sf::st_make_grid(Points_spring, cellsize = cell.size)
    grid_i_spring = sf::st_intersects(Points_spring, grid_spring)
    grid_spring = sf::st_sf(grid_spring, Abundance = tapply(Points_spring$Sum, 
                                                      INDEX = factor(as.numeric(grid_i_spring),
                                                                     levels = 1:length(grid_spring)), 
                                                      FUN = mean, na.rm = TRUE))
    grid_spring <- st_transform(grid_spring, "EPSG:4326")
    
    max.D_spring <- max(grid_spring$Abundance, na.rm=T)
    max.D_spring <- ceiling(max.D_spring)
    
    if(max.D_spring > max.D_fall)  {max = max.D_spring; min = max.D_fall 
                                    maxname= 'Spring'; minname = 'Fall'}
    if(max.D_fall   > max.D_spring){max = max.D_fall; min = max.D_spring
                                    maxname = 'Fall'; minname = 'Spring'}
    
    # Plot fall
    fall <- ggplot()+
      geom_sf(data=grid_fall, aes(fill=Abundance, col=Abundance)) +
      scale_color_viridis_c(limits=c(1, max.D_fall),
                            breaks=c(1, max.D_fall),
                            labels=c('Low', 'High'),
                            na.value = 'transparent',
                            option=('rocket'),
                            direction = -1,
                            alpha = 0.8) +
      scale_fill_viridis_c(limits=c(1, max.D_fall),
                           breaks=c(1, max.D_fall),
                           labels=c('Low',  'High'),
                           na.value = 'transparent',
                           option=('rocket'),
                           direction = -1,
                           alpha = 0.8) +

      new_scale_color() + new_scale_fill() +
      
      geom_sf_pattern(data=test,
                      fill=NA, 
                      lwd=0.05,
                      pattern_density=0.01,
                      pattern_spacing=0.01,
                      pattern_size=0.09,
                      pattern_alpha=0.5,
                      alpha = 0.1,
                      aes(color=Substrate,
                          pattern_color=Substrate,
                          pattern_angle=Substrate)) +
      
      scale_color_manual(values=c('darkblue', 'darkgreen')) +
      scale_pattern_color_manual(values=c('darkblue', 'darkgreen')) +
      scale_pattern_angle_manual(values=c(30, 330)) +

      geom_sf(data=area.outline, fill=NA, col='black', pch=19, cex=0.5)+
      geom_sf(data=coast, fill='gray')+
      
      geom_sf(data=wea, fill=NA, col='red') +
      
      coord_sf(xlim=c(-76, -66),
               ylim=c(36.5,45),
               crs="EPSG:4326")+
      #labs(title=paste0(names(D.list)[i], " distribution ", Year.sub$Year[1])) +
      theme(legend.position = c(0.85,0.22),
            legend.background = element_rect(fill='white', linewidth = 0)) +
      guides(fill = guide_legend(ncol = 2),
             color= guide_legend(ncol=2)) +
      theme(legend.key.size = unit(0.2, 'in'))

    fall <- fall + ggtitle('Fall')
    fall <- fall + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    fall
    
    # Plot spring
    spring <- ggplot()+
      geom_sf(data=grid_spring, aes(fill=Abundance, col=Abundance)) +
      scale_color_viridis_c(limits=c(1, max.D_spring),
                            breaks=c(1, max.D_spring),
                            labels=c('Low', 'High'),
                            na.value = 'transparent',
                            option=('rocket'),
                            direction = -1,
                            alpha = 0.8) +
      scale_fill_viridis_c(limits=c(1, max.D_spring),
                           breaks=c(1, max.D_spring),
                           labels=c('Low',  'High'),
                           na.value = 'transparent',
                           option=('rocket'),
                           direction = -1,
                           alpha = 0.8) +
      
      new_scale_color() + new_scale_fill() +
      
      geom_sf_pattern(data=test,
                      fill=NA, 
                      lwd=0.05,
                      pattern_density=0.01,
                      pattern_spacing=0.01,
                      pattern_size=0.09,
                      pattern_alpha=0.5,
                      alpha = 0.1,
                      aes(color=Substrate,
                          pattern_color=Substrate,
                          pattern_angle=Substrate)) +
      
      scale_color_manual(values=c('darkblue', 'darkgreen')) +
      scale_pattern_color_manual(values=c('darkblue', 'darkgreen')) +
      scale_pattern_angle_manual(values=c(30, 330)) +
      
      geom_sf(data=area.outline, fill=NA, col='black', pch=19, cex=0.5)+
      geom_sf(data=coast, fill='gray')+
      
      geom_sf(data=wea, fill=NA, col='red') +
      
      coord_sf(xlim=c(-76, -66),
               ylim=c(36.5,45),
               crs="EPSG:4326")+
      #labs(title=paste0(names(D.list)[i], " distribution ", Year.sub$Year[1])) +
      theme(legend.position = c(0.85,0.22),
            legend.background = element_rect(fill='white', linewidth = 0)) +
      guides(fill = guide_legend(ncol = 2),
             color= guide_legend(ncol=2)) +
      theme(legend.key.size = unit(0.2, 'in'))
    
    spring <- spring + ggtitle('Spring')
    spring <- spring + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    spring
    
    # Arrange to plot
    both <- ggarrange(spring, fall, nrow=1) + bgcolor('white')
    both <- annotate_figure(both, top=text_grob(paste0(names(D.list)[i],
                                                        ' size class (',
                                                       size_dist[i],
                                                       'cm)\n',
                                                       'Total abundance 2012 - 2021'),
                                                color='black',
                                                #face='bold',
                                                size=16,
                                                vjust=1.85))
    both <- annotate_figure(both, bottom = text_grob(paste0('Maximum in-cell density ',
                                                     round(max/min),
                                                     ' times higher in ',
                                                     maxname, 
                                                     ' than ',
                                                     minname),
                            color='black', size= 12, vjust = -5))
    
    
    # Save
    ggsave(both, 
           filename = 
             paste0(here(),
                    '/Plot_output_wea/',
                    names(D.list)[i], "_total_spatialdensity_wea_last10.png"),
           device="png",
           width = 11, height = 8.5, units='in'
           )
  
}
