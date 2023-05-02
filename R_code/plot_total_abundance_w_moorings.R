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

# Load data
load(here("VAST_runs/refine_effort/refine_effort.Rdata"))

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                legend.background = element_blank(),
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=16, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Pull vector of years
years <- year.labs

# Remake map list locally for recreating plots
mdl <- make_map_info(Region = settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)


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

# Category 4
D_gt.4 <- fit$Report$D_gct[,4,] # drop the category
dimnames(D_gt.4) <- list(cell=1:nrow(D_gt.4), year=years)
D_gt.4 <- D_gt.4 %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D.4 <- merge(D_gt.4, mdl$PlotDF, by.x='cell', by.y='x2i')
#D.4 <- separate(D.4, Year, into = c("Year", "Season"), sep = " (?=[^ ]+$)")
D.4$Cat <- 4

# Rebind to new shape
D <- rbind(D.1, D.2, D.3, D.4)

# Adjust data to log abundance, strip units
D$D <- strip_units(D$D)
D$logD <- log(D$D)

#Rebind to list
D.list <- split(D, f=D$Cat)
names(D.list) <- c("Small", "Medium", "Large", "Unknown size")

# Set CRS 
projargs <- fit$extrapolation_list$projargs
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS(projargs)

moorings <- read.csv(here('Furey_2023_MooringLocations.csv'))
moorings <- st_as_sf(moorings, coords=c('Longitude.Dec', 'Latitude.Dec'))
moorings <- subset(moorings, Region !='NH')
st_crs(moorings) <- 'EPSG:4326'

# Outer loop: Categories
for(i in c(1:3)){
  Cat.sub <- D.list[[i]]
  Year <-  Cat.sub$Year[1]
  
  # Sum for total abundance
  Cat.sub$loc <- paste0(Cat.sub$Lat, Cat.sub$Lon)
  Tot.ab <- Cat.sub[Cat.sub$Year == Year,]
  Sum <- Cat.sub %>% 
    group_by(loc) %>% 
    summarise(Sum=sum(D))
  Sum <- Sum %>% as.data.frame()
  
  # Merge
  tot <- merge(Tot.ab, Sum, by=c('loc'))
  tot <- dplyr::select(tot, Lat, Lon, Sum)
    
    loc_g <- cbind(tot$Lon, 
                   tot$Lat)
    
    n_cells <- dim(loc_g)[1]
    
    Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                             data = data.frame(y = tot$Sum), 
                                             proj4string = CRS_orig)
    
    Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
    
    Points_proj = sp::spTransform(Points_orig, CRS_proj)
    
    cell.size = mean(diff(Points_proj@bbox[1, ]), 
                     diff(Points_proj@bbox[2,]))/floor(sqrt(n_cells))
    Points_sf = sf::st_as_sf(Points_proj)
    grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
    grid_i = sf::st_intersects(Points_sf, grid)
    grid = sf::st_sf(grid, y = tapply(Points_sf$y, 
                                      INDEX = factor(as.numeric(grid_i),
                                                     levels = 1:length(grid)), 
                                      FUN = mean, na.rm = TRUE))
    grid <- st_transform(grid, "EPSG:4326")
    coast <- ecodata::coast
    coast <- st_transform(coast, "EPSG:4326")
    
    box <- c(xmin=-71.5, ymin= 42.5,
             xmax= -69 , ymax=44.5)
    
    grid <- st_crop(grid, box)
    
    max.D <- max(grid$y, na.rm=T)
    max.D <- ceiling(max.D)
    
    # Plot
    p <- ggplot()+
      geom_sf(data=grid, aes(fill=y, col=y)) +
      geom_sf(data=moorings, fill=NA, col='black', pch=19, cex=0.5)+
      scale_color_viridis_c(limits=c(0.01, max.D),
                            breaks=c(100, max.D),
                            labels=c('Low', 'High'),
                            na.value = 'transparent',
                            option=('rocket'),
                            direction = -1,
                            alpha = 0.8) +

      scale_fill_viridis_c(limits=c(0.01, max.D),
                           breaks=c(100, max.D),
                           labels=c('Low',  'High'),
                            na.value = 'transparent',
                            option=('rocket'), 
                           direction = -1, 
                           alpha = 0.8) +
      
      geom_sf(data=coast, fill='gray')+
      coord_sf(xlim=c(-70.5, -69.5),
               ylim=c(43.0,44.0),
               crs="EPSG:4326")+
      #labs(title=paste0(names(D.list)[i], " distribution ", Year.sub$Year[1])) +
      theme(legend.position = c(0.85,0.2),
            legend.background = element_rect(fill='white', linewidth = 0))
    
    p$labels$fill <- "Abund."
    p$labels$colour <- "Abund."
    p <- p + ggtitle(paste0(names(D.list)[i],
                            ' size class\n',
                            'Total abundance 1982 - 2022'))
    p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    p
    
    # Save
    ggsave(p, 
           filename = 
             paste0(here(),
                    '/Plot_output/',
                    names(D.list)[i], "_close_total_spatialdensity_moorings.png"),
           device="png"
           )
  
}

# List files in given directory
files  <- list.files(paste0("C:/Users/klankowicz/Desktop/VAST_examples/",
                            "Mapping3/abund_wind/Small"))
setwd("C:/Users/klankowicz/Desktop/VAST_examples/Mapping3/abund_wind/")
files <- paste0(getwd(), "/Small/", files)

# Set GIF save location
setwd("C:/Users/klankowicz/Desktop/VAST_examples/Mapping3/GIFs/")

# Convert PNGs to GIF
library(gifski)
gifski(files, "Small_abund_wind.gif", delay = 0.25,
       loop=TRUE, width=800, height=600)
