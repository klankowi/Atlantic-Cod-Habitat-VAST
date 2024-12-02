rm(list=ls())

library(VAST)
library(sf)
library(here)
library(tidyverse)

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

# Strata
strat <- st_read(here('Data/GIS/BTS_strata_updated.shp'))

# Whole area
whar <- st_read(here('Data/GIS/cod_region_wgs.shp'))

# Overlap
strat <- st_make_valid(st_transform(strat, st_crs(whar)))
whar <- st_make_valid(whar)
strat <- st_intersection(strat, whar)
rm(whar)
strat <- st_transform(strat, crs="EPSG:26919")

# Load VAST
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))
rm(list=setdiff(ls(), c('strat', 'fit')))

# Annual spatial density per BTS strata
# Extract Data
Y_gt = fit$Report$D_gct[,1,]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
panel_labels = fit$year_labels
file_name = "density"
working_dir = here('VAST_runs/large/Overall_BC/ALL')
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
Big_Data2 <- vector("list", length = length(panel_labels))

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
regions <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

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
  #grid$TS <- panel_labels[tI]
  
  grid2 <- st_intersection(grid, strat)
  grid2 <- grid2 %>% 
    group_by(STRATA_1) %>% 
    summarise(y = mean(y, na.rm=T)) %>% 
    mutate(TS = panel_labels[tI])
  
  #Big_Data[[tI]] <- grid
  Big_Data2[[tI]] <- grid2
  
}

# allData <- do.call("rbind", Big_Data)
# allData <- allData %>% 
#   separate(TS, into=c('Year', 'Season')) %>% 
#   mutate(Year = as.numeric(Year)) %>% 
#   mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

allData2 <- do.call("rbind", Big_Data2)
allData2 <- allData2 %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

#allData$y <- log(allData$y)
allData2$y <- log(allData2$y)

ggplot(data=allData2) +
  geom_line(aes(x=Year, y=y, group=STRATA_1), alpha=0.3) +
  facet_wrap(vars(Season))
rm(list=setdiff(ls(), c('allData2', 'strat', 'coast')))

## Load habitat suitability
hsi <- st_read(here('Habitat_Suitability/Large/large_hsi.shp'))
hsi <- st_intersection(hsi, strat)
hsi <- hsi %>% 
  group_by(Year, Season, STRATA_1) %>% 
  summarise(ga = mean(ga, na.rm=T))

logdens <- allData2
rm(allData2)

logdens <- logdens %>% 
  rename(geometry = grid)
st_geometry(logdens) <- 'geometry'

logdens <- st_centroid(logdens)
logdens <- sfheaders::sf_to_df(logdens, fill=T)

logdens <- logdens %>% 
  dplyr::select(STRATA_1, y, Year, Season)

allData <- left_join(hsi, logdens, by=c('Year', 'Season', 'STRATA_1'))

regions <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
regions <- st_make_valid(st_transform(regions, crs=st_crs(allData)))
regions <- regions %>% 
  dplyr::select(-OBJECTID, -Shape_Leng, -Shape_Area)

allData <- st_intersection(allData, regions)
allData <- allData[!is.na(allData$Season),]

ggplot(data=allData) +
  geom_point(aes(x=ga, y=y, col=Year), cex=0.7, alpha=0.5) +
  #geom_smooth(aes(x=ga, y=y), alpha=0.4) +
  ggh4x::facet_grid2(vars(STOCK), vars(Season)) +
  ylab('Log density') + xlab('Geometric mean HSI')

allData <- allData[with(allData, order(
  STOCK, ga
)),]

ggplot(data=allData) +
  geom_point(aes(x=Year, y=ga, col=y), cex=1, alpha=0.5) +
  #geom_smooth(aes(x=ga, y=y), alpha=0.4) +
  ggh4x::facet_grid2(vars(STOCK), vars(Season)) +
  ylab('Log density') + xlab('Geometric mean HSI')

range01 <- function(x, ...){(x - min(x, ..., na.rm=T)) / (max(x, ..., na.rm=T) - min(x, ..., na.rm=T))}

allData$rangega <- range01(allData$ga)

allData.list <- split(allData, f=allData$Season)

for(i in 1:length(allData.list)){
  allData.list[[i]] <- split(allData.list[[i]], 
                             f=allData.list[[i]]$Year)
}

for(i in 1:length(allData.list)){
  for(j in 1:length(allData.list[[i]])){
    allData.list[[i]][[j]] <- split(allData.list[[i]][[j]],
                                               f=allData.list[[i]][[j]]$STOCK)
  }
}

for(i in 1:length(allData.list)){
  for(j in 1:length(allData.list[[i]])){
    for(k in 1:length(allData.list[[i]][[j]])){
      allData.list[[i]][[j]][[k]]$ts.s.ga <- range01(allData.list[[i]][[j]][[k]]$ga)
    }
    
  }
}

for(i in 1:length(allData.list)){
  for(j in 1:length(allData.list[[i]])){
    allData.list[[i]][[j]] <- do.call(rbind, allData.list[[i]][[j]]) 
  }
}

for(i in 1:length(allData.list)){
  allData.list[[i]] <- do.call(rbind, allData.list[[i]])
}

allData <- do.call(rbind, allData.list)

allData$Season <- factor(allData$Season, levels=c('Spring', 'Fall'))

for(i in 1982:2020){
  map <- 
    ggplot(data=allData[allData$Year == i,]) +
      geom_sf(aes(fill=rangega)) +
      scale_fill_viridis_c(option='viridis',
                           limits=c(0,1)) +
      geom_sf(data=regions, fill=NA, col='gray30', alpha=0.4,
              lwd=0.3) +
      facet_wrap(vars(Season)) +
      labs(fill='Geometric\nmean HSI') +
      ggtitle(paste0(i))
  ggsave(plot=map,
         filename=paste0(here('Habitat_Suitability/large/BTS_Suitability_Range/'),
                         '/', i, '.png'))
}


allC <- st_centroid(allData)
allC <- allC %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(-sfg_id, -point_id, -Index) %>% 
  rename(lat = y..1,
         lon = x)
