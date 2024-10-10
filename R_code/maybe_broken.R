rm(list=ls())          

# Load packages
library(here)
library(tidyverse)
library(VAST)
library(sf)

# Install unitless
install_unit(symbol='unitless', def='unitless', name='unitless')

# Load VAST fit data
load(here('VAST_runs/large/overall_bc/All/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))

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
working_dir = here('VAST_runs/large/Overall_BC')
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
regions <- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)
codstox <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
codstox <- st_make_valid(codstox)

for (tI in 1: ncol(Y_gt)) {for (tI in 1: ncol(Y_gt))
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
  grid$TS <- panel_labels[tI]
  
  Big_Data[[tI]] <- grid
  
}

allData <- do.call("rbind", Big_Data)
allData <- allData %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

allData$y <- log(allData$y)

# Plot all years at once
spd <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Fall',], 
          aes(fill=y, col=y)) +
  scale_fill_viridis_c(option='viridis', alpha=0.8, direction=1,
                       na.value = 'transparent') +
  scale_color_viridis_c(option='viridis', alpha=0.8, direction=1,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Fall Spatial density, large Cod')

ggsave(plot=spd, 
       filename='spatial_density_fall.png', 
       height=8.5, width=11, units='in')

spd <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Spring',], 
          aes(fill=y, col=y)) +
  scale_fill_viridis_c(option='viridis', alpha=0.8, direction=1,
                       na.value = 'transparent') +
  scale_color_viridis_c(option='viridis', alpha=0.8, direction=1,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Spring Spatial density, large Cod')

ggsave(plot=spd, 
       filename='spatial_density_spring.png', 
       height=8.5, width=11, units='in')

# Extract high density, compare
allData.spring <- allData[allData$Season == 'Spring',]
allData.fall <- allData[allData$Season == 'Fall',]

allData.fall <- allData.fall %>% 
  group_by(Year) %>% 
  summarise(maxy.fall = max(y, na.rm = TRUE),
            miny.fall = min(y, na.rm=TRUE)) %>% 
  sfheaders::sf_to_df(fill=TRUE) %>% 
  dplyr::select(Year, maxy.fall, miny.fall) %>% 
  unique()

allData.fall <- allData.fall %>% 
  group_by(Year) %>% 
  summarise(maxy.fall = max(y, na.rm = TRUE),
            miny.fall = min(y, na.rm=TRUE)) %>% 
  sfheaders::sf_to_df(fill=TRUE) %>% 
  dplyr::select(Year, maxy.fall, miny.fall) %>% 
  unique()

compData <- merge(allData.fall, allData.fall, by=c("Year"))
compData$hi.ratio <- compData$maxy.fall/ compData$maxy.fall
compData$lo.ratio <- compData$miny.fall/ compData$miny.fall

ggplot(data=compData) +
   geom_line(aes(x=Year, y=hi.ratio)) +
   geom_hline(yintercept = 1, col='red') +
   labs(y='Fall to Fall high density ratio')

ggplot(data=compData) +
  geom_line(aes(x=Year, y=lo.ratio)) +
  geom_hline(yintercept = 1, col='red') +
  labs(y='Fall to Fall low density ratio')

holdData <- Big_Data

for(i in seq(1,79, 2)){
  holdData[[i]]$y.ratio <- 
    (Big_Data[[i]]$y) / 
    (Big_Data[[(i+1)]]$y)
  
  holdData[[i]]$y.ratio <- 
    log(holdData[[i]]$y.ratio)
}

holdData <- holdData[-c(seq(2, 80, 2))]

compData <- do.call("rbind", holdData)
compData <- compData %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Fall', 'Fall')))

spd <- ggplot() +
  geom_sf(data=compData, 
          aes(fill=y.ratio, col=y.ratio)) +
  scale_fill_gradient2(low=("red"), 
                       min="transparent", midpoint=0,
                       high=("blue"),
                       na.value = 'transparent') +
  scale_color_gradient2(low=("red"), 
                        min="transparent", midpoint=0,
                        high=("blue"),
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Spatial density ratio, small Cod')

ggsave(plot=spd, 
       filename='yearly_density_comparison.png', 
       height=8.5, width=11, units='in')

# Annual changes
allData.fall <- allData[allData$Season == 'Fall',]
allData.fall <- allData[allData$Season == 'Fall',]

# Fall
allData.fall <- split(allData.fall, f=allData.fall$Year)

allData.fall[[1]]$y.change <- NA
for(i in 2:length(allData.fall)){
  allData.fall[[i]]$y.change <- 
    allData.fall[[i]]$y / allData.fall[[(i-1)]]$y
}
allData.fall <- do.call(rbind, allData.fall)
allData.fall$y.change <- log(allData.fall$y.change + 0.1)

spd <- ggplot() +
  geom_sf(data=allData.fall, 
          aes(fill=y.change, col=y.change)) +
  scale_fill_gradient2(low=("red"), 
                       mid="transparent", midpoint=0,
                       high=("blue"),
                       na.value = 'transparent') +
  scale_color_gradient2(low=("red"), 
                        mid="transparent", midpoint=0,
                        high=("blue"),
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=codstox, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Yearly change ratio, small Cod')

ggsave(plot=spd, 
       filename='yearly_change_ratio_fall.png', 
       height=8.5, width=11, units='in')

# Fall
allData.fall <- split(allData.fall, f=allData.fall$Year)

allData.fall[[1]]$y.change <- NA
for(i in 2:length(allData.fall)){
  allData.fall[[i]]$y.change <- 
    allData.fall[[i]]$y / allData.fall[[(i-1)]]$y
}
allData.fall <- do.call(rbind, allData.fall)
allData.fall$y.change <- log(allData.fall$y.change + 0.1)

spd <- ggplot() +
  geom_sf(data=allData.fall, 
          aes(fill=y.change, col=y.change)) +
  scale_fill_gradient2(low=("red"), 
                       mid="transparent", midpoint=0,
                       high=("blue"),
                       na.value = 'transparent') +
  scale_color_gradient2(low=("red"), 
                        mid="transparent", midpoint=0,
                        high=("blue"),
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=codstox, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Yearly change ratio, small Cod')

ggsave(plot=spd, 
       filename='yearly_change_ratio_fall.png', 
       height=8.5, width=11, units='in')

# Total spatial density
#Y_gt <- Y_gt[,1:20]
Y_total <- log(strip_units(Y_gt))
Y_total <- as.data.frame(Y_total)
Y_total <- as_tibble(Y_total)
Y_total.fall <- Y_total %>% dplyr::select(contains('Spring'))
Y_total.fall <- rowSums(Y_total.fall)
Y_t.fall <- matrix(Y_total.fall, ncol=1)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y = Y_t.fall[,1]), proj4string = CRS_orig)
Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
Points_proj = sp::spTransform(Points_orig, CRS_proj)
Zlim = range(Y_t.fall, na.rm = TRUE)
xlim = Points_proj@bbox[1, ]
ylim = Points_proj@bbox[2, ]

cell.size = mean(diff(Points_proj@bbox[1, ]), diff(Points_proj@bbox[2, 
]))/floor(sqrt(n_cells))
Points_sf = sf::st_as_sf(Points_proj)
grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
grid_i = sf::st_intersects(Points_sf, grid)
grid = sf::st_sf(grid, y = tapply(Points_sf$y, INDEX = factor(as.numeric(grid_i), 
                                                              levels = 1:length(grid)), FUN = mean, na.rm = TRUE))
grid$Season <- 'Spring'

footprint <- st_read(('C:/Users/klankowicz/Documents/GitHub/Cod-footprint-processing/Jon_July2024/cod_footprint/footv3.shp'), quiet=T)
footprint <- st_make_valid(footprint)
footprint <- st_transform(footprint, st_crs(grid))

tspd.s <- ggplot() +
  geom_sf(data=grid, 
          aes(fill=y, col=y)) +
  scale_fill_viridis_c(option='viridis', alpha=0.8, direction=1,
                       na.value = 'transparent') +
  scale_color_viridis_c(option='viridis', alpha=0.8, direction=1,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +

  geom_sf(data=footprint, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='log(Abund.)', fill='log(Abund.)') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  facet_wrap(vars(Season))

ggsave(plot=tspd.s,
       filename='summed_spatial_density_fall_footprint.png',
       height=4, width=3.5, units = 'in')

Y_total.fall <- Y_total %>% dplyr::select(contains('Fall'))
Y_total.fall <- rowSums(Y_total.fall)
Y_t.fall <- matrix(Y_total.fall, ncol=1)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y = Y_t.fall[,1]), proj4string = CRS_orig)
Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
Points_proj = sp::spTransform(Points_orig, CRS_proj)
Zlim = range(Y_t.fall, na.rm = TRUE)
xlim = Points_proj@bbox[1, ]
ylim = Points_proj@bbox[2, ]

cell.size = mean(diff(Points_proj@bbox[1, ]), diff(Points_proj@bbox[2, 
]))/floor(sqrt(n_cells))
Points_sf = sf::st_as_sf(Points_proj)
grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
grid_i = sf::st_intersects(Points_sf, grid)
grid = sf::st_sf(grid, y = tapply(Points_sf$y, INDEX = factor(as.numeric(grid_i), 
                                                              levels = 1:length(grid)), FUN = mean, na.rm = TRUE))
grid$Season <- 'Fall'
tspd.f <- ggplot() +
  geom_sf(data=grid, 
          aes(fill=y, col=y)) +
  scale_fill_viridis_c(option='viridis', alpha=0.8, direction=1,
                       na.value = 'transparent') +
  scale_color_viridis_c(option='viridis', alpha=0.8, direction=1,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=footprint, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='log(Abund.)', fill='log(Abund.)') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  facet_wrap(vars(Season))

ggsave(plot=tspd.f,
       filename='summed_spatial_density_fall_footprint.png',
       height=4, width=3.5, units = 'in')

# Extract within footprint
footgrid <- st_intersection(allData, footprint)
footcalc <- footgrid %>% 
  group_by(Year, Season, descr) %>% 
  summarise(sum = sum(y))

ggplot(data=footcalc)+#[footcalc$descr %in% c('jeffreys bank', 'northeast georges',
                                           #'stellwagen bank'),]) +
  geom_line(aes(x=Year, y=sum, col=Season)) +
  facet_wrap(vars(descr), scales = 'free')+
  labs(y="Estimate")

ggplot() +
  geom_sf(data=coast, fill='darkgray', col=NA) +
  geom_sf(data=regions[regions$STOCK %in% c('WGOM', 'GBK'),],
          lty=2, fill=NA) +
  geom_sf(data=footprint, aes(fill=descr)) +
  coord_sf(xlim=c(-71, -67),
           ylim=c(41.5, 43.5),
           crs="EPSG:4326") +
  labs(fill='Fine-scale\nfishing area')
