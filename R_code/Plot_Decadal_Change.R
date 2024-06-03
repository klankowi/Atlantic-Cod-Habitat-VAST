rm(list=ls())          

# Load packages
library(here)
library(tidyverse)
library(VAST)
library(sf)
library(egg)

# Install unitless
#install_unit(symbol='unitless', def='unitless', name='unitless')

# Load VAST fit data
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_larcod_allstrat_natsplin_fsON_ALL.RData'))

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

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
regions <- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

rm(list=setdiff(ls(), c('coast', 'CRS_orig', 'CRS_proj',
                        'loc_g', 'regions', 'fun', 'n_cells',
                        'Y_gt')))

#### Calculate decadal change ####
# Spring
Y_total <- (strip_units(Y_gt))
Y_total <- as.data.frame(Y_total)
Y_total <- as_tibble(Y_total)
Y_total.spring <- Y_total %>% dplyr::select(contains('Spring'))
Y_total.spring.80s <- Y_total.spring %>% 
  dplyr::select(contains(c('1982', '1983', '1984', '1985', '1986',
                           '1987', '1988', '1989', '1990', '1991')))
Y_total.spring.90s <- Y_total.spring %>% 
  dplyr::select(contains(c('1992', '1993', '1994', '1995', '1996',
                           '1997', '1998', '1999', '2000', '2001')))
Y_total.spring.00s <- Y_total.spring %>% 
  dplyr::select(contains(c('2002', '2003', '2004', '2005', '2006',
                           '2007', '2008', '2009', '2010', '2011')))
Y_total.spring.10s <- Y_total.spring %>% 
  dplyr::select(contains(c('2012', '2013', '2014', '2015', '2016',
                           '2017', '2018', '2019', '2020', '2021')))

Y_total.spring.80s <- rowSums(Y_total.spring.80s)
Y_t.spring.80s <- matrix(Y_total.spring.80s, ncol=1)
Y_total.spring.90s <- rowSums(Y_total.spring.90s)
Y_t.spring.90s <- matrix(Y_total.spring.90s, ncol=1)
Y_total.spring.00s <- rowSums(Y_total.spring.00s)
Y_t.spring.00s <- matrix(Y_total.spring.00s, ncol=1)
Y_total.spring.10s <- rowSums(Y_total.spring.10s)
Y_t.spring.10s <- matrix(Y_total.spring.10s, ncol=1)

# Spring
spring <- cbind(Y_t.spring.80s, Y_t.spring.90s, Y_t.spring.00s, 
                Y_t.spring.10s)
colnames(spring) <- c('Eighties', 'Nineties', 'Aughties', 'Tens')
spring.change <- matrix(data=NA, nrow=nrow(spring), ncol=ncol(spring)-1)

for(i in 2:ncol(spring)){
  spring.change[,(i-1)] <- 
    ((spring[,(i)] - spring[,(i-1)]) / spring[,(i-1)]) * 100
}

colnames(spring.change) <- c('% Change 80s - 90s',
                             '% Change 90s - 00s',
                             '% Change 00s - 10s')

spring.firstlast <- spring[,'Tens'] - spring[,'Eighties']

spring.change <- cbind(spring.change, spring.firstlast)
colnames(spring.change)[4] <- '% Change 80s - 10s'

# Spring
grid.spring <- data.frame(
  y =NA,
  lat=1, lon=1,
  Decade = NA,
  Season = NA
)
grid.spring <- st_as_sf(grid.spring, coords=c('lon', 'lat'),
                        crs="EPSG:26919")
st_geometry(grid.spring) <- 'grid'

for(i in 1:ncol(spring.change)){
  Zlim = range(spring.change)
  Points_orig = sp::SpatialPointsDataFrame(
    coords = loc_g, 
    data = data.frame(y = spring.change[,i]), 
    proj4string = CRS_orig)
  
  Points_LongLat = sp::spTransform(Points_orig, 
                                   sp::CRS("+proj=longlat"))
  Points_proj = sp::spTransform(Points_orig, 
                                CRS_proj)
  cell.size = mean(diff(Points_proj@bbox[1, ]), 
                   diff(Points_proj@bbox[2, 
                   ]))/floor(sqrt(n_cells))
  Points_sf = sf::st_as_sf(Points_proj)
  grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
  grid_i = sf::st_intersects(Points_sf, grid)
  grid = sf::st_sf(grid, y = tapply(Points_sf$y, 
                                    INDEX = factor(as.numeric(grid_i), 
                                                   levels = 1:length(grid)),
                                    FUN = mean, na.rm = TRUE))
  grid$Decade <- colnames(spring.change)[i]
  grid$Season <- 'Spring'
  
  grid.spring <- rbind(grid.spring, grid)
  
  rm(grid, grid_i, Points_LongLat, Points_orig, Points_proj,
     Points_sf, Zlim)
}

grid.spring <- grid.spring[!is.na(grid.spring$Season),]
grid.spring$Decade <- factor(grid.spring$Decade,
                             levels=c('% Change 80s - 90s',
                                      '% Change 90s - 00s',
                                      '% Change 00s - 10s',
                                      '% Change 80s - 10s'))

rm(list=setdiff(ls(), c('coast', 'CRS_orig', 'CRS_proj',
                        'loc_g', 'regions', 'fun', 'n_cells',
                        'Y_gt', 
                        'grid.spring')))

# Fall
Y_total <- (strip_units(Y_gt))
Y_total <- as.data.frame(Y_total)
Y_total <- as_tibble(Y_total)
Y_total.fall <- Y_total %>% dplyr::select(contains('Fall'))
Y_total.fall.80s <- Y_total.fall %>% 
  dplyr::select(contains(c('1982', '1983', '1984', '1985', '1986',
                           '1987', '1988', '1989', '1990', '1991')))
Y_total.fall.90s <- Y_total.fall %>% 
  dplyr::select(contains(c('1992', '1993', '1994', '1995', '1996',
                           '1997', '1998', '1999', '2000', '2001')))
Y_total.fall.00s <- Y_total.fall %>% 
  dplyr::select(contains(c('2002', '2003', '2004', '2005', '2006',
                           '2007', '2008', '2009', '2010', '2011')))
Y_total.fall.10s <- Y_total.fall %>% 
  dplyr::select(contains(c('2012', '2013', '2014', '2015', '2016',
                           '2017', '2018', '2019', '2020', '2021')))
Y_total.fall.80s <- rowSums(Y_total.fall.80s)
Y_t.fall.80s <- matrix(Y_total.fall.80s, ncol=1)
Y_total.fall.90s <- rowSums(Y_total.fall.90s)
Y_t.fall.90s <- matrix(Y_total.fall.90s, ncol=1)
Y_total.fall.00s <- rowSums(Y_total.fall.00s)
Y_t.fall.00s <- matrix(Y_total.fall.00s, ncol=1)
Y_total.fall.10s <- rowSums(Y_total.fall.10s)
Y_t.fall.10s <- matrix(Y_total.fall.10s, ncol=1)

fall <- cbind(Y_t.fall.80s, Y_t.fall.90s, Y_t.fall.00s, 
              Y_t.fall.10s)
colnames(fall) <- c('Eighties', 'Nineties', 'Aughties', 'Tens')
fall.change <- matrix(data=NA, nrow=nrow(fall), ncol=(ncol(fall)-1))
for(i in 2:ncol(fall)){
  fall.change[,(i-1)] <- 
    ((fall[,(i)] - fall[,(i-1)]) / fall[,(i-1)]) * 100
}
colnames(fall.change) <- c('% Change 80s - 90s',
                           '% Change 90s - 00s',
                           '% Change 00s - 10s')

fall.firstlast <- fall[,'Tens'] - fall[,'Eighties']

fall.change <- cbind(fall.change, fall.firstlast)
colnames(fall.change)[4] <- '% Change 80s - 10s'

# Spring
grid.fall <- data.frame(
  y =NA,
  lat=1, lon=1,
  Decade = NA,
  Season = NA
)
grid.fall <- st_as_sf(grid.fall, coords=c('lon', 'lat'),
                        crs="EPSG:26919")
st_geometry(grid.fall) <- 'grid'

for(i in 1:ncol(fall.change)){
  Zlim = range(fall.change)
  Points_orig = sp::SpatialPointsDataFrame(
    coords = loc_g, 
    data = data.frame(y = fall.change[,i]), 
    proj4string = CRS_orig)
  
  Points_LongLat = sp::spTransform(Points_orig, 
                                   sp::CRS("+proj=longlat"))
  Points_proj = sp::spTransform(Points_orig, 
                                CRS_proj)
  cell.size = mean(diff(Points_proj@bbox[1, ]), 
                   diff(Points_proj@bbox[2, 
                   ]))/floor(sqrt(n_cells))
  Points_sf = sf::st_as_sf(Points_proj)
  grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
  grid_i = sf::st_intersects(Points_sf, grid)
  grid = sf::st_sf(grid, y = tapply(Points_sf$y, 
                                    INDEX = factor(as.numeric(grid_i), 
                                                   levels = 1:length(grid)),
                                    FUN = mean, na.rm = TRUE))
  grid$Decade <- colnames(fall.change)[i]
  grid$Season <- 'Fall'
  
  grid.fall <- rbind(grid.fall, grid)
  
  rm(grid, grid_i, Points_LongLat, Points_orig, Points_proj,
     Points_sf, Zlim)
}

grid.fall <- grid.fall[!is.na(grid.fall$Season),]
grid.fall$Decade <- factor(grid.fall$Decade,
                             levels=c('% Change 80s - 90s',
                                      '% Change 90s - 00s',
                                      '% Change 00s - 10s',
                                      '% Change 80s - 10s'))

rm(list=setdiff(ls(), c('coast', 'CRS_orig', 'CRS_proj',
                        'loc_g', 'regions', 'fun', 'n_cells',
                        'Y_gt', 
                        'grid.spring', 'grid.fall')))

grid <- rbind(grid.spring, grid.fall)
grid$Season <- factor(grid$Season, 
                      levels=c('Spring', 'Fall'))

rm(list=setdiff(ls(), c('coast', 'regions', 'grid')))

#### Plot decadal change ####
s.8 <- ggplot() +
    geom_sf(data=grid[grid$Decade == '% Change 80s - 90s'&
                        grid$Season == 'Spring',], 
            aes(fill=y, col=y)) +
    scale_fill_gradient2(low='#893D8A', high='#04630F',
                         midpoint=0,
                         na.value = 'transparent') +
    scale_color_gradient2(low='#893D8A', high='#04630F',
                          midpoint=0,
                          na.value = 'transparent') +
    geom_sf(data=coast, fill='gray', col='darkgray') +
    
    geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
    
    coord_sf(xlim=c(st_bbox(grid)[1], 
                    st_bbox(grid)[3]),
             ylim=c(st_bbox(grid)[2], 
                    st_bbox(grid)[4])) +
    labs(col=expression(paste(Delta, " Abund (%)")), 
         fill=expression(paste(Delta, " Abund (%)")),
         y='Spring') +
    theme(legend.position = 'inside',
          legend.direction = 'horizontal',
          legend.position.inside = c(0.5, 0.08),
          legend.background = element_rect(fill='transparent',
                                           linewidth =0),
          axis.text.x = element_text(size=8, angle=20),
          axis.text.y = element_text(size=8),
          strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
          legend.title = element_text(size=8),
          legend.text = element_text(size=6))

s.9 <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 90s - 00s'&
                      grid$Season == 'Spring',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

s.0 <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 00s - 10s'&
                      grid$Season == 'Spring',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

f.8 <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 80s - 90s'&
                      grid$Season == 'Fall',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)")),
       y='Fall') +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

f.9 <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 90s - 00s'&
                      grid$Season == 'Fall',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

f.0 <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 00s - 10s'&
                      grid$Season == 'Fall',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

overall <- egg::ggarrange(
          tag_facet(s.8 +
                      theme(axis.ticks.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(r = 1, l=3) ) +
                      facet_wrap(~"% Change 80s - 90s"),
                    tag_pool = "a"), 
          tag_facet(s.9 + 
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(r = 1, l = 1) ) +
                      facet_wrap(~"% Change 90s - 00s"), 
                    tag_pool = "b" ), 
          tag_facet(s.0 + 
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(l = 1, r=3)  ) +
                      facet_wrap(~"% Change 00s - 10s"),
                    tag_pool = "c"),
          
          tag_facet(f.8 +
                      theme(axis.ticks.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(r = 1, l=3) ),
                    tag_pool = "d"), 
          tag_facet(f.9 + 
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(r = 1, l = 1) ), 
                    tag_pool = "e" ), 
          tag_facet(f.0 + 
                      theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = margin(l = 1, r=3)  ),
                    tag_pool = "f"),
          nrow = 2)

ggsave(filename = here('VAST_runs/large/Overall_BC/Decadal_Change.png'),
       overall,
       height = 8.5, width = 11, units='in')

f.all <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 80s - 10s'&
                      grid$Season == 'Fall',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

s.all <- ggplot() +
  geom_sf(data=grid[grid$Decade == '% Change 80s - 10s'&
                      grid$Season == 'Spring',], 
          aes(fill=y, col=y)) +
  scale_fill_gradient2(low='#893D8A', high='#04630F',
                       midpoint=0,
                       na.value = 'transparent') +
  scale_color_gradient2(low='#893D8A', high='#04630F',
                        midpoint=0,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], 
                  st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], 
                  st_bbox(grid)[4])) +
  labs(col=expression(paste(Delta, " Abund (%)")), 
       fill=expression(paste(Delta, " Abund (%)"))) +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.5, 0.08),
        legend.background = element_rect(fill='transparent',
                                         linewidth =0),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))

ts.change <- egg::ggarrange(
  s.all +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()) +
              facet_wrap(vars(Season)),
  f.all + 
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()) +
              facet_wrap(vars(Season)), 
  nrow = 1)

ggsave(filename = here('VAST_runs/large/Overall_BC/TimeSeries_Change.png'),
       ts.change,
       height = 8.5, width = 11, units='in')
