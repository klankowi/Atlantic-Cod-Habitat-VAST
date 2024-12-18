### Plots for the CINAR PPT ###

rm(list=ls())          

# Load packages
library(here)
library(tidyverse)
library(VAST)
library(sf)

# Load VAST fit data
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))

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
regions <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

closed <- st_transform(st_read(here('Data/GIS/closed_areas_wgs.shp')),
                       st_crs(regions), quiet=T)

# Total spatial density
Y_total <- (strip_units(Y_gt))
Y_total <- as.data.frame(Y_total)
Y_total <- as_tibble(Y_total)
Y_total.spring <- Y_total %>% dplyr::select(contains('Spring'))
Y_total.spring <- rowSums(Y_total.spring)
Y_t.spring <- matrix(Y_total.spring, ncol=1)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y = Y_t.spring[,1]), proj4string = CRS_orig)
Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
Points_proj = sp::spTransform(Points_orig, CRS_proj)
Zlim = range(Y_t.spring, na.rm = TRUE)
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
spring <- grid

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
fall <- grid

grid <- rbind(spring, fall)
grid <- grid %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

tsp <- ggplot() +
  geom_sf(data=grid, 
          aes(fill=log(y), col=log(y))) +
  scale_fill_viridis_c(option='viridis', alpha=0.8, direction=1,
                       na.value = 'transparent') +
  scale_color_viridis_c(option='viridis', alpha=0.8, direction=1,
                        na.value = 'transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA, col='black', lwd=0.1) +
  
  geom_sf_pattern(data=closed,
                  pattern = 'stripe',
                  fill    = 'transparent',
                  colour  = NA,
                  pattern_spacing = 0.007,
                  pattern_alpha=0.15,
                  pattern_fill= 'transparent',
                  pattern_color='gray10') +
  
  geom_sf(data=closed, fill=NA, col='gray40',
          lwd=0.2) +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='log(sum(Abund.))', fill='log(sum(Abund.))') +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.3, 0.1),
        legend.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(size=8, angle=20),
        axis.text.y = element_text(size=8),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm"),
                                  size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(vars(Season))+
  ggtitle('Large cod summed abundance, 1982-2021')
tsp
ggsave(plot=tsp,
       here('Presentations/largeSSD.png'))

rm(list=ls())
small <- read.csv(here('VAST_runs/small/Overall_BC/ALL/Index.csv'))
med <- read.csv(here('VAST_runs/medium/Overall_BC/ALL/Index.csv'))
lar <- read.csv(here('VAST_runs/large/Overall_BC/ALL/Index.csv'))

index <- rbind(small, med, lar)

index <- index %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  separate(Category, into=c('Category', 'Drop')) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Size = factor(Category, levels=c('Small', 'Medium', 'Large')),
         Stratum = factor(Stratum, levels=c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM'))) %>% 
  dplyr::select(-Units, -Category, -Drop) %>% 
  rename(Std.Error = Std..Error.for.Estimate) %>% 
  mutate(Estimate = Estimate / 1000000,
         Std.Error = Std.Error / 1000000)

ind <- ggplot(index[index$Stratum == 'ALL',]) +
  geom_point(aes(x=Year, y=Estimate, col=Size), cex=0.5) +
  geom_line(aes(x=Year, y=Estimate, col=Size)) +
  geom_ribbon(aes(x=Year, ymin=Estimate-Std.Error, ymax=Std.Error+Estimate,
                  fill=Size), alpha=0.3) +
  facet_wrap(vars(Season), ncol =1, scales='free') +
  labs(y='Relative Abundance', x='Year') +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        #legend.position.inside = c(0.7, 0.7),
        legend.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm"),
                                  size=10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5))
ggsave(plot=ind,
       here('Presentations/Index.png'))

stox <- st_transform(st_read(here('Data/GIS/codstox.shp')),
                     crs="EPSG:26919")
coast <- st_transform(ecodata::coast, st_crs(stox))

map <- ggplot() +
  geom_sf(data=stox, aes(fill=STOCK), alpha=0.8) +
  geom_sf(data=coast) +
  coord_sf(xlim=c(-77, -65.5),
           ylim=c(36, 45),
           crs="EPSG:4326") +
  labs(fill='Biological\nstock area') +
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.6, 0.08),
        legend.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm"),
                                  size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
ggsave(plot=map,
       here('Presentations/map.png'))  
