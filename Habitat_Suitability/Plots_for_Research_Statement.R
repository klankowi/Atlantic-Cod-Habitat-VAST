rm(list=ls())

library(VAST)
library(tidyverse)
library(here)
library(gridExtra)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=0.5, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.ticks = element_blank(),
                strip.background = element_rect(color = "black",fill = "grey90", 
                                                linewidth = 0.5),
                strip.text=element_text(size = 14),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

#### Static habitat ####
seg <- readRDS(here('Habitat_Suitability/stationary_effects_grid.RDS'))
coast <- st_make_valid(st_transform(ecodata::coast, st_crs(seg)))

seg <- seg %>% 
  filter(bathy<400) %>% 
  dplyr::select(bathy, gravel, geometry) %>% 
  rename("Static Condition 1" = bathy, 
         "Static Condition 2" = gravel) %>% 
  pivot_longer(cols=c('Static Condition 1', 'Static Condition 2'),
               names_to = 'Variable',
               values_to = 'Values') 


b <- ggplot(data=seg[seg$Variable == 'Static Condition 1',]) +
  geom_sf(aes(fill=Values, col=Values)) +
  scale_fill_viridis_c(na.value = 'transparent',
                       direction = -1) +
  scale_color_viridis_c(na.value = 'transparent',
                        direction = -1) +
  facet_wrap(vars(Variable)) +
  
  geom_sf(data=coast, fill='gray90', col='gray60') +
  
  theme(legend.position = 'none') +
  coord_sf(xlim=c(st_bbox(seg)[1], st_bbox(seg)[3]),
           ylim=c(st_bbox(seg)[2], st_bbox(seg)[4]))
b

g <- ggplot(data=seg[seg$Variable == 'Static Condition 2',]) +
  geom_sf(aes(fill=Values, col=Values)) +
  scale_fill_viridis_c(na.value = 'transparent',
                       direction = 1) +
  scale_color_viridis_c(na.value = 'transparent',
                        direction = 1) +
  facet_wrap(vars(Variable)) +
  
  geom_sf(data=coast, fill='gray90', col='gray60') +
  
  theme(legend.position = 'none') +
  coord_sf(xlim=c(st_bbox(seg)[1], st_bbox(seg)[3]),
           ylim=c(st_bbox(seg)[2], st_bbox(seg)[4]))
g

#### Catch data ####
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))
rm(list=setdiff(ls(), c('fit', 'b', 'g', 'r','coast', 'seg')))

surveys <- fit$data_frame
surveys <- surveys[surveys$t_i %in% c(1, 79),]
surveys <- st_as_sf(surveys, coords=c('Lon_i', 'Lat_i'),
                    crs="EPSG:4326")
surveys$t_i[surveys$t_i == 1] <- 'Catch: Time 1'
surveys$t_i[surveys$t_i == 79] <- 'Catch: Time n'

surveys <- st_transform(surveys, st_crs(coast))

s <- ggplot(data=surveys) +
  geom_sf(aes(size=strip_units(b_i))) +

  facet_wrap(vars(t_i)) +
  
  geom_sf(data=coast, fill='gray90', col='gray60') +
  
  theme(legend.position = 'none') +
  coord_sf(xlim=c(st_bbox(seg)[1], st_bbox(seg)[3]),
           ylim=c(st_bbox(seg)[2], st_bbox(seg)[4]))
s

#### Density ####
# Extract Data
Y_gt = fit$Report$D_gct[,1,c(1, 79)]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
panel_labels = fit$year_labels[c(1, 79)]
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

for (tI in 1: ncol(Y_gt)) {
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
alldata <- do.call(rbind, Big_Data)

rm(list=setdiff(ls(), c('fit', 'b', 'g', 'r','s', 'coast', 'seg',
                        'surveys', 'alldata')))

alldata$TS[alldata$TS == '1982 Spring'] <- 'Density: Time 1'
alldata$TS[alldata$TS == '2021 Spring'] <- 'Density: Time n'

alldata <- st_transform(alldata, st_crs(coast))

d <- ggplot(data=alldata) +
  geom_sf(aes(fill=log(y), col=log(y))) +
  scale_fill_viridis_c(na.value = 'transparent',
                       direction = 1) +
  scale_color_viridis_c(na.value = 'transparent',
                        direction = 1) +
  
  facet_wrap(vars(TS)) +
  
  geom_sf(data=coast, fill='gray90', col='gray60') +
  
  theme(legend.position = 'none') +
  coord_sf(xlim=c(st_bbox(seg)[1], st_bbox(seg)[3]),
           ylim=c(st_bbox(seg)[2], st_bbox(seg)[4]))
d

#### Bottom temp ####
load("~/GitHub/Atlantic-Cod-Habitat-VAST/Habitat_Suitability/Seasonal_Avg_Bottom_Temps2.RData")
rm(list=setdiff(ls(), c('fit', 'b', 'g', 'r','coast', 'seg',
                        'surveys', 'alldata', 's', 'd', 
                        'big.spring')))

temp <- rbind(big.spring[[1]],
              #big.spring[[2]],
              big.spring[[39]])

temp$Year[temp$Year == 1982] <- 'Dynamic conditions: Time 1'
temp$Year[temp$Year == 2020] <- 'Dynamic conditions: Time n'

temp <- st_transform(temp, st_crs(seg))

t <- ggplot(data=temp) +
  geom_sf(aes(fill=bt, col=bt)) +
  scale_fill_viridis_c(na.value = 'transparent',
                       direction = 1) +
  scale_color_viridis_c(na.value = 'transparent',
                        direction = 1) +
  
  facet_wrap(vars(Year)) +
  
  geom_sf(data=coast, fill='gray90', col='gray60') +
  
  theme(legend.position = 'none') +
  coord_sf(xlim=c(st_bbox(seg)[1], st_bbox(seg)[3]),
           ylim=c(st_bbox(seg)[2], st_bbox(seg)[4]))
t

rm(list=setdiff(ls(),
                c('b', 'r', 'g', 't', 's', 'd', 't')))

ggsave(plot = b, here('Habitat_Suitability/BathyPlot.png'))
ggsave(plot = g, here('Habitat_Suitability/GravelPlot.png'))
ggsave(plot = s, here('Habitat_Suitability/SurveyPlot.png'))
ggsave(plot = t, here('Habitat_Suitability/TemperaturePlot.png'))
ggsave(plot = d, here('Habitat_Suitability/DensityPlot.png'))
