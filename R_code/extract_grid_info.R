rm(list=ls())

library(VAST)
library(here)
library(tidyverse)
library(sf)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "inside",
                legend.background = element_rect(fill='transparent', 
                                                 colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

load(here('VAST_runs/small/AIC/none/nonecovs_smallcod_wholearea_natsplin_fsOFF.RData'))

# Extract knot locations
knots <- fit$spatial_list$latlon_x
knots <- as.data.frame(knots)

knots <- st_as_sf(knots, coords=c('Lon', 'Lat'), crs="EPSG:4326")
coast <- st_transform(ecodata::coast, st_crs(knots))

ggplot() +
  geom_sf(data=coast, col='darkgray', fill='lightgray') +
  geom_sf(data=knots, cex=0.5) +
  coord_sf(ylim=c(35, 48),
           xlim=c(-76, -65))

# Find knot nn distance
dmat <- nngeo::st_nn(knots,knots, k=2, returnDist = TRUE)
nndist <- do.call(rbind, dmat$dist)
nnid <- do.call(rbind, dmat$nn)
dmat <- cbind(nnid, nndist)
dmat <- as.data.frame(dmat)
colnames(dmat) <- c('Knot', 'NN', 'remove', 'Dist.m')
dmat$remove <- NULL

summary(dmat$Dist.m) / 1000

# Find grid cell dimensinos for extrapolation grid
# Extract Data
Y_gt = fit$Report$D_gct[,1,]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
fun = mean

# Call data
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

# Total spatial density
Y_total <- log(strip_units(Y_gt))
Y_total <- as.data.frame(Y_total)
Y_total <- as_tibble(Y_total)
Y_total.spring <- Y_total %>% dplyr::select(contains('Spring'))
Y_total.spring <- rowSums(Y_total.spring)
Y_t.spring <- matrix(Y_total.spring, ncol=1)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y = Y_t.spring[,1]), 
                                         proj4string = CRS_orig)
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
grid = sf::st_sf(grid, y = tapply(Points_sf$y, 
                                  INDEX = factor(as.numeric(grid_i),
                                                 levels = 1:length(grid)), 
                                  FUN = mean, na.rm = TRUE))

ggplot() + 
  geom_sf(data=grid, aes(fill=y), col=NA) +
  scale_fill_continuous(na.value=NA) +
  theme(legend.position='n')

onecell <- grid[1,]   

onecell <- spatialEco::extract.vertices(onecell)
onecell <- onecell[1:4,]

dmat <- nngeo::st_nn(onecell, onecell, k=4, returnDist = T)
nndist <- do.call(rbind, dmat$dist)
nnid <- do.call(rbind, dmat$nn)
dmat <- cbind(nnid, nndist)
dmat <- as.data.frame(dmat)
colnames(dmat) <- c('Knot', 'NN1', 'NN2', 'NN3', 
                    'Self', 'Dist.NN1', 'Dist.NN2', 
                    'Dist.NN3')

summary(dmat$Dist.NN1) / 1000
