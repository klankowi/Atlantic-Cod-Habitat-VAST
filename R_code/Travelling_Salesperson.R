# "Traveling salesperson" analysis of stratified random selection sites

rm(list=ls())

# Load library
library(tidyverse)
library(sf)
library(here)
library(TSP)
library(udunits2)

# Set seed
set.seed(5423)

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

# Load data
survs <- st_read(here('Data/GIS/HabSurv_Locations.shp'))

# Reachable from Portland
port <- survs[survs$port != 'GLOUCESTER',]
table(port$fs)

# Split by FS
fslist <- split(port, f=port$fs)

# Random selection
for(i in 1:length(fslist)){
  keepers <- sample(1:nrow(fslist[[i]]), size=25, replace = F)
  fslist[[i]] <- fslist[[i]][keepers,]
}
port <- do.call(rbind, fslist)
rownames(port) <- NULL

st_geometry(port) <- 'geometry'

coast <- st_transform(ecodata::coast, st_crs(port))

ggplot() + 
  geom_sf(data=coast, fill='gray', col='darkgray') +
  geom_sf(data=port,
          aes(col=fs),
          lwd=1)+
  coord_sf(xlim=c(-70.75, -69.25),
           ylim=c(42.8, 43.4)) +
  labs(col='Strata') +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=10, angle=20),
        axis.text.y = element_text(size=10))

# Get midpoint
mids <- st_centroid(port)
mids <- mids %>% 
  dplyr::select(-survey, -year)

# Load Portland dock
portland <- data.frame(
  id = "Portland",
  cobbl_p = NA, gravl_p =NA, rock_p=NA, mud_p=NA, sand_p=NA,
  bthy_dp = NA, bttm_ty=NA, depstrt=NA, fs=NA, port='PORTLAND',
  lon = -70.255005, 
  lat = 43.651156)
portland <- st_as_sf(portland, coords=c('lon', 'lat'))
st_crs(portland) <- 'EPSG:4326'

rand.points <- rbind(portland, mids)
row.names(rand.points) <- NULL
head(rand.points)
length(unique(rand.points$id))

ggplot() + 
  geom_sf(data=coast, fill='gray', col='darkgray') +
  geom_sf(data=mids,
          aes(col=fs),
          cex=1)+
  coord_sf(xlim=c(-70.75, -69.25),
           ylim=c(42.8, 43.4)) +
  labs(col='Strata') +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=10, angle=20),
        axis.text.y = element_text(size=10))

# Traveling salesperson problem - what is the most efficient order 
tsp <- ((st_distance(rand.points)))
colnames(tsp) <- rand.points$id
tsp <- TSP(as.dist(tsp))
tour.list <- vector(mode='list', 1000)
for(i in 1:length(tour.list)){
  tour.list[[i]] <- solve_TSP(tsp, method = 'farthest_insertion')
}
length.list <- vector(mode='list', 1000)
for(i in 1:length(length.list)){
  length.list[[i]] <- tour_length(tour.list[[i]])
}
length.list <- do.call(rbind, length.list)
length.list <- as.data.frame(length.list)
length.list$run <- seq(1, nrow(length.list))
shortest.run <- length.list$run[length.list$V1 == min(length.list$V1)]
shortest.run <- as.integer(shortest.run)
# Shortest path of 100 runs
path <- cut_tour(tour.list[[930]], "Portland", exclude_cut = FALSE)
path.lines <- st_sfc(
  st_cast(
    do.call(c, st_geometry(rand.points[c(path),])),
    'LINESTRING'
  ), crs = 4326
)
path.lines <- st_as_sf(path.lines)

path.lines <- nngeo::st_segments(path.lines, progress = TRUE)
path.lines$id <- seq(1:nrow(path.lines))

ggplot() + 
  geom_sf(data=path.lines)
