rm(list=ls())          

# Load packages
library(here)
library(tidyverse)
library(VAST)
library(sf)

# Install unitless
#install_unit(symbol='unitless', def='unitless', name='unitless')

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

load(here('VAST_runs/medium/Overall_BC/All/Overall_BC_medcod_allstrat_natsplin_fsON_ALL.RData'))

# Extract knots
knots <- fit$spatial_list$Kmeans$centers
knotsm <- knots* 1000
knotsm <- as.data.frame(knotsm)

knots <- st_as_sf(knotsm, coords=c('E_km', 'N_km'), crs="EPSG:26919")

knotsdist <- as.matrix(st_distance(knots))
summary(knotsdist)

# Extract FS
knots <- fit$spatial_list$loc_g
knotsm <- knots* 1000
knotsm <- as.data.frame(knotsm)

knots <- st_as_sf(knotsm, coords=c('E_km', 'N_km'), crs="EPSG:26919")

knotsdist <- as.matrix(st_distance(knots))

f <- nngeo::st_nn(knots, knots, k=2)
f <- do.call(rbind, f)
mean(knotsdist[f]) / 1000

rm(list=ls())

# Pull anisotropic distances
load(here('VAST_runs/small/Overall_BC/All/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.RData'))
small <- fit
load(here('VAST_runs/medium/Overall_BC/All/Overall_BC_medcod_allstrat_natsplin_fsON_ALL.RData'))
medium <- fit
load(here('VAST_runs/large/Overall_BC/All/Overall_BC_larcod_allstrat_natsplin_fsON_ALL.RData'))
large <- fit
rm(list=setdiff(ls(), c('small', 'medium', 'large')))

aniso.dist <- data.frame(
  size=c(rep('Small', 2),
         rep('Medium', 2),
         rep('Large', 2)),
  Kappa=c(rep(c(1,2), 3)),
  Aniso=c(small$Report$Range_raw1,
          small$Report$Range_raw2,
          medium$Report$Range_raw1,
          medium$Report$Range_raw2,
          large$Report$Range_raw1,
          large$Report$Range_raw2)
)
