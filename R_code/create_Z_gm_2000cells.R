# Create Z_gm matrix dataset to calculate center of gravity shifts specific to
# each stock area.

# Would typically run the model for the whole spatial domain first.

# Prepare workspace
rm(list=ls())

# Load libraries
library(VAST)
library(tidyverse)
library(sf)
library(here)
library(ggpubr)
library(ggnewscale)
library(patchwork)

# Load functions
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load VAST fit data
load(here('VAST_runs/add_climate_aja3/add_climate_aja3.Rdata'))
rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs")))

# Call cell centers
cell.centers <- fit$data_list$Z_gm
cell.centers <- st_as_sf(as.data.frame(cell.centers),
                         coords=c('E_km', 'N_km'))
st_crs(cell.centers) <- (fit$extrapolation_list$projargs)
cell.centers$ID <- seq(1, nrow(cell.centers))

# Load polygons
stocks <- st_read(here('Data/GIS/cod_all_strata.shp'))
stocks <- st_transform(stocks, st_crs(cell.centers))

# Plot
coast <- ecodata::coast
coast <- st_transform(coast, st_crs(cell.centers))

ggplot() +
  geom_sf(data=stocks,
          alpha=0.2,
          aes(fill=STRATA, col=STRATA)) +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=cell.centers, 
          pch=19, cex=0.4) +
  coord_sf(xlim=c(st_bbox(stocks)[1],
                  st_bbox(stocks)[3]),
           ylim=c(st_bbox(stocks)[2],
                  st_bbox(stocks)[4]))

# Find overlap
cell.overlaps <- st_intersection(cell.centers, stocks)

# WHich ones did we miss
cell.missing <- cell.centers[cell.centers$ID %notin% cell.overlaps$ID,]
ggplot() +
  geom_sf(data=stocks,
          alpha=0.2,col=NA
          #aes(fill=STRATA, col=STRATA)
          ) +
  #geom_sf(data=coast, fill='gray') +
  geom_sf(data=cell.missing[cell.missing$ID > 1000,], 
          aes(col=as.factor(ID)),
          pch=19, cex=1) +
  coord_sf(xlim=c(300,
                  st_bbox(stocks)[3]),
           ylim=c(4510.893,
                  st_bbox(stocks)[4]))
# It's coastal ones in EGOM, WGOM, SNE. Can assign them properly.
cell.missing$STRATA <- c('SNE', 'WGOM', 'EGOM', 'EGOM', 'SNE',
                         'WGOM', 'EGOM', 'EGOM', 'EGOM', 'EGOM', 'WGOM')
cell.total <- rbind(cell.overlaps, cell.missing)
cell.total <- cell.total[with(cell.total, order(ID)),]

ggplot() +
  geom_sf(data=stocks,
          alpha=0.2,col='gray'
          #aes(fill=STRATA, col=STRATA)
  ) +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=cell.total, 
          aes(col=as.factor(STRATA)),
          pch=19, cex=1) +
  coord_sf(xlim=c(st_bbox(stocks)[1],
                  st_bbox(stocks)[3]),
           ylim=c(st_bbox(stocks)[2],
                  st_bbox(stocks)[4]))

# That works, save.
cell.total <- sfheaders::sf_to_df(cell.total, fill=T)
cell.total <- dplyr::select(cell.total,
                            ID, STRATA, x, y)
colnames(cell.total) <- c('ID', 'STRATA', 'E_km', 'N_km')
write.csv(cell.total,
          here('Data/VAST_input/Z_gm_2000cells.csv'))
