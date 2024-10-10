rm(list=ls())

library(VAST)
library(here)
library(tidyverse)
library(sf)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "inside",
                legend.background = element_rect(fill='transparent', colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

cog <- read.csv(here("VAST_runs/small/Overall_BC/ALL/RangeEdges_ALL.csv"))

setwd(here('VAST_runs/small/Overall_BC'))

cog <- cog %>% 
  rename(northing = Northing.Est,
         easting = Easting.Est,
         n.sd = Northing.SD,
         e.sd = Easting.SD) %>% 
  filter(Quantile != 0.5)

coast <- st_transform(ecodata::coast, crs="EPSG:32619")
codstox <- st_transform(
  st_read(here("Data/GIS/codstox.shp"), quiet=T),
  crs="EPSG:32619"
)

# Spatial plots
cog <- cog %>% 
  mutate(easting = easting * 1000,
         northing = northing * 1000)
cog.sf <- st_as_sf(cog, coords=c('easting', 'northing'), crs="EPSG:32619")

cog.lin <- cog.sf %>% 
  group_by(Quantile) %>% 
  filter(strata == 'ALL') %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

# Spring
SD_plotting.spring <- subset(cog.sf, Season =='Spring')

SD_list <- split(SD_plotting.spring, f=SD_plotting.spring$Quantile)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.spring <- st_multilinestring(do.call("rbind", linestrings))
  t.spring <-  nngeo::st_segments(t.spring)
  t.spring <- st_sf(t.spring)
  t.spring$Year <- seq(1982, 2020, 1)
  st_crs(t.spring) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Quantile[1]
  
  SD_list[[i]] <- t.spring
  SD_list[[i]]$Quantile <- mod
  rm(points, n, linestrings, t.spring, mod)
}
yeartracks <- do.call(rbind, SD_list)

# Plot
spring.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=codstox, fill='transparent', lwd=0.25) +
  geom_sf(data=yeartracks, aes(col=Year), pch=19, cex=0.5, lwd=0.7, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Quantile))+
  coord_sf(xlim=c(-76, -65.5),
           ylim=c(36.5, 44.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Spring Range Edges') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=8, angle=20),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=8))

# Fall
SD_plotting.spring <- subset(cog.sf, Season =='Fall')

SD_list <- split(SD_plotting.spring, f=SD_plotting.spring$Quantile)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.spring <- st_multilinestring(do.call("rbind", linestrings))
  t.spring <-  nngeo::st_segments(t.spring)
  t.spring <- st_sf(t.spring)
  t.spring$Year <- seq(1982, 2020, 1)
  st_crs(t.spring) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Quantile[1]
  
  SD_list[[i]] <- t.spring
  SD_list[[i]]$Quantile <- mod
  rm(points, n, linestrings, t.spring, mod)
}
yeartracks <- do.call(rbind, SD_list)

# Plot
fall.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=codstox, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks, aes(col=Year), pch=19, cex=0.5, lwd=0.7, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Quantile)) +
  coord_sf(xlim=c(-76, -65.5),
           ylim=c(36.5, 44.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Fall Range Edges') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=8, angle=20),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=8))

ggsave("Spatial_RE_Spring.png", spring.vis, width=10, height=10)
ggsave("Spatial_RE_Fall.png", fall.vis, width=10, height=10)

