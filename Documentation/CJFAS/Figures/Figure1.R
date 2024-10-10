### Make Figure 1, Map of spatial domain ####
#### CJFAS 2024 ####
rm(list=ls())

library(sf)
library(tidyverse)
library(here)
library(marmap)
library(stars)
library(raster)
library(ggpattern)

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

# Load coast
coast <- st_transform(ecodata::coast,
                 crs="EPSG:4326")

# Load stocks
stox <- st_read(here('Data/GIS/codstox.shp'), 
                quiet=T)
stox <- st_transform(stox,
                     crs="EPSG:4326")

# Load closed areas
closed <- st_read(here("Data/GIS/closed_areas_wgs.shp"), quiet=T)

# Load isobath
bathy <- read_stars(here('Data/Bathymetry/gebco_2023.tif'))
bathy <- as(bathy, 'Raster')
contour <- rasterToContour(bathy, levels=c(-50,-100, -200))
contour <- st_as_sf(contour)
contour$level <- factor(contour$level, levels=c(-50, -100, -200))

# Plot
fig1 <- ggplot() +

  geom_sf(data=contour, aes(col=level), linewidth=0.4) +
  scale_color_manual(values=c('gray80', 'gray70', 'gray60')) +
  
  geom_sf_pattern(data=closed,
                  pattern = 'stripe',
                  fill    = 'transparent',
                  colour  = NA,
                  pattern_spacing = 0.007,
                  pattern_alpha=0.15,
                  pattern_fill= 'transparent',
                  pattern_color='gray10') +

  geom_sf(data=stox, alpha=0.3, stroke=NA,
          aes(fill=STOCK), col=NA) +
  
  geom_sf(data=coast) +
  
  
  labs(fill='Biological\nStock Unit',
       col='Bathymetric\nContour (m)',
       x='',y='') +
  coord_sf(xlim=c(-76, -65),
           ylim=c(36.75, 45)) +
  
  annotate("text", x = -71.5, y = c(42.275, 42.425), label = list("Bay", "Mass."),
           size=2.6) +
  annotate("text", x = -71.2, y = 42.35, label = "---",
           size=2.6) +

  annotate("text", x = -71.9, y = c(41.625, 41.775), label = list("Bay", "Narra."),
           size=2.6) +
  annotate("text", x = -71.6, y = 41.7, label = "---",
           size=2.6) +

  annotate("text", x = -69.3, y = 41.8, label = "CC Bay",
           size=2.6) +
  annotate("text", x = -69.9, y = 41.8, label = "-----",
           size=2.6) +

  annotate("text", x = -72.0, y = c(41.275, 41.125), label = list("Block Is.","Sound"),
           size=2.6, angle=10) +

  annotate("text", x = -66.2, y = c(41.725, 41.875), label = list("Peak", "Northeast"),
           size=2.6, angle=30) +

  annotate("text", x = -67.8, y = c(41.425, 41.575), label = list("Shoal", "Georges"),
           size=2.6) +

  annotate("text", x = -70.2, y = c(42.125, 42.275), label = list("Bank", "Stellwagen"),
           size=2.6, angle=-30) +

  annotate("text", x = -70.1, y = c(42.925, 43.075), label = list("Ledge", "Jeffreys"),
           size=2.6, angle=50) +

  annotate("text", x = -69.0, y = c(42.825, 42.975), label = list("Ledge", "Cashes"),
          size=2.6, angle=10) +

  annotate("text", x = -69.9, y = c(40.925, 41.075), label = list("Shoals", "Nantucket"),
           size=2.6) +
  
  annotate("text", x = -69.5, y = c(42.425, 42.575), label = list("Basin", "Wilkinson"),
           size=2.6) +
  
  annotate("text", x = -67.6, y = c(43.525, 43.675), label = list("Basin", "Jordan"),
           size=2.6) +
 
  annotate("text", x = -67.02, y = 44.7, label = "G. Manan Ch.",
           size=2.6, angle=55) +
  
  
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.byrow = TRUE,
        legend.position.inside = c(0.98,0.14),
        legend.justification.inside = 'right',
        legend.background = element_rect(fill='transparent', 
                                         linetype = 0)) +
  guides(fill=guide_legend(ncol=2))

ggsave(here('Documentation/Figures/Fig 1.pdf'),
       fig1,
       height=18.2, width=18.2, units='cm',
       dpi = 600)

ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 1.pdf",
       fig1,
       height=18.2, width=18.2, units='cm',
       dpi = 600)
