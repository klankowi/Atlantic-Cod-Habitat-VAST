### Make Figure 4, stock-, size-, and season-specific COG maps ####
#### CJFAS 2024 ####
rm(list=ls())

library(tidyverse)
library(here)
library(sf)
library(stars)
library(raster)
library(ggh4x)

# Function to tag facets without removing strips (based on egg::tag_facet)
tag_facet2 <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

# Color palette function
colorpal <- function(pal, n){colorRampPalette(pal)(n)}

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

# Load background data
# Load coast
coast <- st_transform(ecodata::coast,
                      crs="EPSG:32619")

# Load stocks
stox <- st_read(here('Data/GIS/codstox.shp'), 
                quiet=T)
stox <- st_transform(stox,
                     crs="EPSG:32619")

# Load closed areas
closed <- st_read(here("Data/GIS/closed_areas_wgs.shp"), quiet=T)
closed <- st_transform(closed, crs="EPSG:32619")

# Load isobath
bathy <- read_stars(here('Data/Bathymetry/gebco_2023.tif'))
bathy <- as(bathy, 'Raster')
contour <- rasterToContour(bathy, levels=c(-50,-100, -200))
contour <- st_as_sf(contour)
contour$level <- factor(contour$level, levels=c(-50, -100, -200))
contour <- st_transform(contour, crs="EPSG:32619")

# Load all indices from all sizes
sizes <- c('small', 'medium', 'large')
strata <- c('EGOM', 'GBK', 'SNE', 'WGOM')
for(i in 1:length(sizes)){
  cog <- read.csv(paste0(here('VAST_runs'), '/',
                         sizes[i], 
                         '/Overall_BC/ALL/COG_ALL.csv'))
  for(j in 1:length(strata)){
    strat <- read.csv(paste0(here('VAST_runs'), '/',
                             sizes[i], 
                             '/Overall_BC/',
                           strata[j],
                           '/COG_',
                           strata[j],
                           '.csv'))
    cog <- rbind(cog, strat)
    rm(strat)
  }
  cog$Size <- paste0(str_to_sentence(sizes[i]))
  assign(paste0(sizes[i]), cog)
  rm(cog)
}

# Combine to one DF
cog <- rbind(small, medium, large)

# Clean
cog <- cog %>% 
  mutate(Year = as.numeric(Year),
         Size = factor(Size, levels=c('Small', 'Medium', 'Large')),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Stratum = factor(strata, levels=c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM'))) %>% 
  dplyr::select(-units, -strata)

# Base plot
# Spatial plots
cog <- cog %>% 
  filter(Stratum != 'ALL')
cog.sf <- st_as_sf(cog, coords=c('easting', 'northing'), crs="EPSG:32619")

cog.lin <- cog.sf %>% 
  group_by(Size, Season, Stratum) %>% 
  arrange(Year) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

# Spring
SD_plotting <- cog.sf
SD_plotting$ID <- paste0(SD_plotting$Size, ' ',
                         SD_plotting$Season, ' ', 
                         SD_plotting$Stratum)

SD_list <- split(SD_plotting, f=SD_plotting$ID)

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
  t.line <- st_multilinestring(do.call("rbind", linestrings))
  t.line <-  nngeo::st_segments(t.line)
  t.line <- st_sf(t.line)
  t.line$Year <- seq(1982, 2020, 1)
  st_crs(t.line) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$ID[1]
  strat <- strsplit(SD_list[[i]]$ID[1], split=' ')[[1]][3]
  if(strat == 'EGOM'){
    usecol = c('#760C04', '#c2a3a3')
    paluse = colorpal(usecol, nrow(t.line))
  }
  if(strat == 'GBK'){
    usecol = c('#384D05', '#CBDCA7')
    paluse = colorpal(usecol, nrow(t.line))
  }
  if(strat == 'SNE'){
    usecol = c('#0B3147', '#A3DEE1')
    paluse = colorpal(usecol, nrow(t.line))
  }
  if(strat == 'WGOM'){
    usecol = c('#2D0949', '#EBD3FD')
    paluse = colorpal(usecol, nrow(t.line))
  }
  t.line$colors <- paluse
  
  rm(usecol, paluse)
  
  SD_list[[i]] <- t.line
  SD_list[[i]]$ID <- mod
  rm(points, n, linestrings, t.line, mod)
}
yeartracks <- do.call(rbind, SD_list)
rownames(yeartracks) <- NULL
yeartracks <- yeartracks %>% 
  separate(ID, into=c('Size', 'Season', 'Stratum')) %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large')),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Stratum = factor(Stratum, levels=c('EGOM', 'GBK', 'SNE', 'WGOM')))
st_geometry(yeartracks) <- 't.line'

# Plot
fig4 <- ggplot() +

  geom_sf(data=coast) +
  
  geom_sf(data=stox, aes(fill=STOCK), alpha=0) +
  
  geom_sf(data=yeartracks[yeartracks$Season == 'Fall' &
                          yeartracks$Stratum == 'EGOM' &
                          yeartracks$Size == 'Small',],
          aes(col=Year), alpha=0,
          linewidth=0.6)+
  
  scale_color_continuous(
    low='black', high='gray90',
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks,
          aes(col=colors),
          pch=19, cex=0.5, lwd=0.5, alpha=0.7) +
  
  scale_color_identity(guide='none') +

  facet_grid2(c("Season", "Size")) +
  
  labs(x='', y='', fill='Biological\nStock') +
  
  coord_sf(xlim=c(-74, -66),
           ylim=c(39.5, 45),
           crs="EPSG:4326") +
  
  theme(legend.position = 'bottom',
        legend.key = element_blank(),
        axis.text.x=element_text(size=9, angle=30),
        axis.text.y=element_text(size=9, angle=0),
        strip.background = element_rect(fill='lightgray'),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.box.spacing = margin(0,0,0,0)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)),
         color = guide_legend(override.aes = list(alpha = 1)))

# Label facets
fig4 <- tag_facet2(fig4, open='', close='', col='gray30')

# Save local
ggsave(here('Documentation/Figures/Fig 4.pdf'),
      fig4,
      height=(18.2 * 0.75), width=18.2, units='cm',
      dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 4.pdf",
      fig4,
      height=(18.2 * 0.75), width=18.2, units='cm',
      dpi = 600)