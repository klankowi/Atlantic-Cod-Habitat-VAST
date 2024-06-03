# Compare range edge shifts by size and season

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

# Load data
small.ind <- read.csv(here("VAST_runs/small/Overall_BC/With_AMO/All/RangeEdges_All.csv"))
small.ind$Size <- 'Small'

medium.ind <- read.csv(here("VAST_runs/medium/Overall_BC/All/RangeEdges_All.csv"))
medium.ind$Size <- 'Medium'

large.ind <- read.csv(here("VAST_runs/large/Overall_BC/With_Cobble/All/RangeEdges_All.csv"))
large.ind$Size <- 'Large'

ind <- rbind(small.ind, medium.ind, large.ind)

ind <- ind %>% 
  mutate(Quantile = factor(Quantile, levels=c('0.05', '0.5', '0.95')),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Size = factor(Size, levels=c('Small', 'Medium','Large'))) %>% 
  dplyr::select(-strata)

northing <- ggplot() +
  geom_line(data=ind[ind$Quantile == '0.05',],
            aes(x=Year, y=Northing.Est)) +
  geom_ribbon(data=ind[ind$Quantile == '0.05',],
              aes(x=Year, ymin=Northing.Est - Northing.SD,
                  ymax=Northing.Est + Northing.SD),
              alpha=0.3) +
  facet_grid(Size ~ Season) +
  ggtitle("Seasonal Northing by size")

easting <- ggplot() +
  geom_line(data=ind[ind$Quantile == '0.05',],
            aes(x=Year, y=Easting.Est)) +
  geom_ribbon(data=ind[ind$Quantile == '0.05',],
              aes(x=Year, ymin=Easting.Est - Easting.SD,
                  ymax=Easting.Est + Easting.SD),
              alpha=0.3) +
  facet_grid(Size ~ Season) +
  ggtitle('Seasonal Easting by size')

north <- ind %>% 
  filter(Quantile == '0.95') %>% 
  group_by(Size, Season) %>% 
  mutate(change = lm(Northing.Est ~ Year)$coefficient[2]) %>% 
  dplyr::select(Size, Season, change) %>% 
  unique() %>% 
  as.data.frame()

east <- ind %>% 
  filter(Quantile == '0.95') %>% 
  group_by(Size, Season) %>% 
  mutate(change = lm(Easting.Est ~ Year)$coefficient[2]) %>% 
  dplyr::select(Size, Season, change) %>% 
  unique() %>% 
  as.data.frame()