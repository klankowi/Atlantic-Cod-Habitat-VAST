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
small.ao <- read.csv(here('VAST_runs/small/Overall_BC/with_AMO/ALL/AreaOcc.csv'))
small.ao$Size <- 'Small'
medium.ao <- read.csv(here('VAST_runs/medium/Overall_BC/ALL/AreaOcc.csv'))
medium.ao$Size <- 'Medium'
large.ao <- read.csv(here('VAST_runs/large/Overall_BC/with_cobble/ALL/AreaOcc.csv'))
large.ao$Size <- 'Large'

ao <- rbind(small.ao, medium.ao, large.ao)

regions <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
regions <- st_make_valid(regions)
regions$area <- st_area(regions)

regions <- regions %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(STOCK, area) %>% 
  rename(strata = STOCK) %>% 
  unique() %>% 
  mutate(area = strip_units(area))

regions <- rbind(regions, data.frame(strata='ALL', area=sum(regions$area)))
rownames(regions) <- NULL

ao <- left_join(ao, regions, by=c('strata'))

ao$pct <- (ao$area.occ / ao$area) * 100
ao$hi <- ((ao$area.occ + ao$std.err) / ao$area) * 100
ao$lo <- ((ao$area.occ - ao$std.err) / ao$area) * 100

ggplot() + 
  geom_line(data=ao,
            aes(x=Year, y=pct, col=Size)) +
  geom_ribbon(data=ao,
              aes(x=Year, ymin=lo, ymax=hi, fill=Size),
              alpha = 0.2) +
  # geom_smooth(data=ao,
  #             aes(x=Year, y=pct, col=Size, fill=Size),
  #             alpha=0.2, method='lm')
  facet_grid(strata ~ Season,
             scales='free_y') +
  labs(y='% of strata occupied')
