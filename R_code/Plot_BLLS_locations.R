# Plot relative position of BLLS to GMRI

rm(list=ls())

library(sf)
library(tidyverse)

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

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

blls <- read.csv('C:/Users/klankowicz/Downloads/lls_spp_station_length_1421.csv')
bllsorig <- blls
colnames(blls)

blls <- blls %>% 
  dplyr::select(year, season, station, stratum_use, bottom_type,
                decdeg_beglat_set, decdeg_beglon_set,
                decdeg_endlat_set, decdeg_endlon_set,
                cruise) %>% 
  rename(beglat = decdeg_beglat_set,
         beglon = decdeg_beglon_set,
         endlat = decdeg_endlat_set,
         endlon = decdeg_endlon_set) %>% 
  mutate(season = factor(season, levels=c('SPRING', 'FALL')),
         bottom_type = factor(bottom_type, levels=c('SMOOTH', 'ROUGH'))) %>% 
  unique() %>% 
  as.data.frame()

blls <- reshape(blls,
                varying=6:9,
                direction='long',
                v.names=c('lon', 'lat'))

blls <- blls[with(blls, order(cruise, year, season, station)),]
rownames(blls) <- NULL
blls$cond[blls$time == 1] <- 'start'
blls$cond[blls$time == 2] <- 'stop'

blls <- blls %>% 
  dplyr::select(id, cruise, year, season, station, stratum_use,
                lon, lat, cond, bottom_type, -time)

blls.sf <- st_as_sf(blls, coords=c('lon', 'lat'), crs="EPSG:4326")

blls.sf <- blls.sf %>% 
  group_by(id) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

strat <- st_read(here('Data/GIS/BLLS_Strat.shp'), quiet=T)

st_crs(strat) == st_crs(blls.sf)

coast <- st_transform(ecodata::coast, st_crs(blls.sf))
gmri <- data.frame(lon = -70.254979, lat = 43.651172)
gmri <- st_as_sf(gmri, coords=c('lon', 'lat'), crs="EPSG:4326")
         
gmrimap <- ggplot() + 
  geom_sf(data=coast, fill='gray', col='darkgray') +
  geom_sf(data=gmri) +
  geom_sf(data=strat, fill=NA, col='darkgray') +
  geom_sf(data=blls.sf) +
  coord_sf(xlim=c(-71, -67.5),
           ylim=c(42, 43.8)) +
  theme(legend.position = 'none') +
  ggtitle('Relative position of GMRI to 2014-2021 BLLS')

ggsave(plot=gmrimap, 
       filename=here('Data/Survey_Data/BLLS_Map.png'))
st_write(blls.sf,
         here('Data/GIS/BLLS_Locations.shp'))
write.csv(blls, here('Data/Survey_Data/BLLS_Locations.csv'),
          row.names = F)
