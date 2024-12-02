# Find closest 30 BLLS transects to GMRI

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

blls <- read.csv('C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/BTS Reports/possible_tsplocs_2010.csv')
bllsorig <- blls
colnames(blls)

blls <- blls %>% 
  filter(survey == "NEFSC BLLS") %>% 
  dplyr::select(startlon, startlat, endlon, endlat, year, fs) %>% 
  unique() %>% 
  as.data.frame()

blls$id <- seq(1, nrow(blls), by=1)

blls <- blls %>% 
  rename(beglat = startlat,
         beglon = startlon) %>% 
  unique() %>% 
  as.data.frame()

blls <- blls %>% 
  pivot_longer(cols=c('beglat', 'endlat'), names_to='position', values_to = 'lat') %>% 
  pivot_longer(cols=c('beglon', 'endlon'), names_to='position2', values_to = 'lon') %>% 
  mutate(check = paste0(position, ' ', position2)) %>% 
  filter(check == 'beglat beglon' | check == 'endlat endlon') %>% 
  as.data.frame()

blls$cond[blls$check == 'beglat beglon'] <- 'start'
blls$cond[blls$check == 'endlat endlon'] <- 'stop'

blls <- blls %>% 
  dplyr::select(-position, -position2, -check)

blls <- blls[with(blls, order(year, id)),]
rownames(blls) <- NULL

blls <- blls %>% 
  filter(year >=2014)

blls.sf <- st_as_sf(blls, coords=c('lon', 'lat'), crs="EPSG:4326")

blls.sf <- st_transform(blls.sf, "EPSG:26919")

strat <- st_read(
  'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/GIS/BLLS_Strat.shp', quiet=T)

strat <- st_transform(strat, st_crs(blls.sf))

coast <- st_read('C:/Users/klankowicz/Documents/GitHub/CBASS/GIS/us_medium_shoreline_poly.shp', quiet=T)
coast <- st_transform(coast, st_crs(blls.sf))

# Add port
gmri <- data.frame(LAT = 43.650216, 
                   LON = -70.253123)
gmri <- st_as_sf(gmri, coords=c('LON', 'LAT'), crs="EPSG:4326")
gmri$port <- 'PORTLAND'
gmri <- st_transform(gmri, crs="EPSG:26919")

# Find closest 50
rownames(blls.sf) <- NULL
blls.sf$dist <- st_distance(x = blls.sf, y = gmri, by_element=F)

blls <- sfheaders::sf_to_df(blls.sf, fill=T)

blls.dat <- blls %>% 
  dplyr::select(id, dist, fs, year) %>% 
  unique() %>% 
  as.data.frame()

blls.sf <- blls.sf[with(blls.sf, order(dist, decreasing = F)),]
tabfreq <- as.data.frame(table(blls.sf$id[1:56]))
nrow(tabfreq)

# Closest 30
close <- blls.sf[blls.sf$id  %in% 
                 blls.sf$id[1:56],]

close <- close %>% 
  group_by(id) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

close.sf <- left_join(close, blls.dat, by=c('id'))

st_geometry(close.sf) <- 'geometry'

closemap <- ggplot() + 
  geom_sf(data=coast, fill='gray', col='darkgray') +
  geom_sf(data=gmri,col='black') +
  geom_sf(data=strat, fill=NA, col='darkgray') +
  geom_sf(data=close.sf,
          aes(col=fs), lwd=1) +
  facet_wrap(vars(id)) +
  geom_sf_text(data=close.sf, aes(label = id, col=fs), alpha=0.8,
               ) +
  
  coord_sf(xlim=c(-70.28, -69.77),
           ylim=c(43.08, 43.34),
           crs="EPSG:4326") +
  theme(legend.position = 'bottom') +
  labs(col='Strata', x='', y='')
closemap


close.sf <- st_transform(close.sf, crs="EPSG:4326")

close <- sfheaders::sf_to_df(close.sf, fill=T)
close <- close %>% 
  dplyr::select(-sfg_id, -linestring_id, -dist) %>%
  unique() %>% 
  rename(lon=x, lat=y) %>% 
  as.data.frame()

close <- close[with(close, order(id)),]

close$line.end <- rep(c('start', 'stop'), (nrow(close) / 2))

close <- dplyr::select(close, id, line.end, lon, lat, fs, year)

ggsave(plot=closemap,
       filename='C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/Survey_Data/Closest30_unique.pdf',
       width = 11, height = 8.5, units = 'in', dpi = 300)
st_write(close.sf,
         'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/GIS/Closest30_HabSurv.shp',
         delete_layer = TRUE)
write.csv(close,
          'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/Survey_Data/Closest30_HabSurv.csv',
          row.names = F)
