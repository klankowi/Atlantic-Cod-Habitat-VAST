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

blls <- read.csv('C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/BTS Reports/possible_tsplocs_2010.csv')
bllsorig <- blls
colnames(blls)

blls <- blls %>% 
  filter(survey == "NEFSC BLLS") %>% 
  rename(beglat = startlat,
         beglon = startlon) %>% 
  mutate(depstrat = factor(depstrat, levels=c('SHALLOW', 'DEEP')),
         bottom_type = factor(bottom_type, levels=c('SMOOTH', 'ROUGH'))) %>% 
  unique() %>% 
  as.data.frame()

blls <- reshape(blls,
                varying=13:16,
                direction='long',
                v.names=c('lat', 'lon'))
blls <- blls %>% 
  rename(lat = lon,
         lon = lat)

blls <- blls[with(blls, order(year, survey, id)),]
rownames(blls) <- NULL
blls$cond[blls$time == 1] <- 'start'
blls$cond[blls$time == 2] <- 'stop'

blls <- blls %>% 
  dplyr::select(-time) %>% 
  filter(year >=2014)

blls.sf <- st_as_sf(blls, coords=c('lon', 'lat'), crs="EPSG:4326")

blls.sf <- blls.sf %>% 
  group_by(id) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

blls.dat <- blls %>% 
  dplyr::select(-cond, -lat, -lon) %>% 
  unique() %>% 
  as.data.frame()

strat <- st_read(
  'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/GIS/BLLS_Strat.shp', quiet=T)

st_crs(strat) == st_crs(blls.sf)

coast <- st_transform(ecodata::coast, st_crs(blls.sf))

# Add 50 nmi buffer from Gloucester and GMRI (Zach's suggestion)
bufdist <- swfscMisc::convert.distance(50, from='nm', to='km') * 1000

gloucester <- data.frame(LAT = 42.607250, 
                         LON = -70.662152)
gloucester <- st_as_sf(gloucester, coords=c('LON', 'LAT'), crs="EPSG:4326")
gloucester$port <- 'GLOUCESTER'
gloucester <- st_transform(gloucester, crs="EPSG:26919")
gloucester <- st_buffer(gloucester, dist=bufdist)

gmri <- data.frame(LAT = 43.650216, 
                   LON = -70.253123)
gmri <- st_as_sf(gmri, coords=c('LON', 'LAT'), crs="EPSG:4326")
gmri$port <- 'PORTLAND'
gmri <- st_transform(gmri, crs="EPSG:26919")
gmri <- st_buffer(gmri, dist=bufdist)

blls.sf <- st_transform(blls.sf, st_crs(gmri))
sf.gmri <- st_intersection(blls.sf, gmri)
sf.glou <- st_intersection(blls.sf, gloucester)

sf.eith <- sf.glou[sf.glou$id %in% sf.gmri$id,]
sf.gmri <- sf.gmri[sf.gmri$id %notin% sf.glou$id,]
sf.glou <- sf.glou[sf.glou$id %notin% sf.gmri$id,]

sf.eith$port <- 'EITHER'
sf.gmri$port <- 'PORTLAND'
sf.glou$port <- 'GLOUCESTER'

blls.sf <- rbind(sf.eith, sf.gmri, sf.glou)
blls.sf <- st_transform(blls.sf, crs="EPSG:4326")

bases <- rbind(gloucester, gmri)
bases <- st_transform(bases, crs="EPSG:4326")

blls.dat$port <- NULL
blls.sf <- merge(blls.dat, blls.sf, by=c('id'))
st_geometry(blls.sf) <- 'geometry'

gmrimap <- 
ggplot() + 
  geom_sf(data=coast, fill='gray', col='darkgray') +
  geom_sf(data=bases, fill=NA, col='black') +
  geom_sf(data=strat, fill=NA, col='darkgray') +
  geom_sf(data=blls.sf,
          aes(col=fs)) +
  facet_wrap(vars(survey)) +
  coord_sf(xlim=c(-71, -68.75),
           ylim=c(42, 43.5)) +
  theme(legend.position = 'bottom') +
  labs(col='Strata') +
  ggtitle('Relative position of GMRI/Gloucester to 2014-2021 surveys')

blls <- sfheaders::sf_to_df(blls.sf, fill=T)
blls <- blls %>% 
  dplyr::select(-sfg_id, -linestring_id) %>% 
  rename(lon=x, lat=y)

blls <- blls[with(blls, order(id)),]

blls$line.end <- rep(c('start', 'stop'), (nrow(blls) / 2))

ggsave(plot=gmrimap, 
       filename='C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/Survey_Data/Combined_Map.png')
st_write(blls.sf,
         'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/GIS/HabSurv_Locations.shp',
         delete_layer = TRUE)
write.csv(blls, 
          'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/Data/Survey_Data/HabSurv_Locations.csv',
          row.names = F)
