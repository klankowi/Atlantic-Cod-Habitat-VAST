rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)
library(wesanderson)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Shapefiles and centroid calculation
gom <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
gom <- st_make_valid(gom)
coast <- st_transform(ecodata::coast, crs=st_crs(gom))
centroid <- st_centroid(gom[gom$STOCK == 'WGOM',])

#### Center of gravity ####
cog <- read.csv(here("VAST_runs/medium/Overall_BC/COG_WGOM.csv"))

rm(list=setdiff(ls(), c('cog', 'gom', 'centroid', 'coast')))

cog <- cog %>%
  mutate(easting.km = easting / 1000,
         northing.km = northing / 1000,
         e.sd.km = e.sd / 1000,
         n.sd.km = n.sd / 1000) %>%
  dplyr::select(-units)

cog$Model <- 'BC'

ccord <- sfheaders::sf_to_df(centroid, fill=T)
ccord <- ccord %>% 
  dplyr::select(x, y) %>% 
  rename(lon = x,
         lat = y) %>% 
  mutate(easting.km = lon/1000,
         northing.km = lat/1000)

cog$center.lon <- cog$easting.km - ccord$easting.km
cog$center.lat <- cog$northing.km - ccord$northing.km

# Spatial plots
cog.sf <- st_as_sf(cog, coords=c('easting', 'northing'), crs="EPSG:32619")

cog.lin <- cog.sf %>% 
  group_by(Season) %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

wgom.cog <- ggplot() +
  geom_sf(data=coast, col=NA) +
  geom_sf(data=gom, fill=NA) +
  geom_sf(data=cog.lin,
          lwd=0.3) +
  geom_sf(data=cog.sf[cog.sf$Year %in% c(1982, 2021),],
           aes(col=as.factor(Year))) +
  facet_wrap(vars(Season)) +
  labs(col='Year') +
  coord_sf(xlim=c(-71, -67),
           ylim=c(40, 44),
           crs="EPSG:4326")
ggsave(here("VAST_runs/medium/Overall_BC/SF_COG_WGOM.png"),
       wgom.cog,
       width = 11, height=8.5)


rm(list=setdiff(ls(), c('coast', 'gom')))

#### Range edges ####
re <- read.csv(here("VAST_runs/medium/Overall_BC/RangeEdges_WGOM.csv"))
re$Model <- 'BC'

rm(list=setdiff(ls(), c('re', 'gom', 'centroid', 'coast')))

re <- re %>% 
  mutate(Northing.Est=Northing.Est*1000,
         Northing.SD=Northing.SD*1000,
         Easting.Est=Easting.Est*1000,
         Easting.SD=Easting.SD*1000)

# Spatial plots
re.sf <- st_as_sf(re, coords=c('Easting.Est', 'Northing.Est'), 
                  crs="EPSG:32619")

re.lin <- re.sf %>% 
  group_by(Quantile, Season) %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

wgom.re <- ggplot() +
  geom_sf(data=coast, col=NA) +
  geom_sf(data=gom, fill=NA) +
  geom_sf(data=re.lin, 
          aes(col=as.factor(Quantile)),
         lwd=0.3) +
  geom_sf(data=re.sf[re.sf$Year %in% c(1982,2021),],
         aes(pch=as.factor(Year))) +
  labs(color='Quantile', pch='Year') +
  facet_wrap(vars(Season)) +
  coord_sf(xlim=c(-72, -66),
           ylim=c(40, 45),
           crs="EPSG:4326") +
  theme(legend.position = 'bottom')

ggsave(here('VAST_runs/medium/Overall_BC/SF_RangeEdges_WGOM.png'),
       wgom.re,
       width = 11, height = 8.5)
