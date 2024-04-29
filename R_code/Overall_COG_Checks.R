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
cog.wgom <- read.csv(here("VAST_runs/medium/Overall_BC/WGOM/COG_WGOM.csv"))
cog.egom <- read.csv(here('VAST_runs/medium/Overall_BC/EGOM/COG_EGOM.csv'))
cog.gbk <- read.csv(here("VAST_runs/medium/Overall_BC/GBK/COG_GBK.csv"))
cog.sne <- read.csv(here('VAST_runs/medium/Overall_BC/SNE/COG_SNE.csv'))
cog.all <- read.csv(here("VAST_runs/medium/Overall_BC/ALL/COG_ALL.csv"))

cog <- rbind(cog.wgom, cog.egom, cog.gbk, cog.sne, cog.all)

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
  group_by(Season, strata) %>% 
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
           aes(col=as.factor(Year)),
          cex=0.8) +
  lemon::facet_rep_grid(Season ~ strata) + 
  labs(col='Year') +
  coord_sf(xlim=c(-74, -66),
           ylim=c(40, 45),
           crs="EPSG:4326")+
  theme(legend.position = 'right')

ggsave(here("VAST_runs/medium/Overall_BC/SF_COG.png"),
       wgom.cog,
       width = 11, height=8.5) 


rm(list=setdiff(ls(), c('coast', 'gom', 'cog')))

all.cog <- ggplot() + 
  geom_line(data=cog,
            aes(x=Year, y=northing.km, col=strata)) +
  geom_point(data=cog,
             aes(x=Year, y=northing.km, col=strata),
             alpha=0.3, cex=0.4) +
  geom_ribbon(data=cog,
              aes(x=Year, ymin=northing.km - n.sd.km,
                  ymax=northing.km + n.sd.km,
                  fill=strata), alpha=0.4) +
  facet_wrap(vars(Season), ncol=1) +
  theme(legend.position = 'bottom')

ggsave(here("VAST_runs/medium/Overall_BC/COG_AllStrata.png"),
       all.cog,
       width = 11, height=8.5) 

#### Range edges ####
re.wgom <- read.csv(here("VAST_runs/medium/Overall_BC/WGOM/RangeEdges_WGOM.csv"))
re.egom <- read.csv(here('VAST_runs/medium/Overall_BC/EGOM/RangeEdges_EGOM.csv'))
re.gbk <- read.csv(here("VAST_runs/medium/Overall_BC/GBK/RangeEdges_GBK.csv"))
re.sne <- read.csv(here('VAST_runs/medium/Overall_BC/SNE/RangeEdges_SNE.csv'))
re.all <- read.csv(here("VAST_runs/medium/Overall_BC/ALL/RangeEdges_ALL.csv"))

re <- rbind(re.wgom, re.egom, re.gbk, re.sne, re.all)

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
  group_by(strata, Quantile, Season) %>% 
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
  lemon::facet_rep_grid(Season ~ strata) + 
  coord_sf(xlim=c(-72, -66),
           ylim=c(40, 45),
           crs="EPSG:4326") +
  theme(legend.position = 'bottom')

ggsave(here('VAST_runs/medium/Overall_BC/SF_RangeEdges_WGOM.png'),
       wgom.re,
       width = 11, height = 8.5)

# Seasonal indices of abundance
rm(list=setdiff(ls(), c('coast', 'gom')))

ind <- read.csv(here("VAST_runs/medium/Overall_BC/Index.csv"))
ind <- ind %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  dplyr::select(-Units, -Category) %>% 
  rename(Std.Err = Std..Error.for.Estimate,
         Std.Err.ln = Std..Error.for.ln.Estimate.)
load(here('VAST_runs/medium/Overall_BC/Overall_BC_mediumcod_allstrat_natsplin_fsON.RData'))
strats <- data.frame(
  Stratum = unique(ind$Stratum),
  Stock = fit$settings$strata.limits
)

ind <- left_join(ind, strats, by=c('Stratum'))

ind <- ind %>% 
  dplyr::select(-Stratum) %>% 
  rename(Stratum = STRATA)

ggplot() + 
  geom_line(data=ind[ind$Stratum != 'ALL',], 
            aes(x=Year, y=Estimate, col=Stratum)) + 
  geom_point(data=ind[ind$Stratum != 'ALL',], 
            aes(x=Year, y=Estimate, col=Stratum),
            alpha=0.4, cex=0.5) + 
  geom_ribbon(data=ind[ind$Stratum != 'ALL',],
              aes(x=Year, ymin=Estimate - Std.Err, ymax=Estimate+Std.Err,
                  fill=Stratum), alpha=0.2) +
  scale_color_viridis_d(option='turbo') +
  scale_fill_viridis_d(option='turbo') +
  facet_wrap(vars(Season), ncol=1, scales='free_y') +
  theme(legend.position = 'bottom')



