rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)
library(wesanderson)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Shapefiles and centroid calculation
gom <- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
gom <- st_make_valid(gom)
coast <- st_transform(ecodata::coast, crs=st_crs(gom))
centroid <- st_centroid(gom)

#### Center of gravity ####
# Root
root <- here('VAST_runs/medium/SingleCovs')
ending <- '/COG.csv'

# Single Covs
all <- read.csv(paste0(root, '/all', ending))
none <- read.csv(paste0(root, '/none', ending))

amo <- read.csv(paste0(root, '/amo', ending))
bathy <- read.csv(paste0(root, '/bathy', ending))
bottomtemp <- read.csv(paste0(root, '/bottomtemp', ending))
cobble <- read.csv(paste0(root, '/cobble', ending))
gravel <- read.csv(paste0(root, '/gravel', ending))
mud <- read.csv(paste0(root, '/mud', ending))
nao <- read.csv(paste0(root, '/nao', ending))
rugos <- read.csv(paste0(root, '/rugos', ending))
sand <- read.csv(paste0(root, '/sand', ending))

rm(list=setdiff(ls(), c('all', 'none', 
                        'amo', 'bathy', 'bottomtemp',
                        'cobble', 'gravel', 'mud',
                        'nao', 'rugos', 'sand', 
                        'gom', 'coast', 'centroid')))

all$Model <- 'All Covars'
none$Model <- 'No Covars'

amo$Model <- 'AMO Only'
bathy$Model <- 'Bathy Only'
bottomtemp$Model <- 'Bottom Temp Only'
cobble$Model <- 'Cobble Only'
gravel$Model <- 'Gravel Only'
mud$Model <- 'Mud Only'
nao$Model <- 'NAO Only'
rugos$Model <- 'Rugosity Only'
sand$Model <- 'Sand Only'

cog <- rbind(all, none, amo, bathy, bottomtemp,
             cobble, gravel, mud, nao, rugos, sand)

# withgmrf <- read.csv(paste0(root, '/On', ending))
# justcovs <- read.csv(paste0(root, '/Off', ending))
# 
# withgmrf$Model <- 'With GMRF'
# justcovs$Model <- 'Just covariates'
# 
# cog <- rbind(withgmrf, justcovs)

rm(list=setdiff(ls(), c('cog', 'gom', 'centroid', 'coast')))

cog <- cog %>%
  mutate(easting.km = easting / 1000,
         northing.km = northing / 1000,
         e.sd.km = e.sd / 1000,
         n.sd.km = n.sd / 1000) %>%
  dplyr::select(-units)

cog$Model <- factor(cog$Model,
                    levels=c('All Covars', 'No Covars',
                             'AMO Only', 'Bathy Only', 'Bottom Temp Only',
                             'Cobble Only', 'Gravel Only', 'Mud Only',
                             'NAO Only', 'Rugosity Only', 'Sand Only'))
                    # levels=c('With GMRF', 'Just covariates'))

ccord <- sfheaders::sf_to_df(centroid, fill=T)
ccord <- ccord %>% 
  dplyr::select(x, y) %>% 
  rename(lon = x,
         lat = y) %>% 
  mutate(easting.km = lon/1000,
         northing.km = lat/1000)

cog$center.lon <- cog$easting.km - ccord$easting.km
cog$center.lat <- cog$northing.km - ccord$northing.km

coglist <- split(cog, f=cog$Model)

for(i in 2:length(coglist)){
  coglist[[i]]$easting.resid <- (coglist[[i]]$easting.km - coglist[[1]]$easting.km)^2
  coglist[[i]]$northing.resid <- (coglist[[i]]$northing.km - coglist[[1]]$northing.km)^2
}
coglist[[1]]$easting.resid <- 0
coglist[[1]]$northing.resid <- 0

for(i in 1:length(coglist)){
  coglist[[i]]$easting.SSE <- sum(coglist[[i]]$easting.resid)
  coglist[[i]]$northing.SSE <- sum(coglist[[i]]$northing.resid)
}
cog2 <- do.call(rbind, coglist)

cog2 <- dplyr::select(cog2, Model, easting.SSE, northing.SSE)
cog2 <- unique(cog2)
northing.similar <- cog2[with(cog2, order(northing.SSE, decreasing = F)),]
rownames(northing.similar) <- NULL

easting.similar <- cog2[with(cog2, order(easting.SSE, decreasing = F)),]
rownames(easting.similar) <- NULL

n.spring <- ggplot() +
  geom_point(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lat, col=Model),
            cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lat, col=Model),
            lwd=1, alpha=0.7) +
  #scale_color_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('Spring center of gravity northing - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

n.fall <- ggplot() +
  geom_point(data=cog[cog$Season == 'Fall',], 
             aes(x=Year, y=center.lat, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Fall',], 
            aes(x=Year, y=center.lat, col=Model),
            lwd=1, alpha=0.7) +
  #scale_color_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('Fall center of gravity northing - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

e.spring <- ggplot() +
  geom_point(data=cog[cog$Season == 'Spring',], 
             aes(x=Year, y=center.lon, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lon, col=Model),
            lwd=1, alpha=0.7) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from easting centroid (km)') +
  ggtitle('Spring center of gravity easting - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

e.fall <- ggplot() +
  geom_point(data=cog[cog$Season == 'Fall',], 
             aes(x=Year, y=center.lon, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Fall',], 
            aes(x=Year, y=center.lon, col=Model),
            lwd=1, alpha=0.7) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from easting centroid (km)') +
  ggtitle('Fall center of gravity easting - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/cogshifts_northing_spring.png')),
       n.spring, width = 5.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/cogshifts_northing_fall.png')),
       n.fall, width = 5.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/cogshifts_easting_spring.png')),
       e.spring, width = 5.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/cogshifts_easting_fall.png')),
       e.fall, width = 5.5, height=4.25, units='in')

# Spatial plots
cog.sf <- st_as_sf(cog, coords=c('easting', 'northing'), crs="EPSG:32619")

cog.lin <- cog.sf %>% 
  filter(Model=='All Covars') %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

ggplot() +
  geom_sf(data=coast, col=NA) +
  geom_sf(data=gom, fill=NA) +
  geom_sf(data=cog.lin,
          lwd=0.3) +
  geom_sf(data=cog.sf[cog.sf$Model == 'All Covars' & cog.sf$Year %in% c(1982,2021),],
          aes(col=as.factor(Year))) +
  facet_wrap(vars(Season)) +
  coord_sf(xlim=c(-71, -67),
           ylim=c(41.75, 45),
           crs="EPSG:4326")

rm(list=setdiff(ls(), c('coast', 'gom')))

#### Range edges ####
# Root
root <- here('VAST_runs/medium/SingleCovs')
ending <- '/RangeEdges.csv'

# Single Covs
all <- read.csv(paste0(root, '/all', ending))
none <- read.csv(paste0(root, '/none', ending))

amo <- read.csv(paste0(root, '/amo', ending))
bathy <- read.csv(paste0(root, '/bathy', ending))
bottomtemp <- read.csv(paste0(root, '/bottomtemp', ending))
cobble <- read.csv(paste0(root, '/cobble', ending))
gravel <- read.csv(paste0(root, '/gravel', ending))
mud <- read.csv(paste0(root, '/mud', ending))
nao <- read.csv(paste0(root, '/nao', ending))
rugos <- read.csv(paste0(root, '/rugos', ending))
sand <- read.csv(paste0(root, '/sand', ending))

rm(list=setdiff(ls(), c('all', 'none', 
                        'amo', 'bathy', 'bottomtemp',
                        'cobble', 'gravel', 'mud',
                        'nao', 'rugos', 'sand', 
                        'gom', 'coast', 'centroid')))

all$Model <- 'All Covars'
none$Model <- 'No Covars'

amo$Model <- 'AMO Only'
bathy$Model <- 'Bathy Only'
bottomtemp$Model <- 'Bottom Temp Only'
cobble$Model <- 'Cobble Only'
gravel$Model <- 'Gravel Only'
mud$Model <- 'Mud Only'
nao$Model <- 'NAO Only'
rugos$Model <- 'Rugosity Only'
sand$Model <- 'Sand Only'

re <- rbind(all, none, amo, bathy, bottomtemp,
             cobble, gravel, mud, nao, rugos, sand)

rm(list=setdiff(ls(), c('re', 'gom', 'centroid', 'coast')))

reorig <- re
re <- reorig[reorig$Season == 'Spring',]

ggplot() +
  geom_point(data=re[re$Quantile == '0.05',], 
             aes(x=Year, y=Northing.Est, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Quantile == '0.05',], 
            aes(x=Year, y=Northing.Est, col=Model),
            lwd=1, alpha=0.7) +
  geom_ribbon(data=re[re$Quantile == '0.05',], 
            aes(x=Year, ymin=Northing.Est-Northing.SD, 
                ymax=Northing.Est+Northing.SD,fill=Model),
            alpha=0.3) +
  #scale_color_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  #scale_fill_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('5% quantile range edge northing - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(vars(Model))

ggplot() +
  geom_point(data=re[re$Quantile == '0.05',], 
             aes(x=Year, y=Easting.Est, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Quantile == '0.05',], 
            aes(x=Year, y=Easting.Est, col=Model),
            lwd=1, alpha=0.7) +
  geom_ribbon(data=re[re$Quantile == '0.05',], 
              aes(x=Year, ymin=Easting.Est-Easting.SD, 
                  ymax=Easting.Est+Easting.SD,fill=Model),
              alpha=0.3) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  # scale_fill_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from Easting centroid (km)') +
  ggtitle('5% quantile range edge Easting - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(vars(Model))

# Spatial plots
re <- re %>% 
  mutate(Northing.Est=Northing.Est*1000,
         Northing.SD=Northing.SD*1000,
         Easting.Est=Easting.Est*1000,
         Easting.SD=Easting.SD*1000)

re.sf <- st_as_sf(re, coords=c('Easting.Est', 'Northing.Est'), 
                  crs="EPSG:32619")

re.lin <- re.sf %>% 
  #filter(Model=='All Covars',
  #       Quantile == '0.05') %>% 
  group_by(Model, Quantile) %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

spring.re <- ggplot() +
  geom_sf(data=coast, col=NA) +
  geom_sf(data=gom, fill=NA) +
  geom_sf(data=re.lin, aes(col=interaction(Quantile)),
         lwd=0.3) +
  # geom_sf(data=re.sf[re.sf$Year %in% c(1982,2021),],
  #        aes(col=Year)) +
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-76, -65),
           ylim=c(36, 45),
           crs="EPSG:4326") +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5),
        strip.text = element_text(size=8))

re <- reorig[reorig$Season == 'Fall',]

ggplot() +
  geom_point(data=re[re$Quantile == '0.05',], 
             aes(x=Year, y=Northing.Est, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Quantile == '0.05',], 
            aes(x=Year, y=Northing.Est, col=Model),
            lwd=1, alpha=0.7) +
  geom_ribbon(data=re[re$Quantile == '0.05',], 
              aes(x=Year, ymin=Northing.Est-Northing.SD, 
                  ymax=Northing.Est+Northing.SD,fill=Model),
              alpha=0.3) +
  #scale_color_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  #scale_fill_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('5% quantile range edge northing - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(vars(Model))

ggplot() +
  geom_point(data=re[re$Quantile == '0.05',], 
             aes(x=Year, y=Easting.Est, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Quantile == '0.05',], 
            aes(x=Year, y=Easting.Est, col=Model),
            lwd=1, alpha=0.7) +
  geom_ribbon(data=re[re$Quantile == '0.05',], 
              aes(x=Year, ymin=Easting.Est-Easting.SD, 
                  ymax=Easting.Est+Easting.SD,fill=Model),
              alpha=0.3) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  # scale_fill_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from Easting centroid (km)') +
  ggtitle('5% quantile range edge Easting - single covariate models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(vars(Model))

# Spatial plots
re <- re %>% 
  mutate(Northing.Est=Northing.Est*1000,
         Northing.SD=Northing.SD*1000,
         Easting.Est=Easting.Est*1000,
         Easting.SD=Easting.SD*1000)

re.sf <- st_as_sf(re, coords=c('Easting.Est', 'Northing.Est'), 
                  crs="EPSG:32619")

re.lin <- re.sf %>% 
  #filter(Model=='All Covars',
  #       Quantile == '0.05') %>% 
  group_by(Model, Quantile) %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')

fall.re <- ggplot() +
  geom_sf(data=coast, col=NA) +
  geom_sf(data=gom, fill=NA) +
  geom_sf(data=re.lin, aes(col=interaction(Quantile)),
          lwd=0.3) +
  # geom_sf(data=re.sf[re.sf$Year %in% c(1982,2021),],
  #        aes(col=Year)) +
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-76, -65),
           ylim=c(36, 45),
           crs="EPSG:4326") +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5),
        strip.text = element_text(size=8))

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/rangeshifts_spring.png')),
       spring.re, width = 5.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/medium/SingleCovs', 
                            '/rangeshifts_fall.png')),
       fall.re, width = 5.5, height=4.25, units='in')

rm(list=setdiff(ls(), c('coast', 'gom')))