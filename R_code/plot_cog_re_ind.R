rm(list=ls())

library(VAST)
library(here)
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

# Data
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_larcod_allstrat_natsplin_fsON_ALL.RData'))
rm(list=setdiff(ls(), c('fit', '%notin%')))

cog.wgom <- read.csv(here("VAST_runs/large/Overall_BC/WGOM/COG_WGOM.csv"))
cog.egom <- read.csv(here('VAST_runs/large/Overall_BC/EGOM/COG_EGOM.csv'))
cog.gbk <- read.csv(here("VAST_runs/large/Overall_BC/GBK/COG_GBK.csv"))
cog.sne <- read.csv(here('VAST_runs/large/Overall_BC/SNE/COG_SNE.csv'))
cog.all <- read.csv(here("VAST_runs/large/Overall_BC/ALL/COG_ALL.csv"))

cog <- rbind(cog.wgom, cog.egom, cog.gbk, cog.sne, 
             cog.all)

re <- read.csv(here('VAST_runs/large/Overall_BC/ALL/RangeEdges_ALL.csv'))

ind <- read.csv(here('VAST_runs/large/Overall_BC/ALL/Index.csv'))
ind <- ind %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season,
                         levels=c('Spring', 'Fall')),
         Stratum = factor(Stratum, levels=c('ALL','EGOM','GBK','SNE','WGOM')))

arrocc <- read.csv(here('VAST_runs/large/Overall_BC/ALL/AreaOcc.csv'))

cog <- cog %>% 
  mutate(easting = easting / 1000,
         northing = northing / 1000,
         e.sd = e.sd / 1000,
         n.sd = n.sd / 1000) %>% 
  dplyr::select(-units)

arrocc <- arrocc %>% 
  mutate(area.occ = area.occ / 1000,
         std.err = std.err / 1000) %>% 
  dplyr::select(-units)

cogshifts <- cog %>% 
  group_by(strata, Season) %>% 
  mutate(movmod.n = lm(northing ~ Year)$coefficients[2],
         movmod.e = lm(easting ~ Year)$coefficients[2]) %>% 
  dplyr::select(strata, Season, movmod.n, movmod.e) %>% 
  unique() %>% 
  as.data.frame()

cogshifts$dist <- sqrt((cogshifts$movmod.n^2) +
                         (cogshifts$movmod.e^2))

cog.northing <- ggplot() +
  # geom_point(data=cog,
  #            aes(x=Year, y=northing),
  #            cex=0.6) +
  geom_line(data=cog[cog$strata == 'ALL',],
            aes(x=Year, y=northing),
            lwd=0.5) +
  geom_ribbon(data=cog[cog$strata == 'ALL',],
              aes(x=Year, ymin = northing - n.sd,
                  ymax=northing + n.sd),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Northing, Center of Gravity')+
  theme(legend.position='bottom')

cog.easting <- ggplot() +
  # geom_point(data=cog,
  #            aes(x=Year, y=easting),
  #            cex=0.6) +
  geom_line(data=cog[cog$strata == 'ALL',],
            aes(x=Year, y=easting),
            lwd=0.5) +
  geom_ribbon(data=cog[cog$strata == 'ALL',],
              aes(x=Year, ymin = easting - e.sd,
                  ymax=easting + e.sd),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Easting, Center of Gravity')+
  theme(legend.position='bottom')

re$Quantile <- as.factor(re$Quantile)

reshifts <- re %>% 
  group_by(Quantile, Season) %>% 
  mutate(movmod.n = lm(Northing.Est ~ Year)$coefficients[2],
         movmod.e = lm(Easting.Est ~ Year)$coefficients[2]) %>% 
  dplyr::select(Quantile, Season, movmod.n, movmod.e) %>% 
  unique() %>% 
  as.data.frame()

reshifts$dist <- sqrt((reshifts$movmod.n^2) +
                         (reshifts$movmod.e^2))

re.northing <- ggplot() +
  geom_line(data=re,
             aes(x=Year, y=Northing.Est, 
                 group=Quantile,
                 color=Quantile),
            lwd=0.5) +
  # geom_point(data=re,
  #           aes(x=Year, y=Northing.Est, 
  #               color=Quantile),
  #           cex=0.6) +
  geom_ribbon(data=re,
              aes(x=Year, ymin=Northing.Est - Northing.SD,
                  ymax=Northing.Est + Northing.SD,
                  fill=Quantile),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Range Edges, Northing')+
  theme(legend.position = 'bottom')

re.easting <- ggplot() +
  geom_line(data=re,
            aes(x=Year, y=Easting.Est, 
                group=Quantile,
                color=Quantile),
            lwd=0.5) +
  # geom_point(data=re,
  #            aes(x=Year, y=Easting.Est, 
  #                color=Quantile),
  #            cex=0.6) +
  geom_ribbon(data=re,
              aes(x=Year, ymin=Easting.Est - Easting.SD,
                  ymax=Easting.Est + Easting.SD,
                  fill=Quantile),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Range Edges, Easting')+
  theme(legend.position = 'bottom')

ind$Stratum <- factor(ind$Stratum,
                      levels=c('EGOM', 'GBK', 'SNE', 'WGOM', 'ALL'))

rel.abund <- ggplot() +
  geom_line(data=ind[ind$Stratum != 'ALL',],
            aes(x=Year, y=Estimate, col=Stratum)) +
  geom_ribbon(data=ind[ind$Stratum != 'ALL',],
              aes(x=Year, ymin=Estimate-Std.Err.Est,
                  ymax=Estimate + Std.Err.Est,
                  fill=Stratum),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  ggtitle('Relative Abundance') +
  theme(legend.position='bottom')

arrocc$strata <- factor(arrocc$strata,
                        levels = c('EGOM', 'GBK', 'SNE', 'WGOM', 'ALL'))
aoshifts <- arrocc %>% 
  group_by(strata, Season) %>% 
  mutate(ao = lm(area.occ ~ Year)$coefficients[2]) %>% 
  dplyr::select(strata, Season, ao) %>% 
  unique() %>% 
  as.data.frame()

eff.arr.oc <- ggplot() +
  geom_line(data=arrocc[arrocc$strata != 'ALL',],
            aes(x=Year, y=area.occ, col=strata)) +
  geom_ribbon(data=arrocc[arrocc$strata != 'ALL',],
              aes(x=Year, 
                  ymin=area.occ-std.err,
                  ymax=area.occ+std.err,
                  fill=strata),
              alpha=0.2) +
  facet_wrap(vars(Season), ncol=2) +
  ggtitle('Effective Area Occupied') +
  ylab('Area Occupied (km)') +
  theme(legend.position = 'bottom')

# Remove intermediates
rm(list=setdiff(ls(), c('cog.easting', 'cog.northing',
                        're.easting', 're.northing', 
                        'rel.abund', 'eff.arr.oc')))
setwd(here('VAST_runs/large/Overall_BC'))
# Save plots
ggsave('cog_easting.png', cog.easting, width=8, height=4)
ggsave('cog_northing.png', cog.northing, width=8, height=4)

ggsave('re_easting.png', re.easting, width=8, height=4)
ggsave('re_northing.png', re.northing, width=8, height=4)

ggsave('rel_abund.png', rel.abund, width=8, height=4)

ggsave('area_occupied.png', eff.arr.oc, width=8, height=4)
