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
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

cog.c <- read.csv(here('VAST_runs/large/GMRF/OFF/COG_ALL.csv'))
cog.c$Model <- 'No GMRF'
cog.f<- read.csv(here('VAST_runs/large/GMRF/ON/COG_ALL.csv'))
cog.f$Model <- 'GMRF'
cog <- rbind(cog.c, cog.f)

re.c <- read.csv(here('VAST_runs/large/GMRF/OFF/RangeEdges_ALL.csv'))
re.c$Model <- 'No GMRF'
re.f<- read.csv(here('VAST_runs/large/GMRF/ON/RangeEdges_ALL.csv'))
re.f$Model <- 'GMRF'
re <- rbind(re.c, re.f)

ind.c <- read.csv(here('VAST_runs/large/GMRF/OFF/Index.csv'))
ind.c$Model <- 'No GMRF'
ind.f<- read.csv(here('VAST_runs/large/GMRF/ON/Index.csv'))
ind.f$Model <- 'GMRF'
ind <- rbind(ind.c, ind.f)

arrocc.c <- read.csv(here('VAST_runs/large/GMRF/OFF/AreaOcc.csv'))
arrocc.c$Model <- 'No GMRF'
arrocc.f<- read.csv(here('VAST_runs/large/GMRF/ON/AreaOcc.csv'))
arrocc.f$Model <- 'GMRF'
arrocc <- rbind(arrocc.c, arrocc.f)

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

lm(northing ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'GMRF' & cog$Season == 'Spring',])
lm(northing ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'GMRF' & cog$Season == 'Fall',])

lm(northing ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'No GMRF' & cog$Season == 'Spring',])
lm(northing ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'No GMRF' & cog$Season == 'Fall',])


lm(easting ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'GMRF',])
lm(easting ~ Year, data=cog[cog$Year <=2021 & cog$Model == 'No GMRF',])

cog.northing <- ggplot() +
  geom_point(data=cog[cog$Year <=2022,],
             aes(x=Year, y=northing, col=Model),
             cex=1) +
  geom_line(data=cog[cog$Year <=2022,],
            aes(x=Year, y=northing, col=Model),
            lwd=1) +
  geom_ribbon(data=cog[cog$Year <=2022,],
              aes(x=Year, ymin = northing - n.sd,
                  ymax=northing + n.sd, fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season), nrow=2) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Northing, Center of Gravity')+
  theme(legend.position = 'right')

cog.easting <- ggplot() +
  geom_point(data=cog[cog$Year <=2022,],
            aes(x=Year, y=easting, col=Model),
            cex=1) +
  geom_line(data=cog[cog$Year <=2022,],
            aes(x=Year, y=easting, col=Model),
            lwd=1) +
  geom_ribbon(data=cog[cog$Year <=2022,],
              aes(x=Year, ymin = easting - e.sd,
                  ymax=easting + e.sd, fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season), nrow=2) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Easting, Center of Gravity') +
  theme(legend.position = 'right')

re$Quantile <- as.factor(re$Quantile)
lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.05' & re$Model == 'GMRF',])
lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.5'& re$Model == 'GMRF',])
lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.95' & re$Model == 'GMRF',])

lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.05' & re$Model == 'GMRF',])
lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.5'& re$Model == 'GMRF',])
lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.95' & re$Model == 'GMRF',])

lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.05' & re$Model == 'No GMRF',])
lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.5'& re$Model == 'No GMRF',])
lm(Northing.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.95' & re$Model == 'No GMRF',])

lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.05' & re$Model == 'No GMRF',])
lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.5'& re$Model == 'No GMRF',])
lm(Easting.Est ~ Year, data=re[re$Year <=2021 & re$Quantile == '0.95' & re$Model == 'No GMRF',])

re.northing <- ggplot() +
  geom_line(data=re[re$Year <=2022,],
             aes(x=Year, y=Northing.Est, 
                 group=interaction(Quantile, Model),
                 col=Model)) +
  geom_point(data=re[re$Year <=2022,],
            aes(x=Year, y=Northing.Est, 
                pch=Quantile, col=Model)) +
  geom_ribbon(data=re[re$Year <=2022,],
              aes(x=Year, ymin=Northing.Est - Northing.SD,
                  ymax=Northing.Est + Northing.SD,
                  group=interaction(Quantile, Model),
                  fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season)) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Northing')+
  theme(legend.position = 'bottom')


re.easting <- ggplot() +
  geom_line(data=re[re$Year <=2022,],
            aes(x=Year, y=Easting.Est, 
                group=interaction(Quantile, Model),
                col=Model)) +
  geom_point(data=re[re$Year <=2022,],
             aes(x=Year, y=Easting.Est, 
                 pch=Quantile, col=Model)) +
  geom_ribbon(data=re[re$Year <=2022,],
              aes(x=Year, ymin=Easting.Est - Easting.SD,
                  ymax=Easting.Est + Easting.SD,
                  group=interaction(Quantile, Model),
                  fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season)) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Easting')+
  theme(legend.position = 'bottom')

# ind$Stratum[ind$Stratum == 'Stratum_1'] <- 'Combined'
# ind$Stratum[ind$Stratum == 'Stratum_2'] <- 'US'
# ind$Stratum[ind$Stratum == 'Stratum_3'] <- 'Canada'
# 
# ind$Stratum <- factor(ind$Stratum,
#                       levels=c('Canada', 'Combined', 'US'))
ind <- ind %>% 
  separate(Time, c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year))

rel.abund <- ggplot() +
  geom_line(data=ind[ind$Year<=2022,],
            aes(x=Year, y=Estimate, col=Model)) +
  geom_ribbon(data=ind[ind$Year<=2022,],
              aes(x=Year, ymin=Estimate-`Std..Error.for.Estimate`,
                  ymax=Estimate + `Std..Error.for.Estimate`,
                  fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season))+
  ggtitle('Relative Abundance') +
  theme(legend.position='right')

lm(area.occ ~ Year, data=arrocc[arrocc$Year <=2021 & arrocc$Model == 'GMRF',])
lm(area.occ ~ Year, data=arrocc[arrocc$Year <=2021 & arrocc$Model == 'No GMRF',])

eff.arr.oc <- ggplot() +
  geom_line(data=arrocc[arrocc$Year<=2022,],
            aes(x=Year, y=area.occ, col=Model)) +
  geom_ribbon(data=arrocc[arrocc$Year<=2022,],
              aes(x=Year, 
                  ymin=area.occ-std.err,
                  ymax=area.occ+std.err, fill=Model),
              alpha=0.2) +
  facet_wrap(vars(Season)) +
  ggtitle('Effective Area Occupied') +
  ylab('Area Occupied (km)') +
  theme(legend.position = 'right')

# Remove intermediates
rm(list=setdiff(ls(), c('cog.easting', 'cog.northing',
                        're.easting', 're.northing', 
                        'rel.abund', 'eff.arr.oc')))
setwd(here('VAST_runs/large/GMRF/Testing Results'))
# Save plots
ggsave('cog_easting.png', cog.easting, width=11, height=8.5)
ggsave('cog_northing.png', cog.northing, width=11, height=8.5)

ggsave('re_easting.png', re.easting, width=11, height=8.5)
ggsave('re_northing.png', re.northing, width=11, height=8.5)

ggsave('rel_abund.png', rel.abund, width=11, height=8.5)

ggsave('area_occupied.png', eff.arr.oc, width=11, height=8.5)
