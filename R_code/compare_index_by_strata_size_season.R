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
small.ind <- read.csv(here("VAST_runs/small/Overall_BC/with_AMO/All/Index.csv"))
load(here("VAST_runs/small/Overall_BC/with_AMO/All/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.RData"))
small.strat <- fit$settings$strata.limits
small.strat$Stratum <- unique(small.ind$Stratum)
small.ind <- left_join(small.ind, small.strat, by=c('Stratum'))
small.ind$Size <- 'Small'
rm(list=setdiff(ls(), c('small.ind')))

medium.ind <- read.csv(here("VAST_runs/medium/Overall_BC/All/Index.csv"))
load(here("VAST_runs/medium/Overall_BC/All/Overall_BC_medcod_allstrat_natsplin_fsON_ALL.RData"))
medium.strat <- fit$settings$strata.limits
medium.strat$Stratum <- unique(medium.ind$Stratum)
medium.ind <- left_join(medium.ind, medium.strat, by=c('Stratum'))
medium.ind$Size <- 'Medium'
rm(list=setdiff(ls(), c('medium.ind', 'small.ind')))

large.ind <- read.csv(here("VAST_runs/large/Overall_BC/with_cobble/All/Index.csv"))
load(here("VAST_runs/large/Overall_BC/with_cobble/All/Overall_BC_larcod_allstrat_natsplin_fsON_ALL.RData"))
large.strat <- fit$settings$strata.limits
large.strat$Stratum <- unique(large.ind$Stratum)
large.ind <- left_join(large.ind, large.strat, by=c('Stratum'))
large.ind$Size <- 'Large'
rm(list=setdiff(ls(), c('large.ind', 'medium.ind', 'small.ind')))

ind <- rbind(small.ind, medium.ind, large.ind)

ind <- ind %>% 
  dplyr::select(-Category, -Stratum, -Units) %>% 
  rename(Stratum = STRATA,
         Std.Err.Est = Std..Error.for.Estimate,
         Std.Err.Ln = Std..Error.for.ln.Estimate.) %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Stratum = factor(Stratum,
                          levels=c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')),
         Size = factor(Size,
                       levels=c('Small','Medium', 'Large')))

change <- ind %>% 
  group_by(Size, Season, Stratum) %>% 
  mutate(change = lm(Estimate ~ Year)$coefficients[2]) %>% 
  dplyr::select(Size, Season, Stratum, change) %>% 
  unique() %>% 
  as.data.frame()

strat <- change[change$Stratum != 'ALL',]

strat <- strat %>% 
  group_by(Size, Season) %>% 
  mutate(Sum = sum(change)) %>% 
  dplyr::select(Size, Season, Sum) %>% 
  unique() %>% 
  as.data.frame()

indplot <- ggplot() +
  geom_line(data=ind,
              aes(x=Year, y=Estimate, col=Size)) +
  geom_ribbon(data=ind,
              aes(x=Year, ymin=Estimate-Std.Err.Est,
                  ymax=Estimate+Std.Err.Est,
                  fill=Size),
              alpha=0.3) +
  facet_grid(Season ~ Stratum,
             scales='free')
