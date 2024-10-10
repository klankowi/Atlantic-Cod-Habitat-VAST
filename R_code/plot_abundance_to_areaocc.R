rm(list=ls())

library(VAST)
library(tidyverse)
library(here)

sizes <- c('small', 'medium', 'large')

for(i in 1:length(sizes)){
  setwd(here(paste0('VAST_runs/', sizes[i], '/Overall_BC/ALL/')))
  load(paste0('Overall_BC_', sizes[i], 'cod_allstrat_natsplin_fsON_ALL.RData'))
  assign(paste0(sizes[i]), fit)
  rm(covars, extrap_info_aja, scaled.covars, settings, strata_use, surveys,
     survs, vast_extrap_grid, hab_formula, seas.labs, working_dir, year.labs)
  
  avgD <- as.data.frame(fit$Report$mean_D_ctl[,,])
  colnames(avgD) <- fit$settings$strata.limits$STRATA
  
  aroc <- read.csv(paste0('AreaOcc.csv'))
  aroc <- aroc %>% 
    rename(Strata = strata,
           std.err.ao = std.err) %>% 
    mutate(Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
    dplyr::select(-units)
  
  index <- read.csv(paste0('Index.csv'))
  index <- index %>% 
    dplyr::select(-Category, -Units, -Std..Error.for.ln.Estimate.) %>% 
    rename(Strata = Stratum,
           std.err = Std..Error.for.Estimate,
           estimate = Estimate) %>% 
    separate(Time, into=c('Year', 'Season')) %>% 
    mutate(Year = as.numeric(Year),
           Season = factor(Season, levels=c('Spring', 'Fall')))
  
  avgD <- avgD %>% 
    mutate(Time = rownames(avgD)) %>% 
    separate(Time, into=c('Year', 'Season')) %>% 
    mutate(Year = as.numeric(Year),
           Season = factor(Season, levels=c('Spring', 'Fall'))) %>%
    pivot_longer(c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')) %>% 
    rename(Strata = name) %>% 
    as.data.frame()
  
  avgD <- avgD %>% 
    rename(Avg.Ind.KM = value) %>% 
    mutate(Avg.Ind.KM = strip_units(Avg.Ind.KM))
  
  avgD <- left_join(avgD, aroc, by=c('Year', 'Season', 'Strata'))
  avgD <- left_join(avgD, index, by=c('Year', 'Season', 'Strata'))
  avgD$Size = paste0(sizes[i])
  
  assign(paste0('avgD_', sizes[i]), avgD)
  
  rm(list=setdiff(ls(), c('avgD_small', 'avgD_medium', 'avgD_large',
                          'sizes', 'i')))
}

avgD <- rbind(
  avgD_small, avgD_medium, avgD_large
)

rm(list=setdiff(ls(), c('avgD')))

ggplot(data=avgD[avgD$Size == 'large',]) +
  geom_point(aes(x=log(estimate), y=log(area.occ), col=Strata),
             cex=0.7) +
  
  geom_smooth(aes(x=log(estimate), y=log(area.occ), col=Strata,
                  fill=Strata),
              alpha=0.3, method='lm') +
  
  facet_wrap(vars(Season)) +
  
  facet_wrap(vars(Season), scales='free') +
  #xlab('log(Abundance)') + ylab('log(Area)') +
  theme(legend.position = 'bottom')

modres <- avgD %>% 
  group_by(Size, Strata, Season) %>% 
  mutate(slope = summary(lm(log(area.occ) ~ log(estimate)))$coefficients[2, 1],
         pVal = summary(lm(log(area.occ) ~ log(estimate)))$coefficients[2, 4],
         adjR = summary(lm(log(area.occ) ~ log(estimate)))$adj.r.squared,
         sig = ' ') %>% 
  dplyr::select(Size, Strata, Season, slope, pVal, sig, adjR) %>% 
  unique() %>% 
  as.data.frame()
modres$sig[modres$pVal <= 0.05] <- '*'

modres
