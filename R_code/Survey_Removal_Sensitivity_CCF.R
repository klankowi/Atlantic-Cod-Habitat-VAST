rm(list=ls())

library(tidyverse)
library(here)
library(sf)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, 
                                            fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

indall <- read.csv(here('VAST_runs/medium/Overall_BC/ALL/Index.csv'))

indall$Mod <- 'All'

indall <- indall %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  rename(Std.Err.Est = Std..Error.for.Estimate) %>% 
  dplyr::select(-Category, -Units, -Std..Error.for.ln.Estimate.)

mods <- c('ASMFC Shrimp Trawl', 
          'NEFSC BLLS',
          'NEFSC BTS', 
          #'DFO Trawl',
          #'GSO Trawl',
          #'MADMF Industry',
          #'MADMF Trawl',
          #'ME-NH Inshore Trawl',
          #'RIDEM Trawl',
          'SMAST Video Trawl'
          )

for(i in 1:length(mods)){
  indno_mi <- read.csv(here('VAST_runs/sensitivity/medium/', mods[i], 
                            '/Medium_Index.csv'))
  indno_mi$Mod <- mods[i]
  indno_mi <- indno_mi %>% 
    dplyr::select(-Category, -Std.Err.ln.Est)
  
  indall <- rbind(indall, indno_mi)
}

ind <- indall %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  mutate(Mod = factor(Mod, levels=c('All', 'ASMFC Shrimp Trawl', #'DFO Trawl',
                                    #'GSO Trawl', 'MADMF Trawl',
                                    #'MADMF Industry', 'ME-NH Inshore Trawl',
                                    'NEFSC BLLS', 'NEFSC BTS', #'RIDEM Trawl',
                                    'SMAST Video Trawl')),
         Size = as.factor('Medium'),
         Stratum = factor(Stratum, levels=c('ALL', 'EGOM', 'GBK', 
                                            'SNE', 'WGOM'))) %>% 
  mutate(Estimate = Estimate / 1000000,
         Std.Err.Est = Std.Err.Est / 1000000) %>% 
  mutate(High = Estimate + Std.Err.Est,
         Low = Estimate - Std.Err.Est) %>% 
  as.data.frame()
summary(ind)

ggplot(data=ind[ind$Mod != 'NEFSC BTS' & 
                ind$Mod != 'All',]) +
  geom_line(aes(x=Year, y=Estimate, col=Mod)) +
  
  geom_line(data=ind[ind$Mod == 'All',], aes(x=Year, y=Estimate),
            col='black') +
  
  
  facet_wrap(vars(Season, Stratum), scales='free_y', nrow = 2) +
  theme(legend.position = 'right') +
  labs(x='Year', y='Index', col='Model Removed')

ind$NAME <- paste0(indall$Season, '_', indall$Stratum)

indd.list <- split(ind, f=ind$NAME)
for(i in 1:length(indd.list)){
  indd.list[[i]] <- split(indd.list[[i]], indd.list[[i]]$Mod)
}

# Convert to time series
indd.ts <- indd.list
for(i in 1:length(indd.list)){
  for(j in 1:length(indd.list[[i]])){
    indd.ts[[i]][[j]] <- ts(indd.list[[i]][[j]]$Estimate)
  }
}

ct <- indd.ts

for(i in 1:length(indd.ts)){
  
  for(j in 2:length(indd.ts[[i]])){
    # print(
    # ggCcf(indd.ts[[i]][[1]], indd.ts[[i]][[j]], lag.max = 5, type='correlation') +
    #   ggtitle(paste0(names(indd.ts)[i], ' ', names(indd.ts[[i]][j]),
    #           ' compared to model with all surveys')) +
    #   ylim(-1, 1)
    # )
    
    ct[[i]][[j]] <- ccf(indd.ts[[i]][[1]], indd.ts[[i]][[j]], lag.max=5, 
                        type='correlation', plot=F,
                        main = paste0(names(indd.ts)[i], ' ', 
                                      names(indd.ts[[i]][j])))
    ct[[i]][[j]] <- data.frame(
      lag = ct[[i]][[j]]$lag,
      acf = ct[[i]][[j]]$acf,
      NAME = names(ct)[i],
      Mod = names(ct[[i]][j]),
      upperCI = qnorm((1 + 0.95)/2)/sqrt(40),
      lowerCI = -qnorm((1 + 0.95)/2)/sqrt(40)
      
    )
    
  }
}

for(i in 1:length(ct)){
  ct[[i]] <- ct[[i]][-1]
}

for(i in 1:length(ct)){
  ct[[i]] <- do.call(rbind, ct[[i]])
}
ct <- do.call(rbind, ct)

rownames(ct) <- NULL

lag0 <- ct[ct$lag == 0,]

lag0 <- lag0 %>% 
  separate(NAME, sep='_', into=c('Season', 'Stock')) %>% 
  mutate(Stock = factor(Stock, levels=c('ALL', 'EGOM', 'GBK','SNE','WGOM')),
         Season = factor(Season, levels=c('Spring', 'Fall')))

lagplot <- ggplot(data=lag0) +
  geom_point(aes(x=as.factor(lag), y=acf, col=Stock),
             cex=2, alpha=0.8) +
  ggh4x::facet_grid2(Season ~ Mod) +
  geom_hline(yintercept = 0.9, col='red', lty=2, lwd=0.3) +
  geom_hline(yintercept = lag0$upperCI[1], col='blue', lty=2, lwd=0.3) +
  labs(x='Lag', y='Correlation') +
  ylim(0, 1)+
  theme(strip.text.x = element_text(size=5))
lagplot

# ggsave(plot=lagplot,
#        here('VAST_runs/sensitivity/medium/lag0_corr.png'),
#        width = 8.5, height = 3, units='in') 

rm(indall, indno_mi, indd.ts, lag0, lagplot, i, j, mods, ct)

# General
  # Fully expect removal of BTS to fuck all these correlations, that's fine

# Smalls:
  # DFO Trawl removal: Spring GBK, Fall GBK, Fall SNE, Fall EGOM
  # MADMF Industry: Spring EGOM
  # MADMF Inshore: Spring All, Spring SNE, Spring WGOM, Fall All, Fall SNE, Fall WGOM
  # ME-NH: Spring EGOM, Fall EGOM
  # RIDEM: Fall SNE
  # SMAST: Fall SNE, Spring WGOM

# Simple ratios
rats <- indd.list

for(i in 1:length(indd.list)){
  
  for(j in 2:length(indd.list[[i]])){

    rats[[i]][[j]] <- data.frame(
      NAME = names(rats)[i],
      Mod = names(rats[[i]][j]),
      ratio = rats[[i]][[j]]$Estimate / rats[[i]][[1]]$Estimate,
      Year = rats[[i]][[j]]$Year,
      Upper = rats[[i]][[1]]$Estimate + rats[[i]][[1]]$Std.Err.Est,
      Lower = rats[[i]][[1]]$Estimate - rats[[i]][[1]]$Std.Err.Est,
      Estimate = rats[[i]][[j]]$Estimate,
      AbsDif = abs(rats[[i]][[j]]$Estimate - rats[[i]][[1]]$Estimate),
      Mag = abs(rats[[i]][[j]]$Estimate - rats[[i]][[1]]$Estimate)/
        rats[[i]][[1]]$Estimate
    )
  }
}

for(i in 1:length(rats)){
  rats[[i]] <- rats[[i]][-1]
}

for(i in 1:length(rats)){
  rats[[i]] <- do.call(rbind, rats[[i]])
}
rats <- do.call(rbind, rats)

rownames(rats) <- NULL

rats$TimeOut <- 1

for(i in 1:nrow(rats)){
  if(rats$Estimate[i] > rats$Lower[i] &
     rats$Estimate[i] < rats$Upper[i]){
    rats$TimeOut[i] <- 0
  }
}

rats <- rats %>% 
  separate(NAME, sep='_', into=c('Season', 'Stock')) %>% 
  mutate(Stock = factor(Stock, levels=c('ALL', 'EGOM', 'GBK','SNE','WGOM')),
         Season = factor(Season, levels=c('Spring', 'Fall')))

ratgroup <- rats %>% 
  group_by(Season, Stock, Mod) %>% 
  summarise(Mag = sum(Mag[TimeOut ==1]),
            AbsDif = sum(AbsDif[TimeOut == 1]),
            TimeOut = sum(TimeOut)) %>% 
  as.data.frame()

ggplot() +
  geom_ribbon(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                  rats$Stock == 'ALL',],
            aes(x=Year, ymin=Lower,
                ymax= Upper),
            alpha=0.5)  +
  geom_line(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                       rats$Stock == 'ALL',],
           aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                        rats$Stock == 'ALL' & rats$TimeOut == 1,],
            aes(x=Year, y=Estimate), col='red', cex=0.4) +
  facet_wrap(vars(Mod)) +
  ggtitle('Spring, all biological stock areas')

ggplot() +
  geom_ribbon(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                          rats$Stock == 'WGOM',],
              aes(x=Year, ymin=Lower,
                  ymax= Upper),
              alpha=0.5)  +
  geom_line(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                        rats$Stock == 'WGOM',],
            aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                         rats$Stock == 'WGOM' & rats$TimeOut == 1,],
             aes(x=Year, y=Estimate), col='red', cex=0.4) +
  facet_wrap(vars(Mod)) +
  ggtitle('Spring, WGOM biological stock area')

ggplot() +
  geom_ribbon(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                          rats$Stock == 'EGOM',],
              aes(x=Year, ymin=Lower,
                  ymax= Upper),
              alpha=0.5)  +
  geom_line(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                        rats$Stock == 'EGOM',],
            aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                         rats$Stock == 'EGOM' & rats$TimeOut == 1,],
             aes(x=Year, y=Estimate), col='red', cex=0.4) +
  facet_wrap(vars(Mod)) +
  ggtitle('Spring, EGOM biological stock area')

ggplot() +
  geom_ribbon(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                          rats$Stock == 'SNE',],
              aes(x=Year, ymin=Lower,
                  ymax= Upper),
              alpha=0.5)  +
  geom_line(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                        rats$Stock == 'SNE',],
            aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                         rats$Stock == 'SNE' & rats$TimeOut == 1,],
             aes(x=Year, y=Estimate), col='red', cex=0.4) +
  facet_wrap(vars(Mod)) +
  ggtitle('Spring, SNE biological stock area')


ggplot() +
  geom_ribbon(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                          rats$Stock == 'GBK',],
              aes(x=Year, ymin=Lower,
                  ymax= Upper),
              alpha=0.5)  +
  geom_line(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                        rats$Stock == 'GBK',],
            aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Season == 'Spring' & rats$Mod != 'NEFSC BTS' &
                         rats$Stock == 'GBK' & rats$TimeOut == 1,],
             aes(x=Year, y=Estimate), col='red', cex=0.4) +
  facet_wrap(vars(Mod)) +
  ggtitle('Spring, GBK biological stock area')

