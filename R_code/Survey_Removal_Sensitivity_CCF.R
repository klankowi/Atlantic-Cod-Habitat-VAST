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
          'DFO Trawl',
          'GSO Trawl',
          'MADMF Industry',
          'MADMF Trawl',
          'ME-NH Inshore Trawl',
          'RIDEM Trawl',
          'SMAST Video Trawl'
          )

short <- c('ASMFC', 
           'BLLS', 
           'BTS',
           'DFO', 
           'GSO', 
           'MADMF Industry',
           'MADMF Inshore', 
           'ME-NH', 
           'RIDEM',
           'SMAST')

short <- data.frame(cbind(mods, short))
colnames(short) <- c('Mod', 'Short')

for(i in 1:length(mods)){
  indno_mi <- read.csv(here('VAST_runs/sensitivity/medium/', mods[i], 
                            '/medium_Index.csv'))
  indno_mi$Mod <- mods[i]
  indno_mi <- indno_mi %>% 
    dplyr::select(-Category, -Std.Err.ln.Est)
  
  indall <- rbind(indall, indno_mi)
}

ind <- indall %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  mutate(Mod = factor(Mod, levels=c('All', 
                                    'ASMFC Shrimp Trawl', 
                                    'DFO Trawl',
                                    'GSO Trawl', 
                                    'MADMF Trawl',
                                    'MADMF Industry', 
                                    'ME-NH Inshore Trawl',
                                    'NEFSC BLLS', 
                                    'NEFSC BTS', 
                                    'RIDEM Trawl',
                                    'SMAST Video Trawl')),
         Size = as.factor('Medium'),
         Stratum = factor(Stratum, levels=c('ALL', 'EGOM', 'GBK', 
                                            'SNE', 'WGOM'))) %>% 
  mutate(Estimate = Estimate / 1000000,
         Std.Err.Est = Std.Err.Est / 1000000) %>% 
  mutate(High = Estimate + Std.Err.Est,
         Low = Estimate - Std.Err.Est) %>% 
  as.data.frame()

ind <- merge(ind, short, by=c('Mod'), all=T)

summary(ind)

ggplot(data=ind[ind$Mod != 'NEFSC BTS' & 
                ind$Mod != 'All',]) +
  geom_line(aes(x=Year, y=Estimate, col=Short)) +
  
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

lag0 <- merge(lag0, short, by=c("Mod"), all=T)

confint.poly <- data.frame(
  x=c(-1, -1, 6, 6),
  y=c(lag0$lowerCI[1], lag0$upperCI[1],
      lag0$upperCI[1], lag0$lowerCI[1])
)

lag0 <- lag0[lag0$Stock != 'ALL' &
             lag0$Mod != 'NEFSC BTS',]
lag0$Stock <- droplevels(lag0$Stock)

lagplot <- ggplot(data=lag0) +
  geom_segment(aes(x=as.numeric(Stock),
                   xend=as.numeric(Stock),
                   y=0,
                   yend=acf,
                   col=Stock),
               lwd=2) +
  ggh4x::facet_grid2(Season ~ Short) +
  geom_hline(yintercept = 0.8, col='red', lty=2, lwd=0.3) +
  
  geom_polygon(data=confint.poly,
               aes(x=x, y=y), fill='black',
               alpha=0.1) +
  labs(x='Stock', y='Cross-Correlation at Lag 0') +
  coord_cartesian(ylim=c(-0.5, 1.05),
                  xlim=c(0.5, 4.5)) +
  theme(strip.text.x = element_text(size=7),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
lagplot

ggsave(plot=lagplot,
       here('VAST_runs/sensitivity/medium/lag0_corr.png'),
       width = 8.5, height = 3, units='in', dpi=600)

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

rats <- merge(rats, short, by=c('Mod'), all=T)

all <- ind %>% 
  filter(Mod == 'All') %>% 
  rename(Stock = Stratum) %>% 
  mutate(Short = as.factor('All'),
         Upper = Estimate + Std.Err.Est,
         Lower = Estimate - Std.Err.Est) %>% 
  dplyr::select(-NAME) %>% 
  mutate(ratio = 1,
         AbsDif = 0,
         Mag = 0,
         TimeOut = 0) %>% 
  dplyr::select(Mod, Season, Stock, ratio, Year, Upper,
                Lower, Estimate, AbsDif, Mag, TimeOut, Short)

rats <- rbind(all, rats)

ratgroup <- rats %>% 
  group_by(Season, Stock, Short) %>% 
  summarise(Rat = mean(ratio),
            Mag = sum(Mag[TimeOut ==1]),
            AbsDif = sum(AbsDif[TimeOut == 1]),
            TimeOut = sum(TimeOut)) %>% 
  as.data.frame() %>% 
  mutate(Size = 'Medium')

write.csv(ratgroup, here('VAST_runs/sensitivity/medium/medium_Ratios.csv'),
          row.names = F)

rats <- rats[rats$Stock != 'ALL',]

for(i in 1:nrow(short)){
  use <- short$Short[i]
  print(
  ggplot() +
    geom_ribbon(data=rats[rats$Mod == 'All',],
                aes(x=Year, ymin=Lower, ymax=Upper),
                alpha=0.4) +
    geom_line(data=rats[rats$Short == paste0(use) &
                          rats$Stock != 'ALL',],
              aes(x=Year, y=Estimate), col='red', lwd=0.2) +
    geom_point(data=rats[rats$Short == paste0(use) &
                           rats$TimeOut ==1,],
               aes(x=Year, y=Estimate), col='red', cex=0.4) +
    # geom_polygon(data=survuse.poly4,
    #              aes(x=x, y=y), fill='black',
    #              alpha=0.1) +
    ggh4x::facet_grid2(Stock ~ Season, scales = 'free') +
    ggtitle(paste0(use, ' survey removed'))
  )
    
}

ggplot() +
  geom_ribbon(data=rats[rats$Mod == 'All',],
              aes(x=Year, ymin=Lower, ymax=Upper),
              alpha=0.4) +
  geom_line(data=rats[rats$Short == 'DFO',],
            aes(x=Year, y=Estimate), col='red', lwd=0.2) +
  geom_point(data=rats[rats$Short == 'DFO' &
                       rats$TimeOut ==1,],
             aes(x=Year, y=Estimate), col='red', cex=0.4) +
  ggh4x::facet_grid2(Stock ~ Season, scales = 'free')


ggplot(data=ratgroup[ratgroup$Short != 'BTS' & ratgroup$TimeOut>2,]) +
  geom_point(aes(x=as.numeric(as.factor(Stock)), 
                  y=Mag, col=Stock)) +
  ggh4x::facet_grid2(Season ~ Short) +
  #ylim(0, 2) + 
  labs(x='') +
  scale_x_continuous(limits=c(0.5, 5.5),
                     breaks=0) +
  theme(axis.text.x = element_text(size=7)) +
  labs(y='Mean Ratio')

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
  facet_wrap(vars(Short)) +
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


