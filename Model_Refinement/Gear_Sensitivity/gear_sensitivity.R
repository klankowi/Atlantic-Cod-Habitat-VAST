# set up workspace
rm(list=ls())

# Load libraries
library(here)
library(tidyverse)

# Set GGplot auto theme
theme_set(theme(plot.margin = margin(t=0.25, b=0.25, l=0.5, r=0.25, 'cm'),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                legend.background = element_blank(),
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12),
                strip.text = element_text(size=10)))

# Add surveys
surveys <- read.csv(here('data/Dataframes/Survey_Data.csv'))

# What year did BLLS start
min(surveys$YEAR[surveys$SURVEY == 'NEFSC BLLS'])

# What year did Sentinel start
min(surveys$YEAR[surveys$SURVEY == 'Sentinel'])

# Read in VAST-generated index data
noSent <- read.csv(here('VAST_runs/StrataDensCats_7_noSentinel/Index.csv'))
noBLLS <- read.csv(here('VAST_runs/StrataDensCats_7_noBLLS/Index.csv'))
allSur <- read.csv(here('VAST_runs/StrataDensCats_7/Index.csv'))
noJigs <- read.csv(here('VAST_runs/StrataDensCats_7_allTrawl/Index.csv'))

# Combine to list
inlist <- list(noSent, noBLLS, noJigs, allSur)
names(inlist) <- c('Trawls+BLLS,\nno Sentinel', 'Trawls+Sentinel,\nno BLLS','Only Trawls', 'All Surveys')

# Adjust columns
for(i in 1:length(inlist)){
  # Name runs
  inlist[[i]]$Run <- names(inlist)[i]
  
  # Calculate upper CI
  inlist[[i]]$Upper <- 
    inlist[[i]]$Estimate +
    inlist[[i]]$Std..Error.for.Estimate
  
  # Calculate lower CI
  inlist[[i]]$Lower <- 
    inlist[[i]]$Estimate -
    inlist[[i]]$Std..Error.for.Estimate
}

# Recombine to df
ests <- do.call(rbind, inlist)
row.names(ests) <- NULL

# Set year and season
for(i in 1:nrow(ests)){
  ests$Year[i] <- strsplit(ests$Time[i], split = " ")[[1]][1]
  ests$Season[i] <- strsplit(ests$Time[i], split= " ")[[1]][2]
}

# Remove unneeded columns
ests <- dplyr::select(ests,
                      Run, Category, Year, Season,
                      Estimate, Upper, Lower)

# Rename age groups to size-based groups
ests$Category[ests$Category == 'Ages [0-2)'] <- 'Small (< 40cm)'
ests$Category[ests$Category == 'Ages [2-5)'] <- 'Medium (> 40cm & < 70cm)'
ests$Category[ests$Category == 'Ages [5+]'] <- 'Large (> 70cm)'
ests$Category[ests$Category == 'Unknown Ages'] <- 'Unknown sizes'

# Set levels of size group
ests$Category <- factor(ests$Category, 
                        levels=c('Small (< 40cm)',
                                 'Medium (> 40cm & < 70cm)',
                                 'Large (> 70cm)',
                                 'Unknown sizes'))
# set levels of run
ests$Run <- factor(ests$Run,
                   levels=c('All Surveys',
                            'Only Trawls',
                            'Trawls+Sentinel,\nno BLLS',
                            'Trawls+BLLS,\nno Sentinel'))

# Set levels of season
ests$Season <- factor(ests$Season,
                      levels=c('Spring', 'Fall'))

# Adjust year to numeric
ests$Year <- as.numeric(ests$Year)

# Fall index
fa <- ggplot() +
  geom_line(data=ests[ests$Season == 'Fall',], 
            aes(x=Year, y=Estimate, col=Run)) +
  geom_ribbon(data=ests[ests$Season == "Fall",],
              aes(x=Year, ymin=Lower, ymax=Upper,
                  fill=Run), 
              alpha=0.2, linetype=0) +
  facet_wrap(vars(Category), 
             scales="free_y") +
  xlab('Year') + ylab('Abundance') +
  theme(plot.margin = margin(t=0.25, b=0.25, l=0.5, r=0.25, 'cm')) +
  ggtitle('Fall index')
ggsave(plot = fa,
       path = 'C:/Users/klankowicz/Desktop//',
       filename = 'Fall_Abundance_GearDrop.png',
       device='png')

# Spring index
sa <- ggplot() +
  geom_line(data=ests[ests$Season == 'Spring',], 
            aes(x=Year, y=Estimate, col=Run)) +
  geom_ribbon(data=ests[ests$Season == "Spring",],
              aes(x=Year, ymin=Lower, ymax=Upper,
                  fill=Run), 
              alpha=0.2, linetype=0) +
  facet_wrap(vars(Category),
             scales='free_y') +
  xlab('Year') + ylab('Abundance') +
  theme(plot.margin = margin(t=0.25, b=0.25, l=0.5, r=0.25, 'cm')) +
  ggtitle('Spring index')
ggsave(plot = sa,
       path = 'C:/Users/klankowicz/Desktop/',
       filename = 'Spring_Abundance_GearDrop.png',
       device='png')


# Proportional comparison
for(i in 1:nrow(allSur)){
  noBLLS$Dif[i] <- noBLLS$Estimate[i] / allSur$Estimate[i]
  noBLLS$Upper[i] <- noBLLS$Estimate[i] + noBLLS$Std..Error.for.Estimate[i]
  noBLLS$UpDif[i] <- noBLLS$Upper[i] / allSur$Estimate[i]
  noBLLS$Lower[i] <- noBLLS$Estimate[i] - noBLLS$Std..Error.for.Estimate[i]
  noBLLS$LoDif[i] <- noBLLS$Lower[i] / allSur$Estimate[i]
  noBLLS$Run[i] <- 'Trawls+Sentinel,\nno BLLS'
  
  noSent$Dif[i] <- noSent$Estimate[i] / allSur$Estimate[i]
  noSent$Upper[i] <- noSent$Estimate[i] + noSent$Std..Error.for.Estimate[i]
  noSent$UpDif[i] <- noSent$Upper[i] / allSur$Estimate[i]
  noSent$Lower[i] <- noSent$Estimate[i] - noSent$Std..Error.for.Estimate[i]
  noSent$LoDif[i] <- noSent$Lower[i] / allSur$Estimate[i]
  noSent$Run[i] <- 'Trawls+BLLS,\nno Sentinel'
  
  noJigs$Dif[i] <- noJigs$Estimate[i] / allSur$Estimate[i]
  noJigs$Upper[i] <- noJigs$Estimate[i] + noJigs$Std..Error.for.Estimate[i]
  noJigs$UpDif[i] <- noJigs$Upper[i] / allSur$Estimate[i]
  noJigs$Lower[i] <- noJigs$Estimate[i] - noJigs$Std..Error.for.Estimate[i]
  noJigs$LoDif[i] <- noJigs$Lower[i] / allSur$Estimate[i]
  noJigs$Run[i] <- 'Only Trawls'
  
  allSur$Dif[i] <- 1
  allSur$Upper[i] <- allSur$Estimate[i] + allSur$Std..Error.for.Estimate[i]
  allSur$UpDif[i] <- allSur$Upper[i] / allSur$Estimate[i]
  allSur$Lower[i] <- allSur$Estimate[i] - allSur$Std..Error.for.Estimate[i]
  allSur$LoDif[i] <- allSur$Lower[i] / allSur$Estimate[i]
  allSur$Run[i] <- 'All Surveys'
}

# Combine to list
difdf <- rbind(noBLLS, noSent, noJigs, allSur)
row.names(difdf) <- NULL

# Set year and season
for(i in 1:nrow(difdf)){
  difdf$Year[i] <- strsplit(difdf$Time[i], split = " ")[[1]][1]
  difdf$Season[i] <- strsplit(difdf$Time[i], split= " ")[[1]][2]
}

# Remove unneeded columns
difdf <- dplyr::select(difdf,
                      Run, Category, Year, Season,
                      Estimate, Upper, Lower,
                      Dif, UpDif, LoDif)

# Rename age groups to size-based groups
difdf$Category[difdf$Category == 'Ages [0-2)'] <- 'Small (< 40cm)'
difdf$Category[difdf$Category == 'Ages [2-5)'] <- 'Medium (> 40cm & < 70cm)'
difdf$Category[difdf$Category == 'Ages [5+]'] <- 'Large (> 70cm)'
difdf$Category[difdf$Category == 'Unknown Ages'] <- 'Unknown sizes'

# Set levels of size group
difdf$Category <- factor(difdf$Category, 
                        levels=c('Small (< 40cm)',
                                 'Medium (> 40cm & < 70cm)',
                                 'Large (> 70cm)',
                                 'Unknown sizes'))
# set levels of run
difdf$Run <- factor(difdf$Run,
                   levels=c('All Surveys',
                            'Only Trawls',
                            'Trawls+Sentinel,\nno BLLS',
                            'Trawls+BLLS,\nno Sentinel'))

# Set levels of season
difdf$Season <- factor(difdf$Season,
                      levels=c('Spring', 'Fall'))

# Adjust year to numeric
difdf$Year <- as.numeric(difdf$Year)

# Difference over time - Fall
fp <- ggplot() +
  geom_line(data=difdf[difdf$Season == "Fall",],
            aes(x=Year, y=Dif, col=Run)) +
  geom_ribbon(data=difdf[difdf$Season == "Fall",],
              aes(x=Year, ymin=LoDif, ymax=UpDif,
                  fill=Run), 
              alpha=0.2, linetype=0) +
  facet_wrap(vars(Category),
            scales='free_y') +
  ggtitle('Fall index')
ggsave(plot = fp,
       path = 'C:/Users/klankowicz/Desktop/',
       filename = 'Fall_Prop_GearDrop.png',
       device='png')
  

# Difference over time - Spring
sp <- ggplot() +
  geom_line(data=difdf[difdf$Season == "Spring",],
            aes(x=Year, y=Dif, col=Run)) +
  geom_ribbon(data=difdf[difdf$Season == "Spring",],
              aes(x=Year, ymin=LoDif, ymax=UpDif,
                  fill=Run), 
              alpha=0.2, linetype=0) +
  facet_wrap(vars(Category),
             scales='free_y') +
  ggtitle('Spring index')
ggsave(plot = sp,
       path = 'C:/Users/klankowicz/Desktop/',
       filename = 'Spring_Prop_GearDrop.png',
       device='png')
