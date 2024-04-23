# Prepare workspace
rm(list=ls())

# Load libraries
library(VAST)
library(tidyverse)
library(sf)
library(here)
library(ggpubr)
library(ggnewscale)
library(patchwork)

# Load functions
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load data
indd <- read.csv(here("VAST_runs/herring_graham/Index.csv"))

# Name separations
season_names <- c('Winter-Spring Feeding', 'Summer Feeding-Spawning')
strata_names <- c( "All", 'GOMI', 'GOMO', 'GBK', 'SNE')
category_names <- 'Herring'

# Separate year and season
indd <- indd %>% 
  separate(Time, 
           into = c('Year', 'Season'),
           sep = " ")

# Rename strata
indd$Stratum[indd$Stratum == 'Stratum_1'] <- 'All'
indd$Stratum[indd$Stratum == 'Stratum_2'] <- 'GOMI'
indd$Stratum[indd$Stratum == 'Stratum_3'] <- 'GOMO'
indd$Stratum[indd$Stratum == 'Stratum_4'] <- 'GBK'
indd$Stratum[indd$Stratum == 'Stratum_5'] <- 'SNE'

indd$Season[indd$Season == 'Jan-Jun'] <- 'Winter-Spring Feeding'
indd$Season[indd$Season == 'Jul-Dec'] <- 'Summer Feeding-Spawning'

# Convert to factor
indd$Stratum <- factor(indd$Stratum,
                       levels=c("All",'GBK', 'GOMI', 'GOMO', 'SNE'))
indd$Season <- factor(indd$Season,
                      levels=c('Winter-Spring Feeding', 
                               'Summer Feeding-Spawning'))
indd$Category <- 'Herring'

# Multilevel list
for(i in 1:length(category_names)){ # Category
  temp.cat <- subset(indd, Category == paste0(category_names[i]))
  for(j in 1:length(season_names)){ # Season
    temp.sea <- subset(temp.cat, Season == paste0(season_names[[j]]))
    temp.sea <- split(temp.sea, f=temp.sea$Stratum)
    
    for(k in 1:length(strata_names)){
      temp.sea[[k]]$prop <- 
        temp.sea[[k]]$Estimate /
        temp.sea[[1]]$Estimate *
        100
    }
    
    temp.sea <- temp.sea[-1]
    temp.sea <- do.call(rbind, temp.sea)
    head(temp.sea)
    row.names(temp.sea) <- NULL
    temp.sea$Year <- as.numeric(temp.sea$Year)
    
    pcod <- ggplot() +
      geom_line(data=temp.sea,
                aes(x=Year, y=prop, color=Stratum),
                lwd=1) +
      ylab('Proportion') +
      ggtitle(paste0('Proportion of herring stock in each strata, ',
                     temp.sea$Season[1]))
    
    ggsave(pcod,
           filename=paste0(here('VAST_runs'), "/herring_graham/",
                           season_names[j], '.png'))
    
  }
  
}

