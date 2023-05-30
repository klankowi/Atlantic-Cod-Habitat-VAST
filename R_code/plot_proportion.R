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
indd <- read.csv(here("VAST_runs/add_climate_aja4/GBK/Index.csv"))

# Name separations
category_names <- c('Small', 'Medium', 'Large')
season_names <- c('Spring', 'Fall')
strata_names <- c( "EGOM", "GBK", "SNE", "WGOM", "ALL")

# Separate year and season
indd <- indd %>% 
  separate(Time, 
           into = c('Year', 'Season'),
           sep = " ")

# Rename strata
indd$Stratum[indd$Stratum == 'Stratum_1'] <- 'GBK'
indd$Stratum[indd$Stratum == 'Stratum_2'] <- 'EGOM'
indd$Stratum[indd$Stratum == 'Stratum_3'] <- 'SNE'
indd$Stratum[indd$Stratum == 'Stratum_4'] <- 'WGOM'
indd$Stratum[indd$Stratum == 'Stratum_5'] <- 'ALL'

# Convert to factor
indd$Stratum <- factor(indd$Stratum,
                       levels=c('EGOM', 'GBK', 'SNE', 'WGOM', 'ALL'))
indd$Category <- factor(indd$Category,
                        levels=c('Small', 'Medium', 'Large'))
indd$Season <- factor(indd$Season,
                      levels=c('Spring', 'Fall'))

# Multilevel list
for(i in 1:length(category_names)){ # Category
  temp.cat <- subset(indd, Category == paste0(category_names[i]))
  for(j in 1:length(season_names)){ # Season
    temp.sea <- subset(temp.cat, Season == paste0(season_names[[j]]))
    temp.sea <- split(temp.sea, f=temp.sea$Stratum)
    
    for(k in 1:4){
      temp.sea[[k]]$prop <- 
        temp.sea[[k]]$Estimate /
        temp.sea[[5]]$Estimate *
        100
    }
    
    temp.sea <- temp.sea[-5]
    temp.sea <- do.call(rbind, temp.sea)
    head(temp.sea)
    row.names(temp.sea) <- NULL
    temp.sea$Year <- as.numeric(temp.sea$Year)
    
    pcod <- ggplot() +
      geom_line(data=temp.sea,
                aes(x=Year, y=prop, color=Stratum),
                lwd=1) +
      ylab('Proportion') +
      ggtitle(paste0(temp.sea$Category[1],
                     ' cod, ',
                     temp.sea$Season[1]))
    
    ggsave(pcod,
           filename=paste0(here(), "/Plot_output_3/proportion.info.",
                           category_names[i], '.',
                           season_names[j], '.png'))
    
  }
  
}

