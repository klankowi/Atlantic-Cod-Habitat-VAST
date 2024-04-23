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
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load data
indd <- read.csv(here("VAST_runs/small_climate/Index.csv"))

seasons <- data.frame(
  Time = seq(70, 80),
  Season = c('2016 Fall',
             '2017 Spring', '2017 Fall',
             '2018 Spring', '2018 Fall',
             '2019 Spring', '2019 Fall',
             '2020 Spring', '2020 Fall',
             '2021 Spring', '2021 Fall')
)

indd <- merge(indd, seasons)

indd <- indd %>% 
  separate(Season, into=c('Year', 'Season'))

# Clean
indd <- indd %>% 
  filter(Year != 2016) %>% 
  dplyr::select(-Category, -Units, -Std..Error.for.ln.Estimate.) %>% 
  rename(Std.Err = Std..Error.for.Estimate)

# Rename strata
indd$Stratum[indd$Stratum == 'Stratum_1'] <- 'All'
indd$Stratum[indd$Stratum == 'Stratum_2'] <- 'GOM'
indd$Stratum[indd$Stratum == 'Stratum_3'] <- 'GBK'
indd$Stratum[indd$Stratum == 'Stratum_4'] <- 'WEA'

# Convert to factor
indd$Stratum <- factor(indd$Stratum,
                       levels=c('All', 'GOM', 'GBK', 'WEA'))
indd$Year <- as.numeric(indd$Year)

# Plot
prop.wea <- ggplot(data=indd[indd$Stratum != 'All',]) +
  geom_ribbon(aes(x=Year, ymin=Estimate - Std.Err, ymax=Estimate + Std.Err,
                  fill=Stratum), alpha=0.4, col=NA) +
  geom_line(aes(x=Year, y=Estimate, col=Stratum), lwd=1) +
  facet_wrap(vars(Season))+
  ggtitle('Medium cod relative abundance')

# Save proportion
ggsave(prop.wea,
       filename=here("VAST_runs/small_climate/rel_abundance.png"))

# Split to determine proportion
temp.sea <- split(indd, f=indd$Stratum)
    
for(k in 2:length(temp.sea)){
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
temp.sea$Stratum <- factor(temp.sea$Stratum,
                           levels=c('GOM', 'GBK', 'WEA'))

# Plot proportion
pcod <- ggplot() +
  geom_line(data=temp.sea,
            aes(x=Year, y=prop, color=Stratum),
                lwd=1) +
  facet_wrap(vars(Season)) +
  ylab('Proportion') +
  ggtitle('Proportion of estimated Medium Cod abundance per strata')

# Save proportion
ggsave(pcod,
       filename=here("VAST_runs/small_climate/proportion_info.png"))

