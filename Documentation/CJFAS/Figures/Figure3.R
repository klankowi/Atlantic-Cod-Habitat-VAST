### Make Figure 3, size- and season-specific population cog lineplots ####
#### CJFAS 2024 ####
rm(list=ls())

library(tidyverse)
library(here)
library(ggh4x)
library(ggpubr)
library(ggtext)
library(gridExtra)

# Function to tag facets without removing strips (based on egg::tag_facet)
tag_facet2 <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

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

# Load all indices from all sizes
sizes <- c('small', 'medium', 'large')
for(i in 1:length(sizes)){
  cog <- read.csv(paste0(here('VAST_runs'), '/',
                        sizes[i], 
                        '/Overall_BC/ALL/COG_ALL.csv'))
  cog$Size <- paste0(str_to_sentence(sizes[i]))
  cog <- cog %>% 
    mutate(northing = northing / 1000,
           easting = easting / 1000,
           n.sd = n.sd / 1000,
           e.sd = e.sd / 1000)
  assign(paste0(sizes[i]), cog)
  rm(cog)
}

# Combine to one DF
cog <- rbind(small, medium, large)

# Clean
cog <- cog %>% 
  filter(strata == 'ALL') %>% 
  rename(northing.sd = n.sd,
         easting.sd = e.sd) %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large'))) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  dplyr::select(Size, Year, Season, northing, northing.sd,
                easting, easting.sd) %>% 
  pivot_longer(cols = c('northing', 'easting'),
               names_to = 'Direction',
               values_to = 'Estimate') %>% 
  pivot_longer(cols=c('northing.sd', 'easting.sd'),
               names_to = 'Direction2',
               values_to = 'SD') %>% 
  separate(Direction2, into=c('Direction2', 'Trash2'), sep = "[^[:alnum:]]+") %>% 
  filter(Direction == Direction2) %>% 
  dplyr::select(-Direction2,-Trash2) %>% 
  as.data.frame()

# CLean workspace
rm(small, medium, large, sizes, i)

# Base plot, Northing
fig3a <- ggplot(data=cog[cog$Direction == 'northing',]) +
  geom_point(aes(x=Year, y=Estimate), cex=0.5) +
  geom_line(aes(x=Year, y=Estimate)) +
  geom_ribbon(aes(x=Year, ymin = Estimate - SD, 
                  ymax= Estimate + SD), alpha=0.4) +
  labs(y='Northing (km)', x='') +
  scale_x_continuous(limits=c(1982, 2021), expand=c(0.01,0.01)) +
  scale_y_continuous(limits=c(4575, 4825), expand=c(0,0)) +
  facet_grid2(c("Size", "Season"))+
  
  ggtitle('Northing') +
  
  theme(strip.background = element_rect(fill='lightgray'),
        legend.position = 'none',
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.length.x = unit(0, "cm"),
        plot.title = element_markdown(box.color = "black", 
                                      fill = 'lightgray',
                                      linewidth = 0.1,
                                      r = unit(0, "mm"), 
                                      linetype = 1,
                                      padding = unit(c(2,196.5, 2, 196.5), 'pt')))

# Label facets
fig3a <- tag_facet2(fig3a, open='', close='', col='gray30') 
fig3a

# Base plot, Easting
fig3b <- ggplot(data=cog[cog$Direction == 'easting',]) +
  geom_point(aes(x=Year, y=Estimate), cex=0.5) +
  geom_line(aes(x=Year, y=Estimate)) +
  geom_ribbon(aes(x=Year, ymin = Estimate - SD, 
                  ymax= Estimate + SD), alpha=0.4) +
  labs(col = 'Quantile', fill='Quantile',
       y='Easting (km)') +
  scale_x_continuous(limits=c(1982, 2021), expand=c(0.02,0.02)) +
  scale_y_continuous(limits=c(350, 675), expand=c(0,0)) +
  facet_grid2(c("Size", "Season")) +
  
  ggtitle('Easting') +
  
  theme(strip.background = element_rect(fill='lightgray'),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.box.spacing = margin(9,0,0,0),
        strip.text.x = element_blank(),
        plot.title = element_markdown(box.color = "black", 
                                      fill = 'lightgray',
                                      linewidth = 0.1,
                                      r = unit(0, "mm"), 
                                      linetype = 1,
                                      padding = unit(c(2,200, 2, 200), 'pt')))

# Label facets
fig3b <- tag_facet2(fig3b, open='', close='', col='gray30', tag_pool = letters[7:26])
fig3b

# GGarrange to join northings and eastings
fig3 <- ggarrange(fig3a, fig3b, 
                  ncol=1, nrow=2, 
                  common.legend = TRUE, legend="bottom",
                  align='v')

# Save local
ggsave(here('Documentation/Figures/Fig 3.pdf'),
       fig3,
       height=23.7, width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
       fig3,
       height=23.7, width=18.2, units='cm',
       dpi = 600)
