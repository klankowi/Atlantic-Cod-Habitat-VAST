### Make Figure 5, size- and season-specific population range edge lineplots ####
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
  re <- read.csv(paste0(here('VAST_runs'), '/',
                           sizes[i], 
                           '/Overall_BC/ALL_Catchability2/RangeEdges_ALL.csv'))
  re$Size <- paste0(str_to_sentence(sizes[i]))
  assign(paste0(sizes[i]), re)
  rm(re)
}

# Combine to one DF
re <- rbind(small, medium, large)

# Clean
re <- re %>% 
  filter(strata == 'ALL') %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large'))) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Quantile = factor(Quantile, levels=c('0.05', '0.5', '0.95'))) %>% 
  dplyr::select(Size, Year, Season, Quantile, Northing.Est, Northing.SD,
                Easting.Est, Easting.SD) %>% 
  pivot_longer(cols = c('Northing.Est', 'Easting.Est'),
               names_to = 'Direction',
               values_to = 'Estimate') %>% 
  pivot_longer(cols=c('Northing.SD', 'Easting.SD'),
               names_to = 'Direction2',
                values_to = 'SD') %>% 
  separate(Direction, into=c('Direction', 'Trash'), sep = "[^[:alnum:]]+") %>% 
  separate(Direction2, into=c('Direction2', 'Trash2'), sep = "[^[:alnum:]]+") %>% 
  filter(Direction == Direction2) %>% 
  dplyr::select(-Direction2, -Trash, -Trash2) %>% 
  as.data.frame()

# CLean workspace
rm(small, medium, large, sizes, i)

# Base plot, Northing
fig5a <- ggplot(data=re[re$Direction == 'Northing',]) +
  geom_point(aes(x=Year, y=Estimate, col=Quantile), cex=0.5) +
  geom_line(aes(x=Year, y=Estimate, col=Quantile)) +
  geom_ribbon(aes(x=Year, ymin = Estimate - SD, 
                  ymax= Estimate + SD,
                  fill=Quantile), alpha=0.4) +
  labs(col = 'Quantile', fill='Quantile',
       y='Northing (km)', x='') +
  scale_x_continuous(limits=c(1982, 2021), expand=c(0.02,0.02)) +
  scale_y_continuous(limits=c(4350, 5100), expand=c(0,0)) +
  facet_grid2(c("Size", 
                "Season"))+
  
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
fig5a <- tag_facet2(fig5a, open='', close='', col='gray30') 

# Base plot, Easting
fig5b <- ggplot(data=re[re$Direction == 'Easting',]) +
  geom_point(aes(x=Year, y=Estimate, col=Quantile), cex=0.5) +
  geom_line(aes(x=Year, y=Estimate, col=Quantile)) +
  geom_ribbon(aes(x=Year, ymin = Estimate - SD, 
                  ymax= Estimate + SD,
                  fill=Quantile), alpha=0.4) +
  labs(col = 'Quantile', fill='Quantile',
       y='Easting (km)') +
  scale_x_continuous(limits=c(1982, 2021), expand=c(0.02,0.02)) +
  scale_y_continuous(limits=c(100, 900), expand=c(0,0)) +
  facet_grid2(c("Size", 
                "Season")) +
  
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
fig5b <- tag_facet2(fig5b, open='', close='', col='gray30', tag_pool = letters[7:26])


# GGarrange to join northings and eastings
fig5 <- ggarrange(fig5a, fig5b, 
          ncol=1, nrow=2, 
          common.legend = TRUE, legend="bottom",
          align='v')

# Save local
ggsave(here('Documentation/CJFAS/Figures/Fig 5.pdf'),
       fig5,
       height=23.7, width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 5.pdf",
       fig5,
       height=23.7, width=18.2, units='cm',
       dpi = 600)
