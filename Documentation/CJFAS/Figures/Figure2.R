### Make Figure 2, size- and season-specific indices of abundance ####
#### CJFAS 2024 ####
rm(list=ls())

library(tidyverse)
library(here)
library(ggh4x)

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
  index <- read.csv(paste0(here('VAST_runs'), '/',
                           sizes[i], 
                           '/Overall_BC/ALL/Index.csv'))
  assign(paste0(sizes[i]), index)
  rm(index)
}

# Combine to one DF
index <- rbind(small, medium, large)

# Clean
index <- index %>% 
  filter(Stratum != 'ALL') %>% 
  separate(Category, into=c('Size', 'Cod')) %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large'))) %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall')),
         Stratum = factor(Stratum, levels=c('EGOM', 'GBK', 'SNE', 'WGOM'))) %>% 
  rename(Std.Err = Std..Error.for.Estimate,
         Std.Err.ln = Std..Error.for.ln.Estimate.) %>% 
  dplyr::select(Size, Year, Season, Stratum, Estimate, Std.Err, Std.Err.ln)

# Base plot
fig2 <- ggplot(data=index) +
  geom_point(aes(x=Year, y=(Estimate/1000000), col=Stratum), cex=0.5) +
  geom_line(aes(x=Year, y=(Estimate/1000000), col=Stratum)) +
  geom_ribbon(aes(x=Year, ymin = (Estimate/1000000) - (Std.Err/1000000), 
                  ymax= (Estimate/1000000) + (Std.Err/1000000),
                  fill=Stratum), alpha=0.4) +
  labs(col = 'Biological\nStock', fill='Biological\nStock',
       y='Relative Abundance, Millions') +
  scale_x_continuous(limits=c(1982, 2021), expand=c(0.02,0.02)) +
  facet_grid2(c("Size", "Season"), 
              scales = "free", independent = "all") +
  
  theme(strip.background = element_rect(fill='lightgray'),
        legend.box.margin = margin(-10, -10, -10, -10))

# Label facets
fig2 <- tag_facet2(fig2, open='', close='', col='gray30')

# Save local
ggsave(here('Documentation/Figures/Fig 2.pdf'),
       fig2,
       height=(18.2 * 0.75), width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 2.pdf",
       fig2,
       height=(18.2 * 0.75), width=18.2, units='cm',
       dpi = 600)