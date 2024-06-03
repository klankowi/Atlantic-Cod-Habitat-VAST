rm(list=ls())

library(VAST)
library(here)
library(tidyverse)

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

#install_unit(symbol='unitless', def='unitless', name='unitless')

source(here('R_code/utilities/extract_fit_cog.R'))
source(here('R_code/utilities/extract_fit_range_edges.R'))
source(here('R_code/utilities/extract_fit_eff.R'))
source(here('R_code/utilities/extract_fit_index.R'))

# Covariate
covuse <- c('amo',
            #'bathy',
            #'bottomtemp',
            #'cobble',
            #'gravel',
            #'mud',
            'nao',
            'rugos'#,
            #'sand'
            )
#covuse <- c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')

for(i in 1:length(covuse)){
  message(covuse[i])
  # Load VAST fit data
  load(here(paste0('VAST_runs/large/AIC/', covuse[i],
            '/no', covuse[i], '_largecod_wholearea_natsplin_fsOFF',
            #covuse[i], 
            '.RData')))
  
  setwd(here(paste0("VAST_runs/large/AIC/", covuse[i], '/')))
  
  extract_fit_cog(fit=fit)
  extract_fit_eff(fit=fit)
  #extract_fit_index(fit=fit, category_names='Large')
  extract_fit_range_edges(fit=fit,
                          quantvals = c(0.05, 0.5, 0.95))
  
  plot_results(fit, plot_set = c(1, 2, 3, 6, 7, 8, 9,
                                 11, 12, 13, 14, 15, 16, 17, 18,
                                 19, 20, 21))
  
  rm(list=setdiff(ls(), c('covuse', 
                          'extract_fit_cog',
                          'extract_fit_eff',
                          'extract_fit_range_edges', 
                          'extract_fit_index', 
                          'i')))
}
