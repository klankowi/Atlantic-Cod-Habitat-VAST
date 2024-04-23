rm(list=ls())

library(VAST)
library(here)
library(tidyverse)

install_unit(symbol='unitless', def='unitless', name='unitless')

source(here('R_code/utilities/extract_fit_cog.R'))
source(here('R_code/utilities/extract_fit_range_edges.R'))
source(here('R_code/utilities/extract_fit_eff.R'))

# Covariate
covuse <- c('amo', 'bathy', 'bottomtemp', 'cobble', 'gravel',
            'mud', 'nao', 'rugos', 'sand')

for(i in 1:length(covuse)){
  message(covuse[i])
  # Load VAST fit data
  load(here(paste0('VAST_runs/medium/AIC/', covuse[i], '/no',
            covuse[i], '_medcod_wholearea_natsplin_fsOFF.RData')))
  
  setwd(here(paste0("VAST_runs/medium/AIC/", covuse[i])))
  
  extract_fit_cog(fit=fit)
  extract_fit_eff(fit=fit)
  extract_fit_range_edges(fit=fit,
                          quantvals = c(0.05, 0.5, 0.95))
  
  plot_results(fit, plot_set = c(1, 2, 3, 6, 7, 8, 9,
                                 11, 12, 13, 14, 15, 16, 17, 18, 
                                 19, 20, 21))
  
  rm(list=setdiff(ls(), c('covuse', 'extract_fit_cog',
                          'extract_fit_eff',
                          'extract_fit_range_edges', 'i')))
}
