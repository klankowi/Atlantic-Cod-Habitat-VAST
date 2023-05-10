rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)

# load data
load(here("VAST_runs/add_climate_aja3/add_climate_aja3.Rdata"))
#load(here("VAST_runs/refine_effort/refine_effort.RData"))

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Call parameter names
params <- colnames(scaled.covars)
params <- params[params %notin% c('Lon', 'Lat', 'Year')]

# Call category names
catnames <- fit$category_names
ncat = length(catnames)

for(i in 2:length(catnames)){
  vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = fit, 
                                                 params_plot = c(params), 
                                                 params_plot_levels = 100, 
                                                 effects_pad_values = c(), 
                                                 nice_category_names = paste(catnames[i], "Atlantic cod"), 
                                                 out_dir = here('VAST_runs/add_climate_aja3/effects'),
                                                 category_to_use = i,
                                                 ncat = ncat)
  
  # Warnings are...interesting....
  vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, 
                                                  vast_fit = fit, 
                                                  nice_category_names = paste(catnames[i], "Atlantic cod"),  
                                                  out_dir = here('VAST_runs/add_climate_aja3/effects'))
  vast_habcovs_plot
}
