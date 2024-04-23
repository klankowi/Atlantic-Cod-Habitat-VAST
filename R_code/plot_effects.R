rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)

# load data
load(here('VAST_runs/tuna13/usa_first/tuna13_usafirst.RData'))

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

  vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = fit, 
                                                 params_plot = c(params), 
                                                 params_plot_levels = 100, 
                                                 effects_pad_values = c(), 
                                                 nice_category_names = "USA First", 
                                                 out_dir = here('VAST_runs/tuna13/usa_first/'),
                                                 category_to_use = 1,
                                                 ncat = ncat)
  
  vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, 
                                                  vast_fit = fit, 
                                                  nice_category_names = "USA First",  
                                                  out_dir = here('VAST_runs/tuna13/usa_first/'))

