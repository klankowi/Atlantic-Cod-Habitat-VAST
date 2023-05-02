rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)

# load data
load(here("VAST_runs/refine_effort/refine_effort.RData"))

# Load functions
source(here("R_code/utilities/vast_functions.R"))

# Call parameter names
params <- as.character(fit$X1_formula)
params <- str_split(params, " +")
params <- params[[2]]
params <- params[params != "+"]

# Call category names
catnames <- fit$category_names

vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = fit, 
                                               params_plot = c(params), 
                                               params_plot_levels = 100, 
                                               effects_pad_values = c(), 
                                               nice_category_names = paste(catnames[i], "Atlantic cod"), 
                                               out_dir = here::here("", "results/tables"),
                                               category_to_use = i,
                                               ncat = 4)

# Warnings are...interesting....
vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, 
                                                vast_fit = fit, 
                                                nice_category_names = paste(catnames[i], "Atlantic cod"),  
                                                out_dir = here::here("", "results/plots_maps"))
vast_habcovs_plot