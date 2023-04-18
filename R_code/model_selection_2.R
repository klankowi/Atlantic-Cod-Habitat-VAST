# VAST attempt 2 univariate model selection as a script
# modified from 
# https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization
###############################################################################
# Load necessities
rm(list=ls())
gc()

# Load packages
library(here)
library(dplyr)
library(VAST)
library(tidyverse)
library(units)
library(sf)
library(ggcorrplot)
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))
# Add unitless back as possible unit (removed in units package update Mar 2023)
install_unit(symbol='unitless', def='unitless', name='unitless')
# Read in data
# Load data
dat <- read.csv(here('Data/VAST_input/cod_agesep_VASTdata.csv'))
dat$BATHY.DEPTH <- dat$BATHY.DEPTH * -1
head(dat)

# Save sampling data
survs <- dplyr::select(dat,
                       LON, LAT, AREA_SWEPT,
                       TIME, SURVEY, RESPONSE, 
                       AGEGROUP)
survs$RESPONSE <- as_units(survs$RESPONSE, 'counts')
survs$AREA_SWEPT <- as_units(survs$AREA_SWEPT, 'km^2')
#survs$swept <- survs$AREA_SWEPT
survs$SURVEY <- as.numeric(as.factor(survs$SURVEY)) - 1
# vessel    survey
# 0         ASMFC Shrimp Trawl  
# 1         DFO Trawl  
# 2         GSO Trawl  
# 3         MADMF Industry  
# 4         MADMF Inshore Trawl  
# 5         ME-NH Inshore Trawl  
# 6         NEFSC BLLS   
# 7         NEFSC BTS 
# 8         RIDEM Trawl  
# 9         Sentinel   
# 10        SMAST Video Trawl   

#survs$Data_type <- factor(survs$Data_type, levels=c("Count", "Biomass_KG"))

survs$AGEGROUP <- as.numeric(factor(survs$AGEGROUP, levels=c(
  'Unknown',
  'Age0-2', 
  'Age2-5', 
  'Age5+'
)
)) - 1
# Age 2 - 5: 0
# Age 5+   : 1
table(survs$AGE)

names(survs) <- c('Lon', 'Lat', 'swept',
                  'Year', 'vessel',
                  'Response_variable', 
                  'Age')
#table(survs$Year)
str(survs)

# Save covariates
covars <- dplyr::select(dat,
                        LON, LAT, TIME, cobble_P, gravel_P,
                        mud_P, rock_P, sand_P, rugos, BATHY.DEPTH, oisst,
                        h_bt, nao, amo)
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:ncol(covars)])
table(covars$Year)

# Test correlation
# Create correlation matrix
df_cormat <- dplyr::select(covars, BATHY.DEPTH, rugos, sand_P, rock_P, mud_P,
                           gravel_P, cobble_P, oisst, h_bt, nao, amo)
model.matrix(~0+., data=df_cormat) %>%
  cor(use="all.obs", method="spearman") %>%
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=3)
# OISST and h_bt highly correlated. Remove OISST.
# Rock and cobble highly correlated. Remove rock.
covars$rock_P <- NULL; covars$oisst <- NULL
# Create correlation matrix
df_cormat <- dplyr::select(covars, BATHY.DEPTH, rugos, sand_P, mud_P,
                           gravel_P, cobble_P, h_bt, nao, amo)
model.matrix(~0+., data=df_cormat) %>%
  cor(use="all.obs", method="spearman") %>%
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=3)


# Rescale covariates to have mean 0 and SD 1 (author rec)
scaled.covars <- covars[,4:ncol(covars)] %>% 
  mutate(across(where(is.numeric), scale))
scaled.covars <- cbind(covars[,1:3], scaled.covars)
summary(scaled.covars)
scaled.covars <- data.frame(
  Lon         = as.numeric(scaled.covars$Lon),
  Lat         = as.numeric(scaled.covars$Lat),
  Year        = as.numeric(scaled.covars$Year),
  cobble_P    = as.numeric(scaled.covars$cobble_P),
  gravel_P    = as.numeric(scaled.covars$gravel_P),
  mud_P       = as.numeric(scaled.covars$mud_P),
  #rock_P      = as.numeric(scaled.covars$rock_P),
  sand_P      = as.numeric(scaled.covars$sand_P),
  rugos       = as.numeric(scaled.covars$rugos),
  BATHY.DEPTH = as.numeric(scaled.covars$BATHY.DEPTH),
  #oisst       = as.numeric(scaled.covars$oisst),
  h_bt        = as.numeric(scaled.covars$h_bt),
  nao         = as.numeric(scaled.covars$nao),
  amo         = as.numeric(scaled.covars$amo)
)
str(scaled.covars)

#### Make settings ####
user_region <- readRDS(here('Data/VAST_input/user_region_all.rds'))
user_region$STRATA <- 'All'
user_region$Id <- NULL; user_region$row <- NULL
user_region <- user_region[with(user_region, order(Lon, Lat)),]
row.names(user_region) <- NULL
head(user_region)
strata_use <- data.frame('STRATA' = c("All"))

# Remove intermediates
rm(covars, dat, df_cormat)
gc()
###############################################################################
# Describe model selection efforts

# Model selection 1 (spatial, spatio-temporal effects, no covariates) options 

# Field configs
# _alleffectson             FieldConfig default (all IID)
# _noaniso                  FieldConfig default (all IID) and 
#                                use_anistropy = FALSE
# _noomeps2                 FieldConfig 0 for Omega2, Epsilon2
# _noomeps2_noaniso         FieldConfig 0 for Omega2, Epsilon2 and 
#                                use_anistropy = FALSE
# _noomeps2_noeps1          FieldConfig 0 for Omega2, Epsilon2, Epsilon1
# _noomeps2_noeps1_noaniso  FieldConfig 0 for Omega2, Epsilon2, Epsilon1 and 
#                                use_anistropy = FALSE
# _noomeps12                FieldConfig both Omega, Epsilon 0
# _noomeps12_noaniso        FieldConfig both Omega, Epsilon 0 and 
#                                use_anistropy = FALSE

# default configs
FieldConfig = matrix( "IID", ncol=2, nrow=3, 
                      dimnames=list(c("Omega","Epsilon","Beta"),
                                    c("Component_1","Component_2")))

# Rho configs
# not testing alternative RhoConfigs here just noted for completeness
# 0 off (fixed effects)
# 1 independent
# 2 random walk
# 3 constant among years (fixed effect)
# 4 AR1

RhoConfig <- c(
  "Beta1" = 0,      # temporal structure on years (intercepts) 
  "Beta2" = 0, 
  "Epsilon1" = 0,   # temporal structure on spatio-temporal variation
  "Epsilon2" = 0
) 

# Anisotropy
use_anisotropy <- TRUE

OverdispersionConfig	<- c("eta1"=0, "eta2"=0)
# eta0 = no vessel effects
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

# list of data, settings, and directory for output for each option
mod.config <- c("alleffectson", "noaniso", 
                "noomeps2", "noomeps2_noaniso", 
                "noomeps2_noeps1", "noomeps2_noeps1_noaniso",
                "noomeps12", "noomeps12_noaniso")

# Define possible field configurations
FieldConfig1 <- matrix( "IID", ncol=2, nrow=3, 
                        dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")))
FieldConfig2 <- matrix( c("IID","IID","IID",0,0,"IID"), ncol=2, nrow=3, 
                        dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")))
FieldConfig3 <- matrix( c("IID",0,"IID",0,0,"IID"), ncol=2, nrow=3, 
                        dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")))
FieldConfig4 <- matrix( c(0,0,"IID",0,0,"IID"), ncol=2, nrow=3, 
                        dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")))

# Pull field configs into list
mod.FieldConfig <- list(FieldConfig1, FieldConfig1,
                        FieldConfig2, FieldConfig2,
                        FieldConfig3, FieldConfig3,
                        FieldConfig4, FieldConfig4)

# Name list items
names(mod.FieldConfig) <- mod.config

# List possible anisotropy options
mod.use_anistropy <- list(TRUE, FALSE, 
                          TRUE, FALSE,
                          TRUE, FALSE,
                          TRUE, FALSE)
# Name list items
names(mod.use_anistropy) <- mod.config

###############################################################################
# Run  model selection 1

# Subset for model selection 1
# Use only last 10 years of NEFSC BTS data
sel1 <- survs[survs$Year >=61 & survs$vessel == 7,]

# Loop through options
for(i in 1:length(mod.config)) {
  # Define name of model run
  name <- mod.config[i]
  # Set working directory
  working_dir <- here::here(sprintf("Model_Refinement/Model_Structure2/%s", name))
  # Create working directory if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
  }
  # Call model options to be used
  FieldConfig <- mod.FieldConfig[[i]]
  use_anisotropy <- mod.use_anistropy[[i]]
  # Make settings
  settings <- make_settings(n_x = 200,
                            Region = "User",
                            Version = "VAST_v14_0_1",
                            purpose = "index2",
                            bias.correct = FALSE,
                            use_anisotropy = use_anisotropy,
                            FieldConfig = FieldConfig,
                            RhoConfig = RhoConfig, # always default
                            OverdispersionConfig = OverdispersionConfig #default
  )
  
  # Set obsmodel to reflect count response variable measure
  settings$ObsModel[1] <- 4
  settings$ObsModel[2] <- 1
  # Fit model
  fit <- fit_model(
    # Call REML
    Use_REML = TRUE,
    # Call settings
    settings = settings,
    # Call survey data info
    Lat_i = sel1[,'Lat'], 
    Lon_i = sel1[,'Lon'], 
    t_i = sel1[,'Year'],
    b_i = sel1[,'Response_variable'],
    a_i = sel1[,'swept'],
    v_i = sel1[,'vessel'],
    c_iz = sel1[,'Age'],
    # Call spatial
    input_grid = user_region,
    # Set directory
    working_dir = paste0(working_dir, "/"),
    # Tell model to run
    run_model = TRUE)
} # end config loop

###############################################################################
# Compare model selection 1 results
# Set directory 
outdir <- here("mod_selection")
# Call file names
moddirs <- list.dirs(outdir) 
# Remove name of upper level file
moddirs <- moddirs[-1]
# keep folder name
modnames <- list.dirs(outdir, full.names = FALSE)

# function to apply extracting info
getmodinfo <- function(d.name){
  # read settings
  modpath <- stringr::str_split(d.name, "/", simplify = TRUE)
  modname <- modpath[length(modpath)]
  
  settings <- read.table(file.path(d.name, "settings.txt"), comment.char = "",
                         fill = TRUE, header = FALSE)
  
  n_x <- as.numeric(as.character(settings[(which(settings[,1]=="$n_x")+1),2]))
  grid_size_km <- as.numeric(as.character(settings[(
    which(settings[,1]=="$grid_size_km")+1),2]))
  max_cells <- as.numeric(as.character(settings[(
    which(settings[,1]=="$max_cells")+1),2]))
  use_anisotropy <- as.character(settings[(
    which(settings[,1]=="$use_anisotropy")+1),2])
  fine_scale <- as.character(settings[(
    which(settings[,1]=="$fine_scale")+1),2])
  bias.correct <- as.character(settings[(
    which(settings[,1]=="$bias.correct")+1),2])
  
  #FieldConfig
  if(settings[(which(settings[,1]=="$FieldConfig")+1),1]=="Component_1"){
    omega1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+2),2])
    omega2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+3),1])
    epsilon1 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+4),2])
    epsilon2 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+5),1])
    beta1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+6),2])
    beta2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+7),1])
  }
  
  if(settings[(which(settings[,1]=="$FieldConfig")+1),1]=="Omega1"){
    omega1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+3),1])
    omega2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+4),1])
    epsilon1 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+3),2])
    epsilon2 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+4),2])
    beta1 <- "IID"
    beta2 <- "IID"
  }
  
  #RhoConfig
  rho_beta1 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+3),1]))
  rho_beta2 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+3),2]))
  rho_epsilon1 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+4),1]))
  rho_epsilon2 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+4),2]))
  
  # read parameter estimates, object is called parameter_Estimates
  load(file.path(d.name, "parameter_estimates.RData"))
  
  AIC <- parameter_estimates$AIC[1]  
  converged <- parameter_estimates$Convergence_check[1]
  fixedcoeff <- unname(parameter_estimates$number_of_coefficients[2])
  randomcoeff <- unname(parameter_estimates$number_of_coefficients[3])
  
  # return model atributes as a dataframe
  out <- data.frame(modname = modname,
                    n_x = n_x,
                    grid_size_km = grid_size_km,
                    max_cells = max_cells,
                    use_anisotropy = use_anisotropy,
                    fine_scale =  fine_scale,
                    bias.correct = bias.correct,
                    omega1 = omega1,
                    omega2 = omega2,
                    epsilon1 = epsilon1,
                    epsilon2 = epsilon2,
                    beta1 = beta1,
                    beta2 = beta2,
                    rho_epsilon1 = rho_epsilon1,
                    rho_epsilon2 = rho_epsilon2,
                    rho_beta1 = rho_beta1,
                    rho_beta2 = rho_beta2,
                    AIC = AIC,
                    converged = converged,
                    fixedcoeff = fixedcoeff,
                    randomcoeff = randomcoeff
  )
  return(out)
}

# Pull models
modselect <- purrr::map_dfr(moddirs, getmodinfo)

# Build table to compare models
modselect.200 <- modselect %>%
  filter(n_x == 200) %>%
  mutate(converged2 = case_when(str_detect(converged, 
                                           "no evidence") ~ "likely",
                                str_detect(converged, 
                                           "is likely not") ~ "unlikely",
                                TRUE ~ as.character(NA))) %>%
  mutate(deltaAIC = AIC-min(AIC)) %>%
  select(modname, deltaAIC, fixedcoeff,
         randomcoeff, use_anisotropy, 
         omega1, omega2, epsilon1, epsilon2, 
         beta1, beta2, AIC, converged2) %>%
  arrange(AIC)

# Print table
DT::datatable(modselect.200, rownames = FALSE, 
              options= list(pageLength = 25, scrollX = TRUE))

# Evidence suggests best model includes all effects
# Model is marginally better without anisotropy, but effect is less than 1 unit AIC.
# Will retain anisotropy

###############################################################################
# Model selection 2 setup: covariates
# Define covariate combinations
covar.combo <- c('seds', 'rugos', 'bathy', 'bottom', 'nao', 'amo')
covar.combo <- do.call("c", lapply(seq_along(covar.combo), 
                                   function(i) combn(covar.combo, i, FUN = list)))
for(i in 1:length(covar.combo)){
  temp <- covar.combo[[i]]
  use <- length(temp)
  temp2 <- temp[[1]]
  
  if(use > 1){
    for(j in 2:use){
      temp2 <- paste0(temp2, temp[[j]])
    }
  }

  covar.combo[[i]] <- temp2
  #rm(temp, use, temp2)
}
covar.combo <- do.call(rbind, covar.combo)

mod.covar <- c("base", covar.combo[1:63])
length(mod.covar)
# 64 total models


OverdispersionConfig	<- c("eta1"=0, "eta2"=0)
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight
# We have no effects for vessel so this doesn't matter
OverdispersionConfig1 <- c("eta1"=1, "eta2"=0)
OverdispersionConfig2 <- c("eta1"=1, "eta2"=1)

#########################################################
# Run model selection 2
# Subset to 2016-2021 for time
sel2 <- survs[survs$Year >= 71,]
# Scale covariate values
scaled.covars <- scaled.covars[scaled.covars$Year >=71,]

# Define density covariate formulas
covar.combo <- c('seds', 'rugos', 'BATHY.DEPTH', 'h_bt', 'nao', 'amo')
covar.combo <- do.call("c", lapply(seq_along(covar.combo), 
                                   function(i) combn(covar.combo, i, FUN = list)))
for(i in 1:length(covar.combo)){
  temp <- covar.combo[[i]]
  if("seds" %in% temp){
    temp <- c(temp, 'cobble_P', 'gravel_P', 'mud_P', 'sand_P')
  }
  temp <- temp[ !temp == 'seds']
  use <- length(temp)
  temp2 <- paste0( "~ ", temp[[1]])
  
  if(use > 1){
    for(j in 2:use){
      temp2 <- paste0(temp2, " + ", temp[[j]])
    }
  }
  
  covar.combo[[i]] <- as.formula(temp2)
  #rm(temp, use, temp2)
}
base.mod <- as.data.frame(NULL)
mod.Qik <- c(list(base.mod), covar.combo)
names(mod.Qik) <- mod.covar

# Loop through density covariate options
for(i in 2:length(mod.covar)) {
  message(names(mod.Qik)[i])
  # Define name of model
  name <- paste0(mod.covar[i])
  # Name working directory
  working_dir <- here::here(sprintf("Model_Refinement/Covariate_Selection2/%s/", name))
  # Make folder if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
    file.copy(from=here("Model_Refinement", "Covariate_Selection2", "seds",
                        "Kmeans_extrapolation-2000.RData"), 
              to=paste0(working_dir, '/', 
                        'Kmeans_extrapolation-2000.RData'), 
              overwrite = TRUE, 
              recursive = FALSE, 
              copy.mode = TRUE)
    file.copy(from=here("Model_Refinement", "Covariate_Selection2", "seds",
                        "Kmeans_knots-200.RData"), 
              to=paste0(working_dir, '/', 
                        'Kmeans_knots-200.RData'), 
              overwrite = TRUE, 
              recursive = FALSE, 
              copy.mode = TRUE)
  }
  
  if("parameter_estimates.txt" %in% list.files(working_dir)){
    print(paste0('Run already completed for ', names(mod.Qik)[i], ' model'))
    next()
  }
  
  # Set model options
  # winners of model selection 1
  use_anisotropy <- TRUE
  FieldConfig <- FieldConfig1
  OverdispersionConfig <- OverdispersionConfig
  Q_ik <- mod.Qik[[i]]
  # Make settings
  settings <- make_settings( n_x = 200,
                             Region = "User",
                             Version = "VAST_v14_0_1",
                             purpose = "index2",
                             bias.correct = FALSE,
                             use_anisotropy = use_anisotropy,
                             FieldConfig = FieldConfig,
                             RhoConfig = RhoConfig, # always default
                             OverdispersionConfig = OverdispersionConfig
  )
  # Set obsmodel to reflect CPUE response variable measure
  settings$ObsModel[1] <- 4
  settings$ObsModel[2] <- 1
  
  # Fit model
  fit = fit_model(
    # Call settings
    settings = settings,
    # Call survey data info
    Lat_i = sel2[,'Lat'], 
    Lon_i = sel2[,'Lon'], 
    t_i = sel2[,'Year'],
    b_i = sel2[,'Response_variable'],
    a_i = sel2[,'swept'],
    #v_i = sel2[,'vessel'],
    #c_iz = sel2[,'Age'],
    # Call covariate info
    X1_formula = mod.Qik[[i]],
    X2_formula = mod.Qik[[i]],
    covariate_data = scaled.covars,
    # Call spatial
    input_grid=user_region,
    # Set working dir
    working_dir = paste0(working_dir, "/"),
    # Set naming conventions
    #category_names = cat.labs,
    #year_labels = year.labs,
    
    # Tell model to run
    build_model = TRUE,
    run_model = TRUE)
  
  # if(exists("fit") == FALSE){
  #   settings$FieldConfig["Epsilon", "Component_2"] <- 0
  #   
  #   # Fit model
  #   fit = fit_model(
  #     # Call settings
  #     settings = settings,
  #     # Call survey data info
  #     Lat_i = sel2[,'Lat'], 
  #     Lon_i = sel2[,'Lon'], 
  #     t_i = sel2[,'Year'],
  #     b_i = sel2[,'Response_variable'],
  #     a_i = sel2[,'swept'],
  #     #v_i = sel2[,'vessel'],
  #     #c_iz = sel2[,'Age'],
  #     # Call covariate info
  #     X1_formula = mod.Qik[[i]],
  #     X2_formula = mod.Qik[[i]],
  #     covariate_data = scaled.covars,
  #     # Call spatial
  #     input_grid=user_region,
  #     # Set working dir
  #     working_dir = paste0(working_dir, "/"),
  #     # Set naming conventions
  #     #category_names = cat.labs,
  #     #year_labels = year.labs,
  #     
  #     # Tell model to run
  #     build_model = TRUE,
  #     run_model = TRUE)
  #   
  # }
  
  rm(fit, settings)
  
  
} # end covar loop

# Loop through base options (no covars)
for(i in c(1)) {
  # Define name of model
  name <- paste0(mod.covar[i])
  # Name working directory
  working_dir <- here::here(sprintf("Model_Refinement/Covariate_Selection2/%s/", name))
  # Make folder if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
  }
  # Set model options
  # winners of model selection 1
  use_anisotropy <- TRUE
  FieldConfig <- FieldConfig1
  OverdispersionConfig <- OverdispersionConfig
  Q_ik <- mod.Qik[[i]]
  # Make settings
  settings <- make_settings( n_x = 200,
                             Region = "User",
                             Version = "VAST_v14_0_1",
                             purpose = "index2",
                             bias.correct = FALSE,
                             use_anisotropy = use_anisotropy,
                             FieldConfig = FieldConfig,
                             RhoConfig = RhoConfig, # always default
                             OverdispersionConfig = OverdispersionConfig
  )
  # Set obsmodel to reflect CPUE response variable measure
  settings$ObsModel[1] <- 4
  settings$ObsModel[2] <- 1
  
  # Fit model
  fit = fit_model(
    # Call settings
    settings = settings,
    # Call survey data info
    Lat_i = sel2[,'Lat'], 
    Lon_i = sel2[,'Lon'], 
    t_i = sel2[,'Year'],
    b_i = sel2[,'Response_variable'],
    a_i = sel2[,'swept'],
    #v_i = sel2[,'vessel'],
    #c_iz = sel2[,'Age'],
    # Call covariate info
    #X1_formula = mod.Qik[[i]],
    #X2_formula = mod.Qik[[i]],
    #covariate_data = scaled.covars,
    # Call spatial
    input_grid=user_region,
    # Set working dir
    working_dir = paste0(working_dir, "/"),
    # Set naming conventions
    #category_names = cat.labs,
    #year_labels = year.labs,
    
    # Tell model to run
    build_model = TRUE,
    run_model = TRUE)

} # end covar loop
beep(8)
###############################################################################
# Model selection for covariates
# Set folder 
outdir <- here("Model_Refinement/Covariate_Selection2")
# List folders in outer folder
moddirs <- list.dirs(outdir) 
# Remove top level folder
moddirs <- moddirs[-c(1)]
# keep folder name
modnames <- mod.covar

# function to apply extracting info
getmodinfo <- function(d.name){
  # read settings
  modpath <- stringr::str_split(d.name, "/", simplify = TRUE)
  modname <- modpath[length(modpath)]
  
  settings <- read.table(file.path(d.name, "settings.txt"), comment.char = "",
                         fill = TRUE, header = FALSE)
  
  n_x <- as.numeric(as.character(settings[(which(settings[,1]=="$n_x")+1),2]))
  grid_size_km <- as.numeric(as.character(settings[(
    which(settings[,1]=="$grid_size_km")+1),2]))
  max_cells <- as.numeric(as.character(settings[(
    which(settings[,1]=="$max_cells")+1),2]))
  use_anisotropy <- as.character(settings[(
    which(settings[,1]=="$use_anisotropy")+1),2])
  fine_scale <- as.character(settings[(
    which(settings[,1]=="$fine_scale")+1),2])
  bias.correct <- as.character(settings[(
    which(settings[,1]=="$bias.correct")+1),2])
  
  #FieldConfig
  if(settings[(which(settings[,1]=="$FieldConfig")+1),1]=="Component_1"){
    omega1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+2),2])
    omega2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+3),1])
    epsilon1 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+4),2])
    epsilon2 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+5),1])
    beta1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+6),2])
    beta2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+7),1])
  }
  
  if(settings[(which(settings[,1]=="$FieldConfig")+1),1]=="Omega1"){
    omega1 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+3),1])
    omega2 <- as.character(settings[(which(settings[,1]=="$FieldConfig")+4),1])
    epsilon1 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+3),2])
    epsilon2 <- as.character(settings[(
      which(settings[,1]=="$FieldConfig")+4),2])
    beta1 <- "IID"
    beta2 <- "IID"
  }
  
  
  #RhoConfig
  rho_beta1 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+3),1]))
  rho_beta2 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+3),2]))
  rho_epsilon1 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+4),1]))
  rho_epsilon2 <- as.numeric(as.character(settings[(
    which(settings[,1]=="$RhoConfig")+4),2]))
  
  # read parameter estimates, object is called parameter_Estimates
  load(file.path(d.name, "parameter_estimates.RData"))
  
  AIC <- parameter_estimates$AIC[1]  
  converged <- parameter_estimates$Convergence_check[1]
  fixedcoeff <- unname(parameter_estimates$number_of_coefficients[2])
  randomcoeff <- unname(parameter_estimates$number_of_coefficients[3])
  
  
  # return model atributes as a dataframe
  out <- data.frame(modname = modname,
                    n_x = n_x,
                    grid_size_km = grid_size_km,
                    max_cells = max_cells,
                    use_anisotropy = use_anisotropy,
                    fine_scale =  fine_scale,
                    bias.correct = bias.correct,
                    omega1 = omega1,
                    omega2 = omega2,
                    epsilon1 = epsilon1,
                    epsilon2 = epsilon2,
                    beta1 = beta1,
                    beta2 = beta2,
                    rho_epsilon1 = rho_epsilon1,
                    rho_epsilon2 = rho_epsilon2,
                    rho_beta1 = rho_beta1,
                    rho_beta2 = rho_beta2,
                    AIC = AIC,
                    converged = converged,
                    fixedcoeff = fixedcoeff,
                    randomcoeff = randomcoeff
  )
  return(out)
}

# combine into one table for comparison
modselect <- purrr::map_dfr(moddirs, getmodinfo)

# Build table to compare models
modselect.cov <- modselect %>%
  filter(n_x == 200) %>%
  #filter(str_detect(modname, "base|eta|len|_no$")) %>%
  # mutate(converged2 = case_when(str_detect(converged, 
  #                                          "no evidence") ~ "likely",
  #                               str_detect(converged, 
  #                                          "is likely not") ~ "unlikely",
  #                               TRUE ~ as.character(NA))) %>%
  mutate(deltaAIC = AIC-min(AIC)) %>%
  select(modname, deltaAIC, fixedcoeff,
         randomcoeff, use_anisotropy, 
         omega1, omega2, epsilon1, epsilon2, 
         beta1, beta2, AIC) %>%
  arrange(AIC)

# Print table
DT::datatable(modselect.cov, rownames = FALSE, 
              options= list(pageLength = 25, scrollX = TRUE))

# Evidence supports that the best model (lowest AIC) is the one with all covars.