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

# Read in data
load(here('data/RData_Storage/shortcut_VAST_data_final.RData'))

###############################################################################
# Describe model selection efforts

# Model selection 1 (spatial, spatio-temporal effects, no covariates) options 

# Use sel1 dataset -- only NEFSC BTS data

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


# Model selection 2 (covariates) options, FieldConfig default (all IID)
# _base         No covariates
# _gravel       gravel_P covariate
# _cobble       cobble_P covariate
# _mud          mud_P covariate
# _sand         sand_P covariate
# _bathy        bathymetry covariate
# _oisst        sea surface temperature covariate
# _rugos        rugosity covariate
# _seds         gravel, cobble, mud, and sand covariates
# _sedsbath     gravel, cobble, mud, sand, and bathymetry covariates
# _sedssst      gravel, cobble, mud, sand, and sst covariates
# _sedsrug      gravel, cobble, mud, sand, and rugosity covariates
# _sedsbathsst  gravel, cobble, mud, sand, bathymetry, and sst covariates
# _sedsbathrug  gravel, cobble, mud, sand, bathymetry, and rugosity covariates
# _sedsrugsst   gravel, cobble, mud, sand, rugosity, and sst covariates
# _bathsst      bathymetry and sst covariates
# _bathrug      bathymetry and rugosity covariates
# _bathrugsst   bathymetry, rugosity, and sst covariates
# _rugsst       rugosity and sst covariates
# _all          all listed covariates
# _eta1         vessel overdispersion in 1st predictor
# _eta2         vessel overdispersion in 1st and 2nd predictors

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
sel1 <- subset(survs, vessel == 7)

# Loop through options
for(i in 2:length(mod.config)) {
  # Define name of model run
  name <- mod.config[i]
  # Set working directory
  working_dir <- here::here(sprintf("mod_selection/%s", name))
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

# Evidence suggests best model includes all effects and anisotropy.

###############################################################################
# Model selection 2 setup: covariates
# Define covariate combinations

mod.covar <- c("base", 
               "gravel", "cobble", "mud", "sand", "bathy", "sst", "rugos",
               "seds", "sedsbath", "sedsrug", "sedssst", "sedsbathrug", 
                        "sedsbathsst", "sedsrugsst",
               "bathrug", "bathsst", "bathrugsst", 
               "rugsst", "all",
               "eta10", "eta11",
               "sedsbathssteta11", "alleta11")

OverdispersionConfig	<- c("eta1"=0, "eta2"=0)
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

OverdispersionConfig1 <- c("eta1"=1, "eta2"=0)
OverdispersionConfig2 <- c("eta1"=1, "eta2"=1)

mod.eta <- list(OverdispersionConfig, OverdispersionConfig, 
                OverdispersionConfig, OverdispersionConfig, 
                OverdispersionConfig, OverdispersionConfig, 
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig, OverdispersionConfig,
                OverdispersionConfig1, OverdispersionConfig2, 
                OverdispersionConfig2, OverdispersionConfig2)

names(mod.eta) <- mod.covar

#########################################################
# Run model selection 2
# Subset to 10 years of data for time
sel2 <- subset(survs, Year >= 61)
dat <- subset(scaled.covars, Year >=61)

# Define density covariate formulas
Q_ikbase        <-   NULL
Q_ikgravel      <- ~ gravel_P
Q_ikcobble      <- ~ cobble_P
Q_ikmud         <- ~ mud_P
Q_iksand        <- ~ sand_P
Q_ikbathy       <- ~ BATHY.DEPTH
Q_iksst         <- ~ oisst
Q_ikrugos       <- ~ rugos
Q_ikseds        <- ~ gravel_P + cobble_P + mud_P + sand_P
Q_iksedsbath    <- ~ gravel_P + cobble_P + mud_P + sand_P + BATHY.DEPTH
Q_iksedsrug     <- ~ gravel_P + cobble_P + mud_P + sand_P + rugos
Q_iksedssst     <- ~ gravel_P + cobble_P + mud_P + sand_P + oisst
Q_iksedsbathrug <- ~ gravel_P + cobble_P + mud_P + sand_P + BATHY.DEPTH + 
                                rugos
Q_iksedsbathsst <- ~ gravel_P + cobble_P + mud_P + sand_P + BATHY.DEPTH + 
                                oisst
Q_iksedsrugsst  <- ~ gravel_P + cobble_P + mud_P + sand_P + rugos + oisst
Q_ikbathrug     <- ~ BATHY.DEPTH + rugos
Q_ikbathsst     <- ~ BATHY.DEPTH + oisst
Q_ikbathrugsst  <- ~ BATHY.DEPTH + rugos + oisst
Q_ikrugsst      <- ~ rugos + oisst
Q_ikall         <- ~ gravel_P + cobble_P + mud_P + sand_P + BATHY.DEPTH + 
                                rugos + oisst
Q_iksedsbathssteta11 <- ~ gravel_P + cobble_P + mud_P + sand_P + BATHY.DEPTH +
                                oisst
Q_ikalleta11 <- Q_ikall

# Pull formulas into list
mod.Qik <- list(Q_ikbase,
                Q_ikgravel, Q_ikcobble, Q_ikmud, Q_iksand, Q_ikbathy, 
                            Q_iksst, Q_ikrugos,
                Q_ikseds, Q_iksedsbath, Q_iksedsrug, Q_iksedssst,
                          Q_iksedsbathrug, Q_iksedsbathsst, Q_iksedsrugsst,
                Q_ikbathrug, Q_ikbathsst, Q_ikbathrugsst,
                Q_ikrugsst,
                Q_ikall,
                Q_ikbase, Q_ikbase,
                Q_iksedsbathssteta11, Q_ikalleta11)

# Name formula list items
names(mod.Qik) <- mod.covar

# Loop through density covariate options
for(i in 2:20) {
  # Define name of model
  name <- paste0(mod.covar[i])
  # Name working directory
  working_dir <- here::here(sprintf("covar_selection/%s/", name))
  # Make folder if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
  }
  # Set model options
  # winners of model selection 1
  use_anisotropy <- TRUE
  FieldConfig <- FieldConfig1
  OverdispersionConfig <- mod.eta[[i]]
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
    v_i = sel2[,'vessel'],
    c_iz = sel2[,'Age'],
    # Call covariate info
    X1_formula = mod.Qik[[i]],
    covariate_data = scaled.covars,
    # Call spatial 
    input_grid=user_region,
    # Set working dir
    working_dir = paste0(working_dir, "/"),
    # Tell model to run
    run_model = TRUE)
} # end covar loop

# Loop through base options (no covars)
for(i in c(1, 21, 22)) {
  # Define name of model
  name <- paste0(mod.covar[i])
  # Name working directory
  working_dir <- here::here(sprintf("covar_selection/%s/", name))
  # Make folder if it doesn't exist
  if(!dir.exists(working_dir)) {
    dir.create(working_dir)
  }
  # Set model options
  # winners of model selection 1
  use_anisotropy <- TRUE
  FieldConfig <- FieldConfig1
  OverdispersionConfig <- mod.eta[[i]]
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
    v_i = sel2[,'vessel'],
    c_iz = sel2[,'Age'],
    # Call covariate info
    #X1_formula = mod.Qik[[i]],
    #covariate_data = scaled.covars,
    # Call spatial 
    input_grid=user_region,
    # Set working dir
    working_dir = paste0(working_dir, "/"),
    # Tell model to run
    run_model = TRUE)
} # end covar loop
beep(8)
###############################################################################
# Model selection for covariates
# Set folder 
outdir <- here("covar_selection")
# List folders in outer folder
moddirs <- list.dirs(outdir) 
# Remove top level folder
moddirs <- moddirs[-c(1,3)]
# keep folder name
modnames <- c('all', 'alleta11' ,'base', 'bathrug', 'bathrugsst', 'bathsst', 
              'bathy', 
              'cobble', 'eta10', 'eta11', 'gravel', 'mud', 'rugos', 'rugsst',
              'sand', 'seds', 'sedsbath', 'sedsbathrug', 'sedsbathsst', 
              'sedsbathssteta11',
              'sedsrug', 'sedsrugsst', 'sedssst', 'sst')

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
