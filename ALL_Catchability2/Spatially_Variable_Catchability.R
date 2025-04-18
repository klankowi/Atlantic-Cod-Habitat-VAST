#### Cod VAST Modeling ####
#### Medium size only ####
#### No strata divisions ####
#### Fine scale off ####
#### Whole time series ####
#### Purpose: overall model with strata ####


#### Workspace setup ####
# Clear workspace
rm(list=ls())

# Set wd
setwd("~/VAST_Cod")

# Load packages
suppressPackageStartupMessages(library(TMB))
suppressPackageStartupMessages(library(units, quietly=T, verbose=F))
suppressPackageStartupMessages(library(VAST))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse, quietly=T, verbose=F))
suppressPackageStartupMessages(library(sf, quietly=T, verbose=F))
suppressPackageStartupMessages(library(sp, quietly=T, verbose=F))
suppressPackageStartupMessages(library(raster, quietly=T, verbose=F))
suppressPackageStartupMessages(library(splines, quietly=T, verbose=F))
suppressPackageStartupMessages(library(INLAspacetime, quietly=T, verbose=F))

# Add unitless back as possible unit (removed in units package update Mar 2023)
install_unit(symbol='unitless', def='unitless', name='unitless')

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

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Add sample information and covars ####
# Load data
surveys <- read.csv(here("Data/Survey_Data/Bio_Data_Agesep.csv"))
surveys$month <- month(as.POSIXct(surveys$DATE, format='%Y-%m-%d'))

# Only mediums
surveys <- surveys[surveys$AGEGROUP == 'Age2-5',]

# Add seasonality
surveys$season[surveys$month %in% c(1,2,9,10,11,12)] <- 'BFall'
surveys$season[surveys$month %in% c(3,4,5,6,7,8)] <- 'ASpring'

surveys$YEAR[surveys$month %in% c(1,2)] <- 
  surveys$YEAR[surveys$month %in% c(1,2)]-1

surveys$SEASON <- paste0(surveys$YEAR, ' ', surveys$season)

# Fix stupid error
surveys$TIME <- as.numeric(as.factor(paste0(surveys$YEAR, ' ',
                                            surveys$SEASON)))

# Remove points not in model domain (on land)
region_shape<- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)

surveys_sf <- st_as_sf(surveys, coords=c('LON', 'LAT'))
st_crs(surveys_sf) <- 'EPSG:4326'

surveys_sf <- st_intersection(surveys_sf, region_shape)

surveys <- surveys[surveys$HAUL_ID %in% surveys_sf$HAUL_ID,]

# Test for missing values
surveys <- surveys[!is.na(surveys$rugos) &
                     !is.na(surveys$nao) &
                     !is.na(surveys$h_bt),]

#### Finalize sampling data inputs ####
# Save sampling data
survs <- dplyr::select(surveys,
                       LON, LAT, AREA_SWEPT,
                       TIME, SURVEY, AGE_N)

survs$RESPONSE <- as_units(survs$AGE_N, 'counts')
survs$swept <- as_units(survs$AREA_SWEPT, 'km^2')
survs$vessel <- as.numeric(as.factor(survs$SURVEY)) - 1
# vessel    survey
# 0         ASMFC Shrimp Trawl  
# 1         DFO Trawl  
# 2         MADMF Industry  
# 3         MADMF Inshore Trawl  
# 4         ME-NH Inshore Trawl  
# 5         NEFSC BLLS   
# 6         NEFSC BTS 
# 7        SMAST Video Trawl   

survs <- dplyr::select(survs, LON, LAT, TIME, RESPONSE, 
                       vessel, swept)
names(survs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 
                  'vessel', 'swept')

# Save covariates
covars <- dplyr::select(surveys,
                        LON, LAT, TIME, 
                        cobble_P, gravel_P, mud_P, sand_P, 
                        rugos, BATHY.DEPTH, h_bt,
                        nao, amo, AGEGROUP)
# Remove unknowns
covars <- subset(covars, AGEGROUP =='Age2-5')

covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] <- 
  covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] * -1
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:ncol(covars)])

# Rescale covariates to have mean 0 and SD 1 (author rec)
scaled.covars <- covars[,4:ncol(covars)] %>% 
  mutate(across(where(is.numeric), scale))
scaled.covars <- cbind(covars[,1:3], scaled.covars)

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

message('Data inputs made')

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape$OBJECTID <- 'ALL'
colnames(region_shape) <- c("Region", 'geometry')

# Second, get our index area shapefile
# We could just use this same shapefile in the "index_shapes" argument, but to 
# show off the new functions we wrote, we will also want to have a sf multipolygon
# shapefiles with areas defined within this general region
index_areas<- c("WGOM", "EGOM", "GBK", "SNE")

for(i in seq_along(index_areas)){
  index_area_temp<- st_read(paste0(here("Data/GIS"), '/', 
                                   index_areas[i], "_UTM.shp"), quiet=T)
  index_area_temp <- st_transform(index_area_temp, "EPSG:4326")
  index_area_temp <- st_make_valid(index_area_temp)
  index_area_temp <- dplyr::select(index_area_temp, STOCK, geometry)
  colnames(index_area_temp) <- c("Region", 'geometry')
  
  if(i == 1){
    index_area_shapes<- index_area_temp
  } else {
    index_area_shapes<- bind_rows(index_area_shapes, index_area_temp)
  }
}

index_area_shapes <- bind_rows(region_shape, index_area_shapes)
index_area_shapes <- st_make_valid(index_area_shapes)

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("ALL", "WGOM", "GBK", "EGOM", "SNE"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c("scaled.covars", "settings", 'strata_use', 'survs',
                        'vast_extrap_grid', 'covars', 'surveys')))

working_dir <- paste0(here::here("VAST_runs/medium/Overall_BC/ALL_Catchability2/"))
# Create working directory if it doesn't exist
if(!dir.exists(working_dir)) {
  dir.create(working_dir)
}

setwd(working_dir)

source(here('R_code/utilities/vast_function_edits.R'))
use_edited_funs<- TRUE
if (use_edited_funs) {
  source(here::here("R_code/utilities/vast_function_edits.R"))
  assignInNamespace("match_strata_fn", 
                    match_strata_fn, 
                    "FishStatsUtils")
  assignInNamespace("Prepare_User_Extrapolation_Data_Fn", 
                    Prepare_User_Extrapolation_Data_Fn, 
                    "FishStatsUtils")
}

extrap_info_aja <- make_extrapolation_info(Region="User",
                                           strata.limits=strata_use,
                                           input_grid = vast_extrap_grid,
                                           max_cells = 2000
)

message('Grid made')

# Remove intermediates
rm(list=setdiff(ls(), c("scaled.covars", "settings", 'strata_use', 'survs',
                        'vast_extrap_grid', 'extrap_info',
                        'working_dir', 'extrap_info_aja', 'covars',
                        'surveys')))

#### Run model ####
# Set year labels
year.labs <- c(seq(1982, 2021, 1), seq(1982, 2021, 1))
year.labs <- year.labs[order(year.labs)]
seas.labs <- rep(c('Spring', 'Fall'), 40)
year.labs <- paste0(year.labs, " ", seas.labs)

# Make settings
settings <- make_settings(
  n_x = 200,
  purpose = "index2",
  Region = "User",
  fine_scale = TRUE,
  bias.correct = TRUE,
  knot_method = "grid",
  max_cells = 2000,
  strata.limits = strata_use
)

# Change `ObsModel` settings
# Value 1: Delta-lognormal using bias-corrected-mean and log-coeff of var
# Value 2: Poisson-link delta-model" using log-link for numbers density and 
#          log-link for biomass per number
settings$ObsModel[1] <- 4
settings$ObsModel[2] <- 1

# Change `RhoConfig` settings
# Value 1: The spatio-temporal variation structure among time intervals
# is a random effect following a first-order autoregressive process,
# thus estimating variance as fixed effects and a separate first-order
# AR parameter for each category
# Value 2: Same
settings$RhoConfig[3] <- 4
settings$RhoConfig[4] <- 4

# Density covariate model formula
hab_formula<- ~ 
  ns(amo,         df = 3, intercept = FALSE) +
  ns(BATHY.DEPTH, df = 3, intercept = FALSE) +
  ns(h_bt,        df = 3, intercept = FALSE) +
  ns(cobble_P,    df = 1, intercept = FALSE) +
  ns(gravel_P,    df = 3, intercept = FALSE) + 
  ns(mud_P,       df = 3, intercept = FALSE) +
  ns(nao,         df = 3, intercept = FALSE) +
  ns(rugos,       df = 3, intercept = FALSE) + 
  ns(sand_P,      df = 2, intercept = FALSE)

# Catchabiility covariate model formula
surveys$Type <- "0" 
surveys$Type <- ifelse( surveys$SURVEY != "NEFSC BTS", 
                        as.numeric(as.factor(surveys$SURVEY)),
                        surveys$Type ) 
surveys$Type <- as.factor( surveys$Type)
table( surveys$Type, surveys$SURVEY  )
Q1_formula <- ~ Type
Q1config_k <- c(3,3,3,3,3,3,3,3,3) 
Q2config_k <- NULL
Q2_formula <- NULL
catchability_data <- surveys %>% 
  rename(Lon = LON,
         Lat = LAT,
         Year = TIME) %>% 
  dplyr::select(Lon, Lat, Year, Type)

# Turn on anisotropy
settings$use_anisotropy <- TRUE

fit = fit_model( 
  
  # Call settings
  settings = settings,
  
  # Set working directory (intermediate files)
  working_dir = getwd(), 
  
  # Turn on joint precision (need for range edges)
  getJointPrecision = TRUE,
  
  # Call survey data info
  Lat_i = survs[,'Lat'], 
  Lon_i = survs[,'Lon'], 
  t_i = survs[,'Year'],
  b_i = survs[,'Response_variable'],
  a_i = survs[,'swept'],
  v_i = survs[,'vessel'],
  
  # Call covariate info
  X1_formula = hab_formula,
  X2_formula = hab_formula,
  covariate_data = covars,
  
  # Call catchability info
  "Q1_formula" = Q1_formula,
  #"Q2_formula" = Q2_formula,  
  "catchability_data" = catchability_data,
  "Q1config_k" = Q1config_k,
  #"Q2config_k" = Q2config_k,
  
  # Call spatial 
  input_grid=vast_extrap_grid,
  "extrapolation_list" = extrap_info_aja,
  
  # Set naming conventions
  year_labels = year.labs,
  
  # Tell model to build but not run
  build_model = TRUE,
  run_model = TRUE
  
)

#### Save results ####

save.image(paste0('Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL_Catchability.RData'))
