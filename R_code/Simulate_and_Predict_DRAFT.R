# Dummy exercise to simulate capture data and new locations based on
# a fitted VAST model

# In this case, model spatial domain is Western Gulf of Maine
# Modeled species is medium-sized cod
# Temporal domain is 2000-2010
# Strata are arbitrary
# Model structure is similar to final overall cod model, as allowed

rm(list=ls())
set.seed(123)

# Libraries
library(here)
library(sf)
library(tidyverse)
library(VAST)

# Add unitless back as possible unit (removed in units package update Mar 2023)
install_unit(symbol='unitless', def='unitless', name='unitless')

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, 
                                            fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

#### Add sample information and covars ####
# Load data
surveys <- read.csv(here('Data/Prey_Data/species_sep_prey.csv'))

# Remove uncommon gears
surveys <- surveys[surveys$gear %in% c('otter trawl'),]

# Remove points not in model domain (on land)
region_shape <- st_read(here('Data/GIS/NAFO_Continental_Shelf_10kmBuff.shp'), quiet=T)
region_shape <- region_shape %>% 
  filter(#Label != '6B',
    Region == 'US')
region_shape <- as(region_shape, 'Spatial')
region_shape <- st_as_sf(spTransform(region_shape, CRS("+init=epsg:4326")))
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, Region, Label, geometry)

surveys_sf <- st_as_sf(surveys, coords=c('lon', 'lat'))
st_crs(surveys_sf) <- 'EPSG:4326'

surveys_sf <- st_intersection(surveys_sf, region_shape)
surveys <- sfheaders::sf_to_df(surveys_sf, fill=TRUE)
surveys <- surveys %>%
  rename(lon =x) %>%
  rename(lat =y)

#### Finalize sampling data inputs ####
# Save sampling data
survs <- surveys

# Set variables
survs$RESPONSE <- as_units(survs$wt_kgs, 'kg')
survs$swept <- as_units(1, 'km')
#survs$vessel <- as.numeric(as.factor(survs$gear)) - 1
# Gillnet 0, Otter trawl 1, Pair trawl 2
survs$cat <- as.numeric(as.factor(survs$fish)) -1
survs$survmonths <- paste0(survs$year, '-', str_pad(survs$month, 2, 'left', '0'))

survmonths <- c(unique(survs$survmonths), '2020-06', '2020-07',
                '2023-06', '2023-07', '2023-08', '2023-09', '2023-10')
survmonths <- survmonths[order(survmonths)]

survmonths <- as.data.frame(survmonths)
survmonths$time <- seq(1:nrow(survmonths))

survs <- left_join(survs, survmonths, by='survmonths')

# herring 0, mackerel 1, menhaden 2
survs <- dplyr::select(survs, lon, lat, time, RESPONSE, cat, 
                       #vessel, 
                       swept)
names(survs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 
                  'category',
                  #'vessel',
                  'swept')

#### Output success ####
print('Data inputs made successfully')

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape <- st_union(region_shape)
region_shape <- st_sf(data=region_shape)
region_shape$OBJECTID <- 'ALL'
region_shape <- dplyr::select(region_shape, OBJECTID, data)
colnames(region_shape) <- c("Region", 'geometry')
st_geometry(region_shape) <- 'geometry'
index_area_shapes = region_shape

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("ALL"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c("scaled.covars", "settings", 'strata_use', 'survs',
                        'vast_extrap_grid', 'covars', 'surveys')))

working_dir <- paste0(here::here("VAST_runs/medium/DummyStrat/"))
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
                                           max_cells = 1000
)

message('Grid made')

#### Fit original model ####

# Make settings
settings <- make_settings(
  n_x = 100,
  purpose = "index2",
  Region = "User",
  fine_scale = FALSE,
  bias.correct = FALSE,
  knot_method = "grid",
  max_cells = 1000,
  strata.limits = data.frame(STRATA = 'ALL')
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
# Value 2: Copy values from first predictor to second predictor
# settings$RhoConfig[3] <- 0
# settings$RhoConfig[4] <- 0

# Turn on anisotropy
settings$use_anisotropy <- TRUE

# Density covariate model formula
hab_formula<- ~ 
  #bs(BATHY.DEPTH, degree = 3, intercept = FALSE) +
  bs(h_bt,        degree = 3, intercept = FALSE)

fit_strat = fit_model( 
  
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
  
  # Call covariate info
  X1_formula = hab_formula,
  X2_formula = hab_formula,
  covariate_data = covars,
  
  # Call spatial 
  input_grid=vast_extrap_grid,
  "extrapolation_list" = extrap_info_aja,
  
  # Tell model to build but not run
  build_model = TRUE,
  run_model = TRUE
  
)

# Save results
save(fit_strat, file=paste0(working_dir, '/DummyData_Stratified.RData'))
rm(list=setdiff(ls(), c('fit_strat', 'working_dir', 'surveys')))

#### Create new randomly-selected points ####
# Let's try to make proprtion of points in each strata equal to the actual
# proportion of that strata's area to the whole.
# A = 56, B=26, C=18
rats <- c(0.56, 0.26, 0.18)

# Let's keep the same number of samples per year
# Find number of samples in each strata per year
nefsc <- surveys[surveys$SURVEY == 'NEFSC BTS',]
dat.list <- split(nefsc, 
                  f=nefsc$YEAR)
years <- unique(nefsc$YEAR)

fillabledat <- data.frame(
  AGE_N = NA,
  YEAR = NA,
  AREA_SWEPT = NA,
  DUMMYSTRAT = NA,
  LON = NA,
  LAT = NA
)

for(i in 1:length(dat.list)){
  dat.list[[i]] <- split(dat.list[[i]], f=dat.list[[i]]$DUMMYSTRAT)
  n.per <- c(nrow(dat.list[[i]]$A), 
             nrow(dat.list[[i]]$B),
             nrow(dat.list[[i]]$C))
  n.tot <- sum(n.per)
  n.per <- round(rats * n.tot)
  
  yearkeeper <- data.frame(
    AGE_N = rep(NA, sum(n.per)),
    YEAR = rep(years[i], sum(n.per)),
    AREA_SWEPT = rep(0.0384,
                     sum(n.per)),
    DUMMYSTRAT = c(rep('A', n.per[1]),
                   rep('B', n.per[2]),
                   rep('C', n.per[3])),
    LON = rep(NA, sum(n.per)),
    LAT = rep(NA, sum(n.per))
    
  )
  fillabledat <- rbind(fillabledat, yearkeeper)
  
  rm(n.per, n.tot, yearkeeper)
  
}
fillabledat <- fillabledat[!is.na(fillabledat$DUMMYSTRAT),]

A <- st_read(here('Data/GIS/DummyStrat/A_UTM.shp'), quiet=T)
A <- st_transform(A, 'EPSG:4326')
A <- st_make_valid(A)
B <- st_read(here('Data/GIS/DummyStrat/B_UTM.shp'), quiet=T)
B <- st_transform(B, 'EPSG:4326')
B <- st_make_valid(B)
C <- st_read(here('Data/GIS/DummyStrat/C_UTM.shp'), quiet=T)
C <- st_transform(C, 'EPSG:4326')
C <- st_make_valid(C)

A.samps <- st_sample(A, nrow(fillabledat[fillabledat$DUMMYSTRAT == 'A',]))
A.samps <- st_sf(A.samps)
A.samps <- sfheaders::sf_to_df(A.samps, fill = T)
fillabledat$LON[fillabledat$DUMMYSTRAT == 'A'] <- A.samps$x
fillabledat$LAT[fillabledat$DUMMYSTRAT == 'A'] <- A.samps$y

B.samps <- st_sample(B, nrow(fillabledat[fillabledat$DUMMYSTRAT == 'B',]))
B.samps <- st_sf(B.samps)
B.samps <- sfheaders::sf_to_df(B.samps, fill = T)
fillabledat$LON[fillabledat$DUMMYSTRAT == 'B'] <- B.samps$x
fillabledat$LAT[fillabledat$DUMMYSTRAT == 'B'] <- B.samps$y

C.samps <- st_sample(C, nrow(fillabledat[fillabledat$DUMMYSTRAT == 'C',]))
C.samps <- st_sf(C.samps)
C.samps <- sfheaders::sf_to_df(C.samps, fill = T)
fillabledat$LON[fillabledat$DUMMYSTRAT == 'C'] <- C.samps$x
fillabledat$LAT[fillabledat$DUMMYSTRAT == 'C'] <- C.samps$y

fsf <- st_as_sf(fillabledat, coords=c('LON', 'LAT'), crs="EPSG:4326")

rm(list=setdiff(ls(), c('fit_strat', 'surveys', 'fsf', 'working_dir', 'fillabledat')))

#### Find covariate values at new locations ####
# Find most likely times of year for NEFSC BTS
nefsc <- surveys[surveys$SURVEY == 'NEFSC BTS',]
nefsc$DATE <- as.Date(nefsc$DATE)
table(as.Date(paste0('2011-', month(nefsc$DATE), '-', day(nefsc$DATE))))
# Assume 4 tows per day.
# Can be any day April, May / October, November
# Half in spring, half in fall.

springdates <- seq(as.Date("2011-04-01"), as.Date("2011-05-31"), by="days")
springdates <- rep(springdates, 4);springdates <- springdates[order(springdates)]
falldates <- seq(as.Date("2011-10-01"), as.Date("2011-11-30"), by='days')
falldates <- rep(falldates, 4);falldates <- falldates[order(falldates)]
anchor <- st_as_sf(data.frame(x=-70, y=39.9), coords=c('x', 'y'), crs="EPSG:4326")
anchor <- st_transform(anchor, st_crs())

datlist <- split(fillabledat, f=fillabledat$YEAR)
for(i in 1:length(datlist)){
  more <- sample(c('spring', 'fall'), 1)
  n.tot <- nrow(datlist[[i]])
  
  if(more == 'spring'){
    nspring <- ceiling(n.tot/2)
  }
  
  if(more == 'fall'){
    nspring <- floor(n.tot/2)
  }
  datlist[[i]]$SEASON <- NA
  
  springones <- sample(1:nrow(datlist[[i]]), nspring, replace = F)
  
  datlist[[i]]$SEASON[springones] <- 'SPRING'
  datlist[[i]]$SEASON[is.na(datlist[[i]]$SEASON)] <- 'FALL'
  
  datlist[[i]] <- split(datlist[[i]], f=datlist[[i]]$SEASON)
  
  for(j in 1:length(datlist[[i]])){
    datlist[[i]][[j]] <- st_as_sf(datlist[[i]][[j]], coords=c('LON', 'LAT'),
                                  crs="EPSG:4326")
    
    datlist[[i]][[j]]$DIST <- st_distance(
      x= datlist[[i]][[j]],
      y=anchor
    )
    
    datlist[[i]][[j]] <- datlist[[i]][[j]][
      with(datlist[[i]][[j]],
           order(DIST)),
    ]
    
    if(datlist[[i]][[j]]$SEASON[1] == 'SPRING'){
      datlist[[i]][[j]]$DATE <- springdates[1:nrow(datlist[[i]][[j]])]
      datlist[[i]][[j]]$DATE <- as.Date(paste0(datlist[[i]][[j]]$YEAR[1], 
                                               substr(datlist[[i]][[j]]$DATE,
                                                      start=5, stop=10)))
    }
    
    if(datlist[[i]][[j]]$SEASON[1] == 'FALL'){
      datlist[[i]][[j]]$DATE <- falldates[1:nrow(datlist[[i]][[j]])]
      datlist[[i]][[j]]$DATE <- as.Date(paste0(datlist[[i]][[j]]$YEAR[1], 
                                               substr(datlist[[i]][[j]]$DATE,
                                                      start=5, stop=10)))
    }

  }
  datlist[[i]] <- do.call(rbind, datlist[[i]])
  rm(more, n.tot, nspring, springones)
}
datlist <- do.call(rbind, datlist)
datlist <- datlist[with(datlist, order(DATE)),]
rownames(datlist) <- NULL

# Find covariates at these locations
# Add bottom temp
# Sf to df
stations <- sfheaders::sf_to_df(datlist, fill=T)
stations <- stations %>% 
  dplyr::select(-sfg_id, -point_id) %>% 
  rename(LON = x, LAT = y)
stations$matching <- paste0('X',
                            stations$YEAR, 
                            '.',
                            str_pad(lubridate::month(stations$DATE), 2, 'left', '0'),
                            '.',
                            str_pad(lubridate::day(stations$DATE), 2, 'left', '0'))
stations$h_bt <- NA
stations$HAUL_ID <- seq(1:nrow(stations))

# Set years
years <- 2011:2019

# Loop through years
library(raster)
for(i in years){
  message(i)
  brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                             'hubert', 'extended_grd_krig'),
                        '/', i, '.grd'))
  
  dates <- seq(as.Date(paste0(i, "-01-01")),
               as.Date(paste0(i, "-12-31")),by="1 day")
  
  if(length(dates) < nlayers(brick)){
    dates <- c(dates, NA)
  }
  
  names(brick) <- dates
  
  stations.yr <- subset(stations, stations$YEAR == paste0(i))
  stations.yr <- dplyr::select(stations.yr,
                               HAUL_ID, YEAR, matching, LON, LAT)
  
  
  for(q in 1:nrow(stations.yr)){
    print(round(q/ nrow(stations.yr) * 100),1)
    stations.yr$h_bt[q] <- extract(brick[[stations.yr[q,"matching"]]], 
                                   cbind(stations.yr[q,"LON"], 
                                         stations.yr[q,"LAT"]))
  }
  #head(stations.yr)
  
  #stations.yr <- dplyr::select(stations.yr, HAUL_ID, h_bt)
  
  stations$h_bt[stations$HAUL_ID %in% stations.yr$HAUL_ID] <- 
    stations.yr$h_bt
  
  rm(stations.yr, dates, q, brick)
  
}
head(stations)
ns <- stations[is.na(stations$h_bt),]
table(ns$YEAR)

# Cool. Now we can predict using this new data.
rm(list=setdiff(ls(), c('surveys', 'stations', 'fit_strat', 'working_dir')))

#### Predict ####
# Finalize survey data input
survs <- dplyr::select(stations,
                       LON, LAT, AREA_SWEPT,
                       YEAR, AGE_N)

survs$RESPONSE <- as_units(survs$AGE_N, 'counts')
survs$swept <- as_units(survs$AREA_SWEPT, 'km^2')

survs <- dplyr::select(survs, LON, LAT, YEAR, RESPONSE, 
                       swept)
names(survs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 
                  'swept')

# Save covariates
covars <- dplyr::select(stations,
                        LON, LAT, YEAR, 
                        h_bt)
# Remove unknowns
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:ncol(covars)])

message('Data inputs made')

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape <- st_read(here('Data/GIS/WGOM_UTM.shp'), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape$OBJECTID <- 'ALL'
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
colnames(region_shape) <- c("Region", 'geometry')
index_area_shapes = region_shape

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("ALL"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c('surveys', 'stations', 'fit_strat', 'working_dir', 
                        'survs', 'settings', 'strata_use', 'vast_extrap_grid',
                        'covars')))

working_dir <- paste0(here::here("VAST_runs/medium/DummyStrat/"))
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
                                           max_cells = 1000
)

message('Grid made')

#### Fit new predictive model ####
survs$PredTF_i <- 1
covars$PredTF_i <- 1

oldsurvs <- dplyr::select(surveys, LON, LAT, YEAR, AGE_N, AREA_SWEPT)
oldsurvs$PredTF_i <- 0
colnames(oldsurvs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 'swept', 'PredTF_i')
oldsurvs$Response_variable <- as_units(oldsurvs$Response_variable, 'counts')
oldsurvs$swept <- as_units(oldsurvs$swept, 'km^2')

oldcovs <- dplyr::select(surveys, LON, LAT, YEAR, h_bt)
oldcovs$PredTF_i <- 0
colnames(oldcovs) <- c('Lon', 'Lat', 'Year', 'h_bt', 'PredTF_i')

survs <- rbind(oldsurvs, survs)
covars <- rbind(oldcovs, covars)

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape <- st_read(here('Data/GIS/WGOM_UTM.shp'), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape$OBJECTID <- 'ALL'
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
colnames(region_shape) <- c("Region", 'geometry')
index_area_shapes = region_shape

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("ALL"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c('surveys', 'stations', 'fit_strat', 'working_dir', 
                        'survs', 'settings', 'strata_use', 'vast_extrap_grid',
                        'covars')))

working_dir <- paste0(here::here("VAST_runs/medium/DummyStrat/"))
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
                                           max_cells = 1000
)

message('Grid made')

# Make settings
settings <- make_settings(
  n_x = 100,
  purpose = "index2",
  Region = "User",
  fine_scale = FALSE,
  bias.correct = FALSE,
  knot_method = "grid",
  max_cells = 1000,
  strata.limits = data.frame(STRATA = 'ALL')
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
# Value 2: Copy values from first predictor to second predictor
# settings$RhoConfig[3] <- 0
# settings$RhoConfig[4] <- 0

# Turn on anisotropy
settings$use_anisotropy <- TRUE

# Density covariate model formula
hab_formula<- ~ 
  #bs(BATHY.DEPTH, degree = 3, intercept = FALSE) +
  bs(h_bt,        degree = 3, intercept = FALSE)

fit_pred = fit_model( 
  
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
  PredTF_i = survs[,'PredTF_i'],
  
  # Call covariate info
  X1_formula = hab_formula,
  X2_formula = hab_formula,
  covariate_data = covars,
  
  # Call spatial 
  input_grid=vast_extrap_grid,
  "extrapolation_list" = extrap_info_aja,
  
  # Tell model to build but not run
  build_model = TRUE,
  run_model = TRUE
  
)

# Save results
save(fit_pred, file=paste0(working_dir, '/DummyData_Predicted.RData'))
rm(list=setdiff(ls(), c('fit_pred', 'working_dir', 'surveys')))

#### Simulate new data based on predictive model ####
# Simulated points indexed at 3404 to 5574 in data frame
survs <- fit_pred$data_frame[3404:5574,]
covars <- fit_pred$covariate_data[3404:5574,]
covars$PredTF_i <- NULL

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape <- st_read(here('Data/GIS/WGOM_UTM.shp'), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape$OBJECTID <- 'ALL'
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
colnames(region_shape) <- c("Region", 'geometry')
index_area_shapes = region_shape

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c("ALL"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c('surveys', 'stations', 'fit_strat', 'working_dir', 
                        'survs', 'settings', 'strata_use', 'vast_extrap_grid',
                        'covars')))

working_dir <- paste0(here::here("VAST_runs/medium/DummyStrat/"))
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
                                           max_cells = 1000
)

message('Grid made')

# Make settings
settings <- make_settings(
  n_x = 100,
  purpose = "index2",
  Region = "User",
  fine_scale = FALSE,
  bias.correct = FALSE,
  knot_method = "grid",
  max_cells = 1000,
  strata.limits = data.frame(STRATA = 'ALL')
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
# Value 2: Copy values from first predictor to second predictor
# settings$RhoConfig[3] <- 0
# settings$RhoConfig[4] <- 0

# Turn on anisotropy
settings$use_anisotropy <- TRUE

# Density covariate model formula
hab_formula<- ~ 
  #bs(BATHY.DEPTH, degree = 3, intercept = FALSE) +
  bs(h_bt,        degree = 3, intercept = FALSE)

fit_sim = fit_model( 
  
  # Call settings
  settings = settings,
  
  # Set working directory (intermediate files)
  working_dir = getwd(), 
  
  # Turn on joint precision (need for range edges)
  getJointPrecision = TRUE,
  
  # Call survey data info
  Lat_i = survs[,'Lat_i'], 
  Lon_i = survs[,'Lon_i'], 
  t_i = survs[,'t_i'],
  b_i = survs[,'b_i'],
  a_i = survs[,'a_i'],
  
  # Call covariate info
  X1_formula = hab_formula,
  X2_formula = hab_formula,
  covariate_data = covars,
  
  # Call spatial 
  input_grid=vast_extrap_grid,
  "extrapolation_list" = extrap_info_aja,
  
  # Tell model to build but not run
  build_model = TRUE,
  run_model = TRUE
  
)

# Save results
save(fit_sim, file=paste0(working_dir, '/DummyData_Simulated.RData'))
rm(list=setdiff(ls(), c('fit_sim', 'working_dir', 'surveys')))

plot_biomass_index(fit=fit_sim,
                   DirName = paste0(getwd(), '/Sim'),
                   year_labels = fit_sim$year_labels,
                   PlotName = 'Index_Sim')
