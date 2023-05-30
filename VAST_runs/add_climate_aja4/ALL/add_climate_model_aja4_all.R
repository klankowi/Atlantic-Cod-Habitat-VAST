#### Workspace setup ####
# Clear workspace
rm(list=ls())

# Load libraries
# IMPORTANT NOTE: VAST must be running >=V14, will not work with V13.
library(TMB)
library(units)
library(VAST)
library(here)
library(tidyverse)
library(beepr)
library(sf)
library(rgdal)
library(sp)
library(ggcorrplot)
library(splines)  # Used to include basis-splines
library(INLAspacetime) # used for mesh-building
# Add unitless back as possible unit (removed in units package update Mar 2023)
#install_unit(symbol='unitless', def='unitless', name='unitless')
# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
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
surveys <- read.csv(here("Data/VAST_input/cod_agesep_VASTdata.csv"))

#### Finalize sampling data inputs ####
# Save sampling data
survs <- dplyr::select(surveys,
                       LON, LAT, AREA_SWEPT,
                       TIME, SURVEY, RESPONSE, 
                       AGEGROUP)

# Remove unknown size class, makes model not converge
survs <- subset(survs, AGEGROUP != "Unknown")

survs$RESPONSE <- as_units(survs$RESPONSE, 'counts')
survs$swept <- as_units(survs$AREA_SWEPT, 'km^2')
survs$vessel <- as.numeric(as.factor(survs$SURVEY)) - 1
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
# 9        SMAST Video Trawl   

#survs$Data_type <- factor(survs$Data_type, levels=c("Count", "Biomass_KG"))

survs$AGEGROUP <- as.numeric(factor(survs$AGEGROUP, levels=c(
  'Age0-2', 
  'Age2-5', 
  'Age5+'
  )
  )) - 1
# 0-2: 0
# 2-5: 1
# 5+ : 2

survs <- dplyr::select(survs, LON, LAT, TIME, RESPONSE, 
                       AGEGROUP, 
                       vessel, swept)
names(survs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 
                  'Age', 
                  'vessel', 'swept')
str(survs)

# Save covariates
covars <- dplyr::select(surveys,
                        LON, LAT, TIME, 
                        cobble_P, gravel_P, mud_P, sand_P, 
                        rugos, BATHY.DEPTH, h_bt,
                        nao, amo, AGEGROUP)
# Remove unknowns
covars <- subset(covars, AGEGROUP !='Unknown')

covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] <- 
  covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] * -1
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:ncol(covars)])
table(covars$Year)

# Test correlation
# Create correlation matrix
df_cormat <- dplyr::select(covars,                         
                           cobble_P, gravel_P, mud_P, sand_P, 
                           rugos, BATHY.DEPTH, h_bt,
                           nao, amo)
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

#### Make region ####
source(here('R_code/utilities/vast_functions.R'))
# First, we need our region shapefile
region_shape<- st_read(here::here("", "Data/GIS/cod_region_utm.shp"))
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
region_shape$OBJECTID <- 'ALL'
colnames(region_shape) <- c("Region", 'geometry')

# Second, get our index area shapefile
# We could just use this same shapefile in the "index_shapes" argument, but to 
# show off the new functions we wrote, we will also want to have a sf multipolygon
# shapefiles with areas defined within this general region
index_areas<- c("EGOM", "GBK", "SNE", "WGOM")

for(i in seq_along(index_areas)){
  index_area_temp<- st_read(paste0(here::here("", "Data/GIS"), '/', 
                                   index_areas[i], "_UTM.shp"))
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

index_area_shapes <- bind_rows(index_area_shapes, region_shape)

index_area_shapes <- st_make_valid(index_area_shapes)

# Finally, run `vast_make_extrap_grid` function after specifying the strata we want to use
strata_use<- data.frame(STRATA = c( "ALL", "EGOM", "GBK", "SNE", "WGOM"))
vast_extrap_grid<- vast_make_extrap_grid(region_shapefile = region_shape, 
                                         index_shapes = index_area_shapes, 
                                         #strata.limits = strata_use, 
                                         cell_size = 1000)
vast_extrap_grid <- dplyr::select(vast_extrap_grid, Lon, Lat, Area_km2, STRATA)
row.names(vast_extrap_grid) <- NULL

# Remove intermediates
rm(list=setdiff(ls(), c("scaled.covars", "settings", 'strata_use', 'survs',
                        'vast_extrap_grid')))
gc()

working_dir <- here::here("VAST_runs/add_climate_aja4/ALL")
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

extrap_info_aja_all <- make_extrapolation_info(Region="User",
                                            strata.limits=strata_use,
                                            input_grid = vast_extrap_grid,
                                            max_cells = 2000
)
head(extrap_info_aja_all$a_el)
head(extrap_info_aja_all$Data_Extrap)
write_rds(extrap_info_aja_all,
          here('Data/VAST_input/extrap_info_aja_all.RDS'))
#extrap_info_aja_wgom <- readRDS(here('Data/VAST_input/extrap_info_aja_wgom.RDS'))

# Remove intermediates
rm(list=setdiff(ls(), c("scaled.covars", "settings", 'strata_use', 'survs',
                        'vast_extrap_grid', 'extrap_info',
                        'working_dir', 'extrap_info_aja_all')))
gc()

#### Run model ####


# Set year labels
year.labs <- c(seq(1982, 2021, 1), seq(1982, 2021, 1))
year.labs <- year.labs[order(year.labs)]
seas.labs <- rep(c('Spring', 'Fall'), 40)
year.labs <- paste0(year.labs, " ", seas.labs)
cat.labs <- c('Small', 'Medium', 'Large')

# Make settings
settings <- make_settings(
  n_x = 200,
  purpose = "index2",
  Region = "User",
  fine_scale = TRUE,
  bias.correct = FALSE,
  knot_method = "grid",
  max_cells = 2000,
  strata.limits = strata_use
)

getwd()
# Change `ObsModel` settings
# Value 1: Delta-lognormal using bias-corrected-mean and log-coeff of var
# Value 2: Poisson-link delta-model" using log-link for numbers density and 
#          log-link for biomass per number
settings$ObsModel[1] <- 4
settings$ObsModel[2] <- 1

# Model formula
hab_formula<- ~ bs(cobble_P, degree = 3, intercept = FALSE) + 
  bs(gravel_P, degree = 3, intercept = FALSE) + 
  bs(mud_P, degree = 3, intercept = FALSE) +
  bs(sand_P, degree = 3, intercept = FALSE) +
  bs(rugos, degree = 3, intercept = FALSE) + 
  bs(BATHY.DEPTH, degree = 3, intercept = FALSE) +
  bs(h_bt, degree = 3, intercept = FALSE) +
  bs(nao, degree = 3, intercept = FALSE) +
  bs(amo, degree = 3, intercept = FALSE)

fit = fit_model( 

  # Call settings
    settings = settings, 
    #"Version" = settings$Version,
    
  # Call survey data info
    Lat_i = survs[,'Lat'], 
    Lon_i = survs[,'Lon'], 
    t_i = survs[,'Year'],
    b_i = survs[,'Response_variable'],
    a_i = survs[,'swept'],
    v_i = survs[,'vessel'],
    c_iz = survs[,'Age'],
  
  # Call covariate info
    X1_formula = hab_formula,
    X2_formula = hab_formula,
    covariate_data = scaled.covars,
  
  # Call spatial 
    input_grid=vast_extrap_grid,
    "extrapolation_list" = extrap_info_aja_all,
  
  # Set naming conventions
    category_names = cat.labs,
    year_labels = year.labs,
  
  # Tell model to run
    build_model = TRUE,
    run_model = TRUE

    )

  
beep(8)

#### Plot results ####

save.image('add_climate_aja4_all.RData')

plot( fit )
beep(8)
# 
# map_list = make_map_info(Region = fit$settings$Region, 
#                          spatial_list = fit$spatial_list, 
#                          Extrapolation_List = fit$extrapolation_list)
# 
# plot_clusters(fit=fit,
#               var_name = "D_gct",
#               transform_var = log,
#               k=4,
#               method='ward',
#               year_labels = fit$year_labels,
#               category_names = fit$category_names,
#               map_list = map_list,
#               working_dir = "C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/VAST_runs/refine_effort/",
#               file_name = "Class-D_gct",
#               file_name2 = "Class-D_gct-averages",
#               replace_Inf_with_NA = TRUE,
#               size_threshold = 1e+20,
#               col = viridisLite::viridis,
#               yaxis_log = TRUE
#               
# )
