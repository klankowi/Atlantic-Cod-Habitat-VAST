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
surveys <- read.csv(here("Data/Survey_Data/Bio_Data_Agesep.csv"))

# Remove data from 2022 (incomplete)
ex <- subset(surveys, YEAR < 2022)

# Set response to count data
ex$RESPONSE <- ceiling(ex$AGE_N)

# Remove Sentinel survey
ex <- subset(ex, SURVEY !='Sentinel')

# Check that there are no missing responses
nrow(ex[is.na(ex$RESPONSE)==TRUE,])

# Add environmental data
env <- readRDS(here("Data/VAST_input/agg_stn_all_OISST.RDS"))
env <- sfheaders::sf_to_df(env, fill=T)
env <- dplyr::select(env,
                     HAUL_ID, cobble_P, gravel_P, rock_P, mud_P, sand_P,
                     BATHY.DEPTH, oisst)
ex2 <- left_join(ex, env, by="HAUL_ID")

# Add annoying rugosity data
rugos <- readRDS(here("Data/VAST_input/agg_stn_all_OISST_agesep.RDS"))
rugos <- sfheaders::sf_to_df(rugos, fill=T)
rugos <- dplyr::select(rugos, HAUL_ID, rugos)
rugos <- unique(rugos)

ex3 <- left_join(ex2, rugos, by="HAUL_ID")

ex4 <- subset(ex3, !is.na(ex3$rugos))

# Add area swept for Longline survey
# Sets are 1nmi (1.852 km). Soak time is 2 hrs.
# McElroy et al. (2019) states parameters chosen to approximate NEFSC BTS
# So set area swept to equal bts average (0.0384 sq km)
table(ex4$SURVEY[is.na(ex4$AREA_SWEPT)])
ex4$AREA_SWEPT[is.na(ex4$AREA_SWEPT)] <- 0.0384

# Check results
str(ex4)
table(ex4$YEAR)

#### Finalize sampling data inputs ####
ex4 <- subset(ex4, STOCK %in% c('SNE', 'WGOM') )
# Save sampling data
survs <- dplyr::select(ex4,
                       LON, LAT, AREA_SWEPT,
                       TIME, SURVEY, AGE_N, 
                       AGEGROUP)
survs <- subset(survs, AGEGROUP == 'Age0-2')
survs <- dplyr::select(survs, -AGEGROUP)

survs$COD_N <- as_units(survs$AGE_N, 'counts')
survs$swept <- as_units(survs$AREA_SWEPT, 'km^2')
#survs$swept <- survs$AREA_SWEPT
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
# 9         Sentinel   
# 10        SMAST Video Trawl   

#survs$Data_type <- factor(survs$Data_type, levels=c("Count", "Biomass_KG"))

# survs$AGE <- as.numeric(factor(survs$AGEGROUP, levels=c(
#   'Age0-2', 
#   'Age2-5', 
#   'Age5+',
#   'Unknown'
#   )
#   )) - 1
# # Age 2 - 5: 0
# # Age 5+   : 1
# table(survs$AGE)

survs <- dplyr::select(survs, LON, LAT, TIME, COD_N, 
                       #AGE, 
                       vessel, swept)
names(survs) <- c('Lon', 'Lat', 'Year', 'Response_variable', 
                  #'Age', 
                  'vessel', 'swept')
#survs$Response_variable <- as_units(survs$Response_variable, 'counts')
#table(survs$Year)
str(survs)

# Save covariates
covars <- subset(ex4, AGEGROUP == 'Age0-2')
covars <- dplyr::select(covars,
                        LON, LAT, TIME, cobble_P, gravel_P,
                        mud_P, rock_P, sand_P, rugos, BATHY.DEPTH, oisst)
covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] <- 
  covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] * -1
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:11])
table(covars$Year)

# Test correlation
# Create correlation matrix
df_cormat <- dplyr::select(covars, BATHY.DEPTH, rugos, sand_P, #rock_P, 
                           mud_P,
                           gravel_P, cobble_P, oisst)
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
  oisst       = as.numeric(scaled.covars$oisst)
)
str(scaled.covars)

#### Make settings ####
user_region <- readRDS(here('Data/VAST_input/VAST_regions/user_region_egomsneonly.rds'))
user_region$STRATA <- droplevels(user_region$STRATA)
#user_region$STRATA <- 'All'
user_region$Id <- NULL; user_region$row <- NULL
user_region <- user_region[with(user_region, order(Lon, Lat)),]
row.names(user_region) <- NULL
head(user_region)
strata_use <- data.frame('STRATA' = c("ALL"))

# Remove intermediates
rm(covars, ex, ex2, surveys, surveys2, sus, badtab,
   surveys.list, i, `%notin%`, df_cormat, temp, temp.list, f, sacrifices,
   alllocs, badtab2, ex.f, ex5, shit, bads, checktab, j, reponses, responses)
rm(df, env, rugos)
gc()

#### Run model ####
working_dir <- here::here("VAST_runs/refine_effort_onlysmalls")
# Create working directory if it doesn't exist
if(!dir.exists(working_dir)) {
  dir.create(working_dir)
}

setwd(working_dir)

# Set year labels
year.labs <- c(seq(1982, 2021, 1), seq(1982, 2021, 1))
year.labs <- year.labs[order(year.labs)]
seas.labs <- rep(c('Spring', 'Fall'), 40)
year.labs <- paste0(year.labs, " ", seas.labs)
#cat.labs <- c('Small', 'Medium', 'Large', 'Unknown')

# Make settings
settings <- make_settings(
  n_x = 200,
  purpose = "index2",
  Region = "User",
  fine_scale = TRUE,
  bias.correct = FALSE,
  knot_method = "grid",
  max_cells = 5000
)

# Change `ObsModel` settings
# Value 1: Delta-lognormal using bias-corrected-mean and log-coeff of var
# Value 2: Poisson-link delta-model" using log-link for numbers density and 
#          log-link for biomass per number
settings$ObsModel[1] <- 4
settings$ObsModel[2] <- 1

fit = fit_model( 
  # Adjust newton steps
    #newtonsteps = 0,

  # Call settings
    settings = settings, 
    
  # Call survey data info
    Lat_i = survs[,'Lat'], 
    Lon_i = survs[,'Lon'], 
    t_i = survs[,'Year'],
    b_i = survs[,'Response_variable'],
    a_i = survs[,'swept'],
    v_i = survs[,'vessel'],
    #c_iz = survs[,'Age'],
    #e_i = as.numeric(survs[,'e_i'])-1,
  
  # Call catchability info
   #Q1_formula = Q1_formula,
   #catchability_data  = catchability_data,
  
  # Call covariate info
    X1_formula = ~ gravel_P + cobble_P + mud_P + sand_P + # rock_P +
                   BATHY.DEPTH + oisst + rugos,
    X2_formula = ~ gravel_P + cobble_P + mud_P + sand_P + # rock_P +
                   BATHY.DEPTH + oisst + rugos,
    covariate_data = scaled.covars,
  
  # Call spatial 
    input_grid=user_region,
  
  # Set naming conventions
    #category_names = cat.labs,
    year_labels = year.labs,
  
  # Tell model to run
    build_model = TRUE,
    run_model = TRUE
    )

  
beep(8)

#### Plot results ####

save.image('only_smalls.RData')

plot( fit )
beep(8)
