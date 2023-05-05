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
#library(INLAspacetime) # used for mesh-building
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
  'Unknown',
  'Age0-2', 
  'Age2-5', 
  'Age5+'
)
)) - 1
# Unknown: 0
# 0-2: 1
# 2-5: 2
# 5+ : 3


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
                        nao, amo)
covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] <- 
  covars$BATHY.DEPTH[covars$BATHY.DEPTH < 0] * -1
names(covars) <- c('Lon', 'Lat', 'Year', names(covars)[4:ncol(covars)])
table(covars$Year)

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

total <- merge(survs, scaled.covars, by=c('Year', 'Lon', 'Lat'))
table(total$Age)

library(mgcv)

total$Response_variable <- strip_units(total$Response_variable)

gam.1 <- mgcv::gam(Response_variable ~ 
                     s(cobble_P, bs='cs') +
                     s(gravel_P, bs='cs') +
                     s(mud_P, bs='cs') +
                     s(sand_P, bs='cs') +
                     #s(rugos, bs='cs') +
                     s(BATHY.DEPTH, bs='cs') +
                     s(h_bt, bs='cs') +
                     s(nao, bs='cs') +
                     s(amo, bs='cs'),# +
                     #s(Lon*Lat, bs='cs') +
                     #s(Year, bs='cs'),
                          data = total, 
                          family=gaussian(link="identity"),
                          method='REML',
                          select=FALSE)

summary(within.gam1g)
anova(within.gam1g)
gam.check(within.gam1g, k.rep=999)