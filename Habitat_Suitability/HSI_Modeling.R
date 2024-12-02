rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)
library(stars)
library(raster)

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

#### Load model fit data ####
load(here("VAST_runs/medium/Overall_BC/ALL/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.Rdata"))
medium <- fit
# load(here("VAST_runs/small/Overall_BC/ALL/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.Rdata"))
# small <- fit
# load(here("VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.Rdata"))
# large <- fit

# Clean workspace
rm(list=setdiff(ls(), c('small', 'medium', 'large')))

# Set WD for plotting
out_dir = here("VAST_runs/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Extract HSI curves ####
# Name size classes
sizes <- c(#'small', 
  'medium')#, 'large')

# Blank dataframe to fill
keep <- data.frame(
  fit=NA, se=NA, lower=NA, upper=NA, Lin_pred=NA, Value=NA,
  Covariate=NA, Size=NA
)

# Extract data, looping through sizes
for(i in 1:length(sizes)){
  size <- sizes[i]
  if(size == 'small'){fit <- small}
  if(size == 'medium'){fit <- medium}
  if(size == 'large'){fit <- large}
  
  # Call parameter names
  params <- colnames(fit$covariate_data)
  params <- params[params %notin% c('Lon', 'Lat', 'Year')]
  
  # Parameters not used in certain models
  if(size == 'small'){params <- params[params %notin% c('amo')]}
  if(size == 'large'){params <- params[params %notin% c('cobble_P',
                                                        'sand_P', 
                                                        'mud_P')]}
  # Call category names
  catnames <- fit$category_names
  ncat = length(catnames)
  
  # Extract effects
  vast_covariate_effects<- get_vast_covariate_effects(vast_fit = fit, 
                                                      params_plot = c(params), 
                                                      params_plot_levels = 500, 
                                                      effects_pad_values = c(), 
                                                      nice_category_names = paste0(size, ' Cod'),
                                                      out_dir = out_dir,
                                                      category_to_use = 1,
                                                      ncat = ncat)
  
  # Parameters we want to use
  names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")
  
  # Merge non-informative names to informative name descriptions
  betternames=data.frame(
    Covariate = c('cobble_P', 'gravel_P', 'mud_P', 'sand_P', 'rugos', 'BATHY.DEPTH',
                  'h_bt', 'nao', 'amo'),
    Better = c('Cobble', 'Gravel', 'Mud', 'Sand', 'Rugosity', 'Depth (m)', 
               'BT (C)', 'NAO', 'AMO')
  )
  
  # Remove NA values
  vast_cov_eff_l <- vast_covariate_effects %>%
    drop_na(Value)
  
  # Merge to better names
  vast_cov_eff_l <- merge(vast_cov_eff_l, betternames)
  
  # Clean dataframe
  vast_cov_eff_l <- vast_cov_eff_l %>% 
    dplyr::select(-Covariate) %>% 
    rename(Covariate = Better)
  
  # Name size model
  vast_cov_eff_l$Size <- size
  
  # Bind to blank dataframe
  keep <- rbind(keep, vast_cov_eff_l)
  
  # Clean workspace
  rm(size, params, catnames, ncat, vast_covariate_effects, names_stay,
     betternames, vast_cov_eff_l)
  
}
# A raw copy (takes forever to reach this point, don't want to redo)
ejectbutton <- keep

# Remove blank initiation dataframe
keep <- keep[!is.na(keep$Size),]

# Clean
keep <- keep %>% 
  mutate(Size = str_to_sentence(Size)) %>% 
  mutate(Size = factor(Size, levels = c('Small', 'Medium', 'Large')),
         Covariate = factor(Covariate, levels=c('BT (C)',
                                                'Depth (m)',
                                                'Gravel',
                                                'Sand',
                                                'Mud',
                                                'Cobble',
                                                'Rugosity',
                                                'NAO',
                                                'AMO')),
         Lin_pred = factor(Lin_pred, levels=c('X1', 'X2')))
keep <- droplevels(keep)

# Clean workspace
rm(list=setdiff(ls(), c('keep', 'ejectbutton', '%notin%')))

#### Call depth data ####
# Load raster
gebco <- raster(here('Data/Bathymetry/GEBCO_2023.tif'))

# Clip to region
region <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
region <- st_make_valid(region)
region <- region[region$STOCK == 'EGOM',]
depth <- gebco %>% 
  mask(region) %>%
  crop(region)

# Find range of depth in region
depth@data@values <- depth@data@values * -1
summary(depth@data@values)
plot(depth, asp=1)
keepdepth <- depth

val_curv <- keep %>% 
  filter(Lin_pred == 'X1' & Covariate == 'Depth (m)') %>% 
  dplyr::select(Value, fit) %>% 
  rename(x=Value, y=fit)

# Replace raster data values with matched SI values
for(m in 1:length(depth)){
  print(m)
  newval <- depth[m]
  if(is.na(newval)){next()}
  if(newval<=0){
    depth[m] <- NA
    next()}
  newval <- val_curv$y[val_curv$x==newval]
  # If x-value has no exact match, find linear model to n5 points centered on nearest value
  if(length(newval)==0){
    target.index <- which(abs(val_curv$x - depth[m]) == 
                            min(abs(val_curv$x - depth[m])))
    test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                      (target.index + 2),])
    
    newdata <- data.frame(x=depth[m], y=NA)
    predval <- as.numeric(predict.lm(test, newdata))
    newval <- predval
    
    # Remove intermediates
    rm(target.index, test, newdata, predval)
  }
  depth[m] <- newval
}

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

scaled.depth <- depth
scaled.depth@data@values <- range01(scaled.depth@data@values, na.rm=T)
plot(scaled.depth, asp=1)

#### Call temperature data ####
# Load rasterbrick
i <- 2019
brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                           'hubert', 'extended_grd_krig'),
                      '/', i, '.grd'))
rasterlayer <- brick

# Clip to region
region <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
region <- st_make_valid(region)
region <- region[region$STOCK == 'EGOM',]

sflayers <- vector("list", nlayers(brick))
dailymap <- data.frame(
  doy = NA,
  min = NA,
  max = NA
)

for(i in 1:nlayers(brick)){
  print(i)
  
  sflayers[[i]] <- st_as_sf(st_as_stars(brick[[i]]))
  colnames(sflayers[[i]]) <- c('BT', 'geometry')
  st_geometry(sflayers[[i]]) <- 'geometry'
  
  sflayers[[i]] <- st_intersection(sflayers[[i]], region)
  
  internalmap <- data.frame(
    doy = NA,
    min = NA,
    max = NA
  )
  internalmap$doy <- i
  internalmap$min = min(sflayers[[i]]$BT)
  internalmap$max = max(sflayers[[i]]$BT)
  
  dailymap <- rbind(dailymap, internalmap)
  
  rm(internalmap)
}
sfkeep <- sflayers

# Find range of depth in region
min(dailymap$min, na.rm=T); max(dailymap$max, na.rm=T)

val_curv <- keep %>% 
  filter(Lin_pred == 'X1' & Covariate == 'BT (C)') %>% 
  dplyr::select(Value, fit) %>% 
  rename(x=Value, y=fit)

# Replace raster data values with matched SI values
for(i in 1:length(sflayers)){
  message(i)
  for(m in 1:nrow(sflayers[[i]])){
    #print(m)
    newval <- sflayers[[i]]$BT[m]
    if(is.na(newval)){next()}
    if(newval<=min(val_curv$x)){
      sflayers[[i]]$BT[m] <- NA
      next()
      }
    newval <- val_curv$y[val_curv$x==newval]
  # If x-value has no exact match, find linear model to n5 points centered on nearest value
  if(length(newval)==0){
    target.index <- which(abs(val_curv$x - sflayers[[i]]$BT[m]) == 
                            min(abs(val_curv$x - sflayers[[i]]$BT[m])))
    
    if(target.index >=3 & target.index<=498){
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=sflayers[[i]]$BT[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
      next()
    }
    
    if(target.index ==2){
      test <- lm(y ~ x, data=val_curv[(target.index - 1):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=sflayers[[i]]$BT[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
      next()
    }
    
    if(target.index ==1){
      test <- lm(y ~ x, data=val_curv[(target.index):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=sflayers[[i]]$BT[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
      next()
    }
    
    if(target.index ==499){
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 1),])
      
      newdata <- data.frame(x=sflayers[[i]]$BT[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
      next()
    }
    
    if(target.index ==500){
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index),])
      
      newdata <- data.frame(x=sflayers[[i]]$BT[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
      next()
    }
  }
  sflayers[[i]]$BT[m] <- newval
  }
}

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

scaled.BT <- sflayers

plot(scaled.depth, asp=1)
