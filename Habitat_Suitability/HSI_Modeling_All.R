rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)
library(raster)
library(stars)
library(sf)
library(beepr)

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

# load data
load(here("VAST_runs/medium/Overall_BC/ALL/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.Rdata"))
medium <- fit
load(here("VAST_runs/small/Overall_BC/ALL/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.Rdata"))
small <- fit
load(here("VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.Rdata"))
large <- fit

# Clean workspace
rm(list=setdiff(ls(), c('small', 'medium', 'large')))

# Set WD for plotting
out_dir = here("VAST_runs/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Name size classes
sizes <- c('small', 'medium', 'large')

# Blank dataframe to fill
vast.cov.effs <- data.frame(
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
  vast.cov.effs <- rbind(vast.cov.effs, vast_cov_eff_l)
  
  # Clean workspace
  rm(size, params, catnames, ncat, vast_covariate_effects, names_stay,
     betternames, vast_cov_eff_l)
  
}

# Remove blank initiation dataframe
vast.cov.effs <- vast.cov.effs[!is.na(vast.cov.effs$Size),]

# Clean
vast.cov.effs <- vast.cov.effs %>% 
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
                                                'AMO')))

# Save results for presence/absence
vast.cov.effs <- vast.cov.effs %>% 
  filter(Lin_pred == 'X1')

# Clean workspace
rm(list=setdiff(ls(), c('vast.cov.effs', '%notin%')))

#### Call depth data ####
# Load raster
gebco <- raster(here('Data/Bathymetry/GEBCO_2023.tif'))

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(gebco))

depth <- gebco %>% 
  mask(region) %>%
  crop(region)

# Convert so matches curve -- depths are positive values, land is negative
depth@data@values <- depth@data@values * -1

# Remove land values
depth@data@values[depth@data@values <=0.1] <- NA

# Find range of depth in region
summary(depth@data@values)
plot(depth, asp=1)

# Convert to sf
depthsf <- st_as_sf(st_as_stars(depth))
depthsf <- depthsf %>% 
  rename(bathy = gebco_2023) %>% 
  filter(!is.na(bathy))

# Clean workspace
rm(gebco, region, depth)

# Make column for predicted fit to curve
usedepths <- data.frame(
  bathy = unique(depthsf$bathy)
)
usedepths$predfit <- NA

# Run through size classes
sizes <- c('Small', 'Medium', 'Large')

# Blank list to collect suitability of cell values
depth.suit <- vector('list', length=length(sizes))

for(i in 1:length(sizes)){
  size <- sizes[i]
  
  # Extract depth fit curve
  val_curv <- vast.cov.effs %>% 
    filter(Lin_pred == 'X1' & Covariate == 'Depth (m)') %>% 
    filter(Size == paste0(size)) %>% 
    dplyr::select(Value, fit) %>% 
    rename(x=Value, y=fit)
  
  # Extend modeling to values outside seen in surveys
  ggplot(data=val_curv) +
    geom_line(aes(x=x, y=y)) +
    geom_smooth(aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam')
  
  # Model what we have with a GAM, use high knot value to get 100% Dev exp
  gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
  s <- summary(gam)
  message(paste0(round(s$dev.expl * 100, 3), ' % dev exp'))
  
  # Make df with depth values we want to predict to
  newdata <- data.frame(x=c(0.1, seq((max(val_curv$x, na.rm = T) + 1), 
                                     max(depthsf$bathy, na.rm = T), by=1)))
  
  # Predict model fit at these depth values
  pgam <- mgcv::predict.gam(gam, newdata=newdata)
  
  # Append fits to depth values
  newdata$y <- pgam
  
  # Sanity check
  print(
  ggplot() +
    geom_line(data=val_curv,
              aes(x=x, y=y)) +
    geom_smooth(data=val_curv,
                aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam') +
    geom_point(data=newdata, 
               aes(x=x, y=y))
  )
  
  # Bind into single dataframe
  val_curv <- rbind(val_curv, newdata)
  val_curv <- val_curv[with(val_curv, order(x)),]
  rownames(val_curv) <- NULL
  
  class.usedepths <- usedepths

  # Replace depth with curve fit value
  for(m in 1:nrow(class.usedepths)){
    newval <- class.usedepths$bathy[m]
    
    newval <- val_curv$y[val_curv$x==newval]
    
    if(length(newval) == 1){
      class.usedepths$predfit[m] <- newval
    }
    
    # If x-value has no exact match, find linear model to n5 points centered on nearest value
    if(length(newval)==0){
      target.index <- which(abs(val_curv$x - class.usedepths$bathy[m]) == 
                              min(abs(val_curv$x - class.usedepths$bathy[m])))[1]
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=class.usedepths$bathy[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      class.usedepths$predfit[m] <- newval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
    }
    rm(newval)
  }
  
  # Save to suitability list
  colnames(class.usedepths) <- c('bathy', paste0(size))
  
  # Sanity check
  ggplot(data=class.usedepths) + 
    geom_point(aes(x=bathy, y=Small))
  
  depth.suit[[i]] <- class.usedepths

  rm(gam, s, class.usedepths, pgam, val_curv)
  
}

# Clean workspace
rm(usedepths)

# Name list items
names(depth.suit) <- c('Small', 'Medium', 'Large')

# Merge suitability to depth sf (small)
depthsf <- left_join(depthsf, depth.suit[['Small']], by=c('bathy'))
depthsf <- left_join(depthsf, depth.suit[['Medium']], by=c('bathy'))
depthsf <- left_join(depthsf, depth.suit[['Large']], by=c('bathy'))

#### Call Sediment data ####
# Load raster
seds <- read_rds(here('Data/Density_Covariates/Sediment/sediment_interp_grid_hires.RDS'))

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(seds))

seds <- st_intersection(seds, region)
seds <- seds %>% 
  dplyr::select(-Inc, -OBJECTID, -Shape_Leng, -Shape_Area, -area_km, -rock)

# Find range of each sed val
summary(seds$cobble)
summary(seds$gravel)
summary(seds$mud)
summary(seds$sand)

# Name sediments
sedtypes <- c('cobble', 'gravel', 'mud', 'sand')

# Name sizes
sizes <- c('Small', 'Medium', 'Large')

# Blank list for sediments
all.seds.suit <- vector('list', length(sedtypes))

# Loop through sediments
for(q in 1:length(sedtypes)){
  message(sedtypes[q])
  
  # Call sediment type
  seduse <- sedtypes[q]
  
  # Make column for predicted fit to curve
  usesed <- as.vector(seds[,paste0(seduse)])
  usesed <- as.data.frame(usesed[[1]])
  colnames(usesed) <- c(paste0(seduse))
  usesed$predfit <- NA
  
  # Blank list to collect suitability of cell values
  sed.suit <- vector('list', length=length(sizes))
  
  # Run through sizes
  for(i in 1:length(sizes)){
    size <- sizes[i]
    message(paste0(size))
    
    # Extract depth fit curve
    val_curv <- vast.cov.effs %>% 
      filter(Lin_pred == 'X1' & Covariate == paste0(str_to_sentence(seduse))) %>% 
      filter(Size == paste0(size)) %>% 
      dplyr::select(Value, fit) %>% 
      rename(x=Value, y=fit)
    
    if(nrow(val_curv) == 0){
      sed.suit[[i]] <- data.frame(
        No = NA,
        No2 = NA
      )
      
      colnames(sed.suit[[i]]) <- c(paste0(seduse), paste0(size))
      
      next()
    }
    
    # Extend modeling to values outside seen in surveys
    ggplot(data=val_curv) +
      geom_line(aes(x=x, y=y)) +
      geom_smooth(aes(x=x, y=y),
                  formula = y ~ s(x, bs = "cs", k=15),
                  method='gam')
    
    # Model what we have with a GAM, use high knot value to get 100% Dev exp
    gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
    
    # Make df with depth values we want to predict to
    newdata <- data.frame(x=seq(0, min(val_curv$x), by=0.001))
    newdata <- newdata[newdata$x != min(val_curv$x),]
    newdata <- newdata %>% 
      as.data.frame() %>% 
      rename(x = '.')
    
    # Predict model fit at these depth values
    pgam <- mgcv::predict.gam(gam, newdata=newdata)
    
    # Append fits to depth values
    newdata$y <- pgam
    
    # Sanity check
    print(
    ggplot() +
      geom_line(data=val_curv,
                aes(x=x, y=y)) +
      geom_smooth(data=val_curv,
                  aes(x=x, y=y),
                  formula = y ~ s(x, bs = "cs", k=15),
                  method='gam') +
      geom_point(data=newdata, 
                 aes(x=x, y=y))
    )
    
    # Bind into single dataframe
    val_curv <- rbind(val_curv, newdata)
    val_curv <- val_curv[with(val_curv, order(x)),]
    rownames(val_curv) <- NULL
    
    class.usesed <- usesed
    
    # Replace depth with curve fit value
    for(m in 1:nrow(class.usesed)){
      newval <- class.usesed[,paste0(seduse),][m]
      
      newval <- val_curv$y[val_curv$x==newval]
      
      if(length(newval) == 1){
        class.usesed$predfit[m] <- newval
        rm(newval)
        next()
      }
      
      # If x-value has no exact match, find linear model to n5 points centered on nearest value
      if(length(newval)==0){
        target.index <- which(abs(val_curv$x - class.usesed[,paste0(seduse),][m]) == 
                                min(abs(val_curv$x - class.usesed[,paste0(seduse),][m])))[1]
        
        if(target.index >=3 & target.index <=(nrow(val_curv)-2)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index ==2){
          test <- lm(y ~ x, data=val_curv[(target.index - 1):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index ==1){
          test <- lm(y ~ x, data=val_curv[(target.index):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index == (nrow(val_curv)-1)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index + 1),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index == nrow(val_curv)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
      }
      rm(newval)
    }
    
    # Save to suitability list
    colnames(class.usesed) <- c(paste0(seduse), paste0(size))
    
    sed.suit[[i]] <- class.usesed
    
    rm(gam, class.usesed, pgam, val_curv)
    
  }
  
  all.seds.suit[[q]] <- sed.suit
  
}

# Clean workspace
rm(sed.suit, usesed, val_curv,i, m, q, sedtypes, seduse, size, sizes)

# Name list items
names(all.seds.suit) <- c('cobble', 'gravel', 'mud', 'sand')
for(i in 1:length(all.seds.suit)){
  names(all.seds.suit[[i]]) <- c('Small', 'Medium', 'Large')
}

# CLean workspace
rm(region, i)

# Merge to sediments
# Cobble
cobble <- seds %>% dplyr::select(-sand, -mud, -gravel)
cobble <- left_join(cobble, unique(all.seds.suit[['cobble']][['Small']]), 
                    by='cobble')
cobble <- left_join(cobble, unique(all.seds.suit[['cobble']][['Medium']]), 
                    by='cobble')
cobble <- left_join(cobble, unique(all.seds.suit[['cobble']][['Large']]), 
                    by='cobble')

# Gravel
gravel <- seds %>% dplyr::select(-sand, -mud, -cobble)
gravel <- left_join(gravel, unique(all.seds.suit[['gravel']][['Small']]), 
                    by='gravel')
gravel <- left_join(gravel, unique(all.seds.suit[['gravel']][['Medium']]), 
                    by='gravel')
gravel <- left_join(gravel, unique(all.seds.suit[['gravel']][['Large']]), 
                    by='gravel')

# Mud
mud <- seds %>% dplyr::select(-sand, -gravel, -cobble)
mud <- left_join(mud, unique(all.seds.suit[['mud']][['Small']]), by='mud')
mud <- left_join(mud, unique(all.seds.suit[['mud']][['Medium']]), by='mud')
mud <- left_join(mud, unique(all.seds.suit[['mud']][['Large']]), by='mud')

# Sand
sand <- seds %>% dplyr::select(-gravel, -mud, -cobble)
sand <- left_join(sand, unique(all.seds.suit[['sand']][['Small']]), by='sand')
sand <- left_join(sand, unique(all.seds.suit[['sand']][['Medium']]), by='sand')
sand <- left_join(sand, unique(all.seds.suit[['sand']][['Large']]), by='sand')

#### Call rugosity data ####
# Load raster
load(here('Data/Density_Covariates/Rugosity/rast_rugosity.RData'))
rugos <- masked.raster
rm(masked.raster)

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(rugos))

rugos <- rugos %>% 
  mask(region) %>%
  crop(region)

# Find range of depth in region
summary(rugos@data@values)
plot(rugos, asp=1)

# Convert to sf
rugossf <- st_as_sf(st_as_stars(rugos))
rugossf <- rugossf %>% 
  rename(rugosity = band1) %>% 
  filter(!is.na(rugosity))

# Clean workspace
rm(rugos, region)

# Make column for predicted fit to curve
userugos <- data.frame(
  rugosity = unique(round(rugossf$rugosity,4))
)
userugos$predfit <- NA

# Run through size classes
sizes <- c('Small', 'Medium', 'Large')

# Blank list to collect suitability of cell values
rugos.suit <- vector('list', length=length(sizes))

for(i in 1:length(sizes)){
  size <- sizes[i]
  
  # Extract depth fit curve
  val_curv <- vast.cov.effs %>% 
    filter(Lin_pred == 'X1' & Covariate == 'Rugosity') %>% 
    filter(Size == paste0(size)) %>% 
    dplyr::select(Value, fit) %>% 
    rename(x=Value, y=fit)
  
  # Extend modeling to values outside seen in surveys
  ggplot(data=val_curv) +
    geom_line(aes(x=x, y=y)) +
    geom_smooth(aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam')
  
  # Model what we have with a GAM, use high knot value to get 100% Dev exp
  gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
  s <- summary(gam)
  message(paste0(round(s$dev.expl * 100, 3), ' % dev exp'))
  
  # Make df with depth values we want to predict to
  newdata <- data.frame(x=c(seq(min(userugos$rugosity), min(val_curv$x), by=0.001),
                            seq(max(val_curv$x), max(userugos$rugosity), by=0.001))
  )
  
  # Predict model fit at these depth values
  pgam <- mgcv::predict.gam(gam, newdata=newdata)
  
  # Append fits to depth values
  newdata$y <- pgam
  
  # Sanity check
  ggplot() +
    geom_line(data=val_curv,
              aes(x=x, y=y)) +
    geom_smooth(data=val_curv,
                aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam') +
    geom_point(data=newdata, 
               aes(x=x, y=y))
  
  # Bind into single dataframe
  val_curv <- rbind(val_curv, newdata)
  val_curv <- val_curv[with(val_curv, order(x)),]
  rownames(val_curv) <- NULL
  
  class.userugos <- userugos
  
  # Replace depth with curve fit value
  for(m in 1:nrow(class.userugos)){
    newval <- class.userugos$rugosity[m]
    
    newval <- val_curv$y[val_curv$x==newval]
    
    if(length(newval) == 1){
      class.userugos$predfit[m] <- newval
    }
    
    # If x-value has no exact match, find linear model to n5 points centered on nearest value
    if(length(newval)==0){
      target.index <- which(abs(val_curv$x - class.userugos$rugosity[m]) == 
                              min(abs(val_curv$x - class.userugos$rugosity[m])))[1]
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=class.userugos$rugosity[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      class.userugos$predfit[m] <- newval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
    }
    rm(newval)
  }
  
  # Save to suitability list
  colnames(class.userugos) <- c('rugosity', paste0(size))
  
  rugos.suit[[i]] <- class.userugos
  
  rm(gam, s, class.userugos, pgam, val_curv)
  
}

# Clean workspace
rm(userugos)

# Name list items
names(rugos.suit) <- c('Small', 'Medium', 'Large')

# Merge suitability to depth sf (small)
rugossf$OG.rugos <- rugossf$rugosity
rugossf$rugosity <- round(rugossf$rugosity, 4)
rugossf <- left_join(rugossf, rugos.suit[['Small']], by=c('rugosity'))
rugossf <- left_join(rugossf, rugos.suit[['Medium']], by=c('rugosity'))
rugossf <- left_join(rugossf, rugos.suit[['Large']], by=c('rugosity'))

beepr::beep(8)


#### Static habitat suitability ####
# Step one: rasterize all suitability curves to same CRS,proj, and grid
# Load VAST output grid
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))
# Clear workspace
rm(list=setdiff(ls(), c('rugossf', 'cobble', 'sand', 'gravel', 'mud',
                        'depthsf', '%notin%', 'fit')))
# Extract Data
Y_gt = fit$Report$D_gct[,1,]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
panel_labels = fit$year_labels

# Call data
MapSizeRatio = map_list$MapSizeRatio
Y_gt = Y_gt[map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                    0), "x2i"], , drop = FALSE]
n_cells = nrow(Y_gt)

# Call location list
loc_g = map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                0), c("Lon", "Lat")]
# Call projections
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS("EPSG:26919")

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
regions <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y=rep(1, 2000)), 
                                         proj4string = CRS_orig)
Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
Points_proj = sp::spTransform(Points_orig, CRS_proj)
xlim = Points_proj@bbox[1, ]
ylim = Points_proj@bbox[2, ]

cell.size = mean(diff(Points_proj@bbox[1, ]), diff(Points_proj@bbox[2, 
]))/floor(sqrt(n_cells))
Points_sf = sf::st_as_sf(Points_proj)
grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
grid_i = sf::st_intersects(Points_sf, grid)
grid = sf::st_sf(grid, 
                 y = tapply(Points_sf$y, 
                            INDEX = factor(as.numeric(grid_i),
                                           levels = 1:length(grid)), 
                            FUN = mean, na.rm = TRUE))
grid$ID <- seq(1:nrow(grid))

# PLot blank grid
ggplot() + 
  geom_sf(data=grid, col='black', fill='pink') +
  geom_sf(data=coast) +
  coord_sf(xlim=c(-76, -66),
           ylim=c(36, 45),
           crs="EPSG:4326")

# Check grid
grid <- st_transform(grid, crs="EPSG:26919")
grid <- st_make_valid(grid)
grid <- grid %>% 
  st_cast('MULTIPOLYGON') %>% 
  st_cast('POLYGON')

# Depth
depthsf <- depthsf[with(depthsf, order(bathy)),]
rownames(depthsf) <- NULL
depthsf <- st_transform(depthsf, st_crs(grid))
depthsf <- st_make_valid(depthsf)
depthdf <- depthsf %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(bathy, Small, Medium, Large) %>%
  unique() %>% 
  as.data.frame()
  
depthgrid <- st_intersection(depthsf, grid)
depthgrid <- depthgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, bathy) %>% 
  group_by(ID) %>% 
  summarise(bathy =     mean(bathy,     na.rm=T)) %>% 
  mutate(bathy = round(bathy))
depthgrid <- merge(depthgrid, depthdf, by=c('bathy'))
grid <- left_join(grid, depthgrid, by=c('ID'))
beepr::beep()

# Rugosity
rugossf <- rugossf[with(rugossf, order(rugosity)),]
rownames(rugossf) <- NULL
rugossf <- st_transform(rugossf, st_crs(grid))
rugossf <- st_make_valid(rugossf)
rugosdf <- rugossf %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(rugosity, Small, Medium, Large) %>%
  unique() %>% 
  as.data.frame()

rugosgrid <- st_intersection(rugossf, grid)
rugosgrid <- rugosgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, rugosity) %>% 
  group_by(ID) %>% 
  summarise(rugosity =     mean(rugosity,     na.rm=T))


rugosgrid <- merge(rugosgrid, rugosdf, by=c('rugosity'))
grid <- left_join(grid, rugosgrid, by=c('ID'))

# Cobble
cobblesf <- cobble[with(cobble, order(cobble)),]
rownames(cobblesf) <- NULL
cobblesf <- st_transform(cobblesf, st_crs(grid))
cobblesf <- st_make_valid(cobblesf)
cobbledf <- cobblesf %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(cobble, Small, Medium, Large) %>%
  unique() %>% 
  as.data.frame()

cobblegrid <- st_intersection(cobblesf, grid)
rugosgrid <- rugosgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, rugosity) %>% 
  group_by(ID) %>% 
  summarise(rugosity =     mean(rugosity,     na.rm=T)) %>% 
  mutate(bathy = round(bathy))
rugosgrid <- merge(rugosgrid, rugosdf, by=c('bathy'))
grid <- left_join(grid, rugosgrid, by=c('ID'))

#### Old

# Rugosity
rugospt <- st_centroid(rugossf)
rugospt <- st_transform(rugospt, st_crs(grid))
rugosgrid <- st_intersection(rugospt, grid)

rugosgrid <- sfheaders::sf_to_df(rugosgrid, fill=T)
rugosgrid <- rugosgrid %>% 
  group_by(ID) %>% 
  summarise(rugosity  = mean(rugosity,  na.rm=T),
            Sma.Rugos = mean(Small, na.rm=T),
            Med.Rugos = mean(Medium, na.rm=T),
            Lar.Rugos = mean(Large, na.rm=T))

grid <- left_join(grid, rugosgrid, by=c('ID'))

# Cobble
cobblept <- st_centroid(cobble)
cobblept <- st_transform(cobblept, st_crs(grid))
cobblegrid <- st_intersection(cobblept, grid)

cobblegrid <- sfheaders::sf_to_df(cobblegrid, fill=T)
cobblegrid <- cobblegrid %>% 
  group_by(ID) %>% 
  summarise(cobble     = mean(cobble,    na.rm=T),
            Sma.Cobble = mean(Small, na.rm=T),
            Med.Cobble = mean(Medium, na.rm=T),
            Lar.Cobble = mean(Large, na.rm=T))

grid <- left_join(grid, cobblegrid, by=c('ID'))

# Gravel
gravelpt <- st_centroid(gravel)
gravelpt <- st_transform(gravelpt, st_crs(grid))
gravelgrid <- st_intersection(gravelpt, grid)

gravelgrid <- sfheaders::sf_to_df(gravelgrid, fill=T)
gravelgrid <- gravelgrid %>% 
  group_by(ID) %>% 
  summarise(gravel     = mean(gravel,    na.rm=T),
            Sma.Gravel = mean(Small, na.rm=T),
            Med.Gravel = mean(Medium, na.rm=T),
            Lar.Gravel = mean(Large, na.rm=T))

grid <- left_join(grid, gravelgrid, by=c('ID'))

# Sand
sandpt <- st_centroid(sand)
sandpt <- st_transform(sandpt, st_crs(grid))
sandgrid <- st_intersection(sandpt, grid)

sandgrid <- sfheaders::sf_to_df(sandgrid, fill=T)
sandgrid <- sandgrid %>% 
  group_by(ID) %>% 
  summarise(sand     = mean(sand,      na.rm=T),
            Sma.Sand = mean(Small, na.rm=T),
            Med.Sand = mean(Medium, na.rm=T),
            Lar.Sand = mean(Large, na.rm=T))

grid <- left_join(grid, sandgrid, by=c('ID'))

# Mud
mudpt <- st_centroid(mud)
mudpt <- st_transform(mudpt, st_crs(grid))
mudgrid <- st_intersection(mudpt, grid)

mudgrid <- sfheaders::sf_to_df(mudgrid, fill=T)
mudgrid <- mudgrid %>%  
  group_by(ID) %>% 
  summarise(mud     = mean(mud, na.rm=T),
            Sma.Mud = mean(Small, na.rm=T),
            Med.Mud = mean(Medium, na.rm=T),
            Lar.Mud = mean(Large, na.rm=T))

grid <- left_join(grid, mudgrid, by=c('ID'))

saveRDS(grid, 
        file=paste0(here('Habitat_Suitability'), '/wholegrid.RDS'))

# Clean grid for smalls
small.grid <- grid %>% 
  dplyr::select(ID, 
                bathy, rugosity, cobble, gravel, sand, mud,
                Sma.Depth, Sma.Rugos, Sma.Cobble, Sma.Gravel, Sma.Sand, Sma.Mud,
                grid) %>% 
  rename(geometry = grid)

st_geometry(small.grid) <- 'geometry'

# Pivot to long
small.grid <- small.grid %>% 
  pivot_longer(cols=c('Sma.Depth', 'Sma.Rugos', 'Sma.Cobble',
                      'Sma.Gravel', 'Sma.Sand', 'Sma.Mud'),
               names_to = c('Var'),
               values_to = c('Suit'))

# Scale 0-1
range01 <- function(x, ...){(x - min(x, ..., na.rm=T)) / (max(x, ..., na.rm=T) - min(x, ..., na.rm=T))}
small.grid$SI <- range01(small.grid$Suit)

# Pivot back to wide
small.grid <- small.grid %>% 
  pivot_wider(names_from = Var, values_from = SI)

small.grid <- small.grid %>% 
  group_by(ID) %>% 
  summarise(Sma.Depth = mean(Sma.Depth, na.rm=T),
            Sma.Rugos = mean(Sma.Rugos, na.rm=T),
            Sma.Cobble = mean(Sma.Cobble, na.rm=T),
            Sma.Gravel = mean(Sma.Gravel, na.rm=T),
            Sma.Sand = mean(Sma.Sand, na.rm=T),
            Sma.Mud = mean(Sma.Mud, na.rm=T))

# Geometric and arithmetic average HSI for smalls
"geometric.mean" <- function(x,na.rm=TRUE){ 
    exp(mean(log(x),na.rm=na.rm))
}

small.grid$aa <- NA
small.grid$ga <- NA
for(i in 1:nrow(small.grid)){
  small.grid$aa[i] <- mean(c(small.grid$Sma.Depth[i],
                           small.grid$Sma.Rugos[i],
                           small.grid$Sma.Cobble[i],
                           small.grid$Sma.Gravel[i],
                           small.grid$Sma.Sand[i],
                           small.grid$Sma.Mud[i]),
                           na.rm=T)
  
  small.grid$ga[i] <- geometric.mean(c(small.grid$Sma.Depth[i],
                             small.grid$Sma.Rugos[i],
                             small.grid$Sma.Cobble[i],
                             small.grid$Sma.Gravel[i],
                             small.grid$Sma.Sand[i],
                             small.grid$Sma.Mud[i]))
  
  
}

ggplot() +
  geom_sf(data=small.grid, 
          aes(fill=aa, col=aa)) +
  scale_fill_viridis_c(option='viridis',
                       na.value = 'transparent',
                       limits=c(0,1)) +
  scale_color_viridis_c(option='viridis',
                        na.value = 'transparent',
                        limits=c(0,1))

# Clean grid for mediums
medium.grid <- grid %>% 
  dplyr::select(ID, Med.Depth, Med.Rugos, Med.Cobble, Med.Gravel, Med.Sand, Med.Mud,
                grid) %>% 
  rename(geometry = grid)

st_geometry(medium.grid) <- 'geometry'

# Geometric and arithmetic average HSI for mediums
medium.grid$aa <- NA
medium.grid$ga <- NA
for(i in 1:nrow(medium.grid)){
  medium.grid$aa[i] <- mean(c(medium.grid$Med.Depth[i],
                             medium.grid$Med.Rugos[i],
                             medium.grid$Med.Cobble[i],
                             medium.grid$Med.Gravel[i],
                             medium.grid$Med.Sand[i],
                             medium.grid$Med.Mud[i]))
  
  medium.grid$ga[i] <- geometric.mean(c(medium.grid$Med.Depth[i],
                                       medium.grid$Med.Rugos[i],
                                       medium.grid$Med.Cobble[i],
                                       medium.grid$Med.Gravel[i],
                                       medium.grid$Med.Sand[i],
                                       medium.grid$Med.Mud[i]))
  
  
}

ggplot() +
  geom_sf(data=medium.grid, 
          aes(fill=ga, col=ga)) +
  scale_fill_viridis_c(option='viridis',
                       na.value = 'transparent',
                       limits=c(0,1)) +
  scale_color_viridis_c(option='viridis',
                        na.value = 'transparent',
                        limits=c(0,1))

# Clean grid for larges
large.grid <- grid %>% 
  dplyr::select(ID, Lar.Depth, Lar.Rugos, Lar.Gravel,
                grid) %>% 
  rename(geometry = grid)

st_geometry(large.grid) <- 'geometry'

# Geometric and arithmetic average HSI for larges
large.grid$aa <- NA
large.grid$ga <- NA
for(i in 1:nrow(large.grid)){
  large.grid$aa[i] <- mean(c(large.grid$Lar.Depth[i],
                             large.grid$Lar.Rugos[i],
                             large.grid$Lar.Gravel[i]))
  
  large.grid$ga[i] <- geometric.mean(c(large.grid$Lar.Depth[i],
                                       large.grid$Lar.Rugos[i],
                                       large.grid$Lar.Gravel[i]))
  
  
}

ggplot() +
  geom_sf(data=large.grid, 
          aes(fill=ga, col=ga)) +
  scale_fill_viridis_c(option='viridis',
                       na.value = 'transparent',
                       limits=c(0,1)) +
  scale_color_viridis_c(option='viridis',
                        na.value = 'transparent',
                        limits=c(0,1))

# Clean workspace
rm(list=setdiff(ls(), c('small.grid', 'medium.grid', 'large.grid')))

# Add strata
strata <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
strata <- st_transform(strata, st_crs(small.grid))
strata <- strata %>% dplyr::select(-OBJECTID, -Shape_Leng, -Shape_Area)

ggplot() +
  geom_sf(data=large.grid, 
          aes(fill=ga, col=ga)) +
  scale_fill_viridis_c(option='viridis',
                       na.value = 'transparent',
                       limits=c(0,1)) +
  scale_color_viridis_c(option='viridis',
                        na.value = 'transparent',
                        limits=c(0,1)) +
  geom_sf(data=strata, fill=NA, col='yellow')

save.image(file=paste0(here('Habitat_Suitability'), '/SizeGrids.RData'))
