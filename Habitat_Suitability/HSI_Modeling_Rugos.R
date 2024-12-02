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

# Scale 0-1 (inclusive)
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

rugossf$Sma.Scale <- range01(rugossf$Small, na.rm=T)
rugossf$Med.Scale <- range01(rugossf$Medium, na.rm=T)
rugossf$Lar.Scale <- range01(rugossf$Large, na.rm=T)

# Plot for funsies
r.small<-st_rasterize(rugossf %>% dplyr::select(Sma.Scale, geometry))
r.medium<-st_rasterize(rugossf %>% dplyr::select(Med.Scale, geometry))
r.large<-st_rasterize(rugossf %>% dplyr::select(Lar.Scale, geometry))

plot(r.small, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
plot(r.medium, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
plot(r.large, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
rm(r.small, r.medium, r.large)
