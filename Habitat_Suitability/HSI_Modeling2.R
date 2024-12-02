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
  mutate(#Size = factor(Size, levels = c('Small', 'Medium', 'Large')),
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
keep <- keep %>% 
  filter(Lin_pred == 'X1')

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
keepdepth <- depth

# Extract depth fit curve
val_curv <- keep %>% 
  filter(Lin_pred == 'X1' & Covariate == 'Depth (m)') %>% 
  dplyr::select(Value, fit) %>% 
  rename(x=Value, y=fit)

# Extend modeling to values outside seen in surveys
ggplot(data=val_curv) +
  geom_line(aes(x=x, y=y)) +
  geom_smooth(aes(x=x, y=y),
              formula = y ~ s(x, bs = "cs", k=15),
              method='gam')

# Model what we have with a GAM, use high knot value to get 100% Dev exp
gam1 <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
s <- summary(gam1)
message(paste0(round(s$dev.expl * 100, 3), ' % dev exp'))

# Make df with depth values we want to predict to
newdata <- data.frame(x=c(0.1, seq((max(val_curv$x, na.rm = T) + 1), 
                                       max(depth@data@values, na.rm = T), by=1)))

# Predict model fit at these depth values
pgam1 <- mgcv::predict.gam(gam1, newdata=newdata)

# Append fits to depth values
newdata$y <- pgam1

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

# Convert to sf
depthsf <- st_as_sf(st_as_stars(depth))
depthsf <- depthsf %>% 
  rename(bathy = gebco_2023) %>% 
  filter(!is.na(bathy))

# Make column for predicted fit to curve
usedepths <- data.frame(
  bathy = unique(depthsf$bathy)
)
usedepths$predfit <- NA

# Replace depth with curve fit value
for(m in 1:nrow(usedepths)){
  print(m)
  newval <- usedepths$bathy[m]
  
  newval <- val_curv$y[val_curv$x==newval]
  
  if(length(newval) == 1){
    usedepths$predfit[m] <- newval
  }

  # If x-value has no exact match, find linear model to n5 points centered on nearest value
  if(length(newval)==0){
    target.index <- which(abs(val_curv$x - usedepths$bathy[m]) == 
                            min(abs(val_curv$x - usedepths$bathy[m])))
    test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                      (target.index + 2),])
    
    newdata <- data.frame(x=usedepths$bathy[m], y=NA)
    predval <- as.numeric(predict.lm(test, newdata))
    newval <- predval
    
    usedepths$predfit[m] <- newval
    
    # Remove intermediates
    rm(target.index, test, newdata, predval)
  }
}

# Merge to depth sf
test <- left_join(depthsf, usedepths, by=c('bathy'))
