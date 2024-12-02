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

#### Call temperature data ####
# Set years
years <- 1982:2021

big.spring <- vector('list', 40)
big.fall <- vector('list', 40)

# Loop through years
for(i in 1:length(years)){
  year <- years[i]
  print(year)
  
  brick <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                             'hubert', 'extended_grd_krig'),
                        '/', year, '.grd'))
  
  dates <- seq(as.Date(paste0(year, "-01-01")),
               as.Date(paste0(year, "-12-31")),by="1 day")
  
  if(length(dates) < nlayers(brick)){
    dates <- c(dates, NA)
  }
  
  names(brick) <- dates
  
  if(year<2021){
    brick2 <- brick(paste0(here("Data", 'Density_Covariates', 'Bottom_temp',
                                'hubert', 'extended_grd_krig'),
                           '/', (year+1), '.grd'))
    
    dates2 <- seq(as.Date(paste0((year+1), "-01-01")),
                  as.Date(paste0((year+1), "-12-31")),by="1 day")
    
    if(length(dates2) < nlayers(brick2)){
      dates2 <- c(dates2, NA)
    }
    
    names(brick2) <- dates2
    
    springdates <- seq(as.Date(paste0(year, "-03-01")),
                       as.Date(paste0(year, "-08-31")),by="1 day")
    springdates <- gsub(x=springdates, pattern='-', replacement='.')
    for(k in 1:length(springdates)){springdates[k] <- paste0('X', springdates[k])}
    
    springbrick <- brick(brick[[which(names(brick) %in% springdates)]])
    springmean <- stackApply(springbrick, indices=c(rep(1, nlayers(springbrick))),
                             fun=mean)
    
    falldates <- seq(as.Date(paste0(year, "-09-01")),
                     as.Date(paste0((year+1), "-02-28")),by="1 day")
    if(leap_year(year+1)){
      falldates <- c(falldates, paste0(year+1), '-02-29')
    }
    falldates <- gsub(x=falldates, pattern='-', replacement='.')
    for(k in 1:length(falldates)){falldates[k] <- paste0('X', falldates[k])}
    
    fallbrick1 <- brick(brick[[which(names(brick) %in% falldates)]])
    fallbrick2 <- brick(brick[[which(names(brick2) %in% falldates)]])
    fallstack <- stack(fallbrick1, fallbrick2)
    fallbrick <- brick(fallstack)
    
    fallmean <- stackApply(fallbrick, indices=c(rep(1, nlayers(fallbrick))),
                           fun=mean)
    
    big.spring[[i]] <- springmean
    big.fall[[i]] <- fallmean
    
    rm(fallmean, springmean, fallbrick, fallstack, fallbrick1, fallbrick2, falldates,
       springbrick, springdates, brick, brick2, dates, dates2, year)
    next()
  }
  
  if(year == 2021){
    springdates <- seq(as.Date(paste0(year, "-03-01")),
                       as.Date(paste0(year, "-08-31")),by="1 day")
    springdates <- gsub(x=springdates, pattern='-', replacement='.')
    for(k in 1:length(springdates)){springdates[k] <- paste0('X', springdates[k])}
    
    springbrick <- brick(brick[[which(names(brick) %in% springdates)]])
    springmean <- stackApply(springbrick, indices=c(rep(1, nlayers(springbrick))),
                             fun=mean)
    
    falldates <- seq(as.Date(paste0(year, "-09-01")),
                     as.Date(paste0(year, "-12-31")),by="1 day")
    falldates <- gsub(x=falldates, pattern='-', replacement='.')
    for(k in 1:length(falldates)){falldates[k] <- paste0('X', falldates[k])}
    
    fallbrick <- brick(brick[[which(names(brick) %in% falldates)]])
    
    fallmean <- stackApply(fallbrick, indices=c(rep(1, nlayers(fallbrick))),
                           fun=mean)
    
    big.spring[[i]] <- springmean
    big.fall[[i]] <- fallmean
    
    rm(fallmean, springmean, fallbrick, fallstack, fallbrick1, fallbrick2, falldates,
       springbrick, springdates, brick, brick2, dates, dates2, year)
  }
}

names(big.spring) <- years
names(big.fall) <- years

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(gebco))

for(i in 1:length(big.spring)){
  print(years[i])
  big.spring[[i]] <- big.spring[[i]] %>% 
    mask(region) %>% 
    crop(region)
  
  big.fall[[i]] <- big.fall[[i]] %>% 
    mask(region) %>% 
    crop(region)
}

for(i in 1:length(big.fall)){
  fallsf <- st_as_sf(st_as_stars(big.fall[[i]]))
  fallsf <- fallsf %>% 
    rename(bt = index_1) %>% 
    filter(!is.na(falltemp)) %>% 
    mutate(Season = 'Fall',
           Year = years[i])
  
  springsf <- st_as_sf(st_as_stars(big.spring[[i]]))
  springsf <- springsf %>% 
    rename(bt = index_1) %>% 
    filter(!is.na(springtemp)) %>% 
    mutate(Season = 'Spring',
           Year = years[i])
  
  big.fall[[i]] <- fallsf
  big.spring[[i]] <- springsf
  
  rm(fallsf, springsf)
}

# Group for funsies
fall <- do.call(rbind, big.fall)
spring <- do.call(rbind, big.spring)

# Find range of depth in region
summary(fall$bt[fall$Year != 2021])
summary(spring$bt[spring$Year != 2021])

# Clean workspace
rm(spring, fall)

for(w in 1:length(big.spring)){
  print(w)
  # Make column for predicted fit to curve
  usetemps <- data.frame(
    temps = unique(big.spring[[w]]$bt)
  )
  usetemps$predfit <- NA
  
  # Run through size classes
  sizes <- c('Small', 'Medium', 'Large')
  
  # Blank list to collect suitability of cell values
  temp.suit <- vector('list', length=length(sizes))
  
  for(i in 1:length(sizes)){
    size <- sizes[i]
    
    # Extract depth fit curve
    val_curv <- vast.cov.effs %>% 
      filter(Lin_pred == 'X1' & Covariate == 'BT (C)') %>% 
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
    newdata <- data.frame(x=big.spring[[w]]$bt)
    
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
    
    class.usetemps <- usetemps
    
    # Replace depth with curve fit value
    for(m in 1:nrow(class.usetemps)){
      newval <- class.usetemps$temps[m]
      
      newval <- val_curv$y[val_curv$x==newval]
      
      if(length(newval) == 1){
        class.usetemps$predfit[m] <- newval
        rm(newval)
        next()
      }
      
      # If x-value has no exact match, find linear model to n5 points centered on nearest value
      if(length(newval)==0){
        target.index <- which(abs(val_curv$x - class.usetemps$temps[m]) == 
                                min(abs(val_curv$x - class.usetemps$temps[m])))[1]
        test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                          (target.index + 2),])
        
        newdata <- data.frame(x=class.usetemps$temps[m], y=NA)
        predval <- as.numeric(predict.lm(test, newdata))
        newval <- predval
        
        class.usetemps$predfit[m] <- newval
        
        # Remove intermediates
        rm(target.index, test, newdata, predval)
        next()
      }
    }
    
    # Save to suitability list
    colnames(class.usetemps) <- c('bt', paste0(size))
    
    temp.suit[[i]] <- class.usetemps
    
    rm(gam, s, class.usetemps, pgam, val_curv)
  }
  
  names(temp.suit) <- c('Small', 'Medium', 'Large')
  
  # Merge suitability to depth sf (small)
  big.spring[[w]] <- left_join(big.spring[[w]], temp.suit[['Small']], by=c('bt'))
  big.spring[[w]] <- left_join(big.spring[[w]], temp.suit[['Medium']], by=c('bt'))
  big.spring[[w]] <- left_join(big.spring[[w]], temp.suit[['Large']], by=c('bt'))
  
  # Scale 0-1 (inclusive)
  range01 <- function(x, ...){(x - min(x, ..., na.rm=T)) / (max(x, ..., na.rm=T) - min(x, ..., na.rm=T))}
  
  big.spring[[w]]$Sma.Scale <- range01(big.spring[[w]]$Small)
  big.spring[[w]]$Med.Scale <- range01(big.spring[[w]]$Medium)
  big.spring[[w]]$Lar.Scale <- range01(big.spring[[w]]$Large)
  
}

# Group for funsies
spring <- do.call(rbind, big.spring)
spring$Lar.Scale <- range01(spring$Large)

# Plot for funsies
r.small<-st_rasterize(depthsf %>% dplyr::select(Sma.Scale, geometry))
r.medium<-st_rasterize(depthsf %>% dplyr::select(Med.Scale, geometry))
r.large<-st_rasterize(depthsf %>% dplyr::select(Lar.Scale, geometry))

plot(r.small, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
plot(r.medium, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
plot(r.large, asp=1, zlim=c(0,1), breaks=seq(0,1,by=0.05))
rm(r.small, r.medium, r.large)
