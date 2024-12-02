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

