### Post-fit VAST function to extract area occupied values ###

extract_fit_eff <- function(fit) {
  # Create objects needed to plot
  Sdreport = fit$parameter_estimates$SD
  SD = TMB::summary.sdreport(Sdreport)
  TmbData = fit$data_list
  Report = fit$Report
  
  # Name where data are stored in report
  EffectiveName = "effective_area_ctl"
  
  # Pull data
  SD_mean_Z_ctm = array(NA, dim = c(unlist(TmbData[c("n_c", 
                                                     "n_t", 
                                                     "n_l")]), 
  2), 
  dimnames = list(NULL, NULL, NULL,
                  c("Estimate", "Std. Error")))
  # Pull standard error
  SD_mean_Z_ctm[] = SD[which(rownames(SD) == EffectiveName), 
                       c("Estimate", "Std. Error")]
  # Name dimensions      
  names(dim(SD_mean_Z_ctm)) <- c('Category',
                                 'Time',
                                 'Strata',
                                 'ArrOc')
  
  bigarr <- data.frame(
    area.occ=NA,
    std.err=NA,
    Year=NA,
    Season=NA,
    units=NA,
    strata=NA
  )
  
  for(i in 1:dim(SD_mean_Z_ctm)[3]){
    # Keep only US strata, drop category      
    SD_mean_Z_ctm_lil <- SD_mean_Z_ctm[1,,i,]
    
    # COG
    cog <- as.data.frame(SD_mean_Z_ctm_lil[,])
    colnames(cog) <- c('area.occ', 'std.err')
    cog$Year <- fit$year_labels
    cog <- cog %>% 
      separate(Year, into=c('Year', 'Season'))
    cog$Year <- as.numeric(cog$Year)
    cog$Season <- factor(cog$Season, 
                         levels = c('Spring', 'Fall'))
    
    # Convert to meters
    cog$area.occ <- cog$area.occ * 1000
    cog$std.err <- cog$std.err * 1000
    cog$units <- as_units('meter')
    
    # Add strata
    cog$strata <- strata_use$STRATA[i]
    
    # Append
    bigarr <- rbind(bigarr, cog)
    
    rm(cog, SD_mean_Z_ctm_lil)
  }

  # Clean
  bigarr <- bigarr[!is.na(bigarr$strata),]
  
  write.csv(bigarr, row.names = F,
            'AreaOcc.csv')
}