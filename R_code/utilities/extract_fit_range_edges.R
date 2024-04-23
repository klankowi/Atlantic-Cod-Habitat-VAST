### Post-fit VAST function to extract range quantile values ####
extract_fit_range_edges <- function(fit,
                                    quantvals) {
  # Call data needed in standard function
  Obj = fit$tmb_list$Obj
  Report = Obj$report()
  TmbData = Obj$env$data
  Sdreport = fit$parameter_estimates$SD
  SD = TMB::summary.sdreport(Sdreport)
  
  # Call labels
  year_labels = fit$year_labels
  D_name = "D_gct"
  years_to_plot = 1:TmbData$n_t
  strata_names = 1:TmbData$n_l
  category_names = 1:TmbData$n_c
  m_labels = colnames(TmbData$Z_gm)
  
  # Set wd
  working_dir = getwd()
  
  # Set sampling parameters
  quantiles = quantvals
  n_samples = 100
  interval_width = 1
  
  # Set plotting parameters
  width = NULL
  height = NULL
  calculate_relative_to_average = FALSE
  seed = 123456 
  
  # Must provide necessary inputs
  D_gctr = sample_variable(Sdreport = Sdreport, 
                           Obj = Obj,
                           variable_name = D_name, 
                           n_samples = n_samples, 
                           seed = seed)
  if (any(D_gctr == Inf)) 
    stop("`sample_variable` in `plot_range_edge` is producing D=Inf; please use `n_samples=0` to avoid error")
  if (FALSE) {
    D_gctr = sample_variable(Sdreport = Sdreport, Obj = Obj, 
                             variable_name = "L_epsilon2_cf", 
                             n_samples = n_samples, 
                             seed = seed)
  }
  
  # Create arrays
  E_zctm = array(NA, dim = c(length(quantiles), dim(Report[[D_name]])[2:3], 
                             ncol(TmbData$Z_gm)))
  E_zctmr = array(NA, dim = c(length(quantiles), dim(Report[[D_name]])[2:3], 
                              ncol(TmbData$Z_gm), n_samples))
  Mean_cmr = array(NA, dim = c(dim(Report[[D_name]])[2], ncol(TmbData$Z_gm), 
                               n_samples))
  prop_zctm = array(NA, dim = c(dim(Report[[D_name]])[1:3], 
                                ncol(TmbData$Z_gm)))
  prop_zctmr = array(NA, dim = c(dim(Report[[D_name]])[1:3], 
                                 ncol(TmbData$Z_gm), n_samples))
  
  # Sample from distribution 
  for (rI in 0:n_samples) {
    for (mI in 1:ncol(TmbData$Z_gm)) {
      order_g = order(TmbData$Z_gm[, mI], decreasing = FALSE)
      if (rI == 0) 
        prop_zctm[, , , mI] = apply(Report[[D_name]], 
                                    MARGIN = 2:3, FUN = function(vec) {
                                      cumsum(vec[order_g])/sum(vec)
                                    })
      if (rI >= 0) 
        prop_zctmr[, , , mI, rI] = apply(D_gctr[, , , 
                                                rI, drop = FALSE], MARGIN = 2:3, FUN = function(vec) {
                                                  cumsum(vec[order_g])/sum(vec)
                                                })
      for (cI in 1:dim(E_zctm)[2]) {
        if (rI >= 1) {
          if (calculate_relative_to_average == TRUE) {
            Mean_cmr[cI, mI, rI] = weighted.mean(as.vector(TmbData$Z_gm[, 
                                                                        mI] %o% rep(1, dim(Report[[D_name]])[3])), 
                                                 w = as.vector(D_gctr[, cI, , rI]))
          }
          else {
            Mean_cmr[cI, mI, rI] = 0
          }
        }
        for (zI in 1:dim(E_zctm)[1]) {
          for (tI in 1:dim(E_zctm)[3]) {
            if (rI == 0) {
              index_tmp = which.min((prop_zctm[, cI, 
                                               tI, mI] - quantiles[zI])^2)
              E_zctm[zI, cI, tI, mI] = TmbData$Z_gm[order_g[index_tmp], 
                                                    mI]
            }
            if (rI >= 1) {
              index_tmp = which.min((prop_zctmr[, cI, 
                                                tI, mI, rI] - quantiles[zI])^2)
              E_zctmr[zI, cI, tI, mI, rI] = TmbData$Z_gm[order_g[index_tmp], 
                                                         mI] - Mean_cmr[cI, mI, rI]
            }
          }
        }
      }
    }
  }
  
  # Create array of Quantile-Direction-Year with SD
  SE_zctm = apply(E_zctmr, MARGIN = 1:4, FUN = sd)
  Edge_zctm = abind::abind(Estimate = E_zctm, `Std. Error` = SE_zctm, 
                           along = 5)
  dimnames(Edge_zctm)[[1]] = paste0("quantile_", quantiles)
  
  # Remove extraneous
  rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs",
                          "Edge_zctm", 'quantvals')))
  
  # Northing
  mI = 2
  Index_zct = array(Edge_zctm[, , , mI, "Estimate"], dim(Edge_zctm)[1:3])
  sd_Index_zct = array(Edge_zctm[, , , mI, "Std. Error"], 
                       dim(Edge_zctm)[1:3])
  
  est <- aperm(Index_zct, c(2,3,1))
  est <- as.data.frame(est[,,])
  colnames(est) <- quantvals
  est <- as.data.frame(pivot_longer(est, cols=c(colnames(est))))
  yl <- as.data.frame(rep((fit$year_labels),3))
  yl <- yl %>% 
    rename(time = `rep((fit$year_labels), 3)`) %>% 
    separate(time, into=c('Year', 'Season'))
  
  yl$Year <- as.numeric(yl$Year)
  yl$Season <- factor(yl$Season, 
                      levels=c('Spring', 'Fall'))
  yl <- yl[with(yl, order(Year, Season)),]
  
  est$Year <- yl$Year
  est$Season <- yl$Season
  colnames(est)[2] <- 'est'
  
  sdest <- aperm(sd_Index_zct, c(2,3,1))
  sdest <- as.data.frame(sdest[,,])
  colnames(sdest) <- quantvals
  sdest <- as.data.frame(pivot_longer(sdest, cols=c(colnames(sdest))))
  sdest$Year <- yl$Year
  sdest$Season <- yl$Season
  colnames(sdest)[2] <- 'sd'
  
  northing.re <- merge(est, sdest, by=c('name', 'Year', 'Season'))
  colnames(northing.re) <- c('Quantile', 'Year', 'Season',
                             'Northing.Est', 
                             'Northing.SD')
  
  northing.re <- northing.re[with(northing.re,
                             order(Quantile, Year, Season)),]
  
  # Easting 
  mI = 1
  Index_zct = array(Edge_zctm[, , , mI, "Estimate"], dim(Edge_zctm)[1:3])
  sd_Index_zct = array(Edge_zctm[, , , mI, "Std. Error"], 
                       dim(Edge_zctm)[1:3])
  
  est <- aperm(Index_zct, c(2,3,1))
  est <- as.data.frame(est[,,])
  colnames(est) <- quantvals
  est <- as.data.frame(pivot_longer(est, cols=c(colnames(est))))
  est$Year <- yl$Year
  est$Season <- yl$Season
  colnames(est)[2] <- 'est'
  
  sdest <- aperm(sd_Index_zct, c(2,3,1))
  sdest <- as.data.frame(sdest[,,])
  colnames(sdest) <- quantvals
  sdest <- as.data.frame(pivot_longer(sdest, cols=c(colnames(sdest))))
  sdest$Year <- yl$Year
  sdest$Season <- yl$Season
  colnames(sdest)[2] <- 'sd'
  
  easting.re <- merge(est, sdest, by=c('name', 'Year', 'Season'))
  colnames(easting.re) <- c('Quantile', 'Year', 'Season',
                            'Easting.Est', 
                            'Easting.SD')
  
  re <- merge(northing.re, easting.re, by=c('Quantile', 'Year', 'Season'))
  
  re <- re[with(re, order(Quantile, Year, Season)),]
  
  re$strata <- fit$settings$strata.limits$STRATA[1]
  
  write.csv(re, row.names = F,
            paste0('RangeEdges_',
                   fit$settings$strata.limits$STRATA[1],
                   '.csv')
            )
}