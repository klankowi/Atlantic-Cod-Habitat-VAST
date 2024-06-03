extract_fit_index <- function(fit, category_names){
  DirName = getwd()
  year_labels = fit$year_labels
  strata_names = fit$settings$strata.limits$STRATA
  TmbData = fit$data_list
  Sdreport = fit$parameter_estimates$SD
  extrapolation_list = fit$extrapolation_list
  use_biascorr = TRUE
            
  index_name = "Index_ctl"
  log_index_name = "ln_Index_ctl"
            
  fit$Report = amend_output(fit = fit, 
                            year_labels = year_labels, 
                            category_names = category_names, 
                            strata_names = strata_names, 
                            extrapolation_list = extrapolation_list)

  if ("treat_nonencounter_as_zero" %in% 
      names(TmbData$Options_list$Options)) {
    treat_missing_as_zero = 
      TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }
  if (!"treat_nonencounter_as_zero" %in%
      names(TmbData$Options_list$Options)){
    treat_missing_as_zero = FALSE
  }
  
  SD = TMB::summary.sdreport(Sdreport)
  
  par_SE = TMB:::as.list.sdreport(Sdreport, 
                                  what = "Std. Error", 
                                  report = TRUE)
  
  par_hat = TMB:::as.list.sdreport(Sdreport, 
                                   what = "Estimate", 
                                   report = TRUE)
  
  if (use_biascorr == TRUE && "unbiased" %in% names(Sdreport)) {
    par_biascorrect = TMB:::as.list.sdreport(Sdreport, 
                                             what = "Est. (bias.correct)", 
                                             report = TRUE)
    for (int in seq_len(length(par_hat))) {
      par_hat[[int]] = ifelse(is.na(par_biascorrect[[int]]), 
                              par_hat[[int]], par_biascorrect[[int]])
    }
  }
  if (treat_missing_as_zero == TRUE) {
    Num_ctl = abind::adrop(TmbData$Options_list$metadata_ctz[, 
                                                             , 
                                                             "num_notna", 
                                                             drop = FALSE], 
                           drop = 3) %o% rep(1, TmbData$n_l)
    par_hat[[index_name]] = ifelse(Num_ctl == 0, 0, 
                                   par_hat[[index_name]])
    par_SE[[index_name]] = ifelse(Num_ctl == 0, 0, 
                                  par_SE[[index_name]])
    par_hat[[log_index_name]] = ifelse(Num_ctl == 0, 0, 
                                       par_hat[[log_index_name]])
    par_SE[[log_index_name]] = ifelse(Num_ctl == 0, 0, 
                                      par_SE[[log_index_name]])
  }
  for (int in seq_len(length(par_hat))) {
    if (names(par_hat)[int] %in% names(fit$Report)) {
      dimnames(par_SE[[int]]) = dimnames(par_hat[[int]]) = 
        dimnames(fit$Report[[names(par_hat)[int]]])
    }
    if ("units" %in% class(fit$Report[[names(par_hat)[int]]])) 
      units(par_hat[[int]]) = units(fit$Report[[names(par_hat)[int]]])
  }
  
  Array_ctl = par_hat[[index_name]]
  log_Array_ctl = par_SE[[log_index_name]]
  
  Table = cbind(expand.grid(dimnames(par_hat[[index_name]])), 
              Units = make_unit_label(u = units(par_hat[[index_name]]), 
                                      lab = "", parse = FALSE), 
              Estimate = as.vector(par_hat[[index_name]]), 
              Std.Err.Est = as.vector(par_SE[[index_name]]), 
              Std.Err.ln.Est = 
                as.vector(par_SE[[log_index_name]]))
  
  Table <- Table %>% 
    separate(Time, into=c('Year', 'Season')) %>% 
    mutate(Year = as.numeric(Year),
           Season = factor(Season, levels=c("Spring", 'Fall')),
           Stratum = factor(Stratum,
                            levels=c('ALL', 'EGOM', 'GBK',
                                     'SNE', 'WGOM'))) %>% 
    dplyr::select(-Units) %>% 
    as.data.frame()
  
  write.csv(Table, file = file.path(DirName, 
                                    paste0(category_names, 
                                           "_Index.csv")), 
            row.names = FALSE)

  
  p.fall <- ggplot() + 
    geom_point(data=Table[Table$Season == 'Fall' & Table$Stratum != 'ALL',],
               aes(x=Year, y=Estimate, col=Stratum),
               cex=0.4, alpha=0.6) +
    geom_line(data=Table[Table$Season == 'Fall'& Table$Stratum != 'ALL',],
              aes(x=Year, y=Estimate, col=Stratum),
              lwd=0.6, alpha=0.8) +
    geom_ribbon(data=Table[Table$Season == 'Fall'& Table$Stratum != 'ALL',],
                aes(x=Year, ymin=Estimate-Std.Err.Est,
                    ymax=Estimate+Std.Err.Est,
                    fill=Stratum),
                alpha=0.3) +
    labs(y='Relative abundance') 
  
  p.spring <- ggplot() + 
    geom_point(data=Table[Table$Season == 'Spring'& Table$Stratum != 'ALL',],
               aes(x=Year, y=Estimate, col=Stratum),
               cex=0.4, alpha=0.6) +
    geom_line(data=Table[Table$Season == 'Spring'& Table$Stratum != 'ALL',],
              aes(x=Year, y=Estimate, col=Stratum),
              lwd=0.6, alpha=0.8) +
    geom_ribbon(data=Table[Table$Season == 'Spring'& Table$Stratum != 'ALL',],
                aes(x=Year, ymin=Estimate-Std.Err.Est,
                    ymax=Estimate+Std.Err.Est,
                    fill=Stratum),
                alpha=0.3) +
    labs(y='Relative abundance')
  
  gfall <- ggplot_build(p.fall)$layout$panel_params[[1]]$y$breaks
  gspring <- ggplot_build(p.spring)$layout$panel_params[[1]]$y$breaks
  
  log10_floor <- function(x) {
    10^(floor(log10(x)))
  }
  
  if(gfall[2] > gspring[2]){
    scaling.fac <- log10_floor(gspring[2]) 
  }
  if(gfall[2] < gspring[2]){
    scaling.fac <- log10_floor(gfall[2])
  }
  
  p.spring <- p.spring + 
    scale_y_continuous(labels = gspring / scaling.fac)
  
  p.fall <- p.fall +
    scale_y_continuous(labels = gfall / scaling.fac)
  
  p.both <- egg::ggarrange(
    p.spring +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            legend.position = 'n') +
      facet_wrap(vars(Season)),
    p.fall + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      facet_wrap(vars(Season)), 
    nrow = 2)
  
  ggsave(plot=p.both,
         filename=paste0(DirName, '/', 
                         category_names, '_Seasonal_Indices.png'),
         width=11, height=8.5, units = 'in')
  

}