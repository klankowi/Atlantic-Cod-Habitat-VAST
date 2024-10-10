rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)

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

rm(list=setdiff(ls(), c('small', 'medium', 'large')))

# Set WD for plotting
out_dir = here("VAST_runs/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

sizes <- c('small', 'medium', 'large')

keep <- data.frame(
  fit=NA, se=NA, lower=NA, upper=NA, Lin_pred=NA, Value=NA,
  Covariate=NA, Size=NA
)

for(i in 1:length(sizes)){
  size <- sizes[i]
  if(size == 'small'){fit <- small}
  if(size == 'medium'){fit <- medium}
  if(size == 'large'){fit <- large}
  
  # Call parameter names
  params <- colnames(fit$covariate_data)
  params <- params[params %notin% c('Lon', 'Lat', 'Year')]
  
  if(size == 'small'){params <- params[params %notin% c('amo')]}
  if(size == 'large'){params <- params[params %notin% c('cobble_P',
                                                        'sand_P', 
                                                        'mud_P')]}
  #params <- params[params %notin% c('amo')]
  
  
  # Call category names
  catnames <- fit$category_names
  ncat = length(catnames)
  
  vast_covariate_effects<- get_vast_covariate_effects(vast_fit = fit, 
                                                      params_plot = c(params), 
                                                      params_plot_levels = 100, 
                                                      effects_pad_values = c(), 
                                                      nice_category_names = paste0(size, ' Cod'),
                                                      out_dir = out_dir,
                                                      category_to_use = 1,
                                                      ncat = ncat)
  
  names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")
  
  betternames=data.frame(
    Covariate = c('cobble_P', 'gravel_P', 'mud_P', 'sand_P', 'rugos', 'BATHY.DEPTH',
                  'h_bt', 'nao', 'amo'),
    Better = c('Cobble', 'Gravel', 'Mud', 'Sand', 'Rugosity', 'Depth (m)', 
               'Bottom Temp (C)', 'Daily NAO', 'Monthly AMO')
  )
  
  vast_cov_eff_l <- vast_covariate_effects %>%
    drop_na(Value)
  
  vast_cov_eff_l <- merge(vast_cov_eff_l, betternames)
  
  vast_cov_eff_l <- vast_cov_eff_l %>% 
    dplyr::select(-Covariate) %>% 
    rename(Covariate = Better)
  
  vast_cov_eff_l$Size <- size
  
  keep <- rbind(keep, vast_cov_eff_l)
  
  rm(size, params, catnames, ncat, vast_covariate_effects, names_stay,
     betternames, vast_cov_eff_l)
  
}

keep <- keep[!is.na(keep$Size),]

keep <- keep %>% 
  mutate(Lin_pred = factor(Lin_pred, levels=c('X1', 'X2')),
         Covariate = as.factor(Covariate),
         Size = factor(Size, levels=c('small', 'medium', 'large')))

keepx1 <- keep[keep$Lin_pred == 'X2',]
keepx1$Covariate <- droplevels(keepx1$Covariate)
keepx1 <- split(keepx1, f=keepx1$Covariate)
for(i in 1:length(keepx1)){
  keepx1[[i]] <- split(keepx1[[i]], f=keepx1[[i]]$Size)
  keepx1[[i]] <- keepx1[[i]][sapply(keepx1[[i]], nrow)>0]
}

for(i in 1:length(keepx1)){
  for(j in 1:length(keepx1[[i]])){
    keepx1[[i]][[j]]$dif <- NA
    keepx1[[i]][[j]]$pos <- NA
    for(k in 2:nrow(keepx1[[i]][[j]])){
      keepx1[[i]][[j]]$dif[k] <- 
        keepx1[[i]][[j]]$fit[k] - keepx1[[i]][[j]]$fit[(k-1)]
      if(keepx1[[i]][[j]]$dif[k] > 0){keepx1[[i]][[j]]$pos[k] <- 'Positive'}
      if(keepx1[[i]][[j]]$dif[k] < 0){keepx1[[i]][[j]]$pos[k] <- 'Negative'}
        
    }
  }
}


for(i in 1:length(keepx1)){
  for(j in 1:length(keepx1[[i]])){
    message(paste0(names(keepx1[i]), ' ',
                 names(keepx1[[i]][j])))
  print(table(keepx1[[i]][[j]]$pos))
  }
  }
