
rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)

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

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Root
root <- here('VAST_runs/large/SingleCovs/')
ending <- '_largecod_wholearea_natsplin_fsOFF.RData'

# Covariates
covuse <- c('amo', 'bathy', 'bottomtemp', 'cobble', 'gravel',
            'mud', 'nao', 'rugos', 'sand')

# Baseline Mods
load(paste0(root, '/all/allcovs', ending))
all <- fit
load(paste0(root, '/none/nonecovs', ending))
none <- fit

# Singlfitted.values()# Single Covs
for(i in 1:length(covuse)){
  load(paste0(root, '/', covuse[i], 
                     '/only', covuse[i], ending))
  newName <- covuse[i]
  assign(newName, fit, envir=.GlobalEnv)
  rm(fit, newName)
}

rm(list=setdiff(ls(), c('all', 'none', 
                        'amo', 'bathy', 'bottomtemp',
                        'cobble', 'gravel', 'mud',
                        'nao', 'rugos', 'sand')))

aiccomp <- data.frame(
  model=c('amo', 'bathy', 'bottomtemp',
          'cobble', 'gravel', 'mud',
          'nao', 'rugos', 'sand',
          'none', 'all'),
  aic = c(amo$parameter_estimates$AIC,
          bathy$parameter_estimates$AIC,
          bottomtemp$parameter_estimates$AIC,
          cobble$parameter_estimates$AIC,
          gravel$parameter_estimates$AIC,
          mud$parameter_estimates$AIC,
          nao$parameter_estimates$AIC,
          rugos$parameter_estimates$AIC,
          sand$parameter_estimates$AIC,
          none$parameter_estimates$AIC,
          all$parameter_estimates$AIC)
)

aiccomp <- aiccomp[with(aiccomp, order(aic)),]

aiccomp$delta.aic <- NA
for(i in 2:nrow(aiccomp)){
  aiccomp$delta.aic[i] <- 
    aiccomp$aic[i] - aiccomp$aic[(i-1)]
}
aiccomp$delta.aic[1] <- 0
aiccomp

devcomp <- data.frame(
  model=c('amo', 'bathy', 'bottomtemp',
          'cobble', 'gravel', 'mud',
          'nao', 'rugos', 'sand',
          'none', 'all'),
  deviance = c(amo$Report$deviance,
          bathy$Report$deviance,
          bottomtemp$Report$deviance,
          cobble$Report$deviance,
          gravel$Report$deviance,
          mud$Report$deviance,
          nao$Report$deviance,
          rugos$Report$deviance,
          sand$Report$deviance,
          none$Report$deviance,
          all$Report$deviance)
)

devcomp <- devcomp[with(devcomp, order(deviance, decreasing = T)),]

devcomp$pct.dev.exp <- NA
for(i in 1:nrow(devcomp)){
  devcomp$pct.dev.exp[i] <- 
    (1 - devcomp$deviance[i]/none$Report$deviance) * 100
}
devcomp
