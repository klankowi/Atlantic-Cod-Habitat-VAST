
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
root <- here('VAST_runs/medium/GMRF/')
ending <- '_mediumcod_wholearea_natsplin.RData'

# Covariates
covuse <- c('OFF', 'FULL_OFF')

# Singlfitted.values()# Single Covs
for(i in 1:length(covuse)){
  load(paste0(root, '/', covuse[i], 
                     '/GMRF_', covuse[i], ending))
  newName <- covuse[i]
  assign(newName, fit, envir=.GlobalEnv)
  rm(fit, newName)
}

rm(list=setdiff(ls(), c('OFF', 'FULL_OFF')))

aiccomp <- data.frame(
  model=c('OFF', 'FULL_OFF'),
  aic = c(OFF$parameter_estimates$AIC,
          FULL_OFF$parameter_estimates$AIC)
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
  model=c('OFF','FULL_OFF'),
  deviance = c(
          OFF$Report$deviance,
          FULL_OFF$Report$deviance)
)

devcomp <- devcomp[with(devcomp, order(deviance, decreasing = T)),]

devcomp$pct.dev.exp <- NA
for(i in 1:nrow(devcomp)){
  devcomp$pct.dev.exp[i] <- 
    (1 - devcomp$deviance[i]/OFF$Report$deviance) * 100
}
devcomp
