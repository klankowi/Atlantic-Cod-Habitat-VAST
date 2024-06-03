
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
root <- here('VAST_runs/structuring/large')
ending <- '.RData'

load(paste0(root, '/alleffectson/alleffectson', ending))
alleffectson <- fit

load(paste0(root, '/noaniso/noaniso', ending))
noaniso <- fit

load(paste0(root, '/noomeps2/noomeps2', ending))
noomeps2 <- fit

load(paste0(root, '/noomeps2_noaniso/noomeps2_noaniso', ending))
noomeps2_noaniso <- fit

load(paste0(root, '/noomeps2_noeps1/noomeps2_noeps1', ending))
noomeps2_noeps1 <- fit

load(paste0(root, '/noomeps2_noeps1_noaniso/noomeps2_noeps1_noaniso', ending))
noomeps2_noeps1_noaniso <- fit

load(paste0(root, '/noomeps12/noomeps12', ending))
noomeps12 <- fit

load(paste0(root, '/noomeps12_noaniso/noomeps12_noaniso', ending))
noomeps12_noaniso <- fit

rm(list=setdiff(ls(), c('alleffectson', 'noaniso', 
                        'noomeps2', 'noomeps2_noaniso',
                        'noomeps2_noeps1', 'noomeps2_noeps1_noaniso',
                        'noomeps12', 'noomeps12_noaniso'
                        )))

aiccomp <- data.frame(
  model=c('alleffectson', 'noaniso', 
          'noomeps2', 'noomeps2_noaniso',
          'noomeps2_noeps1', 'noomeps2_noeps1_noaniso',
          'noomeps12', 'noomeps12_noaniso'),
  aic = c(alleffectson$parameter_estimates$AIC,
          noaniso$parameter_estimates$AIC,
          noomeps2$parameter_estimates$AIC,
          noomeps2_noaniso$parameter_estimates$AIC,
          noomeps2_noeps1$parameter_estimates$AIC,
          noomeps2_noeps1_noaniso$parameter_estimates$AIC,
          noomeps12$parameter_estimates$AIC,
          noomeps12_noaniso$parameter_estimates$AIC)
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
  model=c('alleffectson', 'noaniso', 
          'noomeps2', 'noomeps2_noaniso',
          'noomeps2_noeps1', 'noomeps2_noeps1_noaniso',
          'noomeps12', 'noomeps12_noaniso'),
  deviance = c(alleffectson$Report$deviance,
                noaniso$Report$deviance,
                noomeps2$Report$deviance,
                noomeps2_noaniso$Report$deviance,
                noomeps2_noeps1$Report$deviance,
                noomeps2_noeps1_noaniso$Report$deviance,
                noomeps12$Report$deviance,
                noomeps12_noaniso$Report$deviance)
)

devcomp <- devcomp[with(devcomp, order(deviance, decreasing = T)),]

devcomp$pct.dev.exp <- NA
for(i in 1:nrow(devcomp)){
  devcomp$pct.dev.exp[i] <- 
    (1 - devcomp$deviance[i]/noomeps12_noaniso$Report$deviance)
}
devcomp
