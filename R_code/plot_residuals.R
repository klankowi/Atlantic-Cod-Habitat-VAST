rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)
library(wesanderson)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "inside",
                legend.background = element_rect(fill='transparent', colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load data
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))
setwd(here("VAST_runs/large/Overall_BC/ALL"))
rm(list=setdiff(ls(), c("fit")))

# Reload model
fit_new = reload_model(x = fit)

# Residual plots
FishStatsUtils::plot_quantile_residuals(fit=fit_new)
