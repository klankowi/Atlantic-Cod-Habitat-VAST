#### Workspace setup ####
# Clear workspace
rm(list=ls())

# Load libraries
# IMPORTANT NOTE: VAST must be running >=V14, will not work with V13.
library(TMB)
library(units)
library(VAST)
library(here)
library(tidyverse)
library(beepr)
library(sf)
library(rgdal)
library(sp)
library(ggcorrplot)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Add sample information and covars ####
# Load data
surveys <- read.csv(here("Data/VAST_input/cod_agesep_VASTdata.csv"))

# Remove data from 2022 (incomplete)
ex <- subset(surveys, YEAR < 2022)

# Set response to count data
ex$RESPONSE <- ceiling(ex$AGE_N)

# Remove Sentinel survey
ex <- subset(ex, SURVEY !='Sentinel')

# Check that there are no missing responses
nrow(ex[is.na(ex$RESPONSE)==TRUE,])

# Add environmental data
env <- readRDS(here("Data/VAST_input/agg_stn_all_OISST.RDS"))
env <- sfheaders::sf_to_df(env, fill=T)
env <- dplyr::select(env,
                     HAUL_ID, cobble_P, gravel_P, rock_P, mud_P, sand_P,
                     BATHY.DEPTH, oisst)
ex2 <- left_join(ex, env, by="HAUL_ID")

# Add annoying rugosity data
rugos <- readRDS(here("Data/VAST_input/agg_stn_all_OISST_agesep.RDS"))
rugos <- sfheaders::sf_to_df(rugos, fill=T)
rugos <- dplyr::select(rugos, HAUL_ID, rugos)
rugos <- unique(rugos)

ex3 <- left_join(ex2, rugos, by="HAUL_ID")

ex4 <- subset(ex3, !is.na(ex3$rugos))

# Add area swept for Longline survey
# Sets are 1nmi (1.852 km). Soak time is 2 hrs.
# McElroy et al. (2019) states parameters chosen to approximate NEFSC BTS
# So set area swept to equal bts average (0.0384 sq km)
table(ex4$SURVEY[is.na(ex4$AREA_SWEPT)])
ex4$AREA_SWEPT[is.na(ex4$AREA_SWEPT)] <- 0.0384

# Check results
str(ex4)
table(ex4$YEAR)

# Add NAO
nao <- read.csv(here('Data/Climate_Indices/daily_NAO.csv'))
nao$DATE <- paste0(nao$year, '-',
                      str_pad(nao$month, 2, "left", "0"), "-",
                      str_pad(nao$day, 2, "left", '0'))
nao <- dplyr::select(nao, DATE, aao_index_cdas)
colnames(nao) <- c('DATE', 'nao')
nao$DATE <- as.POSIXct(nao$DATE, format='%Y-%m-%d')
nao$DATE <- as.Date(nao$DATE)

ex4 <- ex

ex4$DATE <- as.POSIXct(ex4$DATE, 
                       format='%m/%d/%Y')
ex4$DATE <- as.Date(ex4$DATE)
ex4 <- left_join(ex4, nao, by=c("DATE"))
summary(ex4$nao)

# NAO not provided for Apr 30 2003, must remove

ex4 <- ex4 %>% 
  drop_na(nao)

# Merge AMO
amo <- read.csv(here('Data/Climate_Indices/monthly_amo.csv'))
amo$monthno <- match(amo$Month,month.abb)
amo$yrmo <- paste0(amo$Year, '-', 
                   str_pad(amo$monthno, 2, 'left', '0'))
amo <- dplyr::select(amo, yrmo, Value)
colnames(amo) <- c('yrmo', 'amo')
ex4$month <- lubridate::month(ex4$DATE)
ex4$yrmo <- paste0(ex4$YEAR, '-',
                   str_pad(ex4$month, 2, 'left', '0'))
ex4 <- left_join(ex4, amo, by=c('yrmo'))
summary(ex4$amo)

# Merge GSI
gsi <- as.data.frame(ecodata::gsi)
gsi <- separate(gsi, Time, c('year', 'month'))
gsi$month <- str_pad(gsi$month, 2, 'right', '0')
gsi$yrmo <- paste0(gsi$year, '-', gsi$month)
gsi <- dplyr::select(gsi, yrmo, Value)
colnames(gsi) <- c('yrmo', 'gsi')
ex4 <- left_join(ex4, gsi, by=c('yrmo'))
summary(ex4$gsi)

# Remove missing h_bt
ex4 <- ex4 %>% 
  drop_na(h_bt)

# Save
write.csv(ex4, row.names = F,
          here('Data/VAST_input/cod_agesep_VASTdata.csv'))
