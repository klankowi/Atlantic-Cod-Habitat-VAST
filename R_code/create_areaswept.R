# Area swept variable

rm(list=ls())
library(tidyverse)
library(here)

smast <- read.csv('C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Mod_Files_3/Video_Trawl_SMAST_AREASWEPT.csv')

madmf <- read.csv('C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Mod_Files_3/MADMF_Industry_AREASWEPT.csv')

smast <- dplyr::select(smast, HAUL_ID, DATE, AREA_SWEPT)

madmf <- dplyr::select(madmf, HAUL_ID, DATE, AREA_SWEPT)

as <- rbind(smast, madmf)
str(as)
as$DATE <- as.POSIXct(as$DATE,
                      format='%Y-%m-%d')

all <- read.csv(here('Data/Survey_Data/Survey_Data.csv'))
str(all)
all$DATE <- as.POSIXct(all$DATE,
                       format='%m/%d/%Y')

all <- left_join(all, as, by=c('HAUL_ID', 'DATE'))
table(all$SURVEY[is.na(all$AREA_SWEPT)])

all$AREA_SWEPT[all$SURVEY == 'ASMFC Shrimp Trawl'] <- 0.0120
all$AREA_SWEPT[all$SURVEY == 'DFO Trawl'] <- 0.0404
all$AREA_SWEPT[all$SURVEY == 'GSO Trawl'] <- 0.0120
all$AREA_SWEPT[all$SURVEY == 'MADMF Inshore Trawl'] <- 0.0132
all$AREA_SWEPT[all$SURVEY == 'ME-NH Inshore Trawl'] <- 0.0148
all$AREA_SWEPT[all$SURVEY == 'NEFSC BLLS' |
                 all$SURVEY == 'NEFSC BTS'] <- 0.0384
all$AREA_SWEPT[all$SURVEY == 'RIDEM Trawl'] <- 0.0103
all$AREA_SWEPT[all$SURVEY == 'Sentinel'] <- NA

write.csv(all,
          here("Background_Info/area_swept_perhaul.csv"),
          row.names = F)
