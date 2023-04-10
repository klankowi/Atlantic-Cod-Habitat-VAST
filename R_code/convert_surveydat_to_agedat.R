rm(list=ls())

# Libraries and operators
library(tidyverse)
library(here)
'%notin%' <- function(x,y)!('%in%'(x,y))
source(here("R_code/utilities/true_seasons_func.R"))

#### Read in data
# Survey data
surv <- read.csv(here("Data/Survey_Data/Survey_Data.csv"))
surv <- subset(surv, YEAR >=1982)

# Bio data
bioda <- read.csv(here("Data/Survey_Data/Bio_Data.csv"))
bioda <- subset(bioda, YEAR>=1982)

# Append haul_id and date to bio data
surv.sub <- dplyr::select(surv, INDEX_NAME, HAUL_ID, DATE)
bioda <- left_join(bioda, surv.sub, by=c('HAUL_ID'))

# There are several ways to describe fish: age, weight, and length. If/ which
# descriptors are provided varies by survey. There are also, predictably, a 
# few observations that have missing data due to data collector error or other
# at-sea reasons. We are interested in splitting the overall sampled fish pool
# into three age categories: 0-2, 2-5, 5+. We are also interested in preserving
# all data for inclusion in VAST, which can integrate encounter, count, and 
# biomass data into a single index of abundance.

# To achieve this goal, we will combine tow information (`survey`) and bio data
# (`bioda`). The first split will parse out data with weight information per 
# fish.

# Biodata with weight (biomass)
biomass <- bioda[is.na(bioda$WEIGHT)==FALSE,]
biomass$DATE <- as.POSIXct(biomass$DATE, format="%m/%d/%Y %H:%M")
for(i in 1:nrow(biomass)){
  biomass$TRUE_SEASON[i] <- true_seasons(biomass$DATE[i])
}
biomass$SEASON <- NULL
biomass$SEX[biomass$SEX == "UNKNOWN"] <- NA
biomass$SEX <- droplevels(as.factor(biomass$SEX))
biomass$SEX.SEASON <- paste0(biomass$SEX, "_", biomass$TRUE_SEASON)
biomass$SEX.SEASON <- factor(biomass$SEX.SEASON,
                             levels=c('FEMALE_WINTER', 'FEMALE_SPRING', 'FEMALE_SUMMER', 'FEMALE_FALL',
                                      'MALE_WINTER', 'MALE_SPRING', 'MALE_SUMMER', 'MALE_FALL',
                                      'NA_WINTER', 'NA_SPRING', 'NA_SUMMER', 'NA_FALL'))
# Pull out bad values
error999 <- subset(biomass, WEIGHT == 9.999)
biomass <- subset(biomass, WEIGHT != 9.999)
bio.split <- split(biomass, f=biomass$SEX.SEASON)
for(i in 1:length(bio.split)){
  plotsplit <- bio.split[[i]]
  plotsplit <- plotsplit[with(plotsplit, order(LENGTH_CM, WEIGHT)),]
  row.names(plotsplit) <- NULL
  
  ggplot(plotsplit) +
    geom_point(data=plotsplit, aes(x=AGE, y=WEIGHT, col=TRUE_SEASON), alpha=0.9) +
    facet_grid(~SEX.SEASON) +
    geom_smooth(data=plotsplit, aes(x=AGE, y=WEIGHT, col=TRUE_SEASON))
}

ggplot(biomass) +
  geom_point(data=biomass, aes(x=LENGTH_CM, y=WEIGHT, col=TRUE_SEASON), alpha=0.9) +
  facet_grid(~SEX.SEASON) +
  geom_smooth(data=biomass, aes(x=LENGTH_CM, y=WEIGHT, col=TRUE_SEASON))

# Biodata without weight (count)
countda <- bioda[is.na(bioda$WEIGHT)==TRUE,]
countda <- rbind(countda, error999)

# Survey data without bio data (encounter)
encounter <- surv[surv$SURVEY %notin% bioda$SURVEY,]
head(encounter)

# Remove intermediates
rm(surv.sub)

# Pull surveys from biodata
all.surveys <- unique(bioda$SURVEY)
surv <- subset(surv, SURVEY %in% all.surveys)

# Make list
surv.list <- split(surv, f=surv$SURVEY)

# Loop through surveys
for(f in 1:length(unique(bioda$SURVEY))){
  surv.use <- all.surveys[f]
  bioda.use <- subset(bioda, SURVEY == surv.use)
  surv.use.df <- subset(surv, SURVEY == surv.use)
  
  # Age 0-2 cod
  surv.use.df.a0 <- surv.use.df
  surv.use.df.a0$AGE <- 'A0-2'
  for(s in 1:nrow(surv.use.df.a0)){
    surv.use.df.a0$COD_N[s] <- nrow(bioda.use[bioda.use$HAUL_ID == 
                                                surv.use.df$HAUL_ID[s] &
                                                bioda.use$LENGTH_CM <=25,])
    surv.use.df.a0$COD_KG[s] <- nrow(bioda.use[bioda.use$HAUL_ID == 
                                                surv.use.df$HAUL_ID[s] &
                                                bioda.use$LENGTH_CM <=25,])
  }
  
  # Age 2-5 cod
  surv.use.df.a1 <- surv.use.df
  surv.use.df.a1$AGE <- 'A2-5'
  for(s in 1:nrow(surv.use.df.a1)){
    surv.use.df.a1$COD_N[s] <- nrow(bioda.use[bioda.use$HAUL_ID == 
                                            surv.use.df$HAUL_ID[s] &
                                            bioda.use$LENGTH_CM > 25 &
                                            bioda.use$LENGTH_CM <=63,])
  }
  
  # Age 5+ cod
  surv.use.df.a2 <- surv.use.df
  surv.use.df.a2$AGE <- 'A5PLUS'
  for(s in 1:nrow(surv.use.df.a2)){
    surv.use.df.a2$COD_N[s] <- nrow(bioda.use[bioda.use$HAUL_ID == 
                                                surv.use.df$HAUL_ID[s] &
                                                bioda.use$LENGTH_CM >63,])
  }
  
  # Combine
  surv.age.df <- rbind(surv.use.df.a0, surv.use.df.a1, surv.use.df.a2)
  
  surv.age.df <- surv.age.df[with(surv.age.df, order(YEAR, HAUL_ID, AGE)),]
  row.names(surv.age.df) <- NULL
  
  surv.list[[f]] <- surv.age.df
}

surv.final <- do.call(rbind, surv.list)
surv.final <- surv.final[with(surv.final, order(SURVEY, HAUL_ID, AGE)),]
row.names(surv.final) <- NULL
head(surv.final)
write.csv(surv.final, here('Data/Survey_Data/Survey_Data_AgeSep.csv'), row.names=F)
