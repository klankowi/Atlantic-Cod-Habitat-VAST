# There are several ways to describe fish: age, weight, and length. If/ which
# descriptors are provided varies by survey. There are also, predictably, a 
# few observations that have missing data due to data collector error or other
# at-sea reasons. We want to preserve  all data for inclusion in VAST, which can
# integrate encounter, count, and biomass data into a single index of abundance.
# This is possible, but requires some reshaping of the data. We will be able to
# separate and label data that include biomass, data that only include 
# abundance, and data that only reflect presence/ absence.

# We also want to split these data by age groupings, (Ages 0-2, Ages 2-5, 
# Ages 5+). We suspect these age groupings use different habitats, etc. 

# The following script reshapes the survey and bio data to create these branched
# groupings: Biomass data (Ages 0-2, Ages 2-5, Ages 5+), Abundance data (same
# age groupings), and Encounter data (same age groupings).

#### Set up workspace ####
rm(list=ls())

# Libraries and operators
library(tidyverse)
library(FSA)
library(nlstools)
library(car)
library(here)

'%notin%' <- function(x,y)!('%in%'(x,y))

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

#### Read in data ####
# Survey data
surv <- read.csv(here("data/Dataframes/Survey_Data.csv"))
surv <- subset(surv, YEAR >=1982)
surv <- dplyr::select(surv, -SEASON, -TRUE_SEASON)
# Must have spatial information
badsurv <- subset(surv, is.na(surv$LON) | is.na(surv$LAT))
surv <- subset(surv, !is.na(surv$LAT))
surv <- subset(surv, !is.na(surv$LON))
# Assign season
surv$DATE <- as.POSIXct(surv$DATE, format = "%m/%d/%Y %H:%M")
surv$month <- lubridate::month(surv$DATE)
surv$SEASON <- NA
surv$SEASON[surv$month %in% c(3,4,5,6,7,8)] <- 'SPRING'
surv$SEASON[surv$month %in% c(1,2,9,10,11,12)] <- 'FALL'
surv$month <- NULL
# Convert character to factor
facs <- c('INDEX_NAME', 'SURVEY', 'STOCK', 'AREA', 'STRATUM',
          'BOTTOM.TYPE', 'SEASON')
surv[,facs] <- lapply(surv[,facs], factor)  
surv$SEASON <- factor(surv$SEASON, levels = c("SPRING", "FALL"))
str(surv)

# Bio data
bioda <- read.csv(here("data/Dataframes/Bio_Data.csv"))
bioda <- subset(bioda, YEAR>=1982)
bioda <- dplyr::select(bioda, -SEASON, -TRUE_SEASON)

# Append haul_id and date to bio data
surv.sub <- dplyr::select(surv, INDEX_NAME, HAUL_ID, DATE, SEASON)
bioda <- left_join(bioda, surv.sub, by=c('HAUL_ID'))
# Remove biodata from hauls without spatial information
bioda <- bioda[bioda$HAUL_ID %notin% badsurv$HAUL_ID,]

# Combine NA and UKNOWN sex
bioda$SEX[bioda$SEX == "UNKNOWN"] <- NA
bioda$SEX <- droplevels(as.factor(bioda$SEX))
# Combine NA and UNKNOWN maturity stage
bioda$MATURITY_STAGE[bioda$MATURITY_STAGE == "UNKNOWN"] <- NA
# Convert character to factor
facs <- c('SURVEY', 'STOCK', 'SEX', 'MATURITY_STAGE',
          'INDEX_NAME', 'SEASON')
bioda[,facs] <- lapply(bioda[,facs], factor)
str(bioda)

# Remove intermediates
rm(surv.sub, badsurv, facs)

#### von Bertalanffy growth curve ####
# 
# # All fish with recorded biodata have length. Weight and age are more variable.
# # Step one of analysis is to fit a von Bert growth curve to determine age based
# # on length.
# 
# # Split out data with age and length
# bioage <- subset(bioda, is.na(bioda$AGE)==FALSE)
# summary(bioage$AGE); summary(bioage$LENGTH_CM)
# 
# # Define von Bertalanffy growth function
# vbmod <- LENGTH_CM ~ Linf * (1 - exp(-K * (AGE - t0)))
# 
# # Get starting values for each of your parameters using the `vbStarts` function
# # from `FSA`
# starts <- vbStarts(formula = LENGTH_CM ~ AGE, data = bioage)
# 
# # Fit the von Bertalanffy growth function using nonlinear least squares (nls) 
# # optimization
# mymod <- nls(vbmod, data = bioage, start = starts)
# summary(mymod)
# 
# # Get the desired growth function from a list of those that are available in FSA
# vbO <- vbFuns("typical")
# 
# # Fit the model to the data using nls, like we did before
# vb_fit <- nls(LENGTH_CM~vbO(AGE,Linf,K, t0), data=bioage, start=starts)
# 
# # Now, bootstrap the model fitting process
# boot_fit <- nlsBoot(vb_fit)
# 
# # Predict length at age from the model and calculate some bootstrapped CIs
# boot_preds <- data.frame(
#   predict(boot_fit, vbO, t = sort(unique(bioage$AGE))))
# names(boot_preds) <- c("AGE", "fit", "lwr", "upr")
# 
# # Merge CIs to predictions
# age_preds <- merge(bioage, boot_preds, by = "AGE")
# 
# # Plot
# ggplot(age_preds, aes(x = AGE, y = LENGTH_CM)) +
#   geom_jitter(data=age_preds, width = 0.1, alpha = 0.15, size = 2,
#               aes(col=SEX)) +
#   geom_line(aes(y = fit)) +
#   geom_ribbon(
#     aes(x = AGE, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
#   xlab("Age (years)") +
#   ylab("Total length (mm)")
# cuts2 <- round(boot_preds$fit[boot_preds$AGE == 2], 2)
# cuts5 <- round(boot_preds$fit[boot_preds$AGE == 5], 2)
# 
# rm(age_preds, bioage, boot_fit, boot_preds, mymod, starts, vb_fit,
#    vbmod, vbO)

cuts2 <- 39.17
cuts5 <- 70.16

#### Length at weight ####
# age.bio <- subset(bioda, !is.na(bioda$LENGTH_CM) &
#                     !is.na(bioda$WEIGHT) & bioda$WEIGHT != 0)
# 
# # Plot for posterity
# ggplot(data=age.bio, aes(x=LENGTH_CM, y=WEIGHT)) +
#   geom_point()
# # Some of these are biologically implausible and need to be removed.
# 
# # Create flag for implausible values
# age.bio$FLAG <- 0
# 
# # Set flag for implausible values
# age.bio$FLAG[age.bio$WEIGHT == 9.999] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 3.000 & age.bio$LENGTH_CM == 32] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 4.936 & age.bio$LENGTH_CM == 46] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 5.440 & age.bio$LENGTH_CM == 53] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 7.898 & age.bio$LENGTH_CM == 57] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 0.020 & age.bio$LENGTH_CM == 53] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 0.200 & age.bio$LENGTH_CM == 57] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 0.345 & age.bio$LENGTH_CM == 67] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 0.385 & age.bio$LENGTH_CM == 73] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 1.838 & age.bio$LENGTH_CM == 95] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 0.270 & age.bio$LENGTH_CM == 119] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 12.38 & age.bio$LENGTH_CM == 142] <- 1
# age.bio$FLAG[age.bio$WEIGHT == 22.04 & age.bio$LENGTH_CM == 161] <- 1
# 
# # Identify bad biomass values
# bad.biomass <- subset(age.bio, FLAG ==1)
# 
# # Remove bad biomass values
# age.bio <- subset(age.bio, FLAG ==0)
# 
# # Replot
# ggplot(data=age.bio, aes(x=LENGTH_CM, y=WEIGHT, col=SEX)) +
#   geom_point(alpha=0.5)
# # Remaining values look to be reasonable
# 
# # Linear regression of length at weight
# age.bio$logL <- log(age.bio$LENGTH_CM)
# age.bio$logW <- log(age.bio$WEIGHT)
# lm1 <- lm(logW ~ logL, data=age.bio)
# summary(lm1)
# 
# # Find correction factor
# syx <- summary(lm1)$sigma
# cf <- exp((syx^2)/2)
# 
# # Predict weight at length == age 2
# pred2.log <- predict(lm1,data.frame(logL=log(cuts2)),interval="c")
# bias2.pred.orig <- exp(pred2.log)
# pred2.orig <- cf*bias2.pred.orig
# wts2 <- round(pred2.orig[1],2)
# 
# # Predict weight at length == age5
# pred5.log <- predict(lm1,data.frame(logL=log(cuts5)),interval="c")
# bias5.pred.orig <- exp(pred5.log)
# pred5.orig <- cf*bias5.pred.orig
# wts5 <- round(pred5.orig[1],2)
# 
# rm(pred5.orig, bias5.pred.orig, pred5.log, pred2.orig, bias2.pred.orig,
#    pred2.log, lm1, age.bio, bad.biomass, cf, syx)

wts2 <- 0.58
wts5 <- 3.44

#### Assign age by length ####
# Create new age grouping column
bioda$LENAGE <- NA

# Loop through rows
for(i in 1:nrow(bioda)){
  
  # Age 0-2 cutoff
  if(bioda$LENGTH_CM[i] <= cuts2){
    bioda$LENAGE[i] <- 0
  }
  # Ages 2-5
  if(bioda$LENGTH_CM[i] > cuts2 & bioda$LENGTH_CM[i] <=cuts5){
    bioda$LENAGE[i] <- 3
  }
  # Age 5 cutoff
  if(bioda$LENGTH_CM[i] > cuts5){
    bioda$LENAGE[i] <- 6
  }
}

# View results
table(bioda$LENAGE)
nrow(bioda[is.na(bioda$LENAGE)==TRUE,])

# Pull out data that has provided age
bioage <- subset(bioda, is.na(bioda$AGE)==FALSE)
# Reassign age for those data
bioage$LENAGE <- bioage$AGE
# Pull out data that does not have provided age
bioage.noage <- subset(bioda, is.na(bioda$AGE)==TRUE)
# Recombine to make full dataset
bioda <- rbind(bioage, bioage.noage)

# View results
table(bioda$LENAGE)

# Create new age grouping column
bioda$AGEGROUP <- NA

# Loop through rows
for(i in 1:nrow(bioda)){
  # Age 0-2 cutoff
  if(bioda$LENAGE[i] < 2){
    bioda$AGEGROUP[i] <- 'Age0-2'
  }
  # Age 2-5
  if(bioda$LENAGE[i] >= 2 & bioda$LENAGE[i] <5){
    bioda$AGEGROUP[i] <- 'Age2-5'
  }
  # Age 5 cutoff
  if(bioda$LENAGE[i] >= 5){
    bioda$AGEGROUP[i] <- 'Age5+'
  }
}

# View results
table(bioda$AGEGROUP)
table(bioda$AGE, bioda$AGEGROUP)
nrow(bioda[is.na(bioda$AGEGROUP)==TRUE,])
head(bioda)

# Remove intermediates
rm(bioage, bioage.noage, i)

#### Combine biodata and survey data for known age groups ####
# Split by haul
bioda.list <- split(bioda, f=bioda$HAUL_ID)

# Loop through hauls
for(i in
    1:
    #8000:
    length(bioda.list)
    #10000
    ){
  # Pull out biodata associated with haul
  biotemp <- bioda.list[[i]]
  head(biotemp)
  
  # Pull out survey data associated with haul
  survtemp <- surv[surv$HAUL_ID == biotemp$HAUL_ID[1],]
  head(survtemp)
  
  # Create row for each age group
  survtemp.all <- rbind(survtemp, survtemp, survtemp, survtemp)
  survtemp.all$AGEGROUP <- c('Age0-2', 'Age2-5', 'Age5+', 'Unknown')
  head(survtemp.all)
  
  # Assign cod_n
  survtemp.all$AGE_N[survtemp.all$AGEGROUP == "Age0-2"] <- 
    nrow(biotemp[biotemp$AGEGROUP == "Age0-2",])
  survtemp.all$AGE_N[survtemp.all$AGEGROUP == "Age2-5"] <- 
    nrow(biotemp[biotemp$AGEGROUP == "Age2-5",])
  survtemp.all$AGE_N[survtemp.all$AGEGROUP == "Age5+"] <- 
    nrow(biotemp[biotemp$AGEGROUP == "Age5+",])
  survtemp.all$AGE_N[survtemp.all$AGEGROUP == "Unknown"] <- 
    nrow(biotemp[biotemp$AGEGROUP == "Unknown",])
  
  # If survey data indicates more cod than bio data, add extra fish to unknown
  # BUT ONLY EXTRA FISH. Sometimes there are weirdly more reported fish in the
  # biodata than the survey data. Disregard that.
  if(sum(survtemp.all$AGE_N) < survtemp$COD_N){
    survtemp.all$AGE_N[survtemp.all$AGEGROUP == "Unknown"] <- 
      survtemp$COD_N - sum(survtemp.all$AGE_N)
  }
  
  # Assign cod_kg
  survtemp.all$AGE_KG <- NA
  if(length(biotemp$WEIGHT[biotemp$AGEGROUP == "Age0-2"]) > 0){
      survtemp.all$AGE_KG[survtemp.all$AGEGROUP == "Age0-2"] <- 
    sum(biotemp$WEIGHT[biotemp$AGEGROUP == "Age0-2"], na.rm=F)
  }
  if(length(biotemp$WEIGHT[biotemp$AGEGROUP == "Age2-5"]) > 0){
    survtemp.all$AGE_KG[survtemp.all$AGEGROUP == "Age2-5"] <- 
      sum(biotemp$WEIGHT[biotemp$AGEGROUP == "Age2-5"], na.rm=F)
  }
  if(length(biotemp$WEIGHT[biotemp$AGEGROUP == "Age5+"]) > 0){
    survtemp.all$AGE_KG[survtemp.all$AGEGROUP == "Age5+"] <- 
      sum(biotemp$WEIGHT[biotemp$AGEGROUP == "Age5+"], na.rm=F)
  }
  if(length(biotemp$WEIGHT[biotemp$AGEGROUP == "Unknown"]) > 0){
    survtemp.all$AGE_KG[survtemp.all$AGEGROUP == "Unknown"] <- 
      sum(biotemp$WEIGHT[biotemp$AGEGROUP == "Unknown"], na.rm=F)
  }
  head(survtemp.all)
  
  # when only one age group has cod and AGE_KG==NA and COD_KG !=NA, assign 
  # AGE_KG as COD_KG.
  group.nums <- as.vector(survtemp.all$AGE_N)
  length.zeros <- length(group.nums[group.nums == 0])
  if(length.zeros == 3 & is.na(survtemp.all$COD_KG[1])==FALSE){
    survtemp.all$AGE_KG[survtemp.all$AGE_N != 0] <- survtemp.all$COD_KG[1]
    survtemp.all$AGE_KG[survtemp.all$AGE_N == 0] <- 0
  }  
  head(survtemp.all)
  
  # If any NAs remain in AGE_KG where AGE_N == 0, reassign to 0
  survtemp.all$AGE_KG[survtemp.all$AGE_N == 0] <- 0
  
  # Assign data type
  if(nrow(survtemp.all[is.na(survtemp.all$AGE_KG),]) > 0){
    survtemp.all$Data_type <- 'Count'
  }
  if(nrow(survtemp.all[is.na(survtemp.all$AGE_KG),]) == 0){
    survtemp.all$Data_type <- 'Biomass_KG'
  }
  
  head(survtemp.all)
  
  bioda.list[[i]] <- survtemp.all
  
  rm(biotemp, survtemp, survtemp.all, group.nums, length.zeros)
    
}

test <- do.call(rbind, bioda.list)
row.names(test) <- NULL
head(test)
summary(test$AGE_N)
summary(test$AGE_KG)

table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

rm(bioda.list, i)

#### Data with no cod ####
# Pull out surveys that do not have any fish reflected in bio data
survey.encounter <- surv[surv$HAUL_ID %notin% test$HAUL_ID,]
table(survey.encounter$SURVEY)

# Pull out surveys with no recorded fish
enc.0 <- subset(survey.encounter, survey.encounter$COD_N == 0)

# Clone encounter data for each of the four age groups
enc.0.age0 <- enc.0; enc.0.age0$AGEGROUP <- 'Age0-2'
enc.0.age1 <- enc.0; enc.0.age1$AGEGROUP <- 'Age2-5'
enc.0.age2 <- enc.0; enc.0.age2$AGEGROUP <- 'Age5+'
enc.0.age3 <- enc.0; enc.0.age3$AGEGROUP <- 'Unknown'

# Rebind
enc.0 <- rbind(enc.0.age0, enc.0.age1, enc.0.age2, enc.0.age3)

# Create columns in test but not enc
enc.0$AGE_N <- 0; enc.0$AGE_KG <- 0; enc.0$Data_type <- "Biomass_KG"

# Bind cod==0 data into final dataframe
test <- rbind(test, enc.0)
table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

# Remove intermediates
rm(enc.0, enc.0.age0, enc.0.age1, enc.0.age2, enc.0.age3)

#### Assign age by weight (survey data available, bio data not, cod==1)####
unaccounted <- survey.encounter[survey.encounter$HAUL_ID %notin% test$HAUL_ID,]
enc.n1 <- subset(unaccounted, unaccounted$COD_N == 1 & !is.na(unaccounted$COD_KG))

# Create new age grouping column
enc.n1$AGEGROUP <- NA

# Loop through rows
for(i in 1:nrow(enc.n1)){
  
  # Age 0-2 cutoff
  if(enc.n1$COD_KG[i] < wts2){
    enc.n1$AGEGROUP <- 'Age0-2'
  }
  # Ages 2-5
  if(enc.n1$COD_KG[i] >= wts2 & enc.n1$COD_KG[i] <wts5){
    enc.n1$AGEGROUP[i] <- 'Age2-5'
  }
  # Age 5 cutoff
  if(enc.n1$COD_KG[i] >= wts5){
    enc.n1$AGEGROUP[i] <- 'Age5+'
  }
}

table(enc.n1$AGEGROUP)
enc.n1$AGE_N <- enc.n1$COD_N
enc.n1$AGE_KG <- enc.n1$COD_KG
enc.n1$Data_type <- 'Biomass_KG'

# Split into list
enc.n1.list <- split(enc.n1, f=enc.n1$HAUL_ID)
# Name age groups
all.groups <- c('Age0-2', 'Age2-5', 'Age5+', 'Unknown')

# Create rows for missing age groups
for(i in 1:length(enc.n1.list)){
  temp <- enc.n1.list[[i]]
  present.groups <- temp$AGEGROUP

  use.groups <- all.groups[all.groups %notin% present.groups]
  
  blank <- rbind(temp, temp, temp)
  blank$AGE_N <- 0; blank$AGE_KG <- 0
  blank$AGEGROUP <- use.groups
  
  final <- rbind(temp, blank)
  enc.n1.list[[i]] <- final
  
  rm(present.groups, use.groups, blank, final, temp)
}
# Rebind
enc.n1 <- do.call(rbind, enc.n1.list)
enc.n1 <- enc.n1[with(enc.n1, order(HAUL_ID, AGEGROUP)),]
row.names(enc.n1) <- NULL

# Bind into dataframe
test <- rbind(test, enc.n1)
summary(test$AGE_KG)
summary(test$AGE_N)

table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

# Remove intermediates
rm(enc.n1, enc.n1.list, i, all.groups)

#### Leftovers ####
# Surveys that are still not dealt with
unaccounted2 <- survey.encounter[survey.encounter$HAUL_ID %notin% test$HAUL_ID,]

# Surveys that have recorded biomass
enc.unk <- unaccounted2[!is.na(unaccounted2$COD_KG) & unaccounted2$COD_N !=0 & 
                          unaccounted2$COD_N !=1,]

# Biomass of all caught fish is less than fitted weight at age==2
enc.small <- subset(enc.unk, enc.unk$COD_KG <= wts2)

# If the total weight is less than 1 maximum sized Age0-2 fish, all fish
# in the  haul must be Age0-2.
enc.small$AGEGROUP <- 'Age0-2'
enc.small$AGE_N <- enc.small$COD_N
enc.small$AGE_KG <- enc.small$COD_KG
enc.small$Data_type <- 'Biomass_KG'

# Split into list
enc.small.list <- split(enc.small, f=enc.small$HAUL_ID)

# Create rows for missing age groups
for(i in 1:length(enc.small.list)){
  temp <- enc.small.list[[i]]
  
  use.groups <- c('Age2-5', 'Age5+', 'Unknown')
  
  blank <- rbind(temp, temp, temp)
  blank$AGE_N <- 0; blank$AGE_KG <- 0
  blank$AGEGROUP <- use.groups
  
  final <- rbind(temp, blank)
  enc.small.list[[i]] <- final
  
  rm(use.groups, blank, final, temp)
}
# Rebind
enc.small <- do.call(rbind, enc.small.list)
enc.small <- enc.small[with(enc.small, order(HAUL_ID, AGEGROUP)),]
row.names(enc.small) <- NULL

# Bind into final dataframe
test <- rbind(test, enc.small)
summary(test$AGE_N)
summary(test$AGE_KG)
table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

# Remove intermediates
rm(enc.small, enc.small.list, enc.unk, i)

#### Absolute dregs ###
unaccounted3 <- survey.encounter[survey.encounter$HAUL_ID %notin% test$HAUL_ID,]

sum(unaccounted3$COD_N)
# There are 2195 tows, representing 214,141.7 fish, that cannot be assigned 
# age by any method. These will be thrown into the unknown category.

# Some can be grouped in the biomass category. Some are missing biomass, and
# so must be count.

# Start with biomoass
unaccounted.bio <- unaccounted3[!is.na(unaccounted3$COD_KG),]
# This is 994 tows and 204763.4 fish
unaccounted.bio$AGEGROUP <- 'Unknown'
unaccounted.bio$AGE_N <- unaccounted.bio$COD_N
unaccounted.bio$AGE_KG <- unaccounted.bio$COD_KG
unaccounted.bio$Data_type <- 'Biomass_KG'

# Split into list
unaccounted.bio.list <- split(unaccounted.bio, f=unaccounted.bio$HAUL_ID)

# Create rows for missing age groups
for(i in 1:length(unaccounted.bio.list)){
  temp <- unaccounted.bio.list[[i]]
  
  use.groups <- c('Age0-2', 'Age2-5', 'Age5+')
  
  blank <- rbind(temp, temp, temp)
  blank$AGE_N <- 0; blank$AGE_KG <- 0
  blank$AGEGROUP <- use.groups
  
  final <- rbind(blank, temp)
  unaccounted.bio.list[[i]] <- final
  
  rm(use.groups, blank, final, temp)
}
# Rebind
unaccounted.bio <- do.call(rbind, unaccounted.bio.list)
unaccounted.bio <- unaccounted.bio[with(unaccounted.bio, order(HAUL_ID, AGEGROUP)),]
row.names(unaccounted.bio) <- NULL
# Bind into final dataframe
test <- rbind(test, unaccounted.bio)
summary(test$AGE_N)
summary(test$AGE_KG)
table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

# Remove intermediates
rm(unaccounted, unaccounted.bio, unaccounted.bio.list, unaccounted2, i)

# End with Unknown age class, count data
unaccounted.count <- unaccounted3[is.na(unaccounted3$COD_KG),]
# This is 994 tows and 204763.4 fish
unaccounted.count$AGEGROUP <- 'Unknown'
unaccounted.count$AGE_N <- unaccounted.count$COD_N
unaccounted.count$AGE_KG <- unaccounted.count$COD_KG
unaccounted.count$Data_type <- 'Count'

# Split into list
unaccounted.count.list <- split(unaccounted.count, f=unaccounted.count$HAUL_ID)

# Create rows for missing age groups
for(i in 1:length(unaccounted.count.list)){
  temp <- unaccounted.count.list[[i]]
  
  use.groups <- c('Age0-2', 'Age2-5', 'Age5+')
  
  blank <- rbind(temp, temp, temp)
  blank$AGE_N <- 0; blank$AGE_KG <- NA
  blank$AGEGROUP <- use.groups
  
  final <- rbind(blank, temp)
  unaccounted.count.list[[i]] <- final
  
  rm(use.groups, blank, final, temp)
}
# Rebind
unaccounted.count <- do.call(rbind, unaccounted.count.list)
unaccounted.count <- unaccounted.count[with(unaccounted.count, 
                                            order(HAUL_ID, AGEGROUP)),]
row.names(unaccounted.count) <- NULL
# Bind into final dataframe
test <- rbind(test, unaccounted.count)
summary(test$AGE_N)
summary(test$AGE_KG)
table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)

# Remove intermediates
rm(survey.encounter, unaccounted.count, unaccounted.count.list,
   unaccounted3, i)

#### Final testing ####
allhauls <- unique(c(bioda$HAUL_ID, surv$HAUL_ID))
nrow(test[test$HAUL_ID %notin% allhauls,])
sum(test$AGE_N) >= sum(surv$COD_N)
# Remember that biodata had some fish that were not accounted for in surv data

# Remove intermediates
rm(bioda, surv, allhauls, cuts2, cuts5, wts2, wts5)

# Order by survey and date
test <- test[with(test, order(SURVEY, DATE, HAUL_ID, AGEGROUP)),]
row.names(test) <- NULL
head(test)
table(test$Data_type, test$AGEGROUP)
table(test$SURVEY   , test$AGEGROUP)
table(test$SURVEY   , test$Data_type)
summary(test)

# Create Year-season column
test$SEASON <- factor(test$SEASON, levels = c('SPRING', 'FALL'),
                      labels = c('A.SPRING', 'B.FALL'))
test$TIME <- paste0(test$YEAR, test$SEASON)
test$TIME <- as.numeric(as.factor(test$TIME))
table(test$TIME, test$Data_type)

# Set response
test$RESPONSE <- NA
for(i in 1:nrow(test)){
  if(test$Data_type[i] == "Count"){
    test$RESPONSE[i] <- test$AGE_N[i]
  }
  if(test$Data_type[i] == "Biomass_KG"){
    test$RESPONSE[i] <- test$AGE_KG[i]
  }
}
table(test$Data_type)
summary(test$RESPONSE)


# Save.
write.csv(test,
          here("data/Dataframes/Bio_Data_Agesep3.csv"),
          row.names=F)
