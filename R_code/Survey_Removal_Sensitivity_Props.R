rm(list=ls())

library(tidyverse)
library(here)
library(sf)

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

indall <- read.csv(here('VAST_runs/small/Overall_BC/ALL/Index.csv'))
indno_mi <- read.csv(here('VAST_runs/sensitivity/small/MADMF Trawl/Small_Index.csv'))

head(indall)
head(indno_mi)

indall$Mod <- 'All'
indno_mi$Mod <- 'MADMF Trawl'

indno_mi <- indno_mi %>% 
  dplyr::select(-Category, -Std.Err.ln.Est)

indall <- indall %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  rename(Std.Err.Est = Std..Error.for.Estimate) %>% 
  dplyr::select(-Category, -Units, -Std..Error.for.ln.Estimate.)


ind <- rbind(indall, indno_mi)

ind <- ind %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  mutate(Mod = factor(Mod, levels=c('All', 'MADMF Trawl')),
         Size = as.factor('Small'),
         Stratum = factor(Stratum, levels=c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')))
summary(ind)

#### Find areas with data coverage by this survey ####
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Load data
surveys <- read.csv(here("Data/Survey_Data/Bio_Data_Agesep.csv"))
surveys$month <- month(as.POSIXct(surveys$DATE, format='%Y-%m-%d'))

# Only smalls
surveys <- surveys[surveys$AGEGROUP == 'Age0-2',]

# Add seasonality
surveys$season[surveys$month %in% c(1,2,9,10,11,12)] <- 'BFall'
surveys$season[surveys$month %in% c(3,4,5,6,7,8)] <- 'ASpring'

surveys$YEAR[surveys$month %in% c(1,2)] <- 
  surveys$YEAR[surveys$month %in% c(1,2)]-1

surveys$SEASON <- paste0(surveys$YEAR, ' ', surveys$season)

# Fix stupid error
surveys$TIME <- as.numeric(as.factor(paste0(surveys$YEAR, ' ',
                                            surveys$SEASON)))

# Remove points not in model domain (on land)
region_shape<- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)
codstox <- st_transform(st_read('Data/GIS/codstox.shp', quiet=T),
                        crs="EPSG:4326")
codstox <- st_make_valid(codstox)

surveys_sf <- st_as_sf(surveys, coords=c('LON', 'LAT'))
st_crs(surveys_sf) <- 'EPSG:4326'

surveys_sf <- st_intersection(surveys_sf, region_shape)

surveys <- surveys[surveys$HAUL_ID %in% surveys_sf$HAUL_ID,]

# Test for missing values
surveys <- surveys[!is.na(surveys$rugos) &
                     !is.na(surveys$nao) &
                     !is.na(surveys$h_bt),]

# Remove survey
surveys <- surveys[surveys$SURVEY == 'MADMF Inshore Trawl',]

surveys$season[surveys$season == 'BFall'] <- 'Fall'
surveys$season[surveys$season == 'ASpring'] <- 'Spring'

surveys <- surveys %>% 
  dplyr::select(SURVEY, YEAR, DATE, LON, LAT, AGE_N, season) %>% 
  mutate(season = factor(season, levels=c('Spring', 'Fall')))

surveys <- st_as_sf(surveys, coords=c('LON', 'LAT'), crs='EPSG:4326')
surveys <- st_intersection(surveys, codstox)
table(surveys$STOCK, surveys$SURVEY)

surveys <- surveys %>% 
  mutate(NCat = NA)

surveys$NCat[surveys$AGE_N==0] <- 0
surveys$NCat[surveys$AGE_N %in% seq(1, 10, 1)] <- 10
surveys$NCat[surveys$AGE_N %in% seq(11, 50, 1)] <- 50
surveys$NCat[surveys$AGE_N %in% seq(51, 100,1)] <- 100
surveys$NCat[surveys$AGE_N %in% seq(101, 200,1)] <- 200
surveys$NCat[surveys$AGE_N %in% seq(201, 300,1)] <- 300
surveys$NCat[surveys$AGE_N %in% seq(301, 400,1)] <- 400
surveys$NCat[surveys$AGE_N %in% seq(401, 500,1)] <- 500
surveys$NCat[surveys$AGE_N %in% seq(501, 750,1)] <- 750
surveys$NCat[surveys$AGE_N %in% seq(751, 1000,1)] <- 1000
surveys$NCat[surveys$AGE_N %in% seq(1001, 2000,1)] <- 1500

ggplot(data=ind#[ind$Stratum == 'WGOM' & ind$Season == 'Spring',]
       ) +
  geom_point(data=ind[#ind$Stratum == 'WGOM' & ind$Season == 'Spring' &
                        ind$Year %in% surveys$YEAR,],
             aes(x=Year, y=Estimate, col=Mod))+
  geom_line(aes(x=Year, y=Estimate, col=Mod)) +
  geom_ribbon(aes(x=Year, ymin=Estimate-Std.Err.Est,
                  ymax=Estimate + Std.Err.Est, fill=Mod),
              alpha=0.4) +
  ggh4x::facet_grid2(Stratum~Season, scales='free')

ggplot() +
  geom_sf(data=surveys[surveys$YEAR %in% c(2014:2021) &
                         surveys$season =='Spring' &
                         surveys$AGE_N != 0,],
          aes(size = NCat)) +
  geom_sf(data=codstox, aes(fill=STOCK), alpha=0.2) +
  coord_sf(xlim=c(st_bbox(surveys)[1], st_bbox(surveys)[3]), 
           ylim=c(st_bbox(surveys)[2], st_bbox(surveys)[4])) +
  facet_wrap(vars(YEAR), nrow=2)

#### Calc years with sig dif ####
dif <- indno_mi$Estimate / indall$Estimate
comp <- ind[ind$Mod == 'MADMF Trawl',]
comp$dif <- dif

ggplot(data=comp) +
  geom_line(aes(x=Year, y=dif, col=Stratum)) +
  geom_point(aes(x=Year, y=dif, col=Stratum), cex=0.7) +
  facet_wrap(vars(Season))

testcap <- data.frame(
  Year = comp$Year,
  Season = comp$Season,
  Stratum = comp$Stratum,
  ttest = NA
)

for(i in 1:nrow(testcap)){
  testcap$ttest[i] <-      t.test(c(indall$Estimate[i],
                                    indall$Estimate[i] - indall$Std.Err.Est[i],
                                    indall$Estimate[i] + indall$Std.Err.Est[i]), 
                                  c(indno_mi$Estimate[i],
                                    indno_mi$Estimate[i] - indno_mi$Std.Err.Est[i],
                                    indno_mi$Estimate[i] + indno_mi$Std.Err.Est[i]))$p.value
}

ggplot(data=testcap) +
  geom_line(aes(x=Year, y=ttest, col=Stratum)) +
  geom_hline(yintercept = 0.05, col='red', lty=2) +
  facet_wrap(vars(Season))

testcap[testcap$ttest < 0.05,]

long <- indall %>% 
  rename(All.Est = Estimate, All.Std = Std.Err.Est) %>% 
  mutate(Mod.Est = indno_mi$Estimate,
         Mod.Std = indno_mi$Std.Err.Est)

ggplot(data=long) +
  geom_point(aes(x=All.Est/1000000, y=Mod.Est/1000000)) +
  ggh4x::facet_grid2(Season ~ Stratum) +
  geom_abline(slope=1, intercept = 0, col='red', lty=2)

allsplit <- split(indall, f=indall$Season)
for(i in 1:length(allsplit)){
  allsplit[[i]] <- split(allsplit[[i]],
                         f=allsplit[[i]]$Stratum)
}

for(i in 1:length(allsplit)){
  for(j in 1:length(allsplit[[i]])){
    allsplit[[i]][[j]] <- auto.arima(allsplit[[i]][[j]]$Estimate)
  }
}

nosplit <- split(indno_mi, f=indno_mi$Season)
for(i in 1:length(nosplit)){
  nosplit[[i]] <- split(nosplit[[i]],
                         f=nosplit[[i]]$Stratum)
}

for(i in 1:length(nosplit)){
  for(j in 1:length(nosplit[[i]])){
    nosplit[[i]][[j]] <- auto.arima(nosplit[[i]][[j]]$Estimate)
  }
}

ccf(allsplit[[i]][[j]], nosplit[[i]][[j]])
