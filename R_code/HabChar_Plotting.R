rm(list=ls())

# Add libraries
library(tidyverse)
library(here)
library(sf)
library(marmap)
library(raster)

# Set seed
set.seed(123)

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.title = element_text(size=12),
                legend.text = element_text(size=10),
                legend.background = element_blank(),
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=16, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Create extent polygon
poly <- data.frame(
  lon = c(-71, -69),
  lat = c(42, 43.5)
)
poly <- poly %>% 
  st_as_sf(coords = c("lon", "lat"), 
           crs = 'EPSG:4326') %>% 
  st_bbox() %>% 
  st_as_sfc()

# Pull NOAA bathymetry data
Bathy <- getNOAA.bathy(lon1 = -72, lon2 = -66,
                       lat1 = 41, lat2 = 44, 
                       resolution = 1)

# Convert data to matrix
Bathy_Final <- as.matrix(Bathy)
class(Bathy_Final) <- "matrix"

# Reshape for plotting
BathyData <- Bathy_Final %>%
  as.data.frame() %>% 
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric))

# Clip to poly
# BathyData <- BathyData[BathyData$lon <= -69 & BathyData$lon >= -71 &
#                          BathyData$lat >= 42 & BathyData$lat <= 43.5,]

# Remove values on land
# BathyData <- subset(BathyData, value <= 0)
summary(BathyData$value)
#BathyData$value <- BathyData$value * -1

# Set depth breaks and colors
dbreaks <- c(seq(0, 100, 5),
             seq(100, 300, 25))
dbreaks <- dbreaks * -1
dcol <- colorRampPalette(c("#DEF5E5", "#40498E"))
dcol <- dcol(length(dbreaks))
rm(Bathy_Final, Bathy)

# Load sediments
grid <- st_read(here('Data/GIS/Sediment_Krig_1K_Polygons.shp'))

# Reproject
grid <- st_transform(grid, crs="EPSG:4326")

# Clip to polygon
grid <- st_intersection(grid, poly)

# Reshape grid
gridold <- grid
grid <- dplyr::select(grid, cobble_P, gravel_P, mud_P, sand_P, rock_P, 
                      geometry)
colnames(grid) <- c('cobble', 'gravel', 'mud', 'sand', 'rock', 'geometry')

cobble <- grid$cobble
rock <- grid$rock
gravel <- grid$gravel
sand <- grid$sand
mud <- grid$mud

vals <- c(cobble, rock, gravel, sand, mud)
head(vals)

sed.ty <- c(rep("cobble", 23901),
            rep("rock", 23901),
            rep("gravel", 23901),
            rep("sand", 23901),
            rep("mud", 23901))

grid <- rbind(grid, grid, grid, grid, grid)
grid$val <- vals
grid$sed.ty <- sed.ty
grid <- dplyr::select(grid, val, sed.ty, geometry)
row.names(grid) <- NULL

# Load rugosity
rugos <- st_read(paste0("C:/Users/klankowicz/Box/Katie Lankowicz/",
                        "Data_Analysis/Cod/GIS/Rugosity_15as_F.shp"))

# Reproject
rugos <- st_transform(rugos, st_crs(grid))
head(rugos)
rugos <- dplyr::select(rugos, rugosity, COND, geometry)

# Clip to polygon
rugos <- st_intersection(rugos, poly)

# Convert to stars
library(stars)
rugos.stars <- st_as_stars(rugos, name = attr(rugos, "sf_column")) 
rugos.stars$rugosity <- NULL
rugos.stars

# Load surveys
survs <- read.csv(here('data/Dataframes/sciencecenter_wrug.csv'))
survs <- dplyr::select(survs, -FID_Rugosi, -STRATA_1, -pct70, -rugosity,
                       -OID, -FID_scienc)
head(survs)

# Strip to just one age category
survs <- subset(survs, AGEGROUP == 'Age0-2')

# Split to surveys
survs <- subset(survs, SURVEY == 'NEFSC BTS' |
                  SURVEY =='NEFSC BLLS')

# Stip to just years BLLS took place
table(survs$YEAR[survs$SURVEY == "NEFSC BLLS"])
survs <- subset(survs, YEAR %in% seq(2014, 2022))

# Strip to just WGOM
survs <- subset(survs, STOCK == 'WGOM')

# Strip to 9-box footprint
survs <- subset(survs, LON > -71 & LON < -69 &
                  LAT > 42 & LAT < 43.5)

# Check results
table(survs$SURVEY)
table(survs$YEAR)
table(survs$SEASON)

# Create depth strata
survs$BATHY.DEPTH <- survs$BATHY_DEPT * -1
survs$dep.strat[survs$BATHY.DEPTH < 110] <- 'shallow'
survs$dep.strat[survs$BATHY.DEPTH >= 110] <- 'deep'

write.csv(survs, row.names = F, 'possible_tsplocs.csv')

# Convert to sf
splitsurvs <- st_as_sf(survs, coords=c('LON', "LAT"),
                       crs="EPSG:4326")

# View overlaps bt strata
bts <- subset(splitsurvs, SURVEY == "NEFSC BTS")
blls <- subset(splitsurvs, SURVEY == "NEFSC BLLS")
table(bts$dep.strat, bts$COND)
table(blls$dep.strat, blls$COND)

# Load coast
coast <- ecodata::coast
coast <- st_transform(coast, st_crs(splitsurvs))

# Plot basic
ggplot() +
  geom_sf(data=coast) +
  geom_sf(data=splitsurvs, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude')

# Plot against depth
dep <- ggplot() +
  geom_contour_filled(data = BathyData, aes(x = lon, y = lat, z = value),
                      breaks = dbreaks) +  geom_sf(data=coast) +
  scale_fill_manual(values =  dcol) +
  geom_sf(data=splitsurvs, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +  
  ylab('Latitude') +
  xlab('Longitude')
dep <- dep + guides(fill='none')
dep  

# Plot against sediment
ggplot() +
  geom_sf(data=grid, col=NA, 
          aes(fill=val)) +
  scale_fill_viridis_c(option='viridis',
                       begin=0.3, end=1,
                       na.value=NA) +
  geom_sf(data=coast) +
  geom_sf(data=splitsurvs, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude') +
  facet_wrap(vars(sed.ty))

# Plot against rugosity
ggplot() +
  geom_sf(data=rugos, col=NA, 
          aes(fill=COND)) +
  geom_sf(data=coast) +
  geom_sf(data=splitsurvs, cex=1,
          col='black') +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude') +
  facet_wrap(vars(SURVEY))

#### Stratified random selection ####
# GOal: 100 survey points from each of the two surveys (blls, bts)
# Strata: Rough and Smooth (bottom type), Deep (> 110m) and shallow (<110m)
# Therefore four groups (eg, rough/shallow, rough/deep, etc.)
# Need 25 random points from each

# BLLS
blls$rand <- 0
table(blls$COND, blls$dep.strat)
blls$rand[blls$dep.strat == 'deep' & blls$COND == 'ROUGH' &
            blls$HAUL_ID %in%
            sample(blls$HAUL_ID[blls$dep.strat == 'deep' & 
                                  blls$COND == 'ROUGH'],
                   25, replace=F)] <- 1

blls$rand[blls$dep.strat == 'deep' & blls$COND == 'SMOOTH' &
            blls$HAUL_ID %in%
            sample(blls$HAUL_ID[blls$dep.strat == 'deep' & 
                                  blls$COND == 'SMOOTH'],
                   25, replace=F)] <- 1

blls$rand[blls$dep.strat == 'shallow' & blls$COND == 'ROUGH' &
            blls$HAUL_ID %in%
            sample(blls$HAUL_ID[blls$dep.strat == 'shallow' & 
                                  blls$COND == 'ROUGH'],
                   25, replace=F)] <- 1

blls$rand[blls$dep.strat == 'shallow' & blls$COND == 'SMOOTH' &
            blls$HAUL_ID %in%
            sample(blls$HAUL_ID[blls$dep.strat == 'shallow' & 
                                  blls$COND == 'SMOOTH'],
                   25, replace=F)] <- 1

blls.rand <- subset(blls, rand == 1)
table(blls.rand$rand)
table(blls.rand$COND, blls.rand$dep.strat)

# bts
bts$rand <- 0
table(bts$COND, bts$dep.strat)
# Problem: only 23 samples in rough/ shallow. Will use 22.
# Will need to add 1 extra sample to each other category to fill to 100
bts$rand[bts$dep.strat == 'deep' & bts$COND == 'ROUGH' &
            bts$HAUL_ID %in%
            sample(bts$HAUL_ID[bts$dep.strat == 'deep' & 
                                  bts$COND == 'ROUGH'],
                   26, replace=F)] <- 1

bts$rand[bts$dep.strat == 'deep' & bts$COND == 'SMOOTH' &
            bts$HAUL_ID %in%
            sample(bts$HAUL_ID[bts$dep.strat == 'deep' & 
                                  bts$COND == 'SMOOTH'],
                   26, replace=F)] <- 1

bts$rand[bts$dep.strat == 'shallow' & bts$COND == 'ROUGH' &
            bts$HAUL_ID %in%
            sample(bts$HAUL_ID[bts$dep.strat == 'shallow' & 
                                  bts$COND == 'ROUGH'],
                   22, replace=F)] <- 1

bts$rand[bts$dep.strat == 'shallow' & bts$COND == 'SMOOTH' &
            bts$HAUL_ID %in%
            sample(bts$HAUL_ID[bts$dep.strat == 'shallow' & 
                                  bts$COND == 'SMOOTH'],
                   26, replace=F)] <- 1

bts.rand <- subset(bts, rand == 1)
table(bts.rand$rand)
table(bts.rand$COND, bts.rand$dep.strat)

# Separate randomly selected points
selected <- rbind(bts.rand, blls.rand)
selected <- dplyr::select(selected, -rand)
summary(selected)

# Re-plot
# Plot basic
ggplot() +
  geom_sf(data=coast) +
  geom_sf(data=selected, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude')

# Plot against depth
BathyData <- subset(BathyData, value <= 0)
dep <- ggplot() +
  geom_contour_filled(data = BathyData, aes(x = lon, y = lat, z = value),
                      breaks = dbreaks) +  geom_sf(data=coast) +
  scale_fill_manual(values =  dcol) +
  geom_sf(data=selected, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +  
  ylab('Latitude') +
  xlab('Longitude')
dep <- dep + guides(fill='none')
dep  

# Plot against sediment
ggplot() +
  geom_sf(data=grid, col=NA, 
          aes(fill=val)) +
  scale_fill_viridis_c(option='viridis',
                       begin=0.3, end=1,
                       na.value=NA) +
  geom_sf(data=coast) +
  geom_sf(data=selected, cex=1,
          aes(col=SURVEY)) +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude') +
  facet_wrap(vars(sed.ty))

# Plot against rugosity
ggplot() +
  geom_sf(data=rugos, col=NA, 
          aes(fill=COND)) +
  geom_sf(data=coast) +
  geom_sf(data=selected, cex=1,
          col='black') +
  geom_sf(data=poly, fill=NA, col='black') +
  coord_sf(xlim=c(-71, -68.95),
           ylim=c(41.95, 43.5)) +
  ylab('Latitude') +
  xlab('Longitude') +
  facet_wrap(vars(SURVEY))

# Save results
write.csv(selected, row.names = F, 
          'selected_sites.csv')
