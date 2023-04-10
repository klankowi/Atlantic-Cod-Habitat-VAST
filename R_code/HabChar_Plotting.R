rm(list=ls())

# Add libraries
library(tidyverse)
library(here)
library(sf)
library(marmap)
library(raster)

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

# Set depth breaks and colors
dbreaks <- c(seq(0, 20, 2),
             seq(20, 40, 4),
             seq(40, 70, 5),
             seq(70, 100, 10),
             seq(100, 200, 20),
             seq(200, 400, 50),
             seq(400, 1200, 100))
dbreaks <- unique(dbreaks) * -1
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
rugos <- st_transform(rugos, st_crs(grid))
head(rugos)
rugos <- dplyr::select(rugos, rugosity, COND, geometry)

# Clip to polygon
rugos <- st_intersection(rugos, poly)

# Load surveys
survs <- read.csv(here('data/Dataframes/cod_agesep_VASTdata.csv'))

# Strip to just one age category
survs <- subset(survs, AGEGROUP == 'Age0-2')

# Stip to just years BLLS took place
survs <- subset(survs, YEAR %in% seq(2016, 2022))

# Strip to just WGOM
survs <- subset(survs, STOCK == 'WGOM')

# Strip to 9-box footprint
survs <- subset(survs, LON > -71 & LON < -69 &
                  LAT > 42 & LAT < 43.5)

# Check results
table(survs$SURVEY)
table(survs$YEAR)
table(survs$SEASON)

# Split to surveys
survs <- subset(survs, SURVEY == 'NEFSC BTS' |
                  SURVEY =='NEFSC BLLS')

# Convert to sf
splitsurvs <- st_as_sf(survs, coords=c('LON', "LAT"),
                       crs="EPSG:4326")

# Load coast
coast <- ecodata::coast
coast <- st_transform(coast, st_crs(splitbts))

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
