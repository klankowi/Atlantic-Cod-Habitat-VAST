# Average envionmental conditions
# Suitability indices
# All at scale of density grid

rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)
library(raster)
library(stars)
library(sf)
library(beepr)

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

#### Create blank grid ####
# Load grid
load(here('VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.RData'))
# Clear workspace
rm(list=setdiff(ls(), c('fit')))
# Extract Data
Y_gt = fit$Report$D_gct[,1,]
map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
                         spatial_list = fit$spatial_list,
                         Extrapolation_List = fit$extrapolation_list)
panel_labels = fit$year_labels

# Call data
MapSizeRatio = map_list$MapSizeRatio
Y_gt = Y_gt[map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                    0), "x2i"], , drop = FALSE]
n_cells = nrow(Y_gt)

# Call location list
loc_g = map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                0), c("Lon", "Lat")]
# Call projections
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS("EPSG:26919")

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
regions <- st_read(here("Data/GIS/codstox.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                         data = data.frame(y=rep(1, 2000)), 
                                         proj4string = CRS_orig)
Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
Points_proj = sp::spTransform(Points_orig, CRS_proj)
xlim = Points_proj@bbox[1, ]
ylim = Points_proj@bbox[2, ]

cell.size = mean(diff(Points_proj@bbox[1, ]), diff(Points_proj@bbox[2, 
]))/floor(sqrt(n_cells))
Points_sf = sf::st_as_sf(Points_proj)
grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
grid_i = sf::st_intersects(Points_sf, grid)
grid = sf::st_sf(grid, 
                 y = tapply(Points_sf$y, 
                            INDEX = factor(as.numeric(grid_i),
                                           levels = 1:length(grid)), 
                            FUN = mean, na.rm = TRUE))
grid$ID <- seq(1:nrow(grid))

# PLot blank grid
ggplot() + 
  geom_sf(data=grid, col='black', fill='pink') +
  geom_sf(data=coast) +
  coord_sf(xlim=c(-76, -66),
           ylim=c(36, 45),
           crs="EPSG:4326")

# Check grid
grid <- st_transform(grid, crs="EPSG:26919")
grid <- st_make_valid(grid)
grid <- grid %>% 
  st_cast('MULTIPOLYGON') %>% 
  st_cast('POLYGON')

# Clear workspace
rm(list=setdiff(ls(), c('regions', 'coast', 'grid')))

#### Depth within grid cells ####
# Load raster
gebco <- raster(here('Data/Bathymetry/GEBCO_2023.tif'))

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(gebco))

depth <- gebco %>% 
  mask(region) %>%
  crop(region)

# Convert so matches curve -- depths are positive values, land is negative
depth@data@values <- depth@data@values * -1

# Remove land values
depth@data@values[depth@data@values <=0.1] <- NA

# Find range of depth in region
summary(depth@data@values)
plot(depth, asp=1)

# Convert to sf
depthsf <- st_as_sf(st_as_stars(depth))
depthsf <- depthsf %>% 
  rename(bathy = gebco_2023) %>% 
  filter(!is.na(bathy))

# Clean workspace
rm(gebco, region, depth)

# Order by depth
depthsf <- depthsf[with(depthsf, order(bathy)),]
rownames(depthsf) <- NULL

# Set to same projection as grid
depthsf <- st_transform(depthsf, st_crs(grid))
depthsf <- st_make_valid(depthsf)

# Intersect
depthgrid <- st_intersection(depthsf, grid)

# Convert to dataframe
depthgrid <- depthgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, bathy) %>% 
  group_by(ID) %>% 
  # Mean bathymetry in each grid cell
  summarise(bathy =     mean(bathy,     na.rm=T)) %>% 
  # Round to integer
  mutate(bathy = round(bathy))

# Re-merge to grid by ID
grid <- left_join(grid, depthgrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=bathy), color=NA) +
  scale_fill_viridis_c(na.value='transparent',
                       direction = -1)+
  geom_sf(data=regions, fill=NA, col='black')

#### Rugosity within grid cells ####
# Load raster
load(here('Data/Density_Covariates/Rugosity/rast_rugosity.RData'))
rugos <- masked.raster
rm(masked.raster)

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(rugos))

rugos <- rugos %>% 
  mask(region) %>%
  crop(region)

# Find range of depth in region
summary(rugos@data@values)
plot(rugos, asp=1)

# Convert to sf
rugossf <- st_as_sf(st_as_stars(rugos))
rugossf <- rugossf %>% 
  rename(rugosity = band1) %>% 
  filter(!is.na(rugosity))

# Clean workspace
rm(rugos, region)

# Order by rugosity
rugossf <- rugossf[with(rugossf, order(rugosity)),]
rownames(rugossf) <- NULL

# Set to same projection as grid
rugossf <- st_transform(rugossf, st_crs(grid))
rugossf <- st_make_valid(rugossf)

# Intersect
rugosgrid <- st_intersection(rugossf, grid)

# Convert to dataframe
rugosgrid <- rugosgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, rugosity) %>% 
  group_by(ID) %>% 
  # Mean rugosity in each grid cell
  summarise(rugosity =     mean(rugosity,     na.rm=T))

# Re-merge to grid by ID
grid <- left_join(grid, rugosgrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=rugosity), color=NA) +
  scale_fill_viridis_c(na.value='transparent')+
  geom_sf(data=regions, fill=NA, col='black')

#### Sediments within grid cells ####
# Load raster
seds <- read_rds(here('Data/Density_Covariates/Sediment/sediment_interp_grid_hires.RDS'))

# Clip to region
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_make_valid(region)
region <- st_transform(region, st_crs(seds))

seds <- st_intersection(seds, region)

# Cobble
cobblesf <- seds %>% 
  dplyr::select(cobble, geometry)
st_geometry(cobblesf) <- 'geometry'

# Order by cobble
cobblesf <- cobblesf[with(cobblesf, order(cobble)),]
rownames(cobblesf) <- NULL

# Set to same projection as grid
cobblesf <- st_transform(cobblesf, st_crs(grid))
cobblesf <- st_make_valid(cobblesf)

# Intersect
cobblegrid <- st_intersection(cobblesf, grid)

# Convert to dataframe
cobblegrid <- cobblegrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, cobble) %>% 
  group_by(ID) %>% 
  # Mean rugosity in each grid cell
  summarise(cobble =     mean(cobble,     na.rm=T))

# Re-merge to grid by ID
grid <- left_join(grid, cobblegrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=cobble), color=NA) +
  scale_fill_viridis_c(na.value='transparent')+
  geom_sf(data=regions, fill=NA, col='black')

# Gravel
gravelsf <- seds %>% 
  dplyr::select(gravel, geometry)
st_geometry(gravelsf) <- 'geometry'

# Order by gravel
gravelsf <- gravelsf[with(gravelsf, order(gravel)),]
rownames(gravelsf) <- NULL

# Set to same projection as grid
gravelsf <- st_transform(gravelsf, st_crs(grid))
gravelsf <- st_make_valid(gravelsf)

# Intersect
gravelgrid <- st_intersection(gravelsf, grid)

# Convert to dataframe
gravelgrid <- gravelgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, gravel) %>% 
  group_by(ID) %>% 
  # Mean rugosity in each grid cell
  summarise(gravel =     mean(gravel,     na.rm=T))

# Re-merge to grid by ID
grid <- left_join(grid, gravelgrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=gravel), color=NA) +
  scale_fill_viridis_c(na.value='transparent')+
  geom_sf(data=regions, fill=NA, col='black')

# Mud
mudsf <- seds %>% 
  dplyr::select(mud, geometry)
st_geometry(mudsf) <- 'geometry'

# Order by mud
mudsf <- mudsf[with(mudsf, order(mud)),]
rownames(mudsf) <- NULL

# Set to same projection as grid
mudsf <- st_transform(mudsf, st_crs(grid))
mudsf <- st_make_valid(mudsf)

# Intersect
mudgrid <- st_intersection(mudsf, grid)

# Convert to dataframe
mudgrid <- mudgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, mud) %>% 
  group_by(ID) %>% 
  # Mean rugosity in each grid cell
  summarise(mud =     mean(mud,     na.rm=T))

# Re-merge to grid by ID
grid <- left_join(grid, mudgrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=mud), color=NA) +
  scale_fill_viridis_c(na.value='transparent')+
  geom_sf(data=regions, fill=NA, col='black')

# Sand
sandsf <- seds %>% 
  dplyr::select(sand, geometry)
st_geometry(sandsf) <- 'geometry'

# Order by sand
sandsf <- sandsf[with(sandsf, order(sand)),]
rownames(sandsf) <- NULL

# Set to same projection as grid
sandsf <- st_transform(sandsf, st_crs(grid))
sandsf <- st_make_valid(sandsf)

# Intersect
sandgrid <- st_intersection(sandsf, grid)

# Convert to dataframe
sandgrid <- sandgrid %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, sand) %>% 
  group_by(ID) %>% 
  # Mean rugosity in each grid cell
  summarise(sand =     mean(sand,     na.rm=T))

# Re-merge to grid by ID
grid <- left_join(grid, sandgrid, by=c('ID'))
ggplot(data=grid) +
  geom_sf(aes(fill=sand), color=NA) +
  scale_fill_viridis_c(na.value='transparent')+
  geom_sf(data=regions, fill=NA, col='black')

# Remove extraneous objects
rm(list=setdiff(ls(), c('coast', 'grid', 'region', 'regions')))

#### Conditional effect of variables on presence/absence of cod ####
# Load VAST fit data
load(here("VAST_runs/medium/Overall_BC/ALL/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.Rdata"))
medium <- fit
load(here("VAST_runs/small/Overall_BC/ALL/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.Rdata"))
small <- fit
load(here("VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.Rdata"))
large <- fit

# Clean workspace
rm(list=setdiff(ls(), c('small', 'medium', 'large',
                        'coast', 'grid', 'region', 'regions')))

# Set WD for plotting
out_dir = here("VAST_runs/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Name size classes
sizes <- c('small', 'medium', 'large')

# Blank dataframe to fill
vast.cov.effs <- data.frame(
  fit=NA, se=NA, lower=NA, upper=NA, Lin_pred=NA, Value=NA,
  Covariate=NA, Size=NA
)

# Extract data, looping through sizes
for(i in 1:length(sizes)){
  size <- sizes[i]
  if(size == 'small'){fit <- small}
  if(size == 'medium'){fit <- medium}
  if(size == 'large'){fit <- large}
  
  # Call parameter names
  params <- colnames(fit$covariate_data)
  params <- params[params %notin% c('Lon', 'Lat', 'Year')]
  
  # Parameters not used in certain models
  if(size == 'small'){params <- params[params %notin% c('amo')]}
  if(size == 'large'){params <- params[params %notin% c('cobble_P',
                                                        'sand_P', 
                                                        'mud_P')]}
  # Call category names
  catnames <- fit$category_names
  ncat = length(catnames)
  
  # Extract effects
  vast_covariate_effects<- get_vast_covariate_effects(vast_fit = fit, 
                                                      params_plot = c(params), 
                                                      params_plot_levels = 500, 
                                                      effects_pad_values = c(), 
                                                      nice_category_names = paste0(size, ' Cod'),
                                                      out_dir = out_dir,
                                                      category_to_use = 1,
                                                      ncat = ncat)
  
  # Parameters we want to use
  names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")
  
  # Merge non-informative names to informative name descriptions
  betternames=data.frame(
    Covariate = c('cobble_P', 'gravel_P', 'mud_P', 'sand_P', 'rugos', 'BATHY.DEPTH',
                  'h_bt', 'nao', 'amo'),
    Better = c('Cobble', 'Gravel', 'Mud', 'Sand', 'Rugosity', 'Depth (m)', 
               'BT (C)', 'NAO', 'AMO')
  )
  
  # Remove NA values
  vast_cov_eff_l <- vast_covariate_effects %>%
    drop_na(Value)
  
  # Merge to better names
  vast_cov_eff_l <- merge(vast_cov_eff_l, betternames)
  
  # Clean dataframe
  vast_cov_eff_l <- vast_cov_eff_l %>% 
    dplyr::select(-Covariate) %>% 
    rename(Covariate = Better)
  
  # Name size model
  vast_cov_eff_l$Size <- size
  
  # Bind to blank dataframe
  vast.cov.effs <- rbind(vast.cov.effs, vast_cov_eff_l)
  
  # Clean workspace
  rm(size, params, catnames, ncat, vast_covariate_effects, names_stay,
     betternames, vast_cov_eff_l)
  
}

# Remove blank initiation dataframe
vast.cov.effs <- vast.cov.effs[!is.na(vast.cov.effs$Size),]

# Clean
vast.cov.effs <- vast.cov.effs %>% 
  mutate(Size = str_to_sentence(Size)) %>% 
  mutate(Size = factor(Size, levels = c('Small', 'Medium', 'Large')),
         Covariate = factor(Covariate, levels=c('BT (C)',
                                                'Depth (m)',
                                                'Gravel',
                                                'Sand',
                                                'Mud',
                                                'Cobble',
                                                'Rugosity',
                                                'NAO',
                                                'AMO')))

# Save results for presence/absence
vast.cov.effs <- vast.cov.effs %>% 
  filter(Lin_pred == 'X1')

# Clean workspace
rm(list=setdiff(ls(), c('vast.cov.effs', '%notin%',
                        'grid', 'coast', 'region', 'regions')))

#### Predict pres/abs curve fit to depth ####
# Make column for predicted fit to curve
usedepths <- data.frame(
  bathy = unique(grid$bathy)
)
usedepths$predfit <- NA

# Run through size classes
sizes <- c('Small', 'Medium', 'Large')

# Blank list to collect suitability of cell values
depth.suit <- vector('list', length=length(sizes))

for(i in 1:length(sizes)){
  size <- sizes[i]
  
  # Extract depth fit curve
  val_curv <- vast.cov.effs %>% 
    filter(Lin_pred == 'X1' & Covariate == 'Depth (m)') %>% 
    filter(Size == paste0(size)) %>% 
    dplyr::select(Value, fit) %>% 
    rename(x=Value, y=fit)
  
  # Extend modeling to values outside seen in surveys
  ggplot(data=val_curv) +
    geom_line(aes(x=x, y=y)) +
    geom_smooth(aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam')
  
  # Model what we have with a GAM, use high knot value to get 100% Dev exp
  gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
  s <- summary(gam)
  message(paste0(round(s$dev.expl * 100, 3), ' % dev exp'))
  
  # Make df with depth values we want to predict to
  newdata <- data.frame(x=c(0.1, seq((max(val_curv$x, na.rm = T) + 1), 
                                     max(grid$bathy, na.rm = T), by=1)))
  
  # Predict model fit at these depth values
  pgam <- mgcv::predict.gam(gam, newdata=newdata)
  
  # Append fits to depth values
  newdata$y <- pgam
  
  # Sanity check
  print(
    ggplot() +
      geom_line(data=val_curv,
                aes(x=x, y=y)) +
      geom_smooth(data=val_curv,
                  aes(x=x, y=y),
                  formula = y ~ s(x, bs = "cs", k=15),
                  method='gam') +
      geom_point(data=newdata, 
                 aes(x=x, y=y))
  )
  
  # Bind into single dataframe
  val_curv <- rbind(val_curv, newdata)
  val_curv <- val_curv[with(val_curv, order(x)),]
  rownames(val_curv) <- NULL
  
  class.usedepths <- usedepths
  
  # Replace depth with curve fit value
  for(m in 1:nrow(class.usedepths)){
    newval <- class.usedepths$bathy[m]
    
    newval <- val_curv$y[val_curv$x==newval]
    
    if(length(newval) == 1){
      class.usedepths$predfit[m] <- newval
    }
    
    # If x-value has no exact match, find linear model to n5 points centered on nearest value
    if(length(newval)==0){
      target.index <- which(abs(val_curv$x - class.usedepths$bathy[m]) == 
                              min(abs(val_curv$x - class.usedepths$bathy[m])))[1]
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=class.usedepths$bathy[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      class.usedepths$predfit[m] <- newval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
    }
    rm(newval)
  }
  
  # Save to suitability list
  colnames(class.usedepths) <- c('bathy', paste0(paste0(size),'.bathy'))
  
  # Sanity check
  ggplot(data=class.usedepths) + 
    geom_point(aes(x=bathy, y=Small))
  
  depth.suit[[i]] <- class.usedepths
  
  rm(gam, s, class.usedepths, pgam, val_curv)
  
}

# Clean workspace
rm(usedepths)

# Name list items
names(depth.suit) <- c('Small', 'Medium', 'Large')

# Merge depth suitability to grid
grid <- left_join(grid, depth.suit[['Small']], by=c('bathy'))
grid <- left_join(grid, depth.suit[['Medium']], by=c('bathy'))
grid <- left_join(grid, depth.suit[['Large']], by=c('bathy'))

# Clean workspace
rm(list=setdiff(ls(), c('grid', 'region', 'regions', 'coast', 
                        'vast.cov.effs')))

#### Predict pres/abs curve fit to rugosity ####
# Make column for predicted fit to curve
userugos <- data.frame(
  rugosity = unique(grid$rugosity)
)
userugos$predfit <- NA
userugos <- userugos[!is.na(userugos$rugosity),]

# Run through size classes
sizes <- c('Small', 'Medium', 'Large')

# Blank list to collect suitability of cell values
rugos.suit <- vector('list', length=length(sizes))

for(i in 1:length(sizes)){
  size <- sizes[i]
  
  # Extract depth fit curve
  val_curv <- vast.cov.effs %>% 
    filter(Lin_pred == 'X1' & Covariate == 'Rugosity') %>% 
    filter(Size == paste0(size)) %>% 
    dplyr::select(Value, fit) %>% 
    rename(x=Value, y=fit) %>% 
    #mutate(x = round(x, 6)) %>% 
    unique()
  
  # Extend modeling to values outside seen in surveys
  ggplot(data=val_curv) +
    geom_line(aes(x=x, y=y)) +
    geom_smooth(aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam')
  
  # Model what we have with a GAM, use high knot value to get 100% Dev exp
  gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
  s <- summary(gam)
  message(paste0(round(s$dev.expl * 100, 3), ' % dev exp'))
  
  # Make df with depth values we want to predict to
  newdata <- data.frame(x=userugos$rugosity)
  
  # Predict model fit at these depth values
  pgam <- mgcv::predict.gam(gam, newdata=newdata)
  
  # Append fits to depth values
  newdata$y <- pgam
  
  # Sanity check
  ggplot() +
    geom_line(data=val_curv,
              aes(x=x, y=y)) +
    geom_smooth(data=val_curv,
                aes(x=x, y=y),
                formula = y ~ s(x, bs = "cs", k=15),
                method='gam') +
    geom_point(data=newdata, 
               aes(x=x, y=y))
  
  # Bind into single dataframe
  val_curv <- rbind(val_curv, newdata)
  val_curv <- val_curv[with(val_curv, order(x)),]
  rownames(val_curv) <- NULL
  
  class.userugos <- userugos
  
  # Replace depth with curve fit value
  for(m in 1:nrow(class.userugos)){
    newval <- class.userugos$rugosity[m]
    
    newval <- val_curv$y[val_curv$x==newval]
    
    if(length(newval) == 1){
      class.userugos$predfit[m] <- newval
    }
    
    # If x-value has no exact match, find linear model to n5 points centered on nearest value
    if(length(newval)==0){
      target.index <- which(abs(val_curv$x - class.userugos$rugosity[m]) == 
                              min(abs(val_curv$x - class.userugos$rugosity[m])))[1]
      test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                        (target.index + 2),])
      
      newdata <- data.frame(x=class.userugos$rugosity[m], y=NA)
      predval <- as.numeric(predict.lm(test, newdata))
      newval <- predval
      
      class.userugos$predfit[m] <- newval
      
      # Remove intermediates
      rm(target.index, test, newdata, predval)
    }
    rm(newval)
  }
  
  # Save to suitability list
  colnames(class.userugos) <- c('rugosity', paste0(paste0(size), '.rugos'))
  
  rugos.suit[[i]] <- class.userugos
  
  rm(gam, s, class.userugos, pgam, val_curv)
  
}

# Clean workspace
rm(userugos)

# Name list items
names(rugos.suit) <- c('Small', 'Medium', 'Large')

# Merge rugosity suitability to grid
grid <- left_join(grid, rugos.suit[['Small']], by=c('rugosity'))
grid <- left_join(grid, rugos.suit[['Medium']], by=c('rugosity'))
grid <- left_join(grid, rugos.suit[['Large']], by=c('rugosity'))

rm(list=setdiff(ls(), c('coast', 'grid', 'region', 'regions',
                        'vast.cov.effs')))

#### Predict pres/abs curve fit to sediments ####
# Name sediments
sedtypes <- c('cobble', 'gravel', 'mud', 'sand')

# Name sizes
sizes <- c('Small', 'Medium', 'Large')

# Blank list for sediments
all.seds.suit <- vector('list', length(sedtypes))

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Loop through sediments
for(q in 1:length(sedtypes)){
  message(sedtypes[q])
  
  # Call sediment type
  seduse <- sedtypes[q]
  
  # Make column for predicted fit to curve
  usesed <- grid %>% 
    dplyr::select(paste0(seduse)) %>%
    sfheaders::sf_to_df(fill=T) %>% 
    dplyr::select(paste0(seduse))
    
  usesed <- usesed[[1]]
  usesed <- usesed[!is.na(usesed)]
  usesed <- as.data.frame(usesed)
  
  colnames(usesed) <- c(paste0(seduse))
  usesed$predfit <- NA
  
  # Blank list to collect suitability of cell values
  sed.suit <- vector('list', length=length(sizes))
  
  # Run through sizes
  for(i in 1:length(sizes)){
    size <- sizes[i]
    message(paste0(size))
    
    # Extract depth fit curve
    val_curv <- vast.cov.effs %>% 
      filter(Lin_pred == 'X1' & Covariate == paste0(str_to_sentence(seduse))) %>% 
      filter(Size == paste0(size)) %>% 
      dplyr::select(Value, fit) %>% 
      rename(x=Value, y=fit)
    
    if(nrow(val_curv) == 0){
      sed.suit[[i]] <- data.frame(
        No = NA,
        No2 = NA
      )
      
      colnames(sed.suit[[i]]) <- c(paste0(seduse), paste0(size))
      
      next()
    }
    
    # Extend modeling to values outside seen in surveys
    ggplot(data=val_curv) +
      geom_line(aes(x=x, y=y)) +
      geom_smooth(aes(x=x, y=y),
                  formula = y ~ s(x, bs = "cs", k=15),
                  method='gam')
    
    # Model what we have with a GAM, use high knot value to get 100% Dev exp
    gam <- mgcv::gam(y ~ s(x, bs='cs', k=15), data=val_curv)
    
    # Make df with depth values we want to predict to
    newdata <- usesed[,1]
    newdata <- newdata[newdata %notin% val_curv$x]
    
    newdata <- newdata %>% 
      as.data.frame() %>% 
      rename(x = '.')
    
    # Predict model fit at these depth values
    pgam <- mgcv::predict.gam(gam, newdata=newdata)
    
    # Append fits to depth values
    newdata$y <- pgam
    
    # Sanity check
    print(
      ggplot() +
        geom_line(data=val_curv,
                  aes(x=x, y=y)) +
        geom_smooth(data=val_curv,
                    aes(x=x, y=y),
                    formula = y ~ s(x, bs = "cs", k=15),
                    method='gam') +
        geom_point(data=newdata, 
                   aes(x=x, y=y))
    )
    
    # Bind into single dataframe
    val_curv <- rbind(val_curv, newdata)
    val_curv <- unique(val_curv)
    val_curv <- val_curv[with(val_curv, order(x)),]
    rownames(val_curv) <- NULL
    
    class.usesed <- usesed
    
    # Replace depth with curve fit value
    for(m in 1:nrow(class.usesed)){
      newval <- class.usesed[,paste0(seduse),][m]
      
      newval <- val_curv$y[val_curv$x==newval]
      
      if(length(newval) == 1){
        class.usesed$predfit[m] <- newval
        rm(newval)
        next()
      }
      
      # If x-value has no exact match, find linear model to n5 points centered on nearest value
      if(length(newval)==0){
        target.index <- which(abs(val_curv$x - class.usesed[,paste0(seduse),][m]) == 
                                min(abs(val_curv$x - class.usesed[,paste0(seduse),][m])))[1]
        
        if(target.index >=3 & target.index <=(nrow(val_curv)-2)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index ==2){
          test <- lm(y ~ x, data=val_curv[(target.index - 1):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index ==1){
          test <- lm(y ~ x, data=val_curv[(target.index):
                                            (target.index + 2),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index == (nrow(val_curv)-1)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index + 1),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
        if(target.index == nrow(val_curv)){
          test <- lm(y ~ x, data=val_curv[(target.index - 2):
                                            (target.index),])
          
          newdata <- data.frame(x=class.usesed[,paste0(seduse),][m], y=NA)
          predval <- as.numeric(predict.lm(test, newdata))
          newval <- predval
          
          class.usesed$predfit[m] <- newval
          
          # Remove intermediates
          rm(target.index, test, newdata, predval, newval)
          next()
        }
        
      }
      rm(newval)
    }
    
    # Save to suitability list
    class.usesed <- unique(class.usesed)
    colnames(class.usesed) <- c(paste0(seduse), paste0(paste0(size),
                                                       '.', seduse))
    
    sed.suit[[i]] <- class.usesed
    
    rm(gam, class.usesed, pgam, val_curv)
    
  }
  
  all.seds.suit[[q]] <- sed.suit
  
}

# Clean workspace
rm(sed.suit, usesed, val_curv,i, m, q, sedtypes, seduse, size, sizes)

# Name list items
names(all.seds.suit) <- c('cobble', 'gravel', 'mud', 'sand')
for(i in 1:length(all.seds.suit)){
  names(all.seds.suit[[i]]) <- c('Small', 'Medium', 'Large')
}

# Merge to grid
for(q in 1:length(all.seds.suit)){
  for(j in 1:length(all.seds.suit[[q]])){
    grid <- left_join(grid, all.seds.suit[[q]][[j]])
  }
}

# Remove weird column this made
grid$Large <- NULL

# Clean grid
grid <- grid %>% 
  rename(geometry = grid)
st_geometry(grid) <- 'geometry'

# Save grid
saveRDS(grid, here('Habitat_Suitability/stationary_effects_grid.RDS'))
