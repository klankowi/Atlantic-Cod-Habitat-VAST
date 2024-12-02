# Bottom temperature in hotspots

rm(list=ls())

library(here)
library(sf)
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

# Scale function
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

# Load bottom temperature
load(here('Habitat_Suitability/Seasonal_Avg_Bottom_Temps2.RData'))

# Load spatiotemp density
small.dens <- st_read(here('Habitat_Suitability/Small/small_dens.shp'),quiet=T)
med.dens <- st_read(here('Habitat_Suitability/Medium/medium_dens.shp'), quiet=T)
lar.dens <- st_read(here('Habitat_Suitability/Large/large_dens.shp'), quiet=T)

small.dens <- small.dens %>% 
  dplyr::select(-logy, -gi, -p_fldd_) %>% 
  rename(small.density = y,
         small.cf = clssfct,
         small.Class = Class)

med.dens <- med.dens %>% 
  dplyr::select(-logy, -gi, -p_fldd_) %>% 
  rename(med.density = y,
         med.cf = clssfct,
         med.Class = Class) %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(Year, Season, ID, med.density, med.cf, med.Class) %>% 
  unique() %>% as.data.frame()

lar.dens <- lar.dens %>% 
  dplyr::select(-logy, -gi, -p_fldd_) %>% 
  rename(lar.density = y,
         lar.cf = clssfct,
         lar.Class = Class) %>%   
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(Year, Season, ID, lar.density, lar.cf, lar.Class) %>% 
  unique() %>% as.data.frame()

density <- left_join(small.dens, med.dens)
density <- left_join(density, lar.dens)
rm(small.dens, med.dens, lar.dens)

# Load other suitability
static <- read_rds(here("Habitat_Suitability/wholegrid.RDS"))
static <- st_make_valid(static)

static <- static %>% 
  dplyr::select(-y) %>% 
  rename(geometry = grid)
st_geometry(static) <- 'geometry'
static <- st_transform(static, crs="EPSG:26919")

# Call region outline
region <- st_read(here('Data/GIS/cod_region_wgs.shp'),
                  quiet=T)
region <- st_transform(region, crs="EPSG:26919")
region <- region %>% 
  dplyr::select(-Shape_Leng, -Shape_Area, -area_km) %>% 
  rename(holder = OBJECTID)

# Pin static grid to region
static <- st_intersection(static, region)
static$holder <- NULL
static <- static %>% 
  st_cast('MULTIPOLYGON') %>% 
  st_cast('POLYGON')

# Pin density grid to region
density <- st_intersection(density, region)
density$holder <- NULL
density <- density %>% 
  st_cast('MULTIPOLYGON') %>% 
  st_cast('POLYGON')

# Split density grid to fall and spring lists (items yearly)
dens.fall <- density[density$Season == 'Fall',]
dens.fall <- split(dens.fall, f=dens.fall$Year)
dens.fall <- dens.fall[-40]

dens.spring <- density[density$Season == 'Spring',]
dens.spring <- split(dens.spring, f=dens.spring$Year)
dens.spring <- dens.spring[-40]

# Geometric mean function
"geometric.mean" <- function(x,na.rm=FALSE){ 
    exp(mean(log(x),na.rm=TRUE)) 
}

# Scale function
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


# Rasterize BT and static to match density grid
for(i in 1:length(big.spring)){
  print(i)
  
  # Rename bottom temperature variables
  big.spring[[i]] <- big.spring[[i]] %>% 
    rename(Sma.BT = Small,
           Med.BT = Medium,
           Lar.BT = Large,
           )
  
  # Transform to appropriate CRS
  big.spring[[i]] <- st_transform(big.spring[[i]],
                                  st_crs(static))
  
  # Clean
  big.spring[[i]] <- big.spring[[i]] %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  
  # Intersect static grid with bottom temp (for time step in spring)
  test <- st_intersection(static, big.spring[[i]])
  
  # CLean
  test <- test %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  
  # Group into static-grid cells
  test <- test %>% 
    group_by(Year, Season, ID) %>% 
    summarise(bt = mean(bt, na.rm=T),
              Sma.BT = mean(Sma.BT, na.rm=T),
              Med.BT = mean(Med.BT, na.rm=T),
              Lar.BT = mean(Lar.BT, na.rm=T),
              
              bathy = mean(bathy, na.rm=T),
              Sma.Depth = mean(Sma.Depth, na.rm=T),
              Med.Depth = mean(Med.Depth, na.rm=T),
              Lar.Depth = mean(Lar.Depth, na.rm=T),
              
              rugosity = mean(rugosity, na.rm=T),
              Sma.Rugos = mean(Sma.Rugos, na.rm=T),
              Med.Rugos = mean(Med.Rugos, na.rm=T),
              Lar.Rugos = mean(Lar.Rugos, na.rm=T),
              
              cobble = mean(cobble, na.rm=T),
              Sma.Cobble = mean(Sma.Cobble, na.rm=T),
              Med.Cobble = mean(Med.Cobble, na.rm=T),
              
              gravel = mean(gravel, na.rm=T),
              Sma.Gravel = mean(Sma.Gravel, na.rm=T),
              Med.Gravel = mean(Med.Gravel, na.rm=T),
              Lar.Gravel = mean(Lar.Gravel, na.rm=T),
              
              sand = mean(sand, na.rm=T),
              Sma.Sand = mean(Sma.Sand, na.rm=T),
              Med.Sand = mean(Med.Sand, na.rm=T),
              
              mud = mean(mud, na.rm=T),
              Sma.Mud = mean(Sma.Mud, na.rm=T),
              Med.Mud = mean(Med.Mud, na.rm=T)
              
              )%>% 
    sfheaders::sf_to_df(fill=T) %>% 
    dplyr::select(Year, Season, ID, bt, Sma.BT,
                  Med.BT, Lar.BT) %>% 
    unique() %>% 
    as.data.frame()
  
  # Rejoin to original sf grid
  test <- left_join(static, test, by=c('ID'))
  
  # Join with density
  test <- st_intersection(dens.spring[[i]], test)
  
  # Clean
  test$dim <- st_dimension(test)
  test <- test[test$dim == 2,]
  test <- test %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  test <- test %>% 
    dplyr::select(-ID, -ID.1, -dim, -Season.1, -Year.1)
  
  big.spring[[i]] <- test
  
  rm(test)
  
}

for(i in 1:length(big.fall)){
  print(i)
  
  # Rename bottom temperature variables
  big.fall[[i]] <- big.fall[[i]] %>% 
    rename(Sma.BT = Small,
           Med.BT = Medium,
           Lar.BT = Large,
    )
  
  # Transform to appropriate CRS
  big.fall[[i]] <- st_transform(big.fall[[i]],
                                  st_crs(static))
  
  # Clean
  big.fall[[i]] <- big.fall[[i]] %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  
  # Intersect static grid with bottom temp (for time step in fall)
  test <- st_intersection(static, big.fall[[i]])
  
  # CLean
  test <- test %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  
  # Group into static-grid cells
  test <- test %>% 
    group_by(Year, Season, ID) %>% 
    summarise(bt = mean(bt, na.rm=T),
              Sma.BT = mean(Sma.BT, na.rm=T),
              Med.BT = mean(Med.BT, na.rm=T),
              Lar.BT = mean(Lar.BT, na.rm=T),
              
              bathy = mean(bathy, na.rm=T),
              Sma.Depth = mean(Sma.Depth, na.rm=T),
              Med.Depth = mean(Med.Depth, na.rm=T),
              Lar.Depth = mean(Lar.Depth, na.rm=T),
              
              rugosity = mean(rugosity, na.rm=T),
              Sma.Rugos = mean(Sma.Rugos, na.rm=T),
              Med.Rugos = mean(Med.Rugos, na.rm=T),
              Lar.Rugos = mean(Lar.Rugos, na.rm=T),
              
              cobble = mean(cobble, na.rm=T),
              Sma.Cobble = mean(Sma.Cobble, na.rm=T),
              Med.Cobble = mean(Med.Cobble, na.rm=T),
              
              gravel = mean(gravel, na.rm=T),
              Sma.Gravel = mean(Sma.Gravel, na.rm=T),
              Med.Gravel = mean(Med.Gravel, na.rm=T),
              Lar.Gravel = mean(Lar.Gravel, na.rm=T),
              
              sand = mean(sand, na.rm=T),
              Sma.Sand = mean(Sma.Sand, na.rm=T),
              Med.Sand = mean(Med.Sand, na.rm=T),
              
              mud = mean(mud, na.rm=T),
              Sma.Mud = mean(Sma.Mud, na.rm=T),
              Med.Mud = mean(Med.Mud, na.rm=T)
              
    )%>% 
    sfheaders::sf_to_df(fill=T) %>% 
    dplyr::select(Year, Season, ID, bt, Sma.BT,
                  Med.BT, Lar.BT) %>% 
    unique() %>% 
    as.data.frame()
  
  # Rejoin to original sf grid
  test <- left_join(static, test, by=c('ID'))
  
  # Join with density
  test <- st_intersection(dens.fall[[i]], test)
  
  # Clean
  test$dim <- st_dimension(test)
  test <- test[test$dim == 2,]
  test <- test %>% 
    st_cast('MULTIPOLYGON') %>% 
    st_cast('POLYGON')
  test <- test %>% 
    dplyr::select(-ID, -ID.1, -dim, -Season.1, -Year.1)
  
  big.fall[[i]] <- test
  
  rm(test)
  
}

# Rebind
spring <- do.call(rbind, big.spring)
fall <- do.call(rbind, big.fall)

all <- rbind(spring, fall)
rownames(all) <- NULL

st_write(all, 
         here('Habitat_Suitability/Gridded_Density_RawVars.shp'),
         quiet=T)

# Explore
colnames(all) <- tolower(colnames(all))
bathy <- all %>% 
  dplyr::select(bathy, sma.depth, med.depth, lar.depth) %>% 
  pivot_longer(cols=c('sma.depth', 'med.depth', 'lar.depth'),
               names_to = 'size', 
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(bathy, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=bathy) + geom_point(aes(x=bathy, y=response, col=size))

temp <- all %>% 
  dplyr::select(bt, sma.bt, med.bt, lar.bt) %>% 
  pivot_longer(cols = c('sma.bt', 'med.bt', 'lar.bt'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(bt, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=temp) + geom_point(aes(x=bt, y=response, col=size))  

rugos <- all %>% 
  dplyr::select(rugosity, sma.rugos, med.rugos, lar.rugos) %>% 
  pivot_longer(cols = c('sma.rugos', 'med.rugos', 'lar.rugos'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(rugosity, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=rugos) + geom_point(aes(x=rugosity, y=response, col=size))  

cobble <- all %>% 
  dplyr::select(cobble, sma.cobble, med.cobble, lar.cobble) %>% 
  pivot_longer(cols = c('sma.cobble', 'med.cobble', 'lar.cobble'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(cobble, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=cobble) + geom_point(aes(x=cobble, y=response, col=size))  

gravel <- all %>% 
  dplyr::select(gravel, sma.gravel, med.gravel, lar.gravel) %>% 
  pivot_longer(cols = c('sma.gravel', 'med.gravel', 'lar.gravel'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(gravel, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=gravel) + geom_point(aes(x=gravel, y=response, col=size))  

mud <- all %>% 
  dplyr::select(mud, sma.mud, med.mud, lar.mud) %>% 
  pivot_longer(cols = c('sma.mud', 'med.mud', 'lar.mud'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(mud, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=mud) + geom_point(aes(x=mud, y=response, col=size))  

sand <- all %>% 
  dplyr::select(sand, sma.sand, med.sand, lar.sand) %>% 
  pivot_longer(cols = c('sma.sand', 'med.sand', 'lar.sand'),
               names_to = 'size',
               values_to = 'response') %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(sand, size, response) %>% 
  unique() %>% 
  as.data.frame()
ggplot(data=sand) + geom_point(aes(x=sand, y=response, col=size))  


# Geometric and arithmetic average HSI for larges
"geometric.mean" <- 
  function(x,na.rm=TRUE){ 
    exp(mean(log(x),na.rm=na.rm))
    }

#all$aa <- NA
all$ga <- NA
for(i in 1:nrow(all)){
  print(i)
  #all$aa[i] <- mean(c(
  #                           all$Lar.Depth[i],
  #                           all$Lar.Rugos[i],
  #                           all$Lar.Gravel[i],
  #                           all$Lar.BT[i]))
  
  all$ga[i] <- geometric.mean(c(
    all$Lar.Depth[i],
    all$Lar.Rugos[i],
    all$Lar.Gravel[i],
    all$Lar.BT[i]))
  
  
}

all.list <- split(all, f=all$Year)

# Call regions
region <- st_read(here('Data/GIS/cod_region_wgs.shp'), quiet=T)
region <- st_transform(region, st_crs(nbd))

# Split nbd into seasons and years
nbd$Class <- factor(nbd$Class, levels=c('Hot', 'Cold'))
nbdlist <- split(nbd, f=nbd$Year)

# Call persistent spots
persist <- st_read(here("Habitat_Suitability/large_persistence.shp"),
                   quiet=T)
nnames <- data.frame(
  Class = persist$Class,
  Season = persist$Season,
  nghbrhd = persist$nghbrhd,
  Name = c('NE Peak',
           'SNE',
           'EGOM',
           'NE Peak',
           'Stellwagen',
           'SNE',
           'EGOM')
)
persist <- left_join(persist, nnames, by=c('Class', 'Season', 'nghbrhd'))

# Overlap persistent spots with grid
all2 <- st_intersection(all, persist)
all2 <- all2 %>% 
  filter(Season == Season.1) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         Class = factor(Class, levels=c('Hot', 'Cold'))) %>% 
  dplyr::select(-Season.1)

# Annual features
annual <- all2 %>% 
  group_by(Year, Season, Class, Name) %>% 
  summarise(bathy = mean(bathy, na.rm=T),
            rugosity = mean(rugosity, na.rm=T),
            gravel = mean(gravel, na.rm=T),
            bt = mean(bt, na.rm=T),
            Lar.Depth = mean(Lar.Depth, na.rm=T),
            Lar.Rugos = mean(Lar.Rugos, na.rm=T),
            Lar.Gravel = mean(Lar.Gravel, na.rm=T),
            Lar.BT = mean(Lar.BT, na.rm=T),
            ga = mean(ga, na.rm=T))

annual$Name <- as.factor(annual$Name)

# Plot change
ggplot(data=annual) +
  geom_line(aes(x=Year, y=ga, col=Name)) +
  geom_point(aes(x=Year, y=ga, col=Name)) +
  ggh4x::facet_grid2(vars(Season))


# Overlap stock areas with grid
region <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
region <- st_transform(region, st_crs(all2))
region <- st_make_valid(region)
all2 <- st_intersection(all, region)
all2 <- all2 %>% 
  #filter(Season == Season.1) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         STOCK = factor(STOCK, levels=c('EGOM', 'GBK',
                                        'SNE', 'WGOM')))

# Annual features
annual <- all2 %>% 
  group_by(Year, Season, STOCK) %>% 
  summarise(bathy = mean(bathy, na.rm=T),
            rugosity = mean(rugosity, na.rm=T),
            gravel = mean(gravel, na.rm=T),
            bt = mean(bt, na.rm=T),
            Lar.Depth = mean(Lar.Depth, na.rm=T),
            Lar.Rugos = mean(Lar.Rugos, na.rm=T),
            Lar.Gravel = mean(Lar.Gravel, na.rm=T),
            Lar.BT = mean(Lar.BT, na.rm=T),
            ga = mean(ga, na.rm=T))

#annual$STOCK <- as.factor(annual$STOCK)
annual <- annual[!is.na(annual$Season),]

# Plot change
ggplot(data=annual) +
  geom_line(aes(x=Year, y=ga, col=STOCK)) +
  geom_point(aes(x=Year, y=ga, col=STOCK)) +
  ggh4x::facet_grid2(vars(Season), vars(STOCK))

# Call BTS strata
strat <- st_read(here('Data/GIS/BTS_Strata_Updated.shp'),
                 quiet=T)
strat <- st_make_valid(st_transform(strat, st_crs(all)))
wholearea <- st_read(here('Data/GIS/cod_region_wgs.shp'),
                     quiet=T)
wholearea <- st_make_valid(st_transform(wholearea, st_crs(all)))
strat <- st_intersection(strat, wholearea)

strat <- st_intersection(strat, region)

all2 <- st_intersection(all, strat)
all2 <- all2 %>% 
  #filter(Season == Season.1) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         Strata = as.factor(STRATA_1))

# Annual features
annual <- all2 %>% 
  group_by(Year, Season, Strata, STOCK) %>% 
  summarise(bathy = mean(bathy, na.rm=T),
            rugosity = mean(rugosity, na.rm=T),
            gravel = mean(gravel, na.rm=T),
            bt = mean(bt, na.rm=T),
            Lar.Depth = mean(Lar.Depth, na.rm=T),
            Lar.Rugos = mean(Lar.Rugos, na.rm=T),
            Lar.Gravel = mean(Lar.Gravel, na.rm=T),
            Lar.BT = mean(Lar.BT, na.rm=T),
            ga = mean(ga, na.rm=T))

annual <- annual[!is.na(annual$Season),]

annual$Depth <- NA
for(i in 1:nrow(annual)){
  if(annual$bathy[i] <25){annual$Depth[i] = '1'}
  if(annual$bathy[i] >=25 & annual$bathy[i] < 50){annual$Depth[i] = '2'}
  if(annual$bathy[i] >=50 & annual$bathy[i] < 75){annual$Depth[i] = '3'}
  if(annual$bathy[i] >=75 & annual$bathy[i] < 100){annual$Depth[i] = '4'}
  if(annual$bathy[i] >=100 & annual$bathy[i] < 125){annual$Depth[i] = '5'}
  if(annual$bathy[i] >=125 & annual$bathy[i] < 150){annual$Depth[i] = '6'}
  if(annual$bathy[i] >=150 & annual$bathy[i] < 175){annual$Depth[i] = '7'}
  if(annual$bathy[i] >=175 & annual$bathy[i] < 200){annual$Depth[i] = '8'}
  if(annual$bathy[i] >=200 & annual$bathy[i] < 225){annual$Depth[i] = '9'}
  if(annual$bathy[i] >=225 & annual$bathy[i] < 250){annual$Depth[i] = '10'}
}

annual$Depth <- as.numeric(annual$Depth)
annual <- annual[!is.na(annual$STOCK),]

# Plot change
ggplot(data=annual[annual$Season == 'Fall' & annual$Depth ==2,]) +
  geom_line(aes(x=Year, y=Lar.BT, group=Strata)) +
  #geom_point(aes(x=Year, y=bt, col=STRATA_1)) +
  facet_wrap(vars(STOCK)) +
  theme(legend.position = 'n')

# Extract each time step's "best" sections ()
for(i in 1:
    length(all.list)
    ){
  test <- all.list[[i]][all.list[[i]]$ga >=0.5,]
  map <- 
  ggplot(data=test) +
    geom_sf(aes(fill=ga, col=ga)) +
    scale_fill_viridis_c(option='viridis', 
                         limits=c(0.5, 1)) +
    scale_color_viridis_c(option='viridis',
                          limits=c(0.5, 1)) +
    facet_wrap(vars(Season)) +
    ggtitle(paste0(test$Year[1])) +
    geom_sf(data=coast) +
    coord_sf(xlim=c(-140000, 776200),
             ylim=c(4050000, 4965900))
  ggsave(plot=map,
         filename = paste0(here('Habitat_Suitability/Large/Best/'),
                           '/', paste0(test$Year[1]),
                           '.png'))
  rm(test, map)
}

years <- seq(1982, 2020, 1)

for(i in 8:length(all.list)){
  
  all.list[[i]] <- st_intersection(all.list[[i]], region)
  
  pic <- 
    ggplot() +
      geom_sf(data=all.list[[i]],
              aes(fill=aa, col=aa)) +
      scale_fill_viridis_c(option='viridis',
                           limits=c(0,1),
                           na.value = 'transparent') +
      scale_color_viridis_c(option='viridis',
                            limits=c(0,1),
                            guide='none',
                            na.value = 'transparent')+
      ggnewscale::new_scale_color() +
      
      geom_sf(data=nbdlist[[i]],
              fill=NA, lwd=0.6, aes(col=Class)) +
      facet_wrap(vars(Season)) +
      labs(fill='Mean\nHSI') +
      ggtitle(paste0(all.list[[i]]$Year[1], ' Large Cod'))
  
  pic2 <- 
    ggplot() +
    geom_sf(data=all.list[[i]],
            aes(fill=ga, col=ga)) +
    scale_fill_viridis_c(option='viridis',
                         limits=c(0,1),
                         na.value = 'transparent') +
    scale_color_viridis_c(option='viridis',
                          limits=c(0,1),
                          guide='none',
                          na.value = 'transparent')+
    ggnewscale::new_scale_color() +
    
    geom_sf(data=nbdlist[[i]],
            fill=NA, lwd=0.6, aes(col=Class)) +
    facet_wrap(vars(Season)) +
    labs(fill='Geo. Mean\nHSI') +
    ggtitle(paste0(all.list[[i]]$Year[1], ' Large Cod'))
  
  ggsave(plot=pic,
         filename=paste0(here('Habitat_Suitability/Large/AA_HSI/'),
                         '/', years[i], '_aa_hsi.png'),
         width = 11, height = 8.5)
  
  ggsave(plot=pic2,
         filename=paste0(here('Habitat_Suitability/Large/GA_HSI/'),
                         '/', years[i], '_ga_hsi.png'),
         width = 11, height = 8.5)
  
}


# Split nbd into seasons and years
nbdlist <- split(nbd, f=nbd$Season)
for(i in 1:length(nbdlist)){
  nbdlist[[i]] <- split(nbdlist[[i]], f=nbdlist[[i]]$Year)
}

# Copy list
btlist <- nbdlist

# In each timestep, merge hotspots and bottom temp
for(i in 1:length(nbdlist)){
  print(i)
  if(names(nbdlist)[i]=='Fall'){
    for(j in 1:length(nbdlist[[i]])){
      big.fall[[j]] <- st_transform(big.fall[[j]],
                                      st_crs(nbdlist[[i]][[j]]))
      btlist[[i]][[j]] <- st_intersection(big.fall[[j]],
                                          nbdlist[[i]][[j]])
    }
  }
  if(names(nbdlist)[i]=='Spring'){
    for(j in 1:length(nbdlist[[i]])){
      big.spring[[j]] <- st_transform(big.spring[[j]],
                                      st_crs(nbdlist[[i]][[j]]))
      btlist[[i]][[j]] <- st_intersection(big.spring[[j]],
                                          nbdlist[[i]][[j]])
    }
  }
}

# REbind
for(i in 1:length(btlist)){
  btlist[[i]] <- do.call(rbind, btlist[[i]])
}
btlist <- do.call(rbind, btlist)
rownames(btlist) <- NULL

# Clean
bt <- btlist %>% 
  dplyr::select(-Year.1, -Season.1) %>% 
  mutate(Class = factor(Class, levels=c('Hot', 'Cold')),
         Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  filter(!is.na(nghbrhd))

means <- bt %>% 
  group_by(Year, Season, Class, nghbrhd) %>% 
  summarise(meanbt = mean(bt, na.rm=T))

ggplot(data=means[means$Season == 'Fall' & means$Class == 'Hot',]) +
  geom_sf(aes(fill=meanbt)) +
  facet_wrap(vars(Year)) +
  scale_fill_viridis_c(option='viridis')

ggplot(data=means[means$Season == 'Spring' & means$Class == 'Hot',]) +
  geom_point(aes(x=Year, y=meanbt)) +
  geom_smooth(aes(x=Year, y=meanbt), method='gam')

region <- st_read(here('Data/GIS/codstox.shp'), quiet=T)
region <- st_transform(region, st_crs(means))

fall <- do.call(rbind, big.fall)
fall <- st_transform(fall, st_crs(means))
spring <- do.call(rbind, big.spring)
spring <- st_transform(spring, st_crs(means))

all <- rbind(spring, fall)
all <- st_intersection(all, region)
rownames(all) <- NULL

allmeans <- all %>% 
  group_by(Year, Season, STOCK) %>% 
  summarise(meanbt = mean(bt))

ggplot(data=allmeans) +
  geom_line(aes(x=Year, y=meanbt, col=STOCK)) +
  facet_wrap(vars(Season)) +
  labs(y='Mean bottom temperature',
       col='Stock\nArea')

bt.st.means <- st_intersection(bt, region)

bt.st.means2 <- bt.st.means %>% 
  dplyr::select(-OBJECTID, -Shape_Leng, -Shape_Area) %>% 
  filter(!is.na(nghbrhd))

nbd$Class <- factor(nbd$Class, levels=c('Hot', 'Cold'))

years <- seq(1982, 2020, 1)

range01 <- function(x, ...){(x - min(x, ..., na.rm=T)) / (max(x, ..., na.rm=T) - min(x, ..., na.rm=T))}

all$Lar.TTS <- range01(all$Lar.Scale)

for(i in 1:39){
  year <- years[i]
  pic <- 
  
  ggplot() +
    geom_sf(data=all[all$Year == paste0(year),],
            aes(fill=bt, col=bt)) +
    scale_fill_viridis_c('viridis',
                         name='Seasonal Avg\n Bottom temp (C)',
                         limits=c(2, 20)) +
    scale_color_viridis_c('viridis',
                          guide='none',
                          limits=c(2, 20)) +
    ggnewscale::new_scale_color() +
    
    labs(col='Hotspot type') +
    
    geom_sf(data=nbd[nbd$Year == paste0(year),],
            fill=NA, aes(col=Class), lwd=0.6) +
    facet_wrap(vars(Season)) +
    ggtitle(paste0(year))
  
  ggsave(plot=pic,
         filename=paste0(here('Habitat_Suitability/Small/Temp/'),
                         '/', year, '_hotspots_w_overall_suit.png'),
         width = 11, height=8.5)
  
}

# List files in given directory
setwd(here('Habitat_Suitability/Small/Temp'))
files  <- list.files()
files <- paste0(getwd(), "/", files)

# Set GIF save location

# Convert PNGs to GIF
gifski::gifski(files, "small_cod_bottomtemp.gif", loop = FALSE, delay = 0.7)

