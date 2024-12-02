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

# Load bottom temperature
load(here('Habitat_Suitability/Seasonal_Avg_Bottom_Temps2.RData'))
rm(spring)

# Load other suitability
static <- read_rds(here("Habitat_Suitability/wholegrid.RDS"))
static <- static %>% 
  # dplyr::select(bathy, rugosity, gravel, 
  #               Lar.Depth, Lar.Rugos, Lar.Gravel, 
  #               grid) %>% 
  rename(geometry = grid)
st_geometry(static) <- 'geometry'
static$ID <- seq(1, nrow(static), 1)

# Load hotspots
nbd <- st_read(here('Data/GIS/hotspot_neighborhoods.shp'), quiet=T)
# Remove 2021 (no temp data)
nbd <- nbd[nbd$Year != 2021,]

# Load names
# hsn <- read.csv(here('Habitat_Suitability/hotspots.csv'))
# hsn <- hsn %>% 
#   dplyr::select(-y) %>% 
#   rename(nghbrhd = neighborhood)
# nbd <- left_join(nbd, hsn, by=c('Year', 'Season', 'Class', 'nghbrhd'))

# Rasterize BT to match grid of static
for(i in 1:length(big.spring)){
  print(i)
  big.spring[[i]] <- big.spring[[i]] %>% 
    rename(Sma.BT = Small,
           Med.BT = Medium,
           Lar.BT = Large,
           )
  
  big.spring[[i]] <- st_transform(big.spring[[i]],
                                  st_crs(static))
  
  test <- st_intersection(static, big.spring[[i]])
  
  test <- test %>% 
    group_by(ID) %>% 
    summarise(bt = mean(bt, na.rm=T),
              Sma.BT = mean(Sma.BT, na.rm=T),
              Med.BT = mean(Med.BT, na.rm=T),
              Lar.BT = mean(Lar.BT, na.rm=T)) %>% 
    mutate(Year = big.spring[[i]]$Year[1],
           Season = big.spring[[i]]$Season[1]) %>% 
    sfheaders::sf_to_df(fill=T) %>% 
    dplyr::select(-sfg_id, -multipolygon_id, -polygon_id,
                  -linestring_id, -x, -y) %>% 
    unique() %>% 
    as.data.frame()
  
  test <- left_join(static, test, by=c('ID'))
  
  big.spring[[i]] <- test
  
  rm(test)
  
}

# Rasterize BT to match grid of static
for(i in 1:length(big.fall)){
  print(i)
  big.fall[[i]] <- big.fall[[i]] %>% 
    rename(Sma.BT = Small,
           Med.BT = Medium,
           Lar.BT = Large,
    )
  
  big.fall[[i]] <- st_transform(big.fall[[i]],
                                  st_crs(static))
  
  test <- st_intersection(static, big.fall[[i]])
  
  test <- test %>% 
    group_by(ID) %>% 
    summarise(bt = mean(bt, na.rm=T),
      Sma.BT = mean(Sma.BT, na.rm=T),
      Med.BT = mean(Med.BT, na.rm=T),
      Lar.BT = mean(Lar.BT, na.rm=T)) %>% 
    mutate(Year = big.fall[[i]]$Year[1],
           Season = big.fall[[i]]$Season[1]) %>% 
    sfheaders::sf_to_df(fill=T) %>% 
    dplyr::select(-sfg_id, -multipolygon_id, -polygon_id,
                  -linestring_id, -x, -y) %>% 
    unique() %>% 
    as.data.frame()
  
  test <- left_join(static, test, by=c('ID'))
  
  big.fall[[i]] <- test
  
  rm(test)
  
}

# Rebind
spring <- do.call(rbind, big.spring)
fall <- do.call(rbind, big.fall)

all <- rbind(spring, fall)
rownames(all) <- NULL

all <- all[!is.na(all$y),]

# Geometric and arithmetic average HSI for larges
"geometric.mean" <- 
  function(x,na.rm=TRUE){ 
    exp(mean(log(x),na.rm=na.rm))
  }

# Scale function
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

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

