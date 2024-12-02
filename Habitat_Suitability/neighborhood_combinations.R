# Density and HSI hotspot matching

rm(list=ls())

library(here)
library(sf)
library(tidyverse)
library(spdep)
library(sfdep)

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

# Load coast and stock strata
coast <- st_transform(ecodata::coast, crs="EPSG:26919")
regions <- st_transform(st_read(here('Data/GIS/codstox.shp'), quiet=T),
                        crs="EPSG:26919")

# Load density hotspot neighborhoods
dens <- st_read(here('Habitat_Suitability/large/Large_Density_Neighborhoods.shp'),
                quiet=T)
dens <- dens %>% 
  filter(!is.na(Season)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         Class = factor(Class, levels=c('Hot', 'Cold'))) %>%
  rename(density = y,
         Dens.Class = Class)

# Load HSI hotspot values
hsi <- st_read(here('Habitat_Suitability/large/large_hsi_neighborhoods.shp'),
               quiet=T)
hsi <- hsi %>% 
  filter(!is.na(Season)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         Class = factor(Class, levels=c('Hot', 'Cold'))) %>%
  rename(gm.HSI = ga,
         HSI.Class = Class)

years <- seq(1982, 2020, 1)
seasons <- c('Spring', 'Fall')

blank <- st_intersection(hsi[1,], dens[1,])

for(i in seasons){
  hsi.season <- hsi[hsi$Season == i,]
  dens.season <- dens[dens$Season == i,]
  
  for(j in years){
    hsi.year <- hsi.season[hsi.season$Year == j,]
    dens.year <- dens.season[dens.season$Year == j,]
    
    hsi.year <- st_cast(hsi.year, 'POLYGON')
    hsi.year <- unique(hsi.year)
    dens.year <- st_cast(dens.year, 'POLYGON')
    dens.year <- unique(dens.year)
    
    test <- st_intersection(hsi.year, dens.year)
    test <- st_collection_extract(
      test,
      type = c("POLYGON"),
      warn = TRUE
    )
    
    blank <- rbind(blank, test)
    rm(hsi.year, dens.year, test)

  }
  
}

ggplot(data=blank[blank$Year %in% seq(2017, 2021),]) +
  geom_sf(aes(fill=Dens.Class), alpha=0.4, col=NA) +
  geom_sf(aes(fill=HSI.Class), alpha=0.4, col=NA) +
  geom_sf(data=coast, fill='gray90', col='gray80') +
  coord_sf(xlim=c(st_bbox(regions)[1], st_bbox(regions)[3]),
           ylim=c(st_bbox(regions)[2], st_bbox(regions)[4])) +
  ggh4x::facet_grid2(vars(Season), vars(Year))

blank <- blank %>% 
  dplyr::select(-Year.1, -Season.1, -nghbrhd, -nghbrhd.1) %>% 
  rename(area.nb.HSI = area,
         area.nb.dens = area.1) %>% 
  mutate(ID = seq(1, nrow(blank), 1))

# Load gridded HSI
grid.hsi <- st_read(here('Habitat_Suitability/Large/gridded_hsi_suitability.shp'))
grid.hsi <- grid.hsi %>% 
  dplyr::select(-ID, -gi, -clssfct) %>% 
  rename(hsi.p = p_fldd_,
         grid.gm.hsi = ga)

# Load gridded dens
grid.dens <- st_read(here('Habitat_Suitability/Large/gridded_density_hotspots.shp'))
grid.dens <- grid.dens %>% 
  dplyr::select(-ID, -gi, -clssfct, -logdens) %>% 
  rename(dens.p = p_fldd_,
         grid.dens = dens)

# Merge nbs and density
grid.list <- split(blank, f=blank$Season)
for(i in 1:length(grid.list)){
  grid.list[[i]] <- split(grid.list[[i]], f=grid.list[[i]]$Year)
}

for(i in 1:length(grid.list)){
  for(j in 1:length(grid.list[[i]])){
    print(j)
    test <- grid.list[[i]][[j]]
    test.hsi <- grid.hsi[grid.hsi$Season == test$Season[1] &
                         grid.hsi$Year == test$Year[1],]
    test.dens <- grid.dens[grid.dens$Season == test$Season[1] &
                             grid.dens$Year == test$Year[1],]
    
    test <- st_intersection(test, test.hsi)
    test <- st_intersection(test, test.dens)
    
    test <- test %>% 
      dplyr::select(-Year.1, -Year.2, -Season.1, -Season.2)
    
    grid.list[[i]][[j]] <- test
    
  }
}

for(i in 1:length(grid.list)){
  grid.list[[i]] <- do.call(rbind, grid.list[[i]])
}
grid.list <- do.call(rbind, grid.list)

head(grid.list)

grid.list <- st_collection_extract(
  grid.list,
  type = c("POLYGON"),
  warn = TRUE
)
rownames(grid.list) <- NULL

grid.ag <- grid.list %>% 
  group_by(Year, Season, ID) %>% 
  summarise(bathy = mean(bathy),
            rugosty = mean(rugosty),
            gravel = mean(gravel),
            bt = mean(bt),
            Lr_Dpth = mean(Lr_Dpth),
            Lar_Rgs = mean(Lar_Rgs),
            Lr_Grvl = mean(Lr_Grvl),
            Lar_BT = mean(Lar_BT),
            
            gm.HSI = mean(gm.HSI),
            density = mean(density),
            
            grid.gm.hsi = mean(grid.gm.hsi),
            grid.dens = mean(grid.dens),
            
            HSI.Class = HSI.Class[1],
            Dens.Class = Dens.Class[1],
            )

# Split by region
#grid.ag <- st_intersection(grid.ag, regions)
grid.ag$area <- st_area(grid.ag)
grid.ag$area <- strip_units(grid.ag$area) / 1000000


tsarea <- grid.ag %>% 
  group_by(Year, Season, HSI.Class, Dens.Class) %>% 
  summarise(totalarea = sum(area),
            density = mean(density),
            hsi = mean(gm.HSI),
            
            bathy = mean(bathy),
            rugosty = mean(rugosty),
            gravel = mean(gravel),
            bt = mean(bt),
            Lr_Dpth = mean(Lr_Dpth),
            Lar_Rgs = mean(Lar_Rgs),
            Lr_Grvl = mean(Lr_Grvl),
            Lar_BT = mean(Lar_BT),
            
            grid.gm.hsi = mean(grid.gm.hsi),
            grid.dens = mean(grid.dens),)

tsarea.matches <- tsarea[tsarea$Dens.Class == tsarea$HSI.Class,]

ggplot(data=tsarea.matches) +
         geom_line(aes(x=Year, y=Lar_BT, col=Dens.Class)) +
         facet_wrap(vars(Season))

# View
ggplot(data=grid.ag[grid.ag$Dens.Class == 'Hot' &
                    grid.ag$HSI.Class == 'Hot',]) +
  geom_point(aes(x=Year, y=density)) +
  ggh4x::facet_grid2(vars(Season),
                     scales='free_x')

ggplot() +
  ggpattern::geom_sf_pattern(data=grid.ag[grid.ag$Year == 1982 &
                                          grid.ag$HSI.Class == 'Hot',],
                             pattern = 'stripe',
                  fill    = 'transparent',
                  color="#F8766D",
                  pattern_color="#F8766D",
                  pattern_spacing = 0.005,
                  pattern_alpha=1,
                  pattern_fill= 'transparent') +
  ggpattern::geom_sf_pattern(data=grid.ag[grid.ag$Year == 1982 &
                                            grid.ag$HSI.Class == 'Cold',],
                             pattern = 'stripe',
                             fill    = 'transparent',
                             color="#00BFC4",
                             pattern_color="#00BFC4",
                             pattern_spacing = 0.005,
                             pattern_alpha=1,
                             pattern_fill= 'transparent') +
  geom_sf(data=grid.ag[grid.ag$Year == 1982,],
          aes(fill=Dens.Class), alpha=0.4, col=NA) +
  geom_sf(data=regions, fill=NA) +
  facet_wrap(vars(Season))
