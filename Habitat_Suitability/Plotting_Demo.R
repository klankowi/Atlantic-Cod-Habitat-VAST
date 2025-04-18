# Data to identify density and suitability hotspots

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

# Load bottom temperature (list with grid for each year)
load(here('Habitat_Suitability/Seasonal_Avg_Bottom_Temps2.RData'))
rm(big.spring)

# Load other suitability (temporally static grid)
static <- read_rds(here("Habitat_Suitability/wholegrid.RDS"))
static <- static %>% 
  # dplyr::select(bathy, rugosity, gravel, 
  #               Lar.Depth, Lar.Rugos, Lar.Gravel, 
  #               grid) %>% 
  rename(geometry = grid)
st_geometry(static) <- 'geometry'
static$ID <- seq(1, nrow(static), 1)

# Load density data
dens <- st_read(here('Habitat_Suitability/Large/large_dens.shp'), quiet=T)
# Remove 2021 (no temp data)
dens <- dens %>% 
  filter(Season == 'Fall' & Year>=2008)
dens <- dens[dens$Year != 2021,]

# Rasterize BT to match grid of static
for(i in 27:length(big.fall)){
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
fall <- do.call(rbind, big.fall[27:39])

all <- fall
rownames(all) <- NULL

all <- all[!is.na(all$y),]

# Geometric and arithmetic average HSI for larges
"geometric.mean" <- 
  function(x,na.rm=TRUE){ 
    exp(mean(log(x),na.rm=na.rm))
  }

# Scale function
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

all$Lar.AA <- NA
all$Lar.GA <- NA

all <- all %>% 
  dplyr::select(ID, bathy, Lar.Depth,
                rugosity, Lar.Rugos,
                gravel, Lar.Gravel,
                bt, Lar.BT,
                Year, Season, 
                Lar.AA, Lar.GA,
                geometry) %>% 
  filter(Season == 'Fall' & Year >=2008)

for(i in 1:nrow(all)){
  print(i)
  all$Lar.AA[i] <- mean(c(
    (all$Lar.Depth[i]),
    (all$Lar.Rugos[i]),
    (all$Lar.Gravel[i]),
    (all$Lar.BT[i])))
  
  all$Lar.GA[i] <- geometric.mean(c(
    (all$Lar.Depth[i]),
    (all$Lar.Rugos[i]),
    (all$Lar.Gravel[i]),
    (all$Lar.BT[i])))

}

all$Lar.AA <- range01(all$Lar.AA, na.rm=T)
all$Lar.GA <- range01(all$Lar.GA, na.rm=T)

rm(big.fall, fall, static)

#### Getis-Ord Gi* to find suitability hotspots (with GA HSI) ####
allData.list <- split(all, f=all$Year)
# Loop through each time step
for(i in 1:length(allData.list)){
  # Keep track of season
  message(allData.list[[i]]$Season[1])
  #for(j in 1:length(allData.list[[i]])){
    # Keep track of progress
    print(allData.list[[i]]$Year[1])
    
    # Pull subset
    tes_subset <- allData.list[[i]]
    tes_subset <- tes_subset[!is.na(tes_subset$Lar.GA),]
    
    # Identify neighbors with queen contiguity (edge/vertex touching)
    tes_nb <- poly2nb(tes_subset, queen = TRUE)
    
    # Check for empty neighbor sets
    empty_nb <- which(card(tes_nb) == 0)
    
    # If empty neighbor sets exist, remove
    if(length(empty_nb > 0)){
      tes_subset <- tes_data[-empty_nb, ]
    }
    
    # Binary weighting assigns a weight of 1 to all neighboring features 
    # and a weight of 0 to all other features
    tes_w_binary <- nb2listw(tes_nb, style="B")
    
    # Calculate spatial lag of y
    tes_lag <- lag.listw(tes_w_binary, tes_subset$Lar.GA)
    
    # Identify neighbors, create weights, calculate spatial lag
    tes_nbs <- tes_subset |> 
      mutate(
        nb = st_contiguity(geometry),        # neighbors share border/vertex
        wt = st_weights(nb),                 # row-standardized weights
        tes_lag = st_lag(Lar.GA, nb, wt)          # calculate spatial lag of ga
      ) 
    
    # Calculate the Gi using local_g_perm
    tes_hot_spots <- tes_nbs |> 
      mutate(
        Gi = local_g_perm(Lar.GA, nb, wt, nsim = 999)
      ) |> 
      unnest(Gi) 
    
    # Assign gi stat and p-value (p_folded_sim) to sf dataframe
    tes_hot_spots <- tes_hot_spots |> 
      dplyr::select(gi, p_folded_sim, ID) |> 
      mutate(
        # Add a new column called "classification"
        classification = case_when(
          # Classify based on the following criteria:
          gi > 0 & p_folded_sim <= 0.01 ~ "Very suitable",
          gi > 0 & p_folded_sim <= 0.05 ~ "Suitable",
          gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat suitable",
          gi < 0 & p_folded_sim <= 0.01 ~ "Very unsuitable",
          gi < 0 & p_folded_sim <= 0.05 ~ "Unsuitable",
          gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat unsuitable",
          TRUE ~ "Insignificant"
        ),
        # Convert 'classification' into a factor for easier plotting
        classification = factor(
          classification,
          levels = c("Very suitable", "Suitable", "Somewhat suitable",
                     "Insignificant",
                     "Somewhat unsuitable", "Unsuitable", "Very unsuitable")
        )
      )
    
    # Convert to true dataframe, keep only pertinent columns
    tes_hot_spots <- sfheaders::sf_to_df(tes_hot_spots, fill=T)
    tes_hot_spots <- dplyr::select(tes_hot_spots,ID, gi, p_folded_sim, classification)
    tes_hot_spots <- unique(tes_hot_spots)                               
    
    # Merge into spatial features object by ID number
    allData.list[[i]] <- 
      merge(allData.list[[i]], tes_hot_spots, by=c('ID'))
    
    # Remove temporary features
    rm(tes_subset, tes_hot_spots, tes_nbs, tes_nb, tes_w_binary, tes_lag)
    
  }
}
# Rebind
allData2 <- do.call(rbind, allData.list)
rownames(allData2) <- NULL

# Convert classification to factor for easier plotting
allData2$classification <- factor(allData2$classification,
                                  levels = c("Very suitable", "Suitable", "Somewhat suitable",
                                             "Insignificant",
                                             "Somewhat unsuitable", "Unsuitable", "Very unsuitable"))

# Clean and bind all data
rm(all, allData.list, empty_nb)

sten <- allData2 %>% 
  dplyr::select(ID, Year, bathy, rugosity, gravel, bt, Lar.Depth, Lar.Rugos, Lar.Gravel,
                Lar.BT, Lar.GA, classification, geometry) %>% 
  mutate(Suit.Class = classification)
sten$Suit[sten$Suit.Class %in% c('Insignificant')] <- 'Insignificant'
sten$Suit[sten$Suit.Class %in% c('Unsuitable', 'Very unsuitable', 'Somewhat unsuitable')] <- 'Unsuitable'
sten$Suit[sten$Suit.Class %in% c('Suitable', 'Very suitable', 'Somewhat suitable')] <- 'Suitable'

sten$Suit <- factor(sten$Suit, levels = c('Suitable', 'Insignificant', 'Unsuitable'))

dens <- dens %>% 
  dplyr::select(ID, Year, y, logy, clssfct, Class) %>% 
  mutate(Dens.Class = clssfct,
         Dens = Class)

dens$Dens <- factor(dens$Dens, levels = c('Hot', 'Insignificant', 'Cold'))

dens <- sfheaders::sf_to_df(dens, fill=T)

all <- left_join(dplyr::select(sten, -classification, -Suit),
                 unique(dplyr::select(dens, -clssfct, -Class, -x, -sfg_id, -polygon_id, -linestring_id, -y..1)),
                 by = c('ID', 'Year'))

all <- all %>% 
  rename(density = y,
         log.density = logy)

all$Suit[all$Suit.Class %in% c('Insignificant')] <- 'Insignificant'
all$Suit[all$Suit.Class %in% c('Very suitable', 'Suitable', 'Somewhat suitable')] <- 'Suitable'
all$Suit[all$Suit.Class %in% c('Very unsuitable', 'Unsuitable', 'Somewhat unsuitable')] <- 'Unsuitable'

all$Suit <- factor(all$Suit, levels = c('Suitable', 'Insignificant', 'Unsuitable'))

rm(allData2, dens, sten, i)

# Characterize length of time grid cells spend in each condition
test <- all %>% 
  group_by(Dens, Suit) %>% 
  mutate(n_overlaps = lengths(st_within(geometry))) %>% 
  dplyr::select(ID, Year, Dens, Suit, n_overlaps) %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(ID, Year, Dens, Suit, n_overlaps) %>% 
  unique() %>% as.data.frame()

all <- left_join(all, test, by=c('ID', 'Year', 'Dens', 'Suit'))
rm(test)

# Plot duration of time in each state
ggplot() +
  geom_sf(data = all, 
          col=NA, aes(fill=n_overlaps)) +
  scale_fill_viridis_c(limits = c(0,13)) +
  geom_sf(data = ecodata::coast) +
  ggh4x::facet_grid2(Suit ~ Dens) +
  coord_sf(xlim=c(-76, -66),
           ylim=c(35, 45), crs="EPSG:4326") +
  labs(fill = '# Years')

# Visualize density-suitability overlap for selected years
closed <- st_read(here('Data/GIS/closed_areas_wgs.shp'))
closed <- st_transform(closed, st_crs(all))

ggplot() +
  geom_sf(data=all[all$Dens == 'Hot' & all$Year %in% seq(2009, 2020, 1),],
          fill = "#F8766D", color = NA, alpha=0.5, show.legend = T) +
  geom_sf(data = all[all$Suit == 'Suitable' & all$Year %in% seq(2009, 2020, 1),],
          fill = "#00BFC4", color = NA, alpha = 0.5, show.legend = T) +
  geom_sf(data = st_transform(ecodata::coast, st_crs(all)),
          col='gray70', fill='gray90') +
  geom_sf(data = closed, fill=NA, col='black') +
  coord_sf(xlim=c(-76, -66),
           ylim=c(40.5, 45), crs="EPSG:4326") +
  facet_wrap(vars(Year), nrow = 3) +
  ggtitle('Spatial characteristics for large cod in Fall') +
  labs(caption = paste0('Blue cells area areas identified as having significantly',
                        ' high suitability for large cod (based on bottom temperature,',
                        ' rugosity, depth, and gravel content).\nRed cells are areas',
                        ' identified as having significantly high large cod density. Black',
                        ' outlines denote areas with year-round groundfishing closures.'))

test <- all %>% 
  group_by(Year, Dens, Suit) %>% 
  summarise(geometry = st_union(geometry),
            bathy = mean(bathy),
            rugosity = mean(rugosity),
            gravel = mean(gravel),
            bt = mean(bt),
            
            lar.depth = mean(Lar.Depth),
            lar.rugos = mean(Lar.Rugos),
            lar.gravel = mean(Lar.Gravel),
            lar.bt = mean(Lar.BT),
            
            lar.ga = mean(Lar.GA),
            
            density = mean(density),
            log.density = mean(log.density))
cent <- data.frame(
  Year = test$Year,
  Dens = test$Dens,
  Suit = test$Suit
)
for(i in 1:nrow(cent)){
  cent$centroid[i] <- st_centroid(test[i,])$geometry[1]
}
st_geometry(cent) <- 'centroid'

ggplot() +
  geom_sf(data = cent, 
          aes(col = Year)) +
  scale_color_viridis_c() +
  ggh4x::facet_grid2(Suit ~ Dens)

cent <- cent %>% 
  sfheaders::sf_to_df(fill=T)

ggplot() +
  geom_line(data=cent, aes(x=Year, y=x/ 1000)) +
  ggh4x::facet_grid2(Suit ~ Dens)

test$comb <- paste0(test$Dens, ' ', test$Suit)
plodat <- test[test$comb %in% c('Hot Suitable', 'Hot Unsuitable', 
                                'Cold Suitable', 'Cold Unsuitable'),]
plodat$area <- st_area(plodat)
plodat <- plodat %>% 
  sfheaders::sf_to_df(fill=T) %>% 
  dplyr::select(-comb, -sfg_id, -multipolygon_id, -polygon_id, -linestring_id) %>% 
  mutate(area = FishStatsUtils::strip_units(area)) %>% 
  pivot_longer(cols = seq(4, 15, 1))

plodat$name <- factor(plodat$name, 
                      levels=c('bathy', 'gravel', 'rugosity', 'bt',
                               'lar.depth', 'lar.gravel', 'lar.rugos', 'lar.bt',
                               'density', 'area', 'lar.ga'))

ggplot(data=plodat[plodat$name %in% c('bathy', 'gravel', 'rugosity', 'bt',
                                      'density', 'area', 'lar.ga'),]) +
  geom_point(aes(x=Year, y=value, col=Dens, group = Dens),
             cex=0.7, alpha=0.7) +
  geom_line(aes(x=Year, y=value, col=Dens, group = Dens), lwd=1, alpha=0.8) +
  ggh4x::facet_grid2(name ~ Suit, scales = 'free_y') +
  scale_x_continuous(breaks = seq(2008, 2020, 2))

ggplot(data=all[all$Dens == 'Hot',]) +
  geom_sf(aes(fill=log.density), col=NA) +
  scale_fill_viridis_c() +
  geom_sf(data=st_transform(ecodata::coast, st_crs(all))) +
  geom_sf(data=closed, fill=NA, col='black') +
  facet_wrap(vars(Year)) +
  coord_sf(xlim=c(-71, -65),
           ylim=c(41, 44), crs="EPSG:4326") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
