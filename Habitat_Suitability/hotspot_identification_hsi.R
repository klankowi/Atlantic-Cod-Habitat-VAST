# Getis-Ord Gi* analysis of HSI hotspots

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

# Load HSI values
hsi <- st_read(here('Habitat_Suitability/large/large_hsi.shp'),
               quiet=T)
hsi <- hsi %>% 
  filter(!is.na(Season)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

# Split into seasons then years
allData.list <- split(hsi, f=hsi$Season)
for(i in 1:length(allData.list)){
  allData.list[[i]] <- split(allData.list[[i]], f=allData.list[[i]]$Year)
}

#### Getis-Ord Gi* ####
# Loop through each time step
for(i in 1:length(allData.list)){
  # Keep track of season
  message(allData.list[[i]][[1]]$Season[1])
  for(j in 1:length(allData.list[[i]])){
    # Keep track of progress
    print(allData.list[[i]][[j]]$Year[1])
    
    # Pull subset
    tes_subset <- allData.list[[i]][[j]]
    tes_subset <- tes_subset[!is.na(tes_subset$ga),]
    
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
    tes_lag <- lag.listw(tes_w_binary, tes_subset$ga)
    
    # Identify neighbors, create weights, calculate spatial lag
    tes_nbs <- tes_subset |> 
      mutate(
        nb = st_contiguity(geometry),        # neighbors share border/vertex
        wt = st_weights(nb),                 # row-standardized weights
        tes_lag = st_lag(ga, nb, wt)          # calculate spatial lag of ga
      ) 
    
    # Calculate the Gi using local_g_perm
    tes_hot_spots <- tes_nbs |> 
      mutate(
        Gi = local_g_perm(ga, nb, wt, nsim = 999)
      ) |> 
      unnest(Gi) 
    
    # Assign gi stat and p-value (p_folded_sim) to sf dataframe
    tes_hot_spots <- tes_hot_spots |> 
      dplyr::select(gi, p_folded_sim, ID) |> 
      mutate(
        # Add a new column called "classification"
        classification = case_when(
          # Classify based on the following criteria:
          gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
          gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
          gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
          gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
          gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
          gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
          TRUE ~ "Insignificant"
        ),
        # Convert 'classification' into a factor for easier plotting
        classification = factor(
          classification,
          levels = c("Very hot", "Hot", "Somewhat hot",
                     "Insignificant",
                     "Somewhat cold", "Cold", "Very cold")
        )
      )
    
    # Convert to true dataframe, keep only pertinent columns
    tes_hot_spots <- sfheaders::sf_to_df(tes_hot_spots, fill=T)
    tes_hot_spots <- dplyr::select(tes_hot_spots,ID, gi, p_folded_sim, classification)
    tes_hot_spots <- unique(tes_hot_spots)                               
    
    # Merge into spatial features object by ID number
    allData.list[[i]][[j]] <- 
      merge(allData.list[[i]][[j]], tes_hot_spots, by=c('ID'))
    
    # Remove temporary features
    rm(tes_subset, tes_hot_spots, tes_nbs, tes_nb, tes_w_binary, tes_lag)
    
  }
}

# Rebind into true dataframe
allData2 <- allData.list
for(i in 1:length(allData2)){
  allData2[[i]] <- do.call(rbind, allData2[[i]])
}
allData2 <- do.call(rbind, allData2)
rownames(allData2) <- NULL

# Convert classification to factor for easier plotting
allData2$classification <- factor(allData2$classification,
                                  levels = c("Very hot", "Hot", "Somewhat hot",
                                             "Insignificant",
                                             "Somewhat cold", "Cold", "Very cold"))
# Plot and save plots
spring <- ggplot(data=allData2[allData2$Season == 'Spring',]) +
  geom_sf(aes(fill=classification), col=NA) +
  scale_fill_brewer(type = "div", palette = 5) +
  labs(
    fill = "Hot Spot\nClassification"
  ) +
  
  geom_sf(data=coast, fill='gray', col='gray', alpha=0.7) +
  
  coord_sf(xlim=c(-115489.4, 770308.2),
           ylim=c(4060391.8, 4966321.1)) +
  
  facet_wrap(vars(Year)) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=8, angle=35),
        axis.text.y = element_text(size=8))

fall <- ggplot(data=allData2[allData2$Season == 'Fall',]) +
  geom_sf(aes(fill=classification), col=NA) +
  scale_fill_brewer(type = "div", palette = 5) +
  labs(
    fill = "Hot Spot\nClassification"
  ) +
  
  geom_sf(data=coast, fill='gray', col='gray', alpha=0.7) +
  
  coord_sf(xlim=c(-115489.4, 770308.2),
           ylim=c(4060391.8, 4966321.1)) +
  
  facet_wrap(vars(Year)) +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=8, angle=35),
        axis.text.y = element_text(size=8))

ggsave(spring, 
       filename=here('Habitat_Suitability/large/spring_hsi_hotspots.pdf'),
       width = 8.5, height = 11, dpi=300)

ggsave(fall, 
       filename=here('Habitat_Suitability/large/fall_hsi_hotspots.pdf'),
       width = 8.5, height = 11, dpi=300)

st_write(allData2,
         here('Habitat_Suitability/Large/gridded_hsi_suitability.shp'))

# Check temporal trends for each cell
# Larger class groupings (hot vs cold vs insignificant)
allData2$Class[allData2$classification %in% c('Very cold', 'Cold', 
                                              'Somewhat cold'
                                              )] <- 'Cold'
allData2$Class[allData2$classification %in% c('Very hot', 'Hot',
                                              'Somewhat hot'
                                              )] <- 'Hot'
allData2$Class[allData2$classification %in% c('Insignificant'#,
                                              #'Somewhat hot', 
                                              #'Somewhat cold'
                                              )] <- 'Insignificant'

allData2$Class <- factor(allData2$Class, levels=c('Hot', 'Insignificant', 'Cold'))

# Make table of class grouping in spring
springtab <- as.data.frame(table(allData2$ID[allData2$Season == 'Spring'],
                                 allData2$Class[allData2$Season == 'Spring']))
springtab <- springtab %>% 
  rename(ID = Var1,
         Class = Var2) %>% 
  filter(Freq != 0)

springtab <- springtab[with(springtab, order(ID, Class)),]
rownames(springtab) <- NULL

# Find persistence (90% of time series)
persistentcold <- springtab[springtab$Freq >=36 &
                              springtab$Class == 'Cold',]

persistenthot <- springtab[springtab$Freq >=36 &
                             springtab$Class == 'Hot',]

springp <- ggplot() +
  geom_sf(data=allData2[allData2$Year == 1990 & allData2$Season == 'Spring' &
                          allData2$ID %in% persistentcold$ID,],
          fill='blue', col=NA, alpha=0.6) +
  geom_sf(data=allData2[allData2$Year == 1990 & allData2$Season == 'Spring' &
                          allData2$ID %in% persistenthot$ID,],
          fill='red', col=NA, alpha=0.6) +
  geom_sf(data=coast,fill='gray60', col='gray50') +
  geom_sf(data=regions, lwd=0.15, fill=NA, col='gray20') +
  coord_sf(xlim=c(st_bbox(regions)[1], st_bbox(regions)[3]),
           ylim=c(st_bbox(regions)[2], st_bbox(regions)[4])) +
  ggtitle('Spring 90% persistence')

# Save outcomes
hot.spring <- allData2 %>% 
  filter(Year == 1990 &
           Season == 'Spring' &
           ID %in% persistenthot$ID)

cold.spring <- allData2 %>% 
  filter(Year == 1990 &
           Season == 'Spring' &
           ID %in% persistentcold$ID)
# Find neighboring cells
hotneighbors.s <-  poly2nb(hot.spring, queen = TRUE)
colneighbors.s <-  poly2nb(cold.spring, queen=TRUE)
# Find names of "groups" or "neighborhoods"
nc.h <- attr(hotneighbors.s, "ncomp")$comp.id
nc.c <- attr(colneighbors.s, "ncomp")$comp.id
# Identify neighborhood each cell belongs to
hot.spring$neighborhood <- nc.h
cold.spring$neighborhood <- nc.c
# Melt into neighborhoods (must be touching sides or vertices)
hot.spring <- hot.spring %>% 
  group_by(neighborhood) %>% 
  summarise(ga = mean(ga)) %>% 
  mutate(Class = 'Hot',
         Season = 'Spring')
cold.spring <- cold.spring %>% 
  group_by(neighborhood) %>% 
  summarise(ga = mean(ga)) %>% 
  mutate(Class = 'Cold',
         Season = 'Spring')

spring <- rbind(hot.spring, cold.spring)

## Fall
falltab <- as.data.frame(table(allData2$ID[allData2$Season == 'Fall'],
                               allData2$Class[allData2$Season == 'Fall']))
falltab <- falltab %>% 
  rename(ID = Var1,
         Class = Var2) %>% 
  filter(Freq != 0)

falltab <- falltab[with(falltab, order(ID, Class)),]
rownames(falltab) <- NULL

persistentcold <- falltab[falltab$Freq >=36 &
                            falltab$Class == 'Cold',]

persistenthot <- falltab[falltab$Freq >=36 &
                           falltab$Class == 'Hot',]

fallp <- ggplot() +
  geom_sf(data=allData2[allData2$Year == 1990 & allData2$Season == 'Fall' &
                          allData2$ID %in% persistentcold$ID,],
          fill='blue', col=NA, alpha=0.6) +
  geom_sf(data=allData2[allData2$Year == 1990 & allData2$Season == 'Fall' &
                          allData2$ID %in% persistenthot$ID,],
          fill='red', col=NA, alpha=0.6) +
  geom_sf(data=coast,fill='gray60', col='gray50') +
  geom_sf(data=regions, lwd=0.15, fill=NA, col='gray20') +
  coord_sf(xlim=c(st_bbox(regions)[1], st_bbox(regions)[3]),
           ylim=c(st_bbox(regions)[2], st_bbox(regions)[4])) +
  ggtitle('Fall 90% persistence')

# Save outcomes
hot.fall <- allData2 %>% 
  filter(Year == 1990 &
           Season == 'Fall' &
           ID %in% persistenthot$ID)

cold.fall <- allData2 %>% 
  filter(Year == 1990 &
           Season == 'Fall' &
           ID %in% persistentcold$ID)
# Find neighboring cells
hotneighbors.s <-  poly2nb(hot.fall, queen = TRUE)
colneighbors.s <-  poly2nb(cold.fall, queen=TRUE)
# Find names of "groups" or "neighborhoods"
nc.h <- attr(hotneighbors.s, "ncomp")$comp.id
nc.c <- attr(colneighbors.s, "ncomp")$comp.id
# Identify neighborhood each cell belongs to
hot.fall$neighborhood <- nc.h
cold.fall$neighborhood <- nc.c
# Melt into neighborhoods (must be touching sides or vertices)
hot.fall <- hot.fall %>% 
  group_by(neighborhood) %>% 
  summarise(ga = mean(ga)) %>% 
  mutate(Class = 'Hot',
         Season = 'Fall')
cold.fall <- cold.fall %>% 
  group_by(neighborhood) %>% 
  summarise(ga = mean(ga)) %>% 
  mutate(Class = 'Cold',
         Season = 'Fall')

fall <- rbind(hot.fall, cold.fall)

ggsave(springp,
       filename = here('VAST_runs/large/Overall_BC/spring_hsi_persistence.pdf'),
       width = 8.5, height = 11, dpi=300)

ggsave(fallp,
       filename = here('VAST_runs/large/Overall_BC/fall_hsi_persistence.pdf'),
       width = 8.5, height = 11, dpi=300)

persistence <- rbind(spring, fall)
persistence <- persistence %>% 
  dplyr::select(Class, Season, ga, neighborhood, geometry) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')),
         Class = factor(Class, levels=c('Hot', 'Cold')))

st_write(persistence, 
         here('Habitat_Suitability/large_hsi_persistence.shp'),
         append=F)

rm(list=setdiff(ls(), c('coast', 'persistence', 'regions',
                        'allData2')))

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Annual hotspot pull
# Break into timesteps
allData.list <- split(allData2, f=allData2$Season)
for(i in 1:length(allData.list)){
  allData.list[[i]] <- split(allData.list[[i]], f=allData.list[[i]]$Year)
}

# Loop through time steps
for(i in 1:length(allData.list)){
  # Keep track of season
  message(allData.list[[i]][[1]]$Season[1])
  for(j in 1:length(allData.list[[i]])){
    # Keep track of year
    print(allData.list[[i]][[j]]$Year[1])
    
    # Extract year you're working on
    temp <- allData.list[[i]][[j]]
    
    # Extract grids with p<0.05 (sig hot)
    hot <- temp[temp$classification %in% c('Hot', 'Very hot'),]
    # Find neighboring cells
    hotneighbors <-  poly2nb(hot, queen = TRUE)
    # Find names of "groups" or "neighborhoods"
    nc <- attr(hotneighbors, "ncomp")$comp.id
    # Identify neighborhood each cell belongs to
    hot$neighborhood <- nc
    # Melt into neighborhoods (must be touching sides or vertices)
    hot <- hot %>% 
      group_by(neighborhood) %>% 
      summarise(ga = mean(ga)) %>% 
      mutate(Class = 'Hot',
             Year = temp$Year[1],
             Season = temp$Season[1])
    hottab <- as.data.frame(table(nc))
    hottab <- hottab[hottab$Freq == 1,]
    hot <- hot %>% 
      filter(neighborhood %notin% hottab$nc)
    
    # Same thing for cold (p>0.95) cells
    cold <- temp[temp$classification %in% c('Cold', 'Very cold'),]
    coldneighbors <-  poly2nb(cold, queen = TRUE)
    nc <- attr(coldneighbors, "ncomp")$comp.id
    cold$neighborhood <- nc
    cold <- cold %>% 
      group_by(neighborhood) %>% 
      summarise(ga = mean(ga)) %>% 
      mutate(Class = 'Cold',
             Year = temp$Year[1],
             Season = temp$Season[1])
    coltab <- as.data.frame(table(nc))
    coltab <- coltab[coltab$Freq == 1,]
    cold <- cold %>% 
      filter(neighborhood %notin% coltab$nc)
    
    # Bind into single dataframe
    sig <- rbind(hot, cold)
    
    # Paste into list
    allData.list[[i]][[j]] <- sig
    
    # Remove temporary features
    rm(sig, cold, hot, nc, coldneighbors, hotneighbors, temp,
       hottab, coltab)
  }
}

# Rebind into true dataframe
nbd <- allData.list
for(i in 1:length(nbd)){
  nbd[[i]] <- do.call(rbind, nbd[[i]])
}
nbd <- do.call(rbind, nbd)
rownames(nbd) <- NULL

nbd <- nbd %>% 
  dplyr::select(Year, Season, Class, neighborhood, ga, geometry) %>% 
  mutate(Class = factor(Class, levels=c('Hot', 'Cold')))

ggplot(data=nbd[nbd$Year%in% seq(2012, 2021,1),]) +
  geom_sf(aes(fill=Class)) +
  geom_sf(data=regions, fill=NA) +
  ggh4x::facet_grid2(c('Season', 'Year')) +
  ggtitle('Large cod, regime shift')

## Large
# 1982 - 1990 = Failure of WGOM, gain in GBK
# 1991 - 2000 = Slight recovery in WGOM, rapid contraction and stability in GBK
# 2001 - 2012 = Stability in all areas
# 2013 - 2021 = Regime shift

# Area
nbd$area <- st_area(nbd)

annualarea <- nbd %>% 
  group_by(Year, Season, Class) %>% 
  summarise(area = sum(area),
            ga = mean(ga))

ggplot(data=annualarea) + 
  geom_line(aes(x=Year, y=ga, col=Class)) +
  geom_smooth(aes(x=Year, y=ga, col=Class)) +
  geom_vline(xintercept = c(2008,2012), lty=2) +
  facet_wrap(vars(Season))

closed <- st_read(here('Data/GIS/closed_areas_wgs.shp'), quiet = T)
closed <- st_transform(closed, st_crs(nbd))

for(i in 1:10#length(unique(nbd$Year))
){
  year <- unique(nbd$Year)[i]
  use <- nbd[nbd$Year == year,]
  
  print(
    ggplot(data=use) +
      geom_sf(aes(fill=Class)) +
      #geom_sf_text(aes(label=neighborhood), size=10) +
      geom_sf(data=regions, fill=NA) +
      geom_sf(data=closed, fill='gray10', col=NA, alpha=0.2) +
      facet_wrap(vars(Season)) +
      labs(x='', y='') +
      ggtitle(paste0(year)) +
      theme(legend.position = 'none',
            strip.text.x = element_text(size=12))
    
  )
}

nbd$area <- nbd$area / 1000000
nbd$area <- strip_units(nbd$area)
nbd$area <- round(nbd$area, 3)

st_write(nbd, 
         here('Habitat_Suitability/Large/Large_HSI_Neighborhoods.shp'),
         append = F)
