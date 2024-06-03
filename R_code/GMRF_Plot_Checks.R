rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)
library(wesanderson)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "inside",
                legend.background = element_rect(fill='transparent', colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Shapefiles and centroid calculation
gom <- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
gom <- st_make_valid(gom)
coast <- st_transform(ecodata::coast, crs=st_crs(gom))
centroid <- st_centroid(gom)

#### Center of gravity ####
# Root
root <- here('VAST_runs/large/GMRF')
ending <- '/COG_ALL.csv'

# Covariates
covuse <- c('ON', 'OFF')

# Baseline Mods
on <- read.csv(paste0(root, '/ON', ending))
off <- read.csv(paste0(root, '/OFF', ending))
on$Model <- 'GMRF ON'
off$Model <- 'GMRF OFF'

rm(list=setdiff(ls(), c('on', 'off', 'gom', 'coast',
                        'centroid')))

cog <- rbind(on, off)

rm(list=setdiff(ls(), c('cog', 'gom', 'centroid', 'coast')))

cog <- cog %>%
  mutate(easting.km = easting / 1000,
         northing.km = northing / 1000,
         e.sd.km = e.sd / 1000,
         n.sd.km = n.sd / 1000) %>%
  dplyr::select(-units)


cog$Model <- factor(cog$Model,
                    levels=c('GMRF ON', 'GMRF OFF'))

ccord <- sfheaders::sf_to_df(centroid, fill=T)
ccord <- ccord %>% 
  dplyr::select(x, y) %>% 
  rename(lon = x,
         lat = y) %>% 
  mutate(easting.km = lon/1000,
         northing.km = lat/1000)

cog$center.lon <- cog$easting.km - ccord$easting.km
cog$center.lat <- cog$northing.km - ccord$northing.km

coglist <- split(cog, f=cog$Season)
for(i in 1:length(coglist)){
  coglist[[i]] <- split(coglist[[i]], f=coglist[[i]]$Model)
}

for(i in 1:length(coglist)){
  coglist[[i]][[1]]$easting.resid <- 0
  coglist[[i]][[1]]$northing.resid <- 0
  
  for(j in 2:length(coglist[[i]])){
    coglist[[i]][[j]]$easting.resid <- (coglist[[i]][[j]]$easting.km - coglist[[i]][[1]]$easting.km)^2
    coglist[[i]][[j]]$northing.resid <- (coglist[[i]][[j]]$northing.km - coglist[[i]][[1]]$northing.km)^2
  }
}

for(i in 1:length(coglist)){
  for(j in 1:length(coglist[[i]])){
    coglist[[i]][[j]]$easting.SSE <- sum(coglist[[i]][[j]]$easting.resid)
    coglist[[i]][[j]]$northing.SSE <- sum(coglist[[i]][[j]]$northing.resid)
  }
}
cog2 <- do.call(rbind, coglist)
cog2 <- do.call(rbind, cog2)

cog2 <- dplyr::select(cog2, Model, Season, easting.SSE, northing.SSE)
cog2 <- unique(cog2)

# Look at SSE similarity
northing.similar <- cog2[with(cog2, order(northing.SSE, decreasing = F)),]
northing.similar.spring <- northing.similar[northing.similar$Season == 'Spring',]
rownames(northing.similar.spring) <- NULL
northing.similar.spring
northing.similar.fall <- northing.similar[northing.similar$Season == 'Fall',]
rownames(northing.similar.fall) <- NULL
northing.similar.fall

easting.similar <- cog2[with(cog2, order(easting.SSE, decreasing = F)),]
easting.similar.spring <- easting.similar[easting.similar$Season == 'Spring',]
rownames(easting.similar.spring) <- NULL
easting.similar.spring
easting.similar.fall <- easting.similar[easting.similar$Season == 'Fall',]
rownames(easting.similar.fall) <- NULL
easting.similar.fall

n.spring <- ggplot() +
  geom_point(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lat, col=Model),
            cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lat, col=Model),
            lwd=1, alpha=0.7) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('Spring center of gravity northing - GMRF models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

n.fall <- ggplot() +
  geom_point(data=cog[cog$Season == 'Fall',], 
             aes(x=Year, y=center.lat, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Fall',], 
            aes(x=Year, y=center.lat, col=Model),
            lwd=1, alpha=0.7) +
  #scale_color_manual(values=c('black',
  #                            wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from northing centroid (km)') +
  ggtitle('Fall center of gravity northing - GMRF models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

e.spring <- ggplot() +
  geom_point(data=cog[cog$Season == 'Spring',], 
             aes(x=Year, y=center.lon, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Spring',], 
            aes(x=Year, y=center.lon, col=Model),
            lwd=1, alpha=0.7) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from easting centroid (km)') +
  ggtitle('Spring center of gravity easting - GMRF models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

e.fall <- ggplot() +
  geom_point(data=cog[cog$Season == 'Fall',], 
             aes(x=Year, y=center.lon, col=Model),
             cex=1.5, alpha=0.7) +
  geom_line(data=cog[cog$Season == 'Fall',], 
            aes(x=Year, y=center.lon, col=Model),
            lwd=1, alpha=0.7) +
  # scale_color_manual(values=c('black',
  #                             wes_palette('Darjeeling1', n=5))) +
  labs(x='Year', y='Distance from easting centroid (km)') +
  ggtitle('Fall center of gravity easting - GMRF models') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2))

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/cogshifts_northing_spring.png')),
       n.spring, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/cogshifts_northing_fall.png')),
       n.fall, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/cogshifts_easting_spring.png')),
       e.spring, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/cogshifts_easting_fall.png')),
       e.fall, width = 8.5, height=4.25, units='in')

rm(list=setdiff(ls(), c('cog', 'coast', 'gom')))

# Spatial plots
cog.sf <- st_as_sf(cog, coords=c('easting', 'northing'), crs="EPSG:32619")

cog.lin <- cog.sf %>% 
  filter(Model=='All Covars') %>% 
  arrange(Year) %>% 
  #filter(Year != 2022) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast('LINESTRING')
# Spring
SD_plotting.spring <- subset(cog.sf, Season =='Spring')

SD_list <- split(SD_plotting.spring, f=SD_plotting.spring$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.spring <- st_multilinestring(do.call("rbind", linestrings))
  t.spring <-  nngeo::st_segments(t.spring)
  t.spring <- st_sf(t.spring)
  t.spring$Year <- seq(1982, 2020, 1)
  st_crs(t.spring) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.spring
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.spring, mod)
}
yeartracks <- do.call(rbind, SD_list)

# Plot
spring.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks, aes(col=Year), pch=19, cex=0.5, lwd=0.7, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-72, -66),
           ylim=c(41, 43.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Spring Center of Gravity') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=8, angle=20),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_COG_Spring.png"),
       spring.vis,
       width=10, height=10)

# Fall
SD_plotting.fall <- subset(cog.sf, Season =='Fall')

SD_list <- split(SD_plotting.fall, f=SD_plotting.fall$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.fall <- st_multilinestring(do.call("rbind", linestrings))
  t.fall <-  nngeo::st_segments(t.fall)
  t.fall <- st_sf(t.fall)
  t.fall$Year <- seq(1982, 2020, 1)
  st_crs(t.fall) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.fall
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.fall, mod)
}
yeartracks <- do.call(rbind, SD_list)

# Plot
fall.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks, aes(col=Year), pch=19, cex=0.5, lwd=0.7, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-72, -66),
           ylim=c(41, 43.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Fall Center of Gravity') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=8, angle=20),
        axis.text.y=element_text(size=8),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_COG_Fall.png"),
       fall.vis,
       width=10, height=10)

rm(list=setdiff(ls(), c('coast', 'gom')))

#### Range edges ####
# Root
root <- here('VAST_runs/large/GMRF')
ending <- '/RangeEdges_ALL.csv'

# Baseline Mods
on <- read.csv(paste0(root, '/ON/', ending))
off <- read.csv(paste0(root, '/OFF/', ending))
on$Model <- 'GMRF ON'
off$Model <- 'GMRF OFF'

rm(list=setdiff(ls(), c('on', 'off',
                        'gom', 'coast', 'centroid')))

re <- rbind(on, off)

rm(list=setdiff(ls(), c('re', 'gom', 'centroid', 'coast')))

re$Model <- factor(re$Model,
                    levels=c('GMRF ON', 'GMRF OFF'))
re$Quantile <- as.factor(re$Quantile)

# Spaghetti plots
northing.spring <- ggplot() +
  geom_point(data=re[re$Season == 'Spring',], 
             aes(x=Year, y=Northing.Est, col=Model, pch=Quantile),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Season == 'Spring',], 
            aes(x=Year, y=Northing.Est, col=Model, lty=Quantile),
            lwd=1, alpha=0.7) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Spring range edge Northing') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2)) 

easting.spring <- ggplot() +
  geom_point(data=re[re$Season == 'Spring',], 
             aes(x=Year, y=Easting.Est, col=Model, pch=Quantile),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Season == 'Spring',], 
            aes(x=Year, y=Easting.Est, col=Model, lty=Quantile),
            lwd=1, alpha=0.7) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Spring range edge Easting') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2)) 

northing.fall <- ggplot() +
  geom_point(data=re[re$Season == 'Fall',], 
             aes(x=Year, y=Northing.Est, col=Model, pch=Quantile),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Season == 'Fall',], 
            aes(x=Year, y=Northing.Est, col=Model, lty=Quantile),
            lwd=1, alpha=0.7) +
  labs(x='Year', y='Northing (km)') +
  ggtitle('Fall range edge Northing') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2)) 

easting.fall <- ggplot() +
  geom_point(data=re[re$Season == 'Fall',], 
             aes(x=Year, y=Easting.Est, col=Model, pch=Quantile),
             cex=1.5, alpha=0.7) +
  geom_line(data=re[re$Season == 'Fall',], 
            aes(x=Year, y=Easting.Est, col=Model, lty=Quantile),
            lwd=1, alpha=0.7) +
  labs(x='Year', y='Easting (km)') +
  ggtitle('Fall range edge Easting') +
  theme(legend.position = 'bottom',
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  guides(color = guide_legend(nrow = 2)) 

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/re_northing_spring.png')),
       northing.spring, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/re_northing_fall.png')),
       northing.fall, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/re_easting_spring.png')),
       easting.spring, width = 8.5, height=4.25, units='in')

ggsave(filename=paste0(here('VAST_runs/large/GMRF/Testing_Results', 
                            '/re_easting_fall.png')),
       easting.fall, width = 8.5, height=4.25, units='in')

rm(list=setdiff(ls(), c('coast', 'gom', 're')))

# Spatial plots
re <- re %>% 
  mutate(Northing.Est=Northing.Est*1000,
         Northing.SD=Northing.SD*1000,
         Easting.Est=Easting.Est*1000,
         Easting.SD=Easting.SD*1000)

re.sf <- st_as_sf(re, coords=c('Easting.Est', 'Northing.Est'), 
                  crs="EPSG:32619")

#### 5% Quantile ####

# Spring
SD_plotting.spring <- subset(re.sf, Season =='Spring' & Quantile=='0.05')

SD_list <- split(SD_plotting.spring, f=SD_plotting.spring$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.spring <- st_multilinestring(do.call("rbind", linestrings))
  t.spring <-  nngeo::st_segments(t.spring)
  t.spring <- st_sf(t.spring)
  t.spring$Year <- seq(1982, 2020, 1)
  st_crs(t.spring) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.spring
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.spring, mod)
}
yeartracks.spring <- do.call(rbind, SD_list)

# Plot
spring.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks.spring, 
          aes(col=Year), pch=19, cex=0.5, lwd=0.5, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-72, -69),
           ylim=c(40.5, 42.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Spring Southwestern Range Edge (5% Quantile)') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=6, angle=20),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_RangeEdge_SW_Spring.png"),
       spring.vis,
       width=8, height=6)
rm(SD_list, SD_plotting.spring)

# Fall
SD_plotting.fall <- subset(re.sf, Season =='Fall' & Quantile == '0.05')

SD_list <- split(SD_plotting.fall, f=SD_plotting.fall$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.fall <- st_multilinestring(do.call("rbind", linestrings))
  t.fall <-  nngeo::st_segments(t.fall)
  t.fall <- st_sf(t.fall)
  t.fall$Year <- seq(1982, 2020, 1)
  st_crs(t.fall) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.fall
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.fall, mod)
}
yeartracks.fall <- do.call(rbind, SD_list)
rm(SD_list, SD_plotting.fall)

# Plot
fall.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks.fall, 
          aes(col=Year), pch=19, cex=0.5, lwd=0.5, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-72, -69),
           ylim=c(40.5, 42.5),
          crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Fall Southwestern Range Edge (5% Quantile)') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=6, angle=20),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_RangeEdge_SW_Fall.png"),
       fall.vis,
       width=8, height=6)

rm(fall.vis, spring.vis, yeartracks.fall, yeartracks.spring)

#### 95% Quantile ####
# Spring
SD_plotting.spring <- subset(re.sf, Season =='Spring' & Quantile=='0.95')

SD_list <- split(SD_plotting.spring, f=SD_plotting.spring$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.spring <- st_multilinestring(do.call("rbind", linestrings))
  t.spring <-  nngeo::st_segments(t.spring)
  t.spring <- st_sf(t.spring)
  t.spring$Year <- seq(1982, 2020, 1)
  st_crs(t.spring) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.spring
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.spring, mod)
}
yeartracks.spring <- do.call(rbind, SD_list)

# Plot
spring.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks.spring, 
          aes(col=Year), pch=19, cex=0.5, lwd=0.5, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-70, -65),
           ylim=c(42, 45.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Spring Northeastern Range Edge (95% Quantile)') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=6, angle=20),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_RangeEdge_NE_Spring.png"),
       spring.vis,
       width=8, height=6)
rm(SD_list, SD_plotting.spring)

# Fall
SD_plotting.fall <- subset(re.sf, Season =='Fall' & Quantile == '0.95')

SD_list <- split(SD_plotting.fall, f=SD_plotting.fall$Model)

for(i in 1:length(SD_list)){
  points <- st_cast(st_geometry(SD_list[[i]]), "POINT") 
  # Number of total linestrings to be created
  n <- length(points) - 1
  # Build linestrings
  linestrings <- lapply(X = 1:n, FUN = function(x) {
    
    pair <- st_combine(c(points[x], points[x + 1]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  # Split to individual linestrings, associate year
  t.fall <- st_multilinestring(do.call("rbind", linestrings))
  t.fall <-  nngeo::st_segments(t.fall)
  t.fall <- st_sf(t.fall)
  t.fall$Year <- seq(1982, 2020, 1)
  st_crs(t.fall) <- "EPSG:32619"
  
  mod <- SD_list[[i]]$Model[1]
  
  SD_list[[i]] <- t.fall
  SD_list[[i]]$Model <- mod
  rm(points, n, linestrings, t.fall, mod)
}
yeartracks.fall <- do.call(rbind, SD_list)
rm(SD_list, SD_plotting.fall)

# Plot
fall.vis <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=gom, fill='transparent', lwd=0.25) +
  guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
  ggnewscale::new_scale_color() +
  geom_sf(data=yeartracks.fall, 
          aes(col=Year), pch=19, cex=0.5, lwd=0.5, alpha=0.7) +
  scale_color_continuous(
    limits = c(1982,2021), 
    breaks = c(1982,1990, 2000, 2010, 2021),
    labels = c('1982', ' ', ' ', ' ', '2021'),
    guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
  )+
  facet_wrap(vars(Model)) +
  coord_sf(xlim=c(-70, -65),
           ylim=c(42, 45.5),
           crs="EPSG:4326") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle('Fall Northeastern Range Edge (95% Quantile)') +
  theme(legend.position = 'bottom',
        axis.text.x=element_text(size=6, angle=20),
        axis.text.y=element_text(size=6),
        strip.text=element_text(size=8))

ggsave(here("VAST_runs/large/GMRF/Testing_Results/Spatial_RangeEdge_NE_Fall.png"),
       fall.vis,
       width=8, height=6)
