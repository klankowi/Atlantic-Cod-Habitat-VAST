rm(list=ls())

# Libraries
library(VAST)
library(here)
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

# Load data
setwd(here('VAST_runs/medium/Overall_BC/ALL'))
load('Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.RData')
rm(list=setdiff(ls(), c('fit', 'dharmaRes')))
fit <- reload_model(fit)

# Set parameters
n_cells_residuals = ceiling(fit$data_list$n_s) * 4
year_labels = fit$year_labels
years_to_plot = fit$years_to_plot
dharmaRes = summary( fit, what="residuals", working_dir=NA )

# Aggregate quantile residuals
aggregate_pvalues = function( pvec, na.rm=TRUE ){
  chisq = -2*sum(log(pvec), na.rm=na.rm)
  p = 1 - pchisq( chisq, df=2*length(pvec) )
}

# Make plot using plot_variable
PlotDF = cbind( "Lat"=fit$data_frame[,'Lat_i'], 
                "Lon"=fit$data_frame[,'Lon_i'], 
                "x2i"=1:fit$data_list$n_i, 
                "Include"=TRUE)
# map_list = make_map_info(Region = fit$extrapolation_list$Area_km2_x,
#                          spatial_list = fit$spatial_list,
#                          Extrapolation_List = fit$extrapolation_list)

Y_it = matrix(NA, nrow=nrow(fit$data_frame), 
              ncol=fit$data_list$n_t )

Y_it[ cbind(1:fit$data_list$n_i,fit$data_list$t_i+1) ] = 
  dharmaRes$scaledResiduals

Y_it = Y_it[,years_to_plot,drop=FALSE]

zlim = range(Y_it, na.rm=T)
# MapSizeRatio = map_list$MapSizeRatio

# Call location list
loc_g = cbind(PlotDF[,'Lon'], PlotDF[,'Lat'])
colnames(loc_g) <- c('Lon', 'Lat')
# Call projections
CRS_orig = sp::CRS("+proj=longlat")
CRS_proj = sp::CRS("EPSG:26919")

# Make list to append data
Big_Data <- vector("list", length=length(year_labels))

# Call other spatial data
coast <- st_transform(ecodata::coast, CRS_proj)
regions <- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
regions <- st_transform(regions, CRS_proj)
regions <- st_make_valid(regions)

Y_gt = Y_it[PlotDF[which(PlotDF[, "Include"] > 
                                    0), "x2i"], , drop = FALSE]

for (tI in 1: ncol(Y_it)) {
  Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                           data = data.frame(y = Y_gt[, tI]), 
                                           proj4string = CRS_orig)
  Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
  Points_proj = sp::spTransform(Points_orig, CRS_proj)
  Zlim = zlim
  xlim = Points_proj@bbox[1, ]
  ylim = Points_proj@bbox[2, ]
  
  cell.size = mean(diff(Points_proj@bbox[1, ]), 
                   diff(Points_proj@bbox[2, ]))/
    floor(sqrt(n_cells_residuals))
  
  Points_sf = sf::st_as_sf(Points_proj)
  grid = sf::st_make_grid(Points_sf, cellsize = cell.size)
  grid_i = sf::st_intersects(Points_sf, grid)
  grid = sf::st_sf(grid, y = tapply(Points_sf$y, 
                                    INDEX = factor(as.numeric(grid_i), 
                                                   levels = 1:length(grid)), 
                                    FUN = mean, 
                                    na.rm = TRUE))
  grid$TS <- year_labels[tI]
  
  Big_Data[[tI]] <- grid
  
}

allData <- do.call("rbind", Big_Data)

allData <- allData %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

allData$Sig[allData$y <= 0.05] <- 1
allData$Sig[allData$y > 0.05] <- 0
#allData$y <- log(allData$y)

# Plot all years at once
sseries <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Spring',], 
          aes(fill=Sig, col=Sig)) +
  scale_color_gradient2(low="transparent", mid="transparent", high="red", 
                        limits = c(0, 1),
                        midpoint = 0.5,
                        na.value = 'transparent') +
  scale_fill_gradient2(low="transparent", mid="transparent", high="red", 
                       limits = c(0, 1),
                       midpoint=0.5,
                       na.value='transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Spring Residuals, Medium Cod')

fseries <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Fall',], 
          aes(fill=Sig, col=Sig)) +
  scale_color_gradient2(low="transparent", mid="transparent", high="red", 
                        limits = c(0, 1),
                        midpoint = 0.5,
                        na.value = 'transparent') +
  scale_fill_gradient2(low="transparent", mid="transparent", high="red", 
                       limits = c(0, 1),
                       midpoint=0.5,
                       na.value='transparent') +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  facet_wrap(vars(Year)) +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Fall Residuals, Medium Cod')

allYear <- allData %>% 
  mutate(Cell = rep(seq(1, 1088), 80)) %>% 
  group_by(Cell, Season) %>% 
  summarise(Sig = sum(Sig, na.rm = T)) %>% 
  dplyr::select(Cell, Season, grid, Sig)
  
ssum <- ggplot() +
  geom_sf(data=allYear[allYear$Season == 'Spring',], 
          aes(fill=as.factor(Sig), col=as.factor(Sig))) +
  scale_color_manual(values = c('transparent', 'lightblue',
                                'royalblue', 'navy', 'black'),
                     breaks=c('0', '1','2','3','4')) +
  scale_fill_manual(values = c('transparent', 'lightblue',
                                'royalblue', 'navy', 'black'),
                     breaks=c('0', '1','2','3','4')) +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Spring Residuals, Medium Cod')

fsum <- ggplot() +
  geom_sf(data=allYear[allYear$Season == 'Fall',], 
          aes(fill=as.factor(Sig), col=as.factor(Sig))) +
  scale_color_manual(values = c('transparent', 'lightblue',
                                'royalblue', 'navy', 'black'),
                     breaks=c('0', '1','2','3','4')) +
  scale_fill_manual(values = c('transparent', 'lightblue',
                               'royalblue', 'navy', 'black'),
                    breaks=c('0', '1','2','3','4')) +
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.1,col='black') +
  
  coord_sf(xlim=c(st_bbox(grid)[1], st_bbox(grid)[3]),
           ylim=c(st_bbox(grid)[2], st_bbox(grid)[4])) +
  labs(col='Abund.', fill='Abund.') +
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) +
  ggtitle('Fall Residuals, Medium Cod')


ggsave(sseries,
       filename=here('VAST_runs/medium/Overall_BC/ALL/Spring_Mapped_Residuals.png'),
       width = 8.5, height=11, units='in')
ggsave(fseries,
       filename=here('VAST_runs/medium/Overall_BC/ALL/Fall_Mapped_Residuals.png'),
       width = 8.5, height=11, units='in')

ggsave(ssum,
       filename=here('VAST_runs/medium/Overall_BC/ALL/Spring_Total_Residuals.png'),
       width = 8.5, height=11, units='in')
ggsave(fsum,
       filename=here('VAST_runs/medium/Overall_BC/ALL/Fall_Total_Residuals.png'),
       width = 8.5, height=11, units='in')
