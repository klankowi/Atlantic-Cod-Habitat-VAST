rm(list=ls())

library(here)
library(sf)
library(VAST)
library(tidyverse)
library(stars)

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

# Load data
load(here('VAST_runs/medium/Overall_BC/ALL/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.RData'))
setwd(here("VAST_runs/medium/Overall_BC/ALL"))
rm(list=setdiff(ls(), c("fit")))

# Reload model
fit = reload_model(x = fit)

# Call labels
years_to_plot = fit$years_to_plot

# Extract residuals
dharmaRes = summary(fit, what = "residuals", working_dir = NA)

# Extract n_cells for plotting residuals
n_cells = ceiling(fit$data_list$n_s) * 4

# Function to aggregate p values
aggregate_pvalues = function(pvec, na.rm = TRUE) {
  chisq = -2 * sum(log(pvec), na.rm = na.rm)
  p = 1 - pchisq(chisq, df = 2 * length(pvec))
}

# Map data
PlotDF = cbind(Lat = fit$data_frame[, "Lat_i"], 
               Lon = fit$data_frame[, "Lon_i"], 
               x2i = 1:fit$data_list$n_i, 
               Include = TRUE)

# Create matrix of results
Y_it = matrix(NA, nrow = nrow(fit$data_frame), ncol = fit$data_list$n_t)
Y_it[cbind(1:fit$data_list$n_i, fit$data_list$t_i + 1)] = dharmaRes$scaledResiduals
Y_it = Y_it[, years_to_plot, drop = FALSE]

# Color palette for plotting
col_function = colorRampPalette(colors = c("darkblue", "lightblue", 
                                           "white", "pink", "red"))

# Extract data for plotting
Y_gt = Y_it
map_list = list(PlotDF = PlotDF)
Y_gt = Y_gt[map_list$PlotDF[which(map_list$PlotDF[, "Include"] > 
                                    0), "x2i"], , drop = FALSE]
panel_labels = fit$year_labels
file_name = "quantile_residuals_on_map"
working_dir = here('VAST_runs/medium/Overall_BC/ALL')
setwd(working_dir)
fun = aggregate_pvalues
col = col_function

# Call data
MapSizeRatio = map_list$MapSizeRatio

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

Raster_proj = vector("list", length = ncol(Y_gt))
Big_Data <- vector("list", length=length(panel_labels))

for (tI in 1:ncol(Y_gt)) {
  print(tI)
  Points_orig = sp::SpatialPointsDataFrame(coords = loc_g, 
                                           data = data.frame(y = Y_gt[, tI]), 
                                           proj4string = CRS_orig)
  Points_LongLat = sp::spTransform(Points_orig, sp::CRS("+proj=longlat"))
  Points_proj = sp::spTransform(Points_orig, CRS_proj)
  
  Zlim = c(0, 1)
  xlim = Points_proj@bbox[1, ]
  ylim = Points_proj@bbox[2, ]
  
  cell.size = mean(diff(Points_proj@bbox[1, ]), 
                   diff(Points_proj@bbox[2, ]))/floor(sqrt(n_cells))
  Raster_layer = raster::raster(Points_proj, crs = CRS_proj, 
                                nrows = floor(sqrt(n_cells)), ncols = floor(sqrt(n_cells)))
  Raster_proj[[tI]] = raster::rasterize(x = Points_proj@coords, 
                                        y = Raster_layer, 
                                        field = as.numeric(Points_proj@data[, 1]), 
                                        fun = mean)
  
  Big_Data[[tI]] <- st_as_sf(st_as_stars(Raster_proj[[tI]]))
  Big_Data[[tI]]$TS <- panel_labels[[tI]]

}

allData <- do.call("rbind", Big_Data)
allData <- allData %>% 
  separate(TS, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Season = factor(Season, levels=c('Spring', 'Fall')))

allData$color <- NA
lowvals <- colorRampPalette(c('darkred', 'tomato', 'pink'))
lowvals <- lowvals(10)

highvals <- colorRampPalette(c('lightblue', 'cornflowerblue', 'navy'))
highvals <- highvals(10)

for(i in 1:nrow(allData)){
  if(allData$layer[i]>=0.00 & allData$layer[i]<0.01){allData$color[i] <- lowvals[1]}
  if(allData$layer[i]>=0.01 & allData$layer[i]<0.02){allData$color[i] <- lowvals[2]}
  if(allData$layer[i]>=0.02 & allData$layer[i]<0.03){allData$color[i] <- lowvals[3]}
  if(allData$layer[i]>=0.03 & allData$layer[i]<0.04){allData$color[i] <- lowvals[4]}
  if(allData$layer[i]>=0.04 & allData$layer[i]<0.05){allData$color[i] <- lowvals[5]}
  if(allData$layer[i]>=0.05 & allData$layer[i]<0.06){allData$color[i] <- lowvals[6]}
  if(allData$layer[i]>=0.06 & allData$layer[i]<0.07){allData$color[i] <- lowvals[7]}
  if(allData$layer[i]>=0.07 & allData$layer[i]<0.08){allData$color[i] <- lowvals[8]}
  if(allData$layer[i]>=0.08 & allData$layer[i]<0.09){allData$color[i] <- lowvals[9]}
  if(allData$layer[i]>=0.09 & allData$layer[i]<=0.10){allData$color[i] <- lowvals[10]}
  
  if(allData$layer[i]>=0.90 & allData$layer[i]<0.91){allData$color[i] <- highvals[1]}
  if(allData$layer[i]>=0.91 & allData$layer[i]<0.92){allData$color[i] <- highvals[2]}
  if(allData$layer[i]>=0.92 & allData$layer[i]<0.93){allData$color[i] <- highvals[3]}
  if(allData$layer[i]>=0.93 & allData$layer[i]<0.94){allData$color[i] <- highvals[4]}
  if(allData$layer[i]>=0.94 & allData$layer[i]<0.95){allData$color[i] <- highvals[5]}
  if(allData$layer[i]>=0.95 & allData$layer[i]<0.96){allData$color[i] <- highvals[6]}
  if(allData$layer[i]>=0.96 & allData$layer[i]<0.97){allData$color[i] <- highvals[7]}
  if(allData$layer[i]>=0.97 & allData$layer[i]<0.98){allData$color[i] <- highvals[8]}
  if(allData$layer[i]>=0.98 & allData$layer[i]<0.99){allData$color[i] <- highvals[9]}
  if(allData$layer[i]>=0.99 & allData$layer[i]<=1.0){allData$color[i] <- highvals[10]}
  
}


# Plot all years at once
supfig1a <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Spring',], 
          aes(fill=layer, col=layer)) +
  
  scale_color_gradientn(colours=c(lowvals, rep('transparent', 80),
                                  highvals), 
                        limits=c(0,1),
                        breaks=seq(0, 1, by=0.01), 
                        labels=c(0, 
                                 rep('', 9), 0.10,
                                 rep('', 79), 0.90, 
                                 rep('', 9), 1), 
                        na.value = "transparent",
                        guide='none') +
  
  scale_fill_gradientn(colours=c(lowvals, rep('transparent', 80),
                                  highvals), 
                        limits=c(0,1),
                        breaks=seq(0, 1, by=0.01), 
                        labels=c(0, 
                                 rep('', 9), 0.10,
                                 rep('', 79), 0.90, 
                                 rep('', 9), 1), 
                        na.value = "transparent",
                        guide=guide_colorbar(nbin=100)) +
  
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.05,col='gray20') +
  
  coord_sf(xlim=c(-76, -65),
           ylim=c(36.75, 45),
           crs="EPSG:4326") +
  
  facet_wrap(vars(Year)) +
  
  ggtitle('Spring time series') +
  
  labs(color = 'Residuals\np-value', fill='Residuals\np-value') +
  
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) + 
  guides(fill=guide_colorbar(ticks.colour = "transparent"))

supfig1b <- ggplot() +
  geom_sf(data=allData[allData$Season == 'Fall',], 
          aes(fill=layer, col=layer)) +
  
  scale_color_gradientn(colours=c(lowvals, rep('transparent', 80),
                                  highvals), 
                        limits=c(0,1),
                        breaks=seq(0, 1, by=0.01), 
                        labels=c(0, 
                                 rep('', 9), 0.10,
                                 rep('', 79), 0.90, 
                                 rep('', 9), 1), 
                        na.value = "transparent",
                        guide='none') +
  
  scale_fill_gradientn(colours=c(lowvals, rep('transparent', 80),
                                 highvals), 
                       limits=c(0,1),
                       breaks=seq(0, 1, by=0.01), 
                       labels=c(0, 
                                rep('', 9), 0.10,
                                rep('', 79), 0.90, 
                                rep('', 9), 1), 
                       na.value = "transparent",
                       guide=guide_colorbar(nbin=100)) +
  
  geom_sf(data=coast, fill='gray', col='darkgray') +
  
  geom_sf(data=regions, fill=NA,lwd=0.05,col='gray20') +
  
  coord_sf(xlim=c(-76, -65),
           ylim=c(36.75, 45),
           crs="EPSG:4326") +
  
  facet_wrap(vars(Year)) +
  
  ggtitle('Fall time series') +
  
  labs(color = 'Residuals\np-value', fill='Residuals\np-value') +
  
  theme(legend.position = 'right',
        axis.text.x = element_text(size=6, angle=20),
        axis.text.y = element_text(size=6),
        strip.text.x=element_text(margin=margin(0.1,0,0.1,0, "cm")),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8)) + 
  guides(fill=guide_colorbar(ticks.colour = "transparent"))

ggsave(here('Documentation/Figures/Supplementary/Supp Fig 2a.pdf'),
       supfig1a,
       height=8.5, width=11, units='in',
       dpi = 300)

ggsave(here('Documentation/Figures/Supplementary/Supp Fig 2b.pdf'),
       supfig1b,
       height=8.5, width=11, units='in',
       dpi = 300)

ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Supplementary/Supp Fig 2a.pdf",
       supfig1a,
       height=8.5, width=11, units='in',
       dpi = 300)

ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Supplementary/Supp Fig 2b.pdf",
       supfig1b,
       height=8.5, width=11, units='in',
       dpi = 300)
