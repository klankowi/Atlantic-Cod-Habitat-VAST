# Personalized plotting using ggplot for various VAST default plots

# Prepare workspace
rm(list=ls())

# Load libraries
library(VAST)
library(tidyverse)
library(sf)
library(here)
library(ggpubr)
library(ggnewscale)
library(patchwork)

# Load functions
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', size=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load data
load(here('VAST_runs/refine_effort/refine_effort.RData'))
rm(list=setdiff(ls(), c("fit", "%notin%", "year.labs")))

# Create objects needed to plot
Sdreport = fit$parameter_estimates$SD
SD = TMB::summary.sdreport(Sdreport)
TmbData = fit$data_list

# Name where data are stored in report
CogName = "mean_Z_ctm"
EffectiveName = "effective_area_ctl"

# Set labels
category_names = c('Small', 'Medium', 'Large', 'Unknown') 

###############################################################################
####                     Seasonal Center of gravity                        ####
###############################################################################
# Pull data
SD_mean_Z_ctm = array(NA, dim = c(unlist(TmbData[c("n_c", "n_t", "n_m")]), 2), 
                      dimnames = list(NULL, NULL, NULL,
                                      c("Estimate", "Std. Error")))
# Pull standard error
SD_mean_Z_ctm[] = SD[which(rownames(SD) == CogName), 
                           c("Estimate", "Std. Error")]
# Name dimensions      
names(dim(SD_mean_Z_ctm)) <- c('Category', 'Time', 'Dimension', 'Est.Err')
# Separate sizes      
SD_mean_Z_ctm_small <- SD_mean_Z_ctm[1,,,]
SD_mean_Z_ctm_medium <- SD_mean_Z_ctm[2,,,]
SD_mean_Z_ctm_large <- SD_mean_Z_ctm[3,,,]
SD_mean_Z_ctm_unknown <- SD_mean_Z_ctm[4,,,]
# Concatenate sizes to list, name
size.list.cog <- list(SD_mean_Z_ctm_small, SD_mean_Z_ctm_medium,
                      SD_mean_Z_ctm_large, SD_mean_Z_ctm_unknown)
names(size.list.cog) <- c('Small', 'Medium', 'Large', 'Unknown')  
rm(SD_mean_Z_ctm, SD_mean_Z_ctm_large, SD_mean_Z_ctm_unknown,
   SD_mean_Z_ctm_small, SD_mean_Z_ctm_medium)

###############################################################################
####                   Seasonal Effective area occupied                    ####
###############################################################################
# Pull data
SD_effective_area_ctl = 
  SD_log_effective_area_ctl = array(NA, 
                                    dim = c(unlist(
                                      TmbData[c("n_c", 
                                                "n_t", 
                                                "n_l")]),2), 
                                    dimnames = list(NULL, 
                                                    NULL, 
                                                    NULL, 
                                                    c("Estimate", 
                                                      "Std. Error")))
# Pull standard error
SD_effective_area_ctl[] = SD[which(rownames(SD) == 
                                           EffectiveName), 
                                   c("Estimate", "Std. Error")]
# Name dimensions
names(dim(SD_effective_area_ctl)) <- c('Category', 
                                       'Time', 
                                       'Dimension', 
                                       'Est.Err')
# Separate sizes      
SD_effective_area_ctl_small <- SD_effective_area_ctl[1,,,]
SD_effective_area_ctl_medium <- SD_effective_area_ctl[2,,,]
SD_effective_area_ctl_large <- SD_effective_area_ctl[3,,,]
SD_effective_area_ctl_unknown <- SD_effective_area_ctl[4,,,]
# Concatenate sizes to list, name
size.list.eao <- list(SD_effective_area_ctl_small, SD_effective_area_ctl_medium,
                      SD_effective_area_ctl_large, SD_effective_area_ctl_unknown)
names(size.list.eao) <- c('Small', 'Medium', 'Large', 'Unknown')       
rm(SD_effective_area_ctl, SD_effective_area_ctl_large, 
   SD_effective_area_ctl_small, SD_effective_area_ctl_unknown,
   SD_effective_area_ctl_medium)
rm(SD, Sdreport, TmbData, CogName, EffectiveName, SD_log_effective_area_ctl)

###############################################################################
####                       Plot EAO and COG together                       ####
###############################################################################
# Loop through sizes, split into seasons
for(i in 1:length(size.list.cog)){ # Number of sizes
  #for(j in 1:dim(season.list.cog[[i]])[1]){ # Number of categories
    # COG
    cog <- as.data.frame(size.list.cog[[i]][,,])
    colnames(cog) <- c('easting', 'northing', 'e.sd', 'n.sd')
    cog$YearSeas <- year.labs
    cog <- separate(cog, YearSeas, 
                                into = c("Year", "Season"), sep = " (?=[^ ]+$)")
    cog$Year <- as.numeric(cog$Year)
    spring.cog <- subset(cog, Season == 'Spring')
    fall.cog <- subset(cog, Season == 'Fall')
    rm(cog)
    
    # EAO
    eao <- as.data.frame(size.list.eao[[i]][,])
    colnames(eao) <- c('area.occ', 'sd.err')
    eao$YearSeas <- year.labs
    eao <- separate(eao, YearSeas, 
                           into = c("Year", "Season"), sep = " (?=[^ ]+$)")
    eao$Year <- as.numeric(eao$Year)
    spring.eao <- subset(eao, Season == 'Spring')
    fall.eao <- subset(eao, Season == 'Fall')
    rm(eao)
    
    # Both Spring
    SD_plotting.spring <- merge(spring.cog, spring.eao, by=c("Year", "Season"))
    
    # Plot northing
    northing <- ggplot(SD_plotting.spring) +
      geom_line(aes(x=Year, y=northing), col='#00BFC4', lwd=1) +
      geom_ribbon(aes(ymin=northing-n.sd,
                      ymax=northing+n.sd,
                      x=Year),
                  fill=alpha('#00BFC4', 0.2)) +
      ylim(c(4600,4825)) +
      ylab("Northing (km)") +
      xlab("")
    
    # Plot easting
    easting <- ggplot(SD_plotting.spring) +
      geom_line(aes(x=Year, y=easting), col='#00BFC4', lwd=1) +
      geom_ribbon(aes(ymin=easting-e.sd,
                      ymax=easting+e.sd,
                      x=Year),
                  fill=alpha('#00BFC4', 0.2)) +
      ylim(c(300, 700)) +
      ylab("Easting (km)")+
      xlab("")
    
    # Plot effective area occupied
    arr.occ <- ggplot(SD_plotting.spring) +
      geom_line(aes(x=Year, y=area.occ), col='#00BFC4', lwd=1) +
      geom_ribbon(aes(ymin=area.occ-sd.err,
                      ymax=area.occ+sd.err,
                      x=Year),
                  fill=alpha('#00BFC4', 0.2)) +
      ylim(c(-5000, 75000)) +
      ylab(bquote("Area Occupied km "^2))
    
    
    # Arrange to plot
    spring <- ggarrange(northing, easting, arr.occ, nrow=3)
    spring <- annotate_figure(spring, top = text_grob('Spring',
                                                      color='black',
                                                      face='bold',
                                                      size=14,
                                                      vjust=1.4,
                                                      hjust=0.1))
    
    # Both Fall
    SD_plotting.fall <- merge(fall.cog, fall.eao, by=c("Year", "Season"))
    
    # Plot northing
    northing <- ggplot(SD_plotting.fall) +
      geom_line(aes(x=Year, y=northing), col='#F8766D', lwd=1) +
      geom_ribbon(aes(ymin=northing-n.sd,
                      ymax=northing+n.sd,
                      x=Year),
                  fill=alpha('#F8766D', 0.2)) +
      ylim(c(4600,4825)) +
      ylab(" ") +
      xlab("")
    
    # Plot easting
    easting <- ggplot(SD_plotting.fall) +
      geom_line(aes(x=Year, y=easting), col='#F8766D', lwd=1) +
      geom_ribbon(aes(ymin=easting-e.sd,
                      ymax=easting+e.sd,
                      x=Year),
                  fill=alpha('#F8766D', 0.2)) +
      ylim(c(300, 700)) +
      ylab(" ")+
      xlab("")
    
    # Plot effective area occupied
    arr.occ <- ggplot(SD_plotting.fall) +
      geom_line(aes(x=Year, y=area.occ), col='#F8766D', lwd=1) +
      geom_ribbon(aes(ymin=area.occ-sd.err,
                      ymax=area.occ+sd.err,
                      x=Year),
                  fill=alpha('#F8766D', 0.2)) +
      ylim(c(-5000, 75000)) +
      ylab(" ")
    
    
    # Arrange to plot
    fall <- ggarrange(northing, easting, arr.occ, nrow=3)
    fall <- annotate_figure(fall, top = text_grob('Fall',
                                                      color='black',
                                                      face='bold',
                                                      size=14,
                                                      vjust=1.4,
                                                      hjust=0.1))
    
    # Arrange to plot
    both <- ggarrange(spring, fall, nrow=1) + bgcolor('white')
    both <- annotate_figure(both, top=text_grob(paste0(category_names[i],
                                                       ' size class'),
                                                color='black',
                                                face='bold',
                                                size=16,
                                                vjust=2.75))
    #plot(both)
    ggsave(both,
           filename=paste0(here(), "/Plot_Output/location.info.",
                           category_names[i], '.png'),
           width = 10, height = 8, units='in')
}

###############################################################################
####                            Visualize COG                              ####
###############################################################################
# Load spatial information
coast <- ecodata::coast
coast <- st_transform(coast, "EPSG:32619")
stocks <- st_read(here('Data/GIS/codstox.shp'))
stocks <- st_transform(stocks, "EPSG:32619")
new_bb <- st_bbox(stocks)

for(i in 1:length(size.list.cog)){ # Categories
    # COG
    SD_plotting.cog <- size.list.cog[[i]]
    SD_plotting.cog <- as.data.frame(SD_plotting.cog[,,])
    colnames(SD_plotting.cog) <- c('easting', 'northing', 'e.sd', 'n.sd')
    SD_plotting.cog$YearSeas <- year.labs
    SD_plotting.cog$easting <- SD_plotting.cog$easting * 1000
    SD_plotting.cog$northing <- SD_plotting.cog$northing * 1000
    SD_plotting.cog <- separate(SD_plotting.cog, YearSeas, 
                                into = c("Year", "Season"), sep = " (?=[^ ]+$)")
    SD_plotting.cog$Year <- as.numeric(SD_plotting.cog$Year)
    
    # Convert to sf for plotting
    SD_plotting <- st_as_sf(SD_plotting.cog, coords=c("easting", "northing"))
    st_crs(SD_plotting) <- "EPSG:32619"
    
    # Spring
    SD_plotting.spring <- subset(SD_plotting, Season =='Spring')
    points <- st_cast(st_geometry(SD_plotting.spring), "POINT") 
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
    # Plot
    spring <- ggplot() +
      geom_sf(data=coast, fill='gray') +
      geom_sf(data=stocks, aes(col=STOCK), fill='transparent', lwd=0.25) +
      guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
      new_scale_color() +
      geom_sf(data=SD_plotting.spring, aes(col=Year), pch=19, cex=0.5) +
      scale_color_continuous(
        limits = c(1982,2021), 
        breaks = c(1982,1990, 2000, 2010, 2021),
        labels = c('1982', ' ', ' ', ' ', '2021'),
        guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
      )+
      geom_sf(data=t.spring, aes(col=Year)) +
      coord_sf(xlim=c(new_bb[1], new_bb[3]),
               ylim=c(new_bb[2], new_bb[4])) +
      xlab("Longitude") + ylab("Latitude") +
      ggtitle('Spring')
    
    # Fall
    SD_plotting.fall <- subset(SD_plotting, Season =='Fall')
    points <- st_cast(st_geometry(SD_plotting.fall), "POINT") 
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
    # Plot
    fall <- ggplot() +
      geom_sf(data=coast, fill='gray') +
      geom_sf(data=stocks, aes(col=STOCK), fill='transparent', lwd=0.25) +
      guides(col=guide_legend(title="Stock", nrow=2,byrow=TRUE)) +
      new_scale_color() +
      geom_sf(data=SD_plotting.fall, aes(col=Year), pch=19, cex=0.5) +
      scale_color_continuous(
        limits = c(1982,2021), 
        breaks = c(1982,1990, 2000, 2010, 2021),
        labels = c('1982', ' ', ' ', ' ', '2021'),
        guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE)
      )+
      geom_sf(data=t.fall, aes(col=Year)) +
      coord_sf(xlim=c(new_bb[1], new_bb[3]),
               ylim=c(new_bb[2], new_bb[4])) +
      xlab("Longitude") + ylab("Latitude") +
      ggtitle('Fall')
    
    
    combined <- spring + fall & 
      theme(legend.position = "bottom") 
    combined <- combined + plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(category_names[i],
                               " size class"),
      theme = theme(plot.title = element_text(size = 14,
                                              face = 'bold',
                                              hjust=0.5,
                                              vjust=-0.5)))
    # Save
    ggsave(combined,
           filename=paste0(here(), "/Plot_Output/cog.vis.",
                           category_names[i], '.png'),
           width = 10, height = 8, units='in')
}
