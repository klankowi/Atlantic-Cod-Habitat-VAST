rm(list=ls())

library(VAST)
library(sf)
library(here)
library(tidyverse)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "right",
                legend.background = element_rect(fill='transparent', colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Functions
'%notin%' <- function(x,y)!('%in%'(x,y))

# Identify sizes
sizes = c('small', 'medium', 'large')

# Create blank holder dataframes for model results
bigap <- data.frame(
  Strata=NA, Season=NA, Size=NA,
  Slope=NA, pVal = NA, adjR=NA,
  sig=NA
)

bigapy <- data.frame(
  Strata=NA, Year=NA, Season=NA, Area=NA, n.cells=NA, 
  Estimate=NA, Std.Err=NA, Size=NA
)

listcells <- vector(length=1, mode='list')

# Loop through sizes
for(v in 1:length(sizes)){
  # Load data
  load(here(paste0('VAST_runs/', 
            sizes[v],
            '/Overall_BC/ALL/Overall_BC_', 
            sizes[v], 'cod_allstrat_natsplin_fsON_ALL.RData')))
  rm(covars, extrap_info_aja, scaled.covars, settings, strata_use,
     surveys, survs, vast_extrap_grid, hab_formula, seas.labs, 
     working_dir, year.labs)
  assign(paste0(sizes[v]), fit)
  
  # Pull individuals per grid cell  per strata per time step 
  hold <- fit$Report$Index_gctl
  
  # Pull grid cell areas
  gca <- fit$spatial_list$a_gl
  
  # Convert to dataframe
  rownames(gca) <- seq(1:nrow(gca))
  gca <- gca %>% 
    as.data.frame() %>%
    mutate(Site = rownames(gca)) %>% 
    pivot_longer(cols=c('V1', 'V2', 'V3', 'V4', 'V5'))
  
  # Identify strata
  msn <- data.frame(
    name=c('V1', 'V2', 'V3', 'V4', 'V5'),
    strata = fit$settings$strata.limits$STRATA
  )
  
  # Join strata names to grid cell area df
  gca <- left_join(gca, msn, by=c('name'))
  gca <- gca %>% 
    dplyr::select(-name) %>% 
    rename(area = value)
  
  # Call strata
  strats <- fit$settings$strata.limits$STRATA
  
  # Create blank holder dataframe
  mediumapy <- data.frame(
    Strata=NA, Year=NA, Season=NA, Area=NA, n.cells=NA, 
    Estimate=NA, Std.Err=NA, Size=NA
  )
  
  # Loop through strata
  for(q in 1:length(strats)){

    # Strip individuals per cell for that strata
    ipc <- hold[,1,,q]
    
    # Convert to dataframe
    ipc <- ipc %>% 
      as.data.frame() %>% 
      mutate(Site = rownames(ipc)) %>% 
      pivot_longer(cols=colnames(ipc)) %>% 
      #separate(name, into = c('Year', 'Season')) %>% 
      rename(inds = value) %>% 
      as.data.frame()
    
    # Identify size of 90% of population
    ipa <- ipc %>% 
      group_by(name) %>% 
      mutate(nfp = sum(inds) * 0.90) %>% 
      dplyr::select(name, nfp) %>% 
      unique() %>% 
      as.data.frame()
    
    # Blank holder dataframe
    apy <- data.frame(
      Strata = NA,
      Time = NA,
      Area = NA,
      n.cells=NA
    )
    
    # Separate by timestep
    ipc <- split(ipc, f=ipc$name)
    
    # Loop through timesteps
    for(i in 1:length(ipc)){
      # Order most to least fish in cell
      ipc[[i]] <- ipc[[i]][with(ipc[[i]], order(inds,
                                                decreasing = T)),]
      # Loop through cells
      for(j in 1:nrow(ipc[[i]])){
        # If sum of cells 1:j is less than 90% of population, next j
        if(sum(ipc[[i]]$inds[1:j]) < 
           ipa$nfp[ipa$name == ipc[[i]]$name[1]]
        ){
          next()
        }
        # If sum of cells 1:j is more than or equal to, save j
        if(sum(ipc[[i]]$inds[1:j]) >= 
           ipa$nfp[ipa$name == ipc[[i]]$name[1]]
        ){
          #print(j)
          break()
        }
      }
      
      # Save cell names needed to hit 90% of population
      gcu <- ipc[[i]]$Site[1:j]
      
      # Convert cell names to list
      listgcu <- list(gcu)
      
      # Append to list of cells, each item is a time step
      listcells <- c(listcells, listgcu)
      
      # Extract areas of grid cells needed to hit 90% of population
      gcc <- gca[gca$Site %in% gcu &
                   gca$strata == strats[q],]
      
      # Sum area
      au <- sum(gcc$area)
      
      # Append to blank holder dataframe
      apy <- rbind(apy,
                   c(strats[q],
                     fit$year_labels[i],
                     au,
                     j))
      
    }
    
    # Remove starter row
    apy <- apy[!is.na(apy$Strata),]
    
    # Clean
    apy <- apy %>% 
      separate(Time, into=c('Year', 'Season')) %>% 
      mutate(Year = as.numeric(Year),
             Area = as.numeric(Area),
             Season = factor(Season, levels=c('Spring', 'Fall')))
    
    # Load index of abundance
    index <- read.csv(here(paste0("VAST_runs/",
                                  sizes[v], "/Overall_BC/ALL/Index.csv")))
    
    # Clean
    index <- index %>% 
      dplyr::select(-Category, -Units, -Std..Error.for.ln.Estimate.) %>% 
      separate(Time, into=c('Year', 'Season')) %>% 
      filter(Stratum == strats[q]) %>% 
      mutate(Year = as.numeric(Year),
             Season = factor(Season, levels=c("Spring", 'Fall'))) %>% 
      rename(Std.Err = Std..Error.for.Estimate,
             Strata = Stratum)
    
    # Merge area per year and index per year
    apy <- left_join(apy, index, 
                     by=c('Strata', 'Year', 'Season'))
    
    # Plot comparison of area utilized and abundance
    #print(
      # ggplot() +
      #   geom_line(data=apy,
      #             aes(x=Year, y=scale(Area))) +
      #   geom_line(data=apy,
      #             aes(x=Year, y=scale(Estimate)),
      #             col='red') +
      #   scale_y_continuous(
      #     # Features of the first axis
      #     name = "Relative Area Occupied",
      #     # Add a second axis and specify its features
      #     sec.axis = sec_axis(transform~., name="Relative Abundance")
      #   ) +
      #   facet_wrap(vars(Season)) +
      #   ggtitle(paste0(sizes[v], ' ', strats[q]))
    #)
    
    # Plot linear model of log(Estimate) to log(Area)
    #print(
      # ggplot() + 
      #   geom_point(data=apy,
      #              aes(x=log(Estimate), y=log(Area))) +
      #   geom_smooth(data=apy,
      #               aes(x=log(Estimate), y=log(Area)),
      #               method='lm', alpha=0.3) +
      #   facet_wrap(vars(Season)) +
      #   ggtitle(paste0(sizes[v], ' ', strats[q]))
    #)
    
    # Extract model parameters
    apytest <- apy %>% 
      group_by(Season) %>% 
      mutate(Slope = lm((Area) ~ (Estimate))$coefficients[2],
             pVal = summary(lm((Area) ~ (Estimate)))$coefficients[2,4],
             adjR = summary(lm((Area) ~ (Estimate)))$adj.r.squared) %>% 
      dplyr::select(Strata, Season, Slope, pVal, adjR) %>% 
      unique() %>% as.data.frame()
    
    # Make it easier to see significance
    apytest$sig <- ' '
    apytest$sig[apytest$pVal <=0.05] <- '*'
    
    # Append size
    apytest$Size = paste0(sizes[v])
    apy$Size = paste0(sizes[v])
    
    # Append to big dataframes
    bigap <- rbind(bigap, apytest)
    mediumapy <- rbind(mediumapy, apy)

  }
  
  # Append to big dataframes
  bigapy <- rbind(bigapy, mediumapy)
  rm(apytest, gca, gcc, index, ipa, ipc, listgcu,
     msn, au, gcu)
  rm(mediumapy)
}

# Remove starter rows
bigap <- bigap[!is.na(bigap$Strata),]
bigapy <- bigapy[!is.na(bigapy$Strata),]
listcells <- listcells[-c(1)]

names(listcells) <- c(
  paste0('Small ', 'ALL ', fit$year_labels),
  paste0('Small ', 'EGOM ', fit$year_labels),
  paste0('Small ', 'GBK ', fit$year_labels),
  paste0('Small ', 'SNE ', fit$year_labels),
  paste0('Small ', 'WGOM ', fit$year_labels),
  
  paste0('Medium ', 'ALL ', fit$year_labels),
  paste0('Medium ', 'WGOM ', fit$year_labels),
  paste0('Medium ', 'GBK ', fit$year_labels),
  paste0('Medium ', 'EGOM ', fit$year_labels),
  paste0('Medium ', 'SNE ', fit$year_labels),

  paste0('Large ', 'ALL ', fit$year_labels),
  paste0('Large ', 'EGOM ', fit$year_labels),
  paste0('Large ', 'GBK ', fit$year_labels),
  paste0('Large ', 'WGOM ', fit$year_labels),
  paste0('Large ', 'SNE ', fit$year_labels)
)

# Pull centroids of grid cells
locg <- as.data.frame(fit$spatial_list$loc_g)

# Mutate kilometers to meters for plotting
locg <- locg %>% 
  mutate(E_km = E_km * 1000,
         N_km = N_km * 1000)

# Convert to sf
locg.sf <- st_as_sf(locg,
                    coords=c('E_km', 'N_km'), 
                    crs="EPSG:26919")

# Add cell ID values
locg.sf$ID <- seq(1:nrow(locg.sf))

# Load contextual spatial data
coast <- st_transform(ecodata::coast, st_crs(locg.sf))
codstox <- st_transform(
  st_read(here('Data/GIS/codstox.shp'), quiet=T),
  st_crs(locg.sf)
)

# Extract area of cells
locatecells <- fit$extrapolation_list$a_el %>% 
  as.data.frame() %>% 
  rename(ALL = V1, EGOM = V2, GBK = V3, SNE = V4, WGOM = V5) %>% 
  mutate(ALL = strip_units(ALL),
         EGOM = strip_units(EGOM),
         GBK = strip_units(GBK),
         SNE = strip_units(SNE),
         WGOM = strip_units(WGOM),
         ID = seq(1:length(ALL)))

# Only save cells within SNE
instrata <- locatecells[locatecells$ALL != 0,]

indexing <- data.frame(
  index=seq(1, 1200, 1),
  names = names(listcells)
)
indexing <- indexing %>% 
  separate(names, into=c('Size', 'Stock', 'Year', 'Season'))
sizes <- c('Small', 'Medium', 'Large')
stocks <- c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')
seasons <- c('Spring', 'Fall')

# Run through all combos
for(i in 1:length(sizes)){
  size <- sizes[i]
  for(j in 1:length(stocks)){
    stock <- stocks[j]
    for(k in 1:length(seasons)){
      season <- seasons[k]
      
      nums <- indexing[indexing$Size == paste0(size) &
                       indexing$Stock == paste0(stock) &
                       indexing$Season == paste0(season),]$index
      
      for(q in 1:length(nums)){
        
        title <- paste0(indexing$Year[nums[q]], ' ',
                        indexing$Season[nums[q]])
        
        yp <- 
          ggplot() +
          # Coast
          geom_sf(data=coast, fill='lightgray', col=NA) +
          # All grid cells
          geom_sf(data=locg.sf[locg.sf$ID %in% instrata$ID,], 
                  alpha=0.2, cex=0.75) +
          # Stock outlines
          geom_sf(data=codstox, fill=NA, col='black') +
          # Cells that combined contain 95% of population, ordered densest to least dense
          geom_sf(data=locg.sf[locg.sf$ID %in% listcells[[nums[q]]],],
                  col='red', alpha=0.6, cex=0.75) +
          # Zoom
          coord_sf(xlim=c(-76, -66),
                   ylim=c(36.75, 44.5),
                   crs="EPSG:4326") +
          # Year
          ggtitle(paste0(title))
        
        if(season == 'Spring'){char = 's'}
        if(season == 'Fall'){char='f'}
        
        ggsave(plot=yp,
               filename=here(paste0('VAST_runs/',
                                    tolower(paste0(size)),
                                    '/Overall_BC/',
                                    stock,
                                    '/l',
                                    char,
                                    'p/',
                                    title,
                                    '.png')),
               width = 6, height=6.5, units='in')
      }
    }
  }
}

## Animate
setwd(here('VAST_runs/large/Overall_BC/ALL/lsp'))
gifski::gifski(png_files = list.files(here('VAST_runs/large/Overall_BC/ALL/lsp')),
               gif_file = 'spring_kernel.gif', delay=0.2, loop=F,
               width = 800, height = 866)
