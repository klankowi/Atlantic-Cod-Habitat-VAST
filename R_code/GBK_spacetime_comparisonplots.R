rm(list=ls())

# Load libraries
library(tidyverse)
library(sf)
library(here)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
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

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Load survey location data
# Data are publicly available at this repository:
# https://github.com/klankowi/Atlantic-Cod-Habitat-VAST
survs <- read.csv(here('Data/Survey_Data/Survey_Data.csv'))

# Strip to DFO and NEFSC BTS
table(survs$SURVEY)
survs <- survs[survs$SURVEY == 'DFO Trawl' |
               survs$SURVEY == 'NEFSC BTS',]

# Strip to George's Bank
table(survs$STOCK)
survs <- survs[survs$STOCK == 'GBK',]

# Strip to only necessary data
survs <- dplyr::select(survs,
                       SURVEY, SEASON, YEAR,
                       DATE, LAT, LON)

# Data shaping
survs$SURVEY <- as.factor(survs$SURVEY)
survs$SEASON <- as.factor(survs$SEASON)
survs$DATE <- as.POSIXct(survs$DATE,
                                format= "%m/%d/%Y")
survs$DATE <- format(survs$DATE,
                     format = "%Y-%m-%d")
survs$DATE <- as.POSIXct(as.character(survs$DATE),
                         format= "%Y-%m-%d")
survs <- survs[with(survs, order(DATE)),]
survs$SERIES <- paste(survs$SURVEY, survs$SEASON)
survs$SERIES <- as.factor(survs$SERIES)
summary(survs)

# Pull series of years with sampling data
years.used <- unique(survs$YEAR)

# Split data by series
series.list <- split(survs, f=survs$SERIES)

# Loop through series
for(i in 1:length(series.list)){
  # Pull years without data
  years.notused <- years.used[years.used %notin% series.list[[i]]$YEAR]
  # Create temporary dataframe with those years to fill time series
  # If you do not do this, plot will connect points over skipped years,
  # this makes it look like there are samples in those skipped years.
  tempdf <- data.frame(
    SURVEY = series.list[[i]]$SURVEY[1],
    SEASON = series.list[[i]]$SEASON[1],
    YEAR = years.notused,
    DATE = NA,
    LAT = NA,
    LON = NA,
    SERIES = series.list[[i]]$SERIES[1]
  )
  # Merge and re-order
  series.list[[i]] <- rbind(series.list[[i]], tempdf)
  series.list[[i]] <- series.list[[i]][with(series.list[[i]], order(
    YEAR, DATE)
  ),]
  rownames(series.list[[i]]) <- NULL
  # Split results by year to calculate min, max, and median date sampled per year
  series.list[[i]] <- split(series.list[[i]], f=series.list[[i]]$YEAR)
  for(j in 1:length(series.list[[i]])){
    series.list[[i]][[j]]$MINDATE <- min(series.list[[i]][[j]]$DATE)
    series.list[[i]][[j]]$MAXDATE <- max(series.list[[i]][[j]]$DATE)
    series.list[[i]][[j]]$MEDDATE <- median(series.list[[i]][[j]]$DATE)
  }
  # Rebind years to single series
  series.list[[i]] <- do.call(rbind, series.list[[i]])
}
# Rebind series to dataframe, reorder
survs2 <- do.call(rbind, series.list)
rownames(survs2) <- NULL
survs2 <- survs2[with(survs2, order(YEAR, DATE)),]
head(survs2)

# Force min, max, median dates to have same arbitrary year (for plotting)
survs2$MEDDATE <- as.POSIXct(paste0('2000-', 
                                     substr(survs2$MEDDATE, 
                                            start = 6, stop=10)),
                              format='%Y-%m-%d')
survs2$MINDATE <- as.POSIXct(paste0('2000-', 
                                     substr(survs2$MINDATE, 
                                            start = 6, stop=10)),
                              format='%Y-%m-%d')
survs2$MAXDATE <- as.POSIXct(paste0('2000-', 
                                     substr(survs2$MAXDATE, 
                                            start = 6, stop=10)),
                              format='%Y-%m-%d')

# Plot and save time metric
time.met <- ggplot(survs2) +
  geom_line(aes(x=YEAR, y=MEDDATE, col=SERIES)) +
  geom_point(aes(x=YEAR, y=MEDDATE, col=SERIES)) +
  geom_errorbar(aes(x=YEAR, ymin=MINDATE, ymax=MAXDATE, col=SERIES),
                size = 0.2) +
  xlab('Year') + ylab('Survey Dates') +
  scale_y_datetime(limits = as.POSIXct(c("2000-01-01", "2000-12-31")),
                   date_labels = "%b",
                   date_breaks = "2 months")
time.met
ggsave(time.met,
       filename = paste0(here(), "/Plot_output/GBK_time_comparison.png"),
       device='png')  

# Convert to sf for spatial plotting
survs2 <- survs2[!is.na(survs2$LAT),]
survs2 <- survs2[!is.na(survs2$LON),]
survs.sf <- st_as_sf(survs2, coords=c('LON', 'LAT'))
st_crs(survs.sf) <- "EPSG:4326"

# Load coast
coast <- ecodata::coast
coast <- st_transform(coast, st_crs(survs.sf))

# Data shaping
survs.sf$`NEFSC BTS` <- survs.sf$YEAR
survs.sf$`DFO Trawl` <- survs.sf$YEAR

# Remove years prior to 1987 (no temporal overlap)
survs.sf <- subset(survs.sf, YEAR >=1987)

# Facet plot of spring surveys by year
spring.splityear <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=survs.sf[survs.sf$SERIES == 'DFO Trawl SPRING' |
                          survs.sf$SERIES == "NEFSC BTS SPRING",], 
          aes(col=SERIES),
          cex=0.8) +
  coord_sf(xlim=c(st_bbox(survs.sf)[1], st_bbox(survs.sf)[3]),
           ylim=c(st_bbox(survs.sf)[2], st_bbox(survs.sf)[4])) +
  facet_wrap(vars(YEAR))
spring.splityear
ggsave(spring.splityear,
       filename=paste0(here(), "/Plot_output/GBK_SPRING_SPLITYEAR.png"),
       device='png')

# Plot of overlapping spring surveys for all years together
all.space <- ggplot() +
  geom_sf(data=coast, fill='gray') +
  geom_sf(data=survs.sf[survs.sf$SERIES == 'DFO Trawl SPRING',], 
          aes(col=`DFO Trawl`),
          cex=0.8) +
  ggnewscale::new_scale_color() +
  geom_sf(data=survs.sf[survs.sf$SERIES == 'NEFSC BTS SPRING',],
          aes(col=`NEFSC BTS`),
          cex=0.8) +
  scale_color_gradient(low = "#8B0001",
                       high = "#F6BDC0") +
  coord_sf(xlim=c(st_bbox(survs.sf)[1], st_bbox(survs.sf)[3]),
           ylim=c(st_bbox(survs.sf)[2], st_bbox(survs.sf)[4]))
all.space
ggsave(all.space,
       filename=paste0(here(), "/Plot_output/GBK_SPACE_YEARCOMPARISON.png"),
       device="png")
