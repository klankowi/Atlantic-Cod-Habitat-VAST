# function to make rasters into data frame for merge with survey
# needs long df with date split to year, month, day, lat, lon, sst
# crop to NEUS extent
# from https://towardsdatascience.com/transforming-spatial-data-to-tabular-data-in-r-4dab139f311f

raster_to_bt <- function(brick){
  
  sstdf <- as.data.frame(raster::rasterToPoints(brick, spatial = TRUE))
  sstdf <- sstdf %>%
    dplyr::rename(Lon = x,
                  Lat = y) %>%
    tidyr::pivot_longer(cols = starts_with("X"),
                        #names_to = c("year", "month", "day"),
                        names_prefix = "X",
                        #names_sep = "\\.",
                        values_to = "bt") #%>% 
  sstdf$year = i
  sstdf$name = as.numeric(sstdf$name)
  sstdf$month = lubridate::month(as.Date(sstdf$name,
                                         origin = paste0((i-1), "-12-31")))
  sstdf$day = lubridate::day(as.Date(sstdf$name,
                                     origin= paste0((i-1), "-12-31")))
  sstdf$name <- NULL
  
  sstdf <- sstdf[with(sstdf, order(year, month, day)),]
  
  return(sstdf)
}