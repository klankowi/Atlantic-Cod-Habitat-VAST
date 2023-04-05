# function to make rasters into data frame for merge with survey
# needs long df with date split to year, month, day, lat, lon, sst
# crop to NEUS extent
# from https://towardsdatascience.com/transforming-spatial-data-to-tabular-data-in-r-4dab139f311f

raster_to_sstdf <- function(brick,
                            rotate=TRUE){
  
  if(rotate) brick_r <- raster::rotate(brick)
  brick_r <- raster::crop(brick_r, raster::extent(-77,-65,35,45))
  sstdf <- as.data.frame(raster::rasterToPoints(brick_r, spatial = TRUE))
  sstdf <- sstdf %>%
    dplyr::rename(Lon = x,
                  Lat = y) %>%
    tidyr::pivot_longer(cols = starts_with("X"),
                        names_to = c("year", "month", "day"),
                        names_prefix = "X",
                        names_sep = "\\.",
                        values_to = "sst",
    )
  return(sstdf)
}