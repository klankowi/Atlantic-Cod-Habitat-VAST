vast_make_extrap_grid <- function(region_shapefile, index_shapes, cell_size) {
  
  # For debugging
  if (FALSE) {
    region_shapefile<- st_read(here::here("data/supporting/region_shapefile/full_survey_region.shp"))  
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
    cell_size <- NULL
    grid_dim_km = c(10, 10)
    maximum_distance_from_sample = NULL
    
    region_shapefile <- menh_chull
    index_shapes <- menh_regions_out
    cell_size <- 5000
  }
  
  # Transform crs of shapefile to common WGS84 lon/lat format.
  region_wgs84 <- region_shapefile
  
  # Get UTM zone
  lon <- sum(st_bbox(region_wgs84)[c(1, 3)]) / 2
  utm_zone <- floor((lon + 180) / 6) + 1
  
  # Transform to the UTM zone
  crs_utm <- st_crs(paste0("+proj=utm +zone=", utm_zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
  region_utm <- st_transform(region_wgs84, crs = crs_utm)
  
  # Make extrapolation grid with sf
  region_grid <- st_as_sf(st_make_grid(region_utm, cellsize = cell_size, what = "centers"), crs = crs_utm)
  
  # Now get only the points that fall within the shape polygon
  points_keep <- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, region_utm, sparse = FALSE))
  region_grid <- region_grid %>%
    mutate(., "in_poly" = st_intersects(region_grid, region_utm, sparse = FALSE)) %>%
    filter(., in_poly == TRUE)
  
  # Convert back to WGS84 lon/lat, as that is what VAST expects.
  extrap_grid <- region_grid %>%
    st_transform(., crs = 4326)
  
  # Adding in the a strata/region component for stratified abundance. This will depend on index_shapes input.
  if (!is.null(index_shapes)) {
    extrap_grid <- extrap_grid %>%
      st_join(., index_shapes, join = st_within) %>%
      mutate(.,
             "Lon" = as.numeric(st_coordinates(.)[, 1]),
             "Lat" = as.numeric(st_coordinates(.)[, 2])
      ) %>%
      st_drop_geometry() %>%
      dplyr::select(., Lon, Lat, Region) %>%
      mutate(.,
             Area_km2 = ((cell_size / 1000)^2),
             STRATA = factor(Region, levels = index_shapes$Region, labels = index_shapes$Region)
      )
  } else {
    extrap_grid <- extrap_grid %>%
      mutate(.,
             "Lon" = as.numeric(st_coordinates(.)[, 1]),
             "Lat" = as.numeric(st_coordinates(.)[, 2])
      ) %>%
      st_drop_geometry() %>%
      dplyr::select(., Lon, Lat) %>%
      mutate(., Area_km2 = ((cell_size / 1000)^2))
  }
  
  # Return it
  return(extrap_grid)
}