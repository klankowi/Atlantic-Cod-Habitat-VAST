match_strata_fn_aja <- function(points, strata_dataframe, index_shapes) {
  if (FALSE) {
    points <- Tmp
    l <- 1
    strata_dataframe <- strata.limits[l, , drop = FALSE]
    index_shapes <- index_shapes
  }
  if (is.null(index_shapes)) {
    # Default all strata
    match_latitude_TF <- match_longitude_TF <- match_depth_TF <- rep(TRUE, nrow(strata_dataframe))
    if (all(c("south_border", "north_border") %in% names(strata_dataframe))) {
      match_latitude_TF <- as.numeric(x["BEST_LAT_DD"]) > strata_dataframe[, "south_border"] & as.numeric(x["BEST_LAT_DD"]) <= strata_dataframe[, "north_border"]
    }
    if (all(c("west_border", "east_border") %in% names(strata_dataframe))) {
      match_longitude_TF <- as.numeric(x["BEST_LON_DD"]) > strata_dataframe[, "west_border"] & as.numeric(x["BEST_LON_DD"]) <= strata_dataframe[, "east_border"]
    }
    if (all(c("shallow_border", "deep_border") %in% names(strata_dataframe))) {
      match_depth_TF <- as.numeric(x["BEST_DEPTH_M"]) > strata_dataframe[, "shallow_border"] & as.numeric(x["BEST_DEPTH_M"]) <= strata_dataframe[, "deep_border"]
    }
    # Return stuff
    Char <- as.character(strata_dataframe[match_latitude_TF & match_longitude_TF & match_depth_TF, "STRATA"])
    return(ifelse(length(Char) == 0, NA, Char))
  }
  
  # Andrew edit...
  if (!is.null(index_shapes)) {
    Tmp_sf <- data.frame(points) %>%
      st_as_sf(., coords = c("BEST_LON_DD", "BEST_LAT_DD"), crs = st_crs(index_shapes), remove = FALSE)
    match_shape <- Tmp_sf %>%
      st_join(., index_shapes, join = st_within) %>%
      mutate(., "Row_ID" = seq(from = 1, to = nrow(.))) %>%
      st_drop_geometry() %>%
      dplyr::select(., Region) %>%
      as.vector()
    return(match_shape)
  }
}