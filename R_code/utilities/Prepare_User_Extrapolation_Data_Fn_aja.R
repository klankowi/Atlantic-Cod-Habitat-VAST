Prepare_User_Extrapolation_Data_Fn_aja <- function(input_grid, strata.limits = NULL, projargs = NA, zone = NA, flip_around_dateline = TRUE, index_shapes, ...) {
  if (FALSE) {
    # Run make_extrapolation_info_aja first...
    strata.limits <- strata.limits
    input_grid <- input_grid
    projargs <- projargs
    zone <- zone
    flip_around_dateline <- flip_around_dateline
    index_shapes <- index_shapes
  }
  if (is.null(strata.limits)) {
    strata.limits <- data.frame(STRATA = "All_areas")
  }
  message("Using strata ", strata.limits)
  Data_Extrap <- input_grid
  Area_km2_x <- Data_Extrap[, "Area_km2"]
  Tmp <- cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
  if ("Depth" %in% colnames(Data_Extrap)) {
    Tmp <- cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
  }
  a_el <- as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
  for (l in 1:ncol(a_el)) {
    a_el[, l] <- match_strata_fn_aja(points = Tmp, strata_dataframe = strata.limits[l, , drop = FALSE], index_shapes = index_shapes[index_shapes$Region == as.character(strata.limits[l, , drop = FALSE]), ])
    a_el[, l] <- ifelse(is.na(a_el[, l]), 0, Area_km2_x)
  }
  tmpUTM <- project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
  Data_Extrap <- cbind(Data_Extrap, Include = 1)
  if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
    Data_Extrap[, c("E_km", "N_km")] <- tmpUTM[, c("X", "Y")]
  } else {
    Data_Extrap <- cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
  }
  Return <- list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, Area_km2_x = Area_km2_x)
  return(Return)
}