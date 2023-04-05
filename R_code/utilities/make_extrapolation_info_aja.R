make_extrapolation_info_aja <- function(Region, projargs = NA, zone = NA, strata.limits = data.frame(STRATA = "All_areas"), create_strata_per_region = FALSE, max_cells = NULL, input_grid = NULL, observations_LL = NULL, grid_dim_km = c(2, 2), maximum_distance_from_sample = NULL, grid_in_UTM = TRUE, grid_dim_LL = c(0.1, 0.1), region = c("south_coast", "west_coast"), strata_to_use = c("SOG", "WCVI", "QCS", "HS", "WCHG"), epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[1], survey = "Chatham_rise", surveyname = "propInWCGBTS", flip_around_dateline, nstart = 100, area_tolerance = 0.05, backwards_compatible_kmeans = FALSE, DirPath = paste0(getwd(), "/"), index_shapes, ...) {
  if (FALSE) {
    # First run fit_model_aja...
    Region <- settings$Region
    projargs <- NA
    zone <- settings$zone
    strata.limits <- settings$strata.limits
    create_strata_per_region <- FALSE
    max_cells <- settings$max_cells
    input_grid <- input_grid
    observations_LL <- NULL
    grid_dim_km <- settings$grid_size_km
    maximum_distance_from_sample <- NULL
    index_shapes <- index_shapes
  }
  if (is.null(max_cells)) {
    max_cells <- Inf
  }
  for (rI in seq_along(Region)) {
    Extrapolation_List <- NULL
    if (tolower(Region[rI]) == "user") {
      if (is.null(input_grid)) {
        stop("Because you're using a user-supplied region, please provide 'input_grid' input")
      }
      if (!(all(c("Lat", "Lon", "Area_km2") %in% colnames(input_grid)))) {
        stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
      }
      if (missing(flip_around_dateline)) {
        flip_around_dateline <- FALSE
      }
      Extrapolation_List <- Prepare_User_Extrapolation_Data_Fn_aja(strata.limits = strata.limits, input_grid = input_grid, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, index_shapes = index_shapes, ...)
    }
    if (is.null(Extrapolation_List)) {
      if (is.null(observations_LL)) {
        stop("Because you're using a new Region[rI], please provide 'observations_LL' input with columns named `Lat` and `Lon`")
      }
      if (missing(flip_around_dateline)) {
        flip_around_dateline <- FALSE
      }
      Extrapolation_List <- Prepare_Other_Extrapolation_Data_Fn(strata.limits = strata.limits, observations_LL = observations_LL, grid_dim_km = grid_dim_km, maximum_distance_from_sample = maximum_distance_from_sample, grid_in_UTM = grid_in_UTM, grid_dim_LL = grid_dim_LL, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (rI == 1) {
      Return <- Extrapolation_List
    } else {
      Return <- combine_extrapolation_info(Return, Extrapolation_List, create_strata_per_region = create_strata_per_region)
    }
  }
  if (max_cells < nrow(Return$Data_Extrap)) {
    message("# Reducing extrapolation-grid from ", nrow(Return$Data_Extrap), " to ", max_cells, " cells for Region(s): ", paste(Region, collapse = ", "))
    loc_orig <- Return$Data_Extrap[, c("E_km", "N_km")]
    loc_orig <- loc_orig[which(Return$Area_km2_x > 0), ]
    Kmeans <- make_kmeans(n_x = max_cells, loc_orig = loc_orig, nstart = nstart, randomseed = 1, iter.max = 1000, DirPath = DirPath, Save_Results = TRUE, kmeans_purpose = "extrapolation", backwards_compatible_kmeans = backwards_compatible_kmeans)
    Kmeans[["cluster"]] <- RANN::nn2(data = Kmeans[["centers"]], query = Return$Data_Extrap[, c("E_km", "N_km")], k = 1)$nn.idx[, 1]
    aggregate_vector <- function(values_x, index_x, max_index, FUN = sum) {
      tapply(values_x, INDEX = factor(index_x, levels = 1:max_index), FUN = FUN)
    }
    a_el <- matrix(NA, nrow = max_cells, ncol = ncol(Return$a_el))
    for (lI in 1:ncol(Return$a_el)) {
      a_el[, lI] <- aggregate_vector(values_x = Return$a_el[, lI], index_x = Kmeans$cluster, max_index = max_cells)
    }
    Area_km2_x <- aggregate_vector(values_x = Return$Area_km2_x, index_x = Kmeans$cluster, max_index = max_cells)
    Include <- aggregate_vector(values_x = Return$Data_Extrap[, "Include"], index_x = Kmeans$cluster, max_index = max_cells, FUN = function(vec) {
      any(vec > 0)
    })
    lonlat_g <- project_coordinates(X = Kmeans$centers[, "E_km"], Y = Kmeans$centers[, "N_km"], projargs = "+proj=longlat +ellps=WGS84", origargs = Return$projargs)
    Data_Extrap <- cbind(Lon = lonlat_g[, 1], Lat = lonlat_g[, 2], Include = Include, Kmeans$centers)
    Return <- list(a_el = a_el, Data_Extrap = Data_Extrap, zone = Return$zone, projargs = Return$projargs, flip_around_dateline = Return$flip_around_dateline, Area_km2_x = Area_km2_x)
  }
  if (length(Region) > 1 & create_strata_per_region == TRUE) {
    Return$a_el <- cbind(Total = rowSums(Return$a_el), Return$a_el)
  }
  class(Return) <- "make_extrapolation_info"
  return(Return)
}