vast_mesh_to_sf <- function(vast_fit, crs_transform = "+proj=longlat +datum=WGS84 +no_defs") {
  if (FALSE) {
    tar_load(vast_fit)
    crs_transform <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  require(sp) || stop("Install sp, else thine code shan't work for thee")
  require(sf) || stop("Install sf or this code will be a mess")
  require(INLA) || stop("You need the R-INLA package for this, note that it's not crantastic...
                        install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  
  # Get the extrapolation mesh information from the vast_fitted object
  mesh <- vast_fit$spatial_list$MeshList$anisotropic_mesh
  mesh["crs"] <- vast_fit$extrapolation_list$projargs
  
  # Grab the CRS if it exists, NA is fine (NULL spits a warning, but is also fine)
  crs <- sp::CRS(mesh$crs)
  
  # Make sure the CRS isn't a geocentric one, which is won't be if yo look up geocentric..
  # isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  isgeocentric <- inla.crs_is_geocent(mesh$crs)
  # Look up geo-centric coordinate systems, nothing we'll need to worry about, but stop if so
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp and sf' don't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before calling inla.mesh2sf()."
    ))
  }
  # This pulls out from the mesh the triangles as polygons, this was the piece I couldn't figure out.
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(
      lapply(
        1:nrow(mesh$graph$tv),
        function(x) {
          tv <- mesh$graph$tv[x, , drop = TRUE]
          Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)], 1:2, drop = FALSE])), ID = x)
        }
      ),
      proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  
  # This one is easy, just grab the vertices (points)
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  # Make these sf objects
  triangles <- st_as_sf(triangles)
  vertices <- st_as_sf(vertices)
  
  # Transform?
  if (!is.null(crs_transform)) {
    triangles <- st_transform(triangles, crs = crs_transform)
    vertices <- st_transform(vertices, crs = crs_transform)
  }
  
  # Add your output list.
  return_sf <- list(triangles = triangles, vertices = vertices)
  return(return_sf)
}