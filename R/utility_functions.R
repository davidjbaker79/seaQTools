#' Make hexagonal grid for national grids
#'
#' Analysis is conducted on hexagonal grids and this function creates the grid
#' with a specified area.
#'
#' @param region_sf An sf object defining the spatial region.
#' @param cell_size A numeric value specifying the size of the grid cell (m)
#' @param bufferDist distance from coast (m)
#'
#' @import sf
#'
#' @return An sf object with 'id' used as a unique cell id.
makeGrid <- function(region_sf, cell_size, bufferDist) {
  region_bf <-
    st_buffer(region_sf, dist = bufferDist)
  modGrid <-
    st_make_grid(
      region_bf,
      cellsize = cell_size,
      offset = c(0, 0),
      square = TRUE
    )
  modGrid <- st_as_sf(modGrid)
  modGrid$id <- as.numeric(rownames(modGrid))
  modGrid <- st_difference(modGrid, region_sf)
  modGrid$area <- as.numeric(round(st_area(modGrid) / 10 ^ 6, 2))
  modGrid <- cbind(modGrid, st_coordinates(st_centroid(modGrid)))

}
