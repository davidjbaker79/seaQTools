#' Prepare environmental data for model building
#'
#' Function crops out viewable areas for each site, including gridded
#' environmental variables, and returns an
#'
#' @param x An sf object of the major site locations
#' @param y An sf object of gridded environmental data
#' @param rsf An sf object defining the geographic region land area
#' @param rbf An sf object defining the geographic region ocean area
#' @param grd An sf grid object
#' @param h A numeric giving the height of observation
#'
#' @import sf
#'
#' @return An sf object
#'
#' @export
site.enviro.setup <- function(x, y, rsf, grd, h = 30) {
  rbf <- st_buffer(rsf, dist = 15000)
  ms <- lapply(1:nrow(x), function(i) {
    print(i)
    fv_grd_i <- extract.field.view(
      n = i,
      sites_major_sf = x,
      height = h,
      region_sf = rsf,
      region_bf = rbf,
      modGrid = grd
    )
    yi <- st_intersects(y, fv_grd_i)
    if (any(lengths(yi) > 0)) {
      yi <- y[lengths(yi) > 0, ]
      yi <- prepare.enviro.view.area(yi)
      yi <- cbind(yi, st_coordinates(x[i, ]))
      yi$viewArea <- as.numeric(st_area(fv_grd_i) / 10 ^ 6)
      yi$site_major <- x$site_major[i]
      yi
    }
  })
  ms_i <- do.call(rbind, ms)
}
