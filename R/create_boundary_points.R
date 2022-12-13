#' Create points around coastal boundary
#'
#'
#' @param rsf An sf object defining the geographic region land area
#' @param grd An sf grid object
#'
#' @import sf
#' @importFrom units set_units
#'
#' @return An sf object
#'
#' @export
create.boundary.points <- function(rsf, grd) {

  # Cast polygon of coast to linestring
  ls_sf <- st_cast(st_cast(rsf, "POLYGON"), "LINESTRING")

  # Sample points on line
  lp <- st_line_sample(ls_sf,
                       density = units::set_units(0.5, 1/km), type= "regular")
  lp <- st_cast(lp, "POINT")
  df <-  as.data.frame(st_coordinates(lp))
  df$site_major <- 1:nrow(df)
  df <- df[df$X < max(df$X), ]
  df <- st_as_sf(df, coords = c("X", "Y"), crs = st_crs(grd))
  df

}
