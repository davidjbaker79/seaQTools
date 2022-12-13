#' Determine the viewable area from each site and extract enviro data
#'
#' @param n A numeric giving row of sites_major_sf.
#' @param sites_major_sf A dataframe giving coordinates of each major site.
#' @param height A numeric height above sea level.
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#'
#' @import sf
#' @importFrom geosphere horizon
#' @import concaveman
#'
#' @return A data.frame.
#'
#' @export
extract.field.view <- function(n,
                               sites_major_sf,
                               height = 30,
                               region_sf,
                               region_bf,
                               modGrid) {
  #- Use height of site to calculate the horizon distance
  hz <- horizon(height)

  #- Subset to a single site
  p1 <- sites_major_sf[n, ]
  p1$lon <- st_coordinates(p1)[, 1]
  p1$lat <- st_coordinates(p1)[, 2]

  #- Crop out a single site
  p1_b <- st_buffer(p1, hz)

  #- Crop coast
  st_agr(p1_b) <- "constant"
  pc <- st_difference(p1_b, region_sf)
  st_agr(pc) <- "constant"
  pc <-  st_cast(pc, "POLYGON")
  pc <- pc[which(st_area(pc) == max(st_area(pc))), ]

  #- Points on buffer edge
  poly_points <- st_segmentize(p1_b, dfMaxLength = 1000)
  poly_points <- as.data.frame(st_coordinates(poly_points))[, 1:2]
  poly_points <- st_as_sf(poly_points,
                          coords = c("X", "Y"),
                          crs = st_crs(p1_b))

  #- Extract coordinates to points on buffer edge
  df <- st_coordinates(poly_points)

  #- Create lines between site and edge
  rows <- split(df, seq(nrow(df)))
  lines <- lapply(rows, function(row) {
    lmat <-
      matrix(unlist(c(row, st_coordinates(p1))), ncol = 2, byrow = TRUE)
    st_linestring(lmat)
  })
  lines <- st_sfc(lines,
                  crs = st_crs(p1_b))
  lines_sf <- st_sf('geometry' = lines)

  #- crop coast
  lseg <- st_difference(lines_sf, region_sf)
  lseg <- st_cast(lseg, "MULTILINESTRING")
  lseg <- st_cast(lseg, "LINESTRING")
  int <- lengths(st_intersects(lseg, st_buffer(p1, 1000))) > 0
  lseg <- lseg[int, ]

  #- Convert linestring to polygon
  fieldview <- concaveman(lseg)

}
