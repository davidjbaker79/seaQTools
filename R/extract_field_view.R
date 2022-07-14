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
#' @importFrom dplyr left_join
#' @importFrom dplyr bind_rows
#' @import concaveman
#' @importFrom Hmisc wtd.var
#'
#' @return A data.frame.
extract.field.view <- function(n,
                               sites_major_sf,
                               height = 30,
                               region_sf,
                               modGrid,
                               env_dat) {
  #- Region buffer
  region_bf <- st_buffer(region_sf, dist = 30000)

  #- Use height of site to calculate the horizon dist
  hz <- horizon(height)

  #- Subset to a single site
  p1 <- sites_major_sf[n, ]
  p1$lon <- st_coordinates(p1)[, 1]
  p1$lat <- st_coordinates(p1)[, 2]

  #- Crop out a single site
  p1_b <- st_buffer(p1, hz)

  #- Crop coast
  pc <- st_difference(p1_b, region_sf)
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

  #-- Extract environmental data i = 2
  zones <- c(1000, 2000, 5000, 10000)
  zone_nam <- c("0_1km", "1_2km", "2_5km", "5_10km")
  fv_grd_id_i <- lapply(1:4, function(i) {
    if (i == 1) {
      buff_i <- st_buffer(p1, dist = zones[i])
    } else {
      buff_1 <- st_buffer(p1, dist = zones[i - 1])
      buff_2 <- st_buffer(p1, dist = zones[i])
      buff_i <- st_difference(buff_2, buff_1)
    }
    zone_i <- st_intersection(buff_i, fieldview)
    grd_z <- st_intersection(modGrid, zone_i)
    grd_z$area_fv <-
      as.numeric(round(st_area(grd_z) / 10 ^ 6, 2))
    grd_z$dist_fv <-
      as.numeric(round(st_distance(p1, grd_z, by_element = T)))
    grd_z$position_dist_code <- zone_nam[i]
    grd_z$viewArea <- as.numeric(st_area(zone_i) / 10 ^ 6)
    grd_z <-
      grd_z[, c("id",
                "site_major",
                "area_fv",
                "dist_fv",
                "position_dist_code",
                "viewArea")]
    grd_z <- st_drop_geometry(grd_z)
    grd_z <- left_join(grd_z, env_dat)
  })
  fv_grd_id <- do.call(bind_rows, fv_grd_id_i)

}
