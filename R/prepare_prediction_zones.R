#' Prepare grid/zones for making predictions from SDM
#'
#' This function prepares the grids to intersect the environmental data with an
#' original resolution of 1km and to summarise at a coarser resolution,
#' suitable for making predictions, i.e. that better match the resolution of
#' the fitted model. The grids are also clipped at the distance bands
#' (i.e. 0_1km, 1_2km ...). The environmental variables can then be aggregated
#' to the coarse grids using dplyr functions post hoc, and weights equal to the
#' amount of the 1km cell that overlaps the coarser grid can be applied.
#'
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param coarse_res A numeric giving the horizontal resolution of coarse grid.
#'
#' @import sf
#'
#' @return A data.frame.
#'
#' @export
prepare.pred.zones <-
  function(modGrid,
           region_sf,
           coarse_res = 5000) {
    #- Coarse grid
    modGrid_L <- makeGrid(region_sf, coarse_res, 15000)
    names(modGrid_L)[which(names(modGrid_L) == "id")] <- "id_zone"

    #-- Extract environmental data i = 2
    zones <- c(1000, 2000, 5000, 10000)
    zone_nam <- c("0_1km", "1_2km", "2_5km", "5_10km")
    fv_grd_id_i <- lapply(1:4, function(i) {
      if(i == 1) {
        buff_i <- st_buffer(region_sf, dist = zones[i])
      } else {
        buff_1 <- st_buffer(region_sf, dist = zones[i-1])
        buff_2 <- st_buffer(region_sf, dist = zones[i])
        buff_i <- st_difference(buff_2,buff_1)
      }
      zone_i <- st_intersection(modGrid, buff_i)
      grd_z <- st_intersection(modGrid_L, zone_i)
      grd_z$area_fv <-
        as.numeric(round(st_area(grd_z) / 10 ^ 6, 2))
      grd_z$position_dist_code <- zone_nam[i]
      grd_z <-
        grd_z[, c("id",
                  "id_zone",
                  "X",
                  "Y",
                  "area_fv",
                  "position_dist_code")]
      grd_z <- st_drop_geometry(grd_z)
      grd_z <- grd_z[grd_z$area_fv > 0,]
    })
    fv_grd_id <- do.call(bind_rows, fv_grd_id_i)

  }
