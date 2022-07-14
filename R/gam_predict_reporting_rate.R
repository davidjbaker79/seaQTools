#' Predict reporting rate to distance bands / time periods
#'
#' @param m A model object produced by mgcv::gam.
#' @param x A data.frame with environmental data for predictions.
#' @param region_sf
#' @param position_dist_zone
#' @param coarse_res
#'
#' @import sf
#' @importFrom dplyr left_join
#'
#' @return A data.frame with spatial predictions.
gam.predict.rr <-
  function(m,
           x,
           region_sf,
           position_dist_zone,
           coarse_res = 5000) {
    za <- lapply(position_dist_zone, function(i) {
      #- Filter to zone
      ez <- x[x$position_dist_code == i,]

      #- Coarse grid
      if (i == "0_1km") {
        d_m <- 1000
      } else if (i == "1_2km") {
        d_m <- 2000
      } else if (i == "2_5km") {
        d_m <- 5000
      } else if (i == "5_10km") {
        d_m <- 10000
      }
      region_bf <- st_buffer(region_sf, dist = d_m)
      modGrid_L <- makeGrid(region_sf, coarse_res, 15000)
      names(modGrid_L)[1] <- "id_zone"

      #- Create zones
      b1 <- st_buffer(region_sf, dist = 1000)
      b2 <- st_buffer(region_sf, dist = 2000)
      b5 <- st_buffer(region_sf, dist = 5000)
      b10 <- st_buffer(region_sf, dist = 10000)

      #- Grid zones
      if (i == "0_1km") {
        zo <- b1
      } else if (i == "1_2km") {
        zo <- st_difference(b2, b1)
      } else if (i == "2_5km") {
        zo <- st_difference(b5, b2)
      } else if (i == "5_10km") {
        zo <- st_difference(b10, b5)
      }
      zo <-  st_intersection(modGrid_L,   zo)
      modGrid_L <- st_intersection(modGrid_L, zo)

      # Predict model to regional data
      fit <-
        predict(m, newdata = ez, type = "response")
      pred <- cbind.data.frame(ez, fit)
      m_pred_sf <- left_join(modGrid_L, pred)
      m_pred_sf <- na.omit(m_pred_sf)

    })
    za <- do.call(rbind, za)
  }
