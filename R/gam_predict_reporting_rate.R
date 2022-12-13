#' Predict reporting rate
#'
#' Predict reporting rate to new data.
#'
#' @param m A model object.
#' @param nw A data.frame with new data for predictions
#'
#' @import mgcv
#' @import sf
#'
#' @return A data.frame.
#'
#' @export
gam.predict.rr <- function(m, nw){

  x <-
    cbind(nw,
          fit = predict(m, newdata= nw, type = "response"))
  x <- st_as_sf(x, coords = c("X", "Y"), crs = "epsg:27700")
  x <- cbind(x, st_coordinates(x))

}
