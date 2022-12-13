#' Prepared gridded environmental data in viewable area.
#'
#' Summarise data in viewable area.
#'
#' @param x A sf object of environmental data (i.e. 1km gridded data).
#'
#' @import data.table
#' @importFrom stats IQR
#'
#' @return A data.frame.
#'
#' @export
prepare.enviro.view.area <- function(x) {

  z <- setDT(x)[, .(
    sst_m = mean(sst, na.rm = TRUE),
    sst_iqr = IQR(sst, na.rm = TRUE),
    sst_min = min(sst),
    sst_max = max(sst),
    npp_m = mean(npp, na.rm = TRUE),
    npp_sd = IQR(npp, na.rm = TRUE),
    npp_min = min(npp),
    npp_max = max(npp),
    tri_m = mean(TRI, na.rm = TRUE),
    tri_sd = IQR(TRI, na.rm = TRUE),
    tri_min = min(TRI),
    tri_max = max(TRI),
    depth_m = mean(depth, na.rm = TRUE),
    depth_sd = IQR(depth, na.rm = TRUE),
    depth_max = min(depth),
    depth_min = max(depth)
  ), by = c("Year")]
  z <- na.omit(z)

}

