#' Extract and grid VGPM data from Oregon State University
#'
#' Extract and grid VGPM data from
#' http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.vgpm.m.chl.m.sst.php
#'
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#' @param vgpm_path The file path to the vgpm data.
#' @param yr A numeric giving the year of data to extract.
#' @param mth A numeric giving the month (1-12) of data to extract.
#'
#' @import terra
#' @import sf
#' @import data.table
#'
#' @return A data.frame with monthly sst per grid cell.
#'
#' @export
vgpm2grid <- function(region_sf, modGrid, vgpm_path, mth, yr) {

  #- Load vgpm data
  vgpm_dat <- terra::rast(vgpm_path)
  ext(vgpm_dat) <- c(-180, 180, -90, 90)
  vgpm_dat[vgpm_dat < 0] <- NA

  #- vector of months mth_v <- rep(seq(1,12),2)
  if(mth %in% 1:9) mth <- paste0("0", mth)

  #- Buffer and transform region
  region_bf <- st_buffer(region_sf, dist = 30000)
  region_bf <- st_transform(region_bf, crs = "epsg:4326")
  vgpm_grd_c <- crop(vgpm_dat, region_bf)
  vgpm_grd_t <-  project(vgpm_grd_c, "epsg:27700")

  #- Dissaggrate so that missing values can be interoplate later
  vgpm_grd_t <- disagg(vgpm_grd_t, 25)

  #- Rasterize the sw_bf data to the same grid as the gebco data
  region_bf_t <-  st_transform(region_bf, crs = "epsg:27700")
  vgpm_bf_mask <- rasterize(vect(region_bf_t), vgpm_grd_t)

  #- Smooth values
  vgpm_bf <- focal(vgpm_grd_t, w = 55, fun = "mean",  na.rm=TRUE)
  names(vgpm_bf) <- "npp"

  #- Now mask out the data beyond Xkm offshore
  vgpm_bf <- mask(vgpm_bf, vgpm_bf_mask, inverse = FALSE)

  #- Now mask out land
  sw_mask <- rasterize(vect(region_sf), vgpm_grd_t)
  vgpm_bf <- mask(vgpm_bf, sw_mask, inverse = TRUE)

  #- Grid data
  npp_i <- extract(vgpm_bf, vect(modGrid))
  npp_i <- as.data.table(npp_i)
  npp_i <- npp_i[, mean(.data$npp, na.rm=TRUE), by = list("ID")]
  npp_i <- cbind(modGrid, npp_i)
  npp_i <- st_drop_geometry(npp_i)
  npp_i <- npp_i[,c("id", "area", "X", "Y", "V1")]
  names(npp_i)[5] <- "npp"
  npp_i <- npp_i[!is.nan(npp_i$npp),]
  npp_i$mth <- mth
  npp_i$yr <- yr
  npp_i <- npp_i[npp_i$X <= sf::st_bbox(region_sf)$xmax,]


}
