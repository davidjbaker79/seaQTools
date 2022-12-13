#' Extract and grid Modis Sea Surface Temperature from MODIS
#'
#' Download and extract data from
#' http://oceandata.sci.gsfc.nasa.gov/opendap/hyrax/MODIST/L3SMI/ to a supplied
#' grid.
#'
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#' @param yr A numeric giving the year of data to extract.
#' @param mth A numeric giving the month (1-12) of data to extract.
#'
#' @importFrom lubridate days_in_month
#' @importFrom ncdf4 nc_open
#' @importFrom ncdf4 ncvar_get
#' @import terra
#' @import sf
#' @import data.table
#'
#' @return A data.frame with monthly sst per grid cell.
#'
#' @export
modis_sst <- function(region_sf, modGrid, yr, mth) {

  #- vector of months
  if(mth %in% 1:9) mth <- paste0("0", mth)
  # end date
  date <- as.Date(paste0(yr, "-",mth,"-01"), "%Y-%m-%d")
  days_in_mth <- as.character(lubridate::days_in_month(date))

  #- Open url
  url_root <-
    paste0(
      "http://oceandata.sci.gsfc.nasa.gov/opendap/hyrax/MODIST/L3SMI/",
      yr,
      "/",
      mth,
      "01/"
    )
  url_i <-
    paste0("TERRA_MODIS.",
           yr,
           mth,
           "01_",
           yr,
           mth,
           days_in_mth,
           ".L3m.MO.SST.sst.4km.nc")
  url <- paste0(url_root, url_i)
  modis <- nc_open(url)

  #- Extract SST
  lonm <- ncvar_get(modis, 'lon')
  latm <- ncvar_get(modis, 'lat')
  sst <- ncvar_get(modis, "sst")
  sst <- t(as.matrix(sst))
  sst_grd <- rast(sst,
                  crs = "epsg:4326",
                  extent = c(-180,180,-90,90))

  #- Buffer and transform region
  region_bf <- st_buffer(region_sf, dist = 30000)
  region_bf <-  st_transform(region_bf, crs = "epsg:4326")
  sst_grd_c <- crop(sst_grd, region_bf)
  sst_grd_t <- project(sst_grd_c, "epsg:27700")
  sst_grd_t <- disagg(sst_grd_t, 10)

  #- Rasterize the sw_bf data to the same grid as the gebco data
  region_bf_t <-  st_transform(region_bf, crs = "epsg:27700")
  sst_bf_mask <- rasterize(vect(region_bf_t), sst_grd_t)

  #- Fill missing values
  sst_grd_t <- focal(sst_grd_t, fun = "mean", w = 21,
                     na.rm=TRUE)
  names(sst_grd_t) <- "sst"

  #- Now mask out the data beyond Xkm offshore
  sst_bf <- mask(sst_grd_t, sst_bf_mask, inverse = FALSE)

  #- Now mask out land
  sw_mask <- rasterize(vect(region_sf), sst_grd_t)
  sst_bf <- mask(sst_bf, sw_mask, inverse = TRUE)

  #- Grid data
  sst_i <- extract(sst_bf, vect(modGrid))
  sst_i <- as.data.table(sst_i)
  sst_i <- sst_i[, mean(sst, na.rm=TRUE), by = list("ID")]
  sst_i <- cbind(modGrid, sst_i)
  sst_i <- st_drop_geometry(sst_i)
  sst_i <- sst_i[,c("id", "area", "X", "Y", "V1")]
  names(sst_i)[5] <- "sst"
  sst_i <- sst_i[!is.nan(sst_i$sst),]
  sst_i$mth <- mth
  sst_i$yr <- yr
  sst_i <- sst_i[sst_i$X <= sf::st_bbox(region_sf)$xmax,]

}

