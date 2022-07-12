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
#' @import tidyverse
#' @import sf
#' @import terra
#' @import ncdf4
#' @import httr
#'
#' @return A data.frame with monthly sst per grid cell.
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
  region_bf_4326 <- region_sf %>%
    st_buffer(dist = 30000) %>%
    st_transform(crs = "epsg:4326")
  sst_grd_c <- crop(sst_grd, region_bf_4326)
  sst_grd_t <-  terra::project(sst_grd_c, "epsg:27700")
  sst_grd_t <- disagg(sst_grd_t, 10)

  #- Rasterize the sw_bf data to the same grid as the gebco data
  region_bf_t <-  st_transform(region_bf_4326, crs = "epsg:27700")
  sst_bf_mask <- rasterize(vect(region_bf_t), sst_grd_t)

  #- Fill missing values
  sst_grd_t <- focal(sst_grd_t, fun = "mean", w = 21,
                     na.rm=TRUE)

  #- Now mask out the data beyond Xkm offshore
  sst_bf <- terra::mask(sst_grd_t, sst_bf_mask, inverse = FALSE)

  #- Now mask out land
  sw_mask <- rasterize(vect(region_sf), sst_grd_t)
  sst_bf <- terra::mask(sst_bf, sw_mask, inverse = TRUE)

  #- Grid data
  sw_SST_Grd_i <-
    terra::extract(sst_bf, vect(modGrid)) %>%
    rename(sst = "focal_mean") %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    group_by(ID) %>%
    summarise(sst = mean(sst, na.rm = TRUE)) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    bind_cols(modGrid, .) %>%
    st_drop_geometry() %>%
    select(id, area, X, Y, sst) %>%
    filter(!is.nan(sst)) %>%
    mutate(mth = mth,
           yr = yr) %>%
    filter(X <= sf::st_bbox(region_sf)$xmax)

}

#' Extract and grid VGPM data from Oregon State University
#'
#' Extract and grid VGPM data from
#' http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.vgpm.m.chl.m.sst.php
#'
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#' @param yr A numeric giving the year of data to extract.
#' @param mth A numeric giving the month (1-12) of data to extract.
#'
#' @import tidyverse
#' @import sf
#' @import terra
#' @import ncdf4
#' @import httr
#'
#' @return A data.frame with monthly sst per grid cell.
vgpm2grid <- function(region_sf, modGrid, vgpm_path, mth, yr) {

  #- Load vgpm data
  vgpm_dat <- terra::rast(vgpm_path)
  ext(vgpm_dat) <- c(-180, 180, -90, 90)
  vgpm_dat[vgpm_dat < 0] <- NA

  #- vector of months mth_v <- rep(seq(1,12),2)
  if(mth %in% 1:9) mth <- paste0("0", mth)

  #- Buffer and transform region
  region_bf_4326 <- region_sf %>%
    st_buffer(dist = 30000) %>%
    st_transform(crs = "epsg:4326")
  vgpm_grd_c <- crop(vgpm_dat, region_bf_4326)
  vgpm_grd_t <-  terra::project(vgpm_grd_c, "epsg:27700")

  #- Dissaggrate so that missing values can be interoplate later
  vgpm_grd_t <- disagg(vgpm_grd_t, 25)

  #- Rasterize the sw_bf data to the same grid as the gebco data
  region_bf_t <-  st_transform(region_bf_4326, crs = "epsg:27700")
  vgpm_bf_mask <- rasterize(vect(region_bf_t), vgpm_grd_t)

  #- Smooth values
  vgpm_bf <- focal(vgpm_grd_t, w = 55, fun = "mean",  na.rm=TRUE)

  #- Now mask out the data beyond Xkm offshore
  vgpm_bf <- mask(vgpm_bf, vgpm_bf_mask, inverse = FALSE)

  #- Now mask out land
  sw_mask <- rasterize(vect(region_sf), vgpm_grd_t)
  vgpm_bf <- mask(vgpm_bf, sw_mask, inverse = TRUE)

  #- Grid data
  sw_vgpm_Grd_i <-
    terra::extract(vgpm_bf, vect(modGrid)) %>%
    rename(vgpm = 2) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    group_by(ID) %>%
    summarise(npp = mean(vgpm, na.rm = TRUE)) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    bind_cols(modGrid, .) %>%
    st_drop_geometry() %>%
    select(id, area, X, Y, npp) %>%
    filter(!is.nan(npp)) %>%
    mutate(mth = mth,
           yr = yr) %>%
    filter(X <= sf::st_bbox(region_sf)$xmax)

}
