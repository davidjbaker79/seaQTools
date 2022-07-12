#' Associating point locations sites with allocated sites
#'
#' @param point_survey_sites A data.frame containing Seaquest sites with
#' original (as recorded) site coordinates.
#' @param allocated_sites_sf An sf object containing allocated site coordinates.
#' @param region_sf An sf object for the region coastal boundary.
#'
#' @import tidyverse
#' @import sf
#' @import terra
#'
#' @export data.frame with sites_major associated with original sites.
assoc.major.sites <-
  function(point_survey_sites,
           allocated_sites_sf,
           region_sf) {
    #--- Buffer region
    region_bf <- st_buffer(region_sf, dist = 30000)

    #--- Point survey sites
    uni_sites_sf <- point_survey_sites %>%
      select(Survey_ID, site_code, Site_latitude, Site_longitude) %>%
      distinct() %>% # This removes duplicate rows
      st_as_sf(coords = c("Site_longitude", "Site_latitude"),
               crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
      st_transform(crs = st_crs(region_sf)) %>%
      st_crop(st_bbox(sw_bf)) %>%
      mutate(lon = st_coordinates(.)[, 1],
             lat = st_coordinates(.)[, 2])

    #-- Nearest major allocated sites
    dist.mat   <-
      st_distance(uni_sites_sf, allocated_sites_sf, by_element = FALSE)
    uni_sites_sf$site_major <-
      sapply(1:nrow(uni_sites_sf), function(i)
        which.min(dist.mat[i,]))
    uni_sites_sf <- st_drop_geometry(uni_sites_sf)

  }

#' Prepare SDM Seaquest data
#'
#' @param foc_species Name of focal species.
#' @param sq_dat A data.frame containing Seaquest sites with original (as
#' recorded) site coordinates.
#' @param dist_zone Name of distance zone
#' @param enviro_dat A data.frame with environmental data gridded.
#'
#' foc_species = "GREY_SEAL"
#' sq_dat = sq_dat_summer
#' enviro_dat = enviro_dist_bands
prepare.sdm.sq.data <-
  function(foc_species,
           sq_dat,
           enviro_dat,
           max_dist_code) {
    #- Expand out all Surveys to include all zones
    #- This is important for calculating detection
    surv_id <- unique(sq_dat$Survey_ID)
    zones <- c("0_1km" , "1_2km",  "2_5km",  "5_10km")
    surv_zones <- expand.grid(surv_id, zones) %>%
      rename(Survey_ID = "Var1",
             position_dist_code = "Var2") %>%
      arrange(position_dist_code)
    surv_date <- sq_dat %>%
      select(Survey_ID, site_major, lon, lat, Month, Year) %>%
      distinct()
    surv_full <- full_join(surv_date, surv_zones) %>%
      mutate(Year = as.numeric(Year))

    #- P/A for focal species
    sq_dat_i <- sq_dat %>%
      filter(!is.na(position_dist_code),
             species_code == foc_species) %>%
      #mutate(Y = if_else(species_code == foc_species, 1, 0)) %>%
      select(Survey_ID,
             site_major,
             lon,
             lat,
             Year,
             position_dist_code) %>%
      mutate(Y = 1,
             Year = as.numeric(Year))

    #' Fill out the data.frame with each time period a survey took place and
    #' each position_dist_code then calculate the number of surveys per time
    #' period where the target species was detected out of the number of surveys
    sq_dat_z <- left_join(surv_full, sq_dat_i) %>%
      replace_na(list(Y = 0)) %>%
      select(Survey_ID,
             site_major,
             Year,
             position_dist_code,
             Y)  %>%
      group_by(site_major, Year, position_dist_code) %>%
      summarise(nObs = sum(Y),
                nTot = length(Y))

    #- Sites X environment
    sp_env <- left_join(sq_dat_z, enviro_dat) %>%
      na.omit()

  }

#' Prepare SDM Seaquest data for projection
#'
#' @param enviro_dat A data.frame with environmental data gridded.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param coarse_res A numeric giving the horizontal resolution of coarse grid.
#'
#' @import tidyverse
#' @import sf
#'
#' @return A data.frame.
#' modGrid = modGrid1km; enviro_dat = enviro_2020_summer; region_sf = sw_bd
prepare.pred.zones <-
  function(enviro_dat,
           modGrid,
           region_sf,
           coarse_res = 5000) {
    #- Coarse grid
    modGrid_L <- makeGrid(region_sf, coarse_res, 15000) %>%
      rename(id_zone = "id")

    #- Create zones
    buff_1km <- st_buffer(region_sf, dist = 1000)
    buff_2km <- st_buffer(region_sf, dist = 2000)
    buff_5km <- st_buffer(region_sf, dist = 5000)
    buff_10km <- st_buffer(region_sf, dist = 10000)

    #- Grid zones
    zone_1km <- st_intersection(modGrid_L, buff_1km)
    zone_2km <- buff_2km %>%
      st_difference(buff_1km) %>%
      st_intersection(modGrid_L, .)
    zone_5km <- buff_5km %>%
      st_difference(buff_2km) %>%
      st_intersection(modGrid_L, .)
    zone_10km <- buff_10km %>%
      st_difference(buff_5km) %>%
      st_intersection(modGrid_L, .)

    #- 1km zone
    fv_grd_id_0_1km <- st_intersection(modGrid, zone_1km) %>%
      mutate(area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
             position_dist_code = "0_1km") %>%
      select(id, id_zone, X, Y, area_fv, position_dist_code) %>%
      filter(area_fv > 0) %>%
      st_drop_geometry()
    #- 1-2km
    fv_grd_id_1_2km <- st_intersection(modGrid, zone_2km) %>%
      mutate(area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
             position_dist_code = "1_2km") %>%
      select(id, id_zone, X, Y, area_fv, position_dist_code) %>%
      filter(area_fv > 0) %>%
      st_drop_geometry()
    #- 2-5km
    fv_grd_id_2_5km <- st_intersection(modGrid, zone_5km) %>%
      mutate(area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
             position_dist_code = "2_5km") %>%
      select(id, id_zone, X, Y, area_fv, position_dist_code) %>%
      filter(area_fv > 0) %>%
      st_drop_geometry()
    #- 5_10km
    fv_grd_id_5_10km <- st_intersection(modGrid, zone_10km) %>%
      mutate(area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
             position_dist_code = "5_10km") %>%
      select(id, id_zone, X, Y, area_fv, position_dist_code) %>%
      filter(area_fv > 0) %>%
      st_drop_geometry()

    #- Combine
    fv_grd_id <- do.call(
      bind_rows,
      list(
        fv_grd_id_0_1km,
        fv_grd_id_1_2km,
        fv_grd_id_2_5km,
        fv_grd_id_5_10km
      )
    )

    #- Join with environmental data
    env_out <- left_join(fv_grd_id, enviro_dat) %>%
      mutate(position_dist_code =
               factor(
                 position_dist_code,
                 levels = c("0_1km", "1_2km", "2_5km", "5_10km")
               )) %>%
      #filter(!is.na(Year)) %>%
      group_by(id_zone,
               X,
               Y,
               position_dist_code,
               Year) %>%
      summarise(
        sst_m = weighted.mean(sst, area_fv, na.rm = TRUE),
        sst_min = min(sst, na.rm = TRUE),
        sst_max = max(sst, na.rm = TRUE),
        npp_m = weighted.mean(npp, area_fv, na.rm = TRUE),
        npp_min = min(npp, na.rm = TRUE),
        npp_max = max(npp, na.rm = TRUE),
        tri_m = weighted.mean(TRI, area_fv, na.rm = TRUE),
        tri_min = min(TRI, na.rm = TRUE),
        tri_max = max(TRI, na.rm = TRUE),
        depth_m = weighted.mean(depth, area_fv, na.rm = TRUE),
        depth_min = min(depth, na.rm = TRUE),
        depth_max = max(depth, na.rm = TRUE)
      ) %>%
      rename(
        lon = "X",
        lat = "Y"
      ) %>%
      na.omit()

  }
