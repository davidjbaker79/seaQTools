#' Determine the viewable area from each site and extract enviro data
#'
#' @param n A numeric giving row of sites_major_sf.
#' @param sites_major_sf A dataframe giving coordinates of each major site.
#' @param height A numeric height above sea level.
#' @param region_sf An sf object defining the spatial coastal region of study.
#' @param modGrid An sf object defining the grid for analysis created using the
#' makeGrid function.
#'
#' @import tidyverse
#' @import sf
#' @import concaveman
#' @import Hmisc
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
  p1 <- sites_major_sf[n,] %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2])

  #- Crop out a single site
  p1_b <- st_buffer(p1, hz)

  #- crop coast
  pc <- st_difference(p1_b, region_sf) %>%
    st_cast("POLYGON")
  pc <- pc[which(st_area(pc) == max(st_area(pc))),]

  #- Points on buffer edge
  poly_points <- st_segmentize(p1_b, dfMaxLength = 1000) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    select(X, Y) %>%
    st_as_sf(coords = c(X, Y),
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
  # ggplot() +
  #   geom_sf(data = lines_sf) +
  #   geom_sf(data = poly_points, size = 1, colour = "red")

  #- crop coast
  lseg <- st_difference(lines_sf, region_sf)
  lseg <- lseg %>%
    st_cast("MULTILINESTRING") %>%
    st_cast("LINESTRING")
  int <- lengths(st_intersects(lseg, st_buffer(p1, 500))) > 0
  lseg <- lseg[int,]
  # ggplot() +
  #   geom_sf(data = lseg) +
  #   geom_sf(data = poly_points,
  #           size = 1,
  #           colour = "red") +
  #   geom_sf(data = region_sf)

  #- Convert linestring to polygon
  fieldview <- concaveman(lseg)
  #ggplot() + geom_sf(data = fieldview)

  #- Add distance stratification
  buff_1km <- st_buffer(p1, dist = 1000)
  buff_2km <- st_buffer(p1, dist = 2000)
  buff_5km <- st_buffer(p1, dist = 5000)
  buff_10km <- st_buffer(p1, dist = 10000)

  #- Create zones
  zone_1km <- st_intersection(buff_1km, fieldview)
  zone_2km <- buff_2km %>%
    st_difference(buff_1km) %>%
    st_intersection(fieldview)
  zone_5km <- buff_5km %>%
    st_difference(buff_2km) %>%
    st_intersection(fieldview)
  zone_10km <- buff_10km %>%
    st_difference(buff_5km) %>%
    st_intersection(fieldview)

  #-- Extract environmental data
  #- 1km zone
  fv_grd_id_0_1km <- st_intersection(modGrid, zone_1km) %>%
    mutate(
      area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
      dist_fv = as.numeric(round(st_distance(p1, ., by_element = T))),
      position_dist_code = "0_1km",
      viewArea =  as.numeric(st_area(zone_1km) / 10 ^ 6)
    ) %>%
    select(id,
           area,
           site_major,
           area_fv,
           dist_fv,
           position_dist_code,
           viewArea) %>%
    st_drop_geometry() %>%
    left_join(env_dat)
  #- 1-2km
  fv_grd_id_1_2km <- st_intersection(modGrid, zone_2km) %>%
    mutate(
      area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
      dist_fv = as.numeric(round(st_distance(p1, ., by_element = T))),
      position_dist_code = "1_2km",
      viewArea =  as.numeric(st_area(zone_2km) / 10 ^ 6)
    ) %>%
    select(id,
           area,
           site_major,
           area_fv,
           dist_fv,
           position_dist_code,
           viewArea) %>%
    st_drop_geometry() %>%
    left_join(env_dat)
  #- 2-5km
  fv_grd_id_2_5km <- st_intersection(modGrid, zone_5km) %>%
    mutate(
      area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
      dist_fv = as.numeric(round(st_distance(p1, ., by_element = T))),
      position_dist_code = "2_5km",
      viewArea =  as.numeric(st_area(zone_5km) / 10 ^ 6)
    ) %>%
    select(id,
           area,
           site_major,
           area_fv,
           dist_fv,
           position_dist_code,
           viewArea) %>%
    st_drop_geometry() %>%
    left_join(env_dat)
  #- 5-10km
  fv_grd_id_5_10km <- st_intersection(modGrid, zone_10km) %>%
    mutate(
      area_fv = as.numeric(round(st_area(.) / 10 ^ 6, 2)),
      dist_fv = as.numeric(round(st_distance(p1, ., by_element = T))),
      position_dist_code = "5_10km",
      viewArea =  as.numeric(st_area(zone_10km) / 10 ^ 6)
    ) %>%
    select(id,
           area,
           site_major,
           area_fv,
           dist_fv,
           position_dist_code,
           viewArea) %>%
    st_drop_geometry() %>%
    left_join(env_dat)

  #- Combine and processes environmental variables at major site and distance band
  fv_grd_id <- do.call(
    bind_rows,
    list(
      fv_grd_id_0_1km,
      fv_grd_id_1_2km,
      fv_grd_id_2_5km,
      fv_grd_id_5_10km
    )
  ) %>%
    mutate(position_dist_code =
             factor(
               position_dist_code,
               levels = c("0_1km", "1_2km", "2_5km", "5_10km")
             )) %>%
    group_by(site_major,
             lon,
             lat,
             viewArea,
             position_dist_code,
             Year) %>%
    summarise(
      sst_m = weighted.mean(sst, area_fv, na.rm = TRUE),
      sst_sd = sqrt(wtd.var(sst, area_fv, na.rm = TRUE)),
      sst_min = min(sst, na.rm = TRUE),
      sst_max = max(sst, na.rm = TRUE),
      npp_m = weighted.mean(npp, area_fv, na.rm = TRUE),
      npp_sd = sqrt(wtd.var(npp, area_fv, na.rm = TRUE)),
      npp_min = min(npp, na.rm = TRUE),
      npp_max = max(npp, na.rm = TRUE),
      tri_m = weighted.mean(TRI, area_fv, na.rm = TRUE),
      tri_sd = sqrt(wtd.var(TRI, area_fv, na.rm = TRUE)),
      tri_min = min(TRI, na.rm = TRUE),
      tri_max = max(TRI, na.rm = TRUE),
      depth_m = weighted.mean(depth, area_fv, na.rm = TRUE),
      depth_sd = sqrt(wtd.var(depth, area_fv, na.rm = TRUE)),
      depth_min = min(depth, na.rm = TRUE),
      depth_max = max(depth, na.rm = TRUE)
    )

}
