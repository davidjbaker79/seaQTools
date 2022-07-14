#' Associating point locations sites with allocated sites
#'
#' @param x A sf object containing Seaquest sites with original (as recorded)
#' site coordinates.
#' @param y An sf object containing allocated site coordinates.
#' @param region_sf An sf object for the region coastal boundary.
#'
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join
#' @import sf
#'
#' @return A data.frame with sites_major associated with original sites.
assoc.major.sites <-
  function(x,
           y,
           region_sf) {
    #--- Buffer region
    region_bf <- st_buffer(region_sf, dist = 30000)

    #-- Nearest major allocated sites
    dist.mat   <-
      st_distance(x, y, by_element = FALSE)
    x$site_major <-
      sapply(1:nrow(x), function(i)
        which.min(dist.mat[i,]))
    x <- st_drop_geometry(x)

    #-- Output
    y <- st_drop_geometry(y)
    out <- left_join(x, y)

  }
