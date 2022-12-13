#' Associating point locations sites with allocated sites
#'
#' These sites have been manually aggregated where it looked like multiple
#' designations had been given for the same site.
#'
#' @param x A sf object containing Seaquest sites with original (as recorded)
#' site coordinates.
#' @param y An sf object containing allocated site coordinates.
#'
#' @import sf
#'
#' @return A data.frame with sites_major associated with original sites.
#'
#' @export
assoc.major.sites <- function(x, y) {

  #-- Nearest major allocated sites
  cl <- list()
  for(i in seq_len(nrow(x))){
    cl[[i]] <- y[which.min(st_distance(y, x[i,])),]
  }
  cl <- do.call(rbind, cl)
  cl <- st_drop_geometry(cl)

  #-- Output
  x <- cbind(x, st_coordinates(x))
  x <- st_drop_geometry(x)
  x <- cbind(x, cl)
  names(x) <- c("lon", "lat", "site_major")
  x
  }
