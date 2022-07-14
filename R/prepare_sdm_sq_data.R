#' Prepare Seaquest data for modelling
#'
#' @param x A data.frame containing Seaquest sites with original (as
#' recorded) site coordinates
#' @param y A data.frame with environmental data gridded.
#' @param speciesName Name of focal species.
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr distinct
#' @importFrom dplyr full_join
#' @importFrom dplyr left_join
#' @importFrom plyr ddply
#' @import sf
#' @importFrom rlang .data
#'
#' @return A data.frame.
#'
#' @export
prepare.sdm.sq.data <- function(x, y, speciesName) {
  sid <- unique(x$Survey_ID)
  z <- c("0_1km" , "1_2km",  "2_5km",  "5_10km")
  sz <- expand.grid(sid, z)
  names(sz) <- c("Survey_ID", "position_dist_code")
  sz <- arrange(sz, .data$position_dist_code)
  sdat <-
    sdat[, c("Survey_ID", "site_major", "lon", "lat", "Month", "Year")]
  sdat <- distinct(sdat)
  sdat <- full_join(sdat, sz)
  sdat$Year <- as.numeric(sdat$Year)
  sdat <- sdat[!is.na(sdat$site_major), ]

  #- P/A for focal species
  x_i <- x[!is.na(x$position_dist_code),]
  x_i <- x_i[x_i$species_code == speciesName,]
  x_i <-
    x_i[, c("Survey_ID",
            "site_major",
            "lon",
            "lat",
            "Year",
            "position_dist_code")]
  x_i$Y <- 1
  x_i$Year = as.numeric(x_i$Year)
  x_z <- left_join(sdat, x_i)
  x_z$Y[is.na(x_z$Y)] <- 0
  x_z <-
    x_z[, c("Survey_ID", "site_major", "Year", "position_dist_code", "Y")]
  x_z <-
    ddply(
      x_z,
      c("site_major", "Year", "position_dist_code"),
      .data$summarise,
      nObs = sum(.data$Y),
      nTot = length(.data$Y)
    )
  sp_env <- left_join(x_z, y)
  sp_env <- na.omit(sp_env)

}
