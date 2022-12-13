#' Prepare Seaquest data for modelling without distance zones
#'
#' @param x A data.frame containing Seaquest sites with original (as
#' recorded) site coordinates
#' @param y A data.frame with environmental data gridded.
#' @param a A data.frame with allocated sites.
#' @param speciesName Name of focal species.
#'
#' @import lubridate
#' @import sf
#'
#' @return A data.frame.
#'
#' @export
prepare.mod.data <- function(x, y, a, rsf, speciesName) {

  # -- Set up species data
  x <- x[, c(
    "Survey_ID",
    "lon",
    "lat",
    "Survey_date",
    "species_code",
    "ST_visibility_code",
    "ST_sea_state_code"
  )]
  x$Survey_date <- lubridate::ymd(x$Survey_date)
  x$Month <- format(x$Survey_date, "%m")
  x$Year <- format(x$Survey_date, "%Y")

  # -- Associate allocated sites with original sites
  xs <- x[,c("lon", "lat")]
  xs <- xs[!duplicated(xs),]
  xs <- st_as_sf(xs, coords = c("lon", "lat"), crs = st_crs(rsf))
  xs <- assoc.major.sites(xs, a)

  # -- Merge with sq data
  x <- merge(x, xs, by = c("lon", "lat"), all.x = TRUE)

  # -- Prepare species data
  sp <- x[, c("Survey_ID", "site_major", "lon", "lat",  "Month", "Year")]
  sp <- sp[!duplicated(sp),]
  sp$Year <- as.numeric(sp$Year)
  sp <- sp[!is.na(sp$site_major), ]

  # -- P/A for focal species
  xi <- x[x$species_code == speciesName,]
  xi <-  xi[!duplicated(xi),]
  xi <- xi[, c("Survey_ID", "site_major", "Year")]
  xi$Y <- 1
  xi$Year <- as.numeric(xi$Year)

  # -- Join with environmental data
  x_z <- merge(sp, xi, by = c("Survey_ID", "site_major", "Year"), all = TRUE)
  x_z$Y[is.na(x_z$Y)] <- 0
  x_z <-
    x_z[, c("Survey_ID", "site_major", "Month", "Year", "Y")]
  x_zs <- aggregate(Y ~ site_major + Year,
                    data = x_z,
                    FUN = function(x) c(nObs = sum(x), nTot = length(x) ))
  x_zs <- do.call(data.frame, x_zs)
  sp_env <- merge(x_zs, y, by = c("site_major", "Year"))
  sp_env <- na.omit(sp_env)

}
