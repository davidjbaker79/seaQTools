#' Prepare ad hoc data for modelling
#'
#' @param x A data.frame containing ad hoc occurrence data (as
#' recorded) site coordinates
#' @param y A data.frame with environmental data gridded.
#' @param speciesName Name of focal species.
#' @param tgt_group_spp Name of species to be included in the target group. Defaults to all species.
#' @param start_date String giving the start date in the format "YYYY/MM/DD".
#' @param end_date String giving the end date in the format "YYYY/MM/DD".
#' @param mths Numeric vector for the months included in the analysis
#' @param n_bkgrd Numeric for the number of background point per year.
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom scales rescale
#' @importFrom gtools mixedsort
#' @importFrom rlang .data
#'
#' @return A data.frame.
prepare.adhoc.sdm.data <-
  function(x,
           y,
           speciesName,
           tgt_group_spp,
           start_date,
           end_date,
           mths,
           n_bkgrd) {
    #- Associated with 5km grid
    ah_z <- left_join(x, y)
    ah_z <- ah_z[(!is.na(ah_z$lon_z) & !is.na(ah_z$lat_z)),]

    #- Date ranges
    ah_z <-
      ah_z[ah_z$Survey_date >= as.Date(start_date) &
             ah_z$Survey_date <= as.Date(end_date), ]

    #- Remove dead records
    ah_z <- ah_z[ah_z$RecordType != "Dead",]

    #- Filter months
    if (!is.null(months)) {
      ah_z <- ah_z[as.numeric(ah_z$Month) %in% mths,]
    }
    #- Select columns
    ah_z <-
      ah_z[, c("id_zone",
               "lon_z",
               "lat_z",
               "SpeciesVenacular",
               "Year",
               "position_dist_code")]
    #ggplot(ah_f) + geom_point(aes(x= lon_z, y = lat_z))

    #- Filter focal species
    ah_f <- ah_z[ah_z$SpeciesVenacular == speciesName, ]
    ah_f$det <- 1

    #---- Create background sample
    #- Target group data
    if (tgt_group_spp == "all") {
      tgt_group <- unique(ah_z$SpeciesVenacular)
      tgt_group <- tgt_group[!(tgt_group %in% focal_species)]
    }
    tg_dat <-
      ah_z[ah_z$SpeciesVenacular %in% tgt_group, ]

    #- Calculate number of records per 1km grid cell
    tg_dat$nRec <- 1
    sppDat <- aggregate(
      nRec  ~
        id_zone + lon_z + lat_z + Year + position_dist_code,
      data = tg_dat,
      FUN = function(x)
        length(x)
    )

    #- Unique years
    Yrs <- mixedsort(unique(sppDat$Year))

    #- Target group background by month/year i = 2011
    st_bkgd <- lapply(Yrs, function(i) {
      sppDat_i <- filter(sppDat, Year == i)
      sppDat_i$v = scales::rescale(sppDat_i$nRec, to = c(0.01, 1))
      #- Select background sample
      bkgd_i <- sppDat_i[sample(
        1:nrow(sppDat_i),
        size = n_bkgrd,
        replace = TRUE,
        prob = sppDat_i$v
      ), ]
      bkgd_i$Year <- i
      bkgd_i$det <- 0
      bkgd_i
    })
    st_bkgd <- do.call(rbind, st_bkgd)
    st_bkgd <-
      st_bkgd[, c("id_zone",
                  "lon_z",
                  "lat_z",
                  "Year",
                  "position_dist_code",
                  "det")]
    #---- Join with ad hoc occurrence data
    mod_dat <- bind_rows(ah_f, st_bkgd)

    #--- Merge back in environment data
    mod_dat <- left_join(mod_dat, y)

    #-- need to remove NA for model
    mod_dat <- na.omit(mod_dat)


  }

