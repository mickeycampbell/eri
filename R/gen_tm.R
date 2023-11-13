#' Generate Transition Matrix
#'
#' @description
#' Generate a transition matrix using slope, elevation, and vegetation type data using the [gdistance::transition()] function.
#'
#' @details
#' * Travel rates, in meters, are calculated first through the use of the [slope_rate()] function, which estimates travel rate as a function of the terrain slope in the direction of movement.
#' * Then, all slopes greater than 45 degrees, irrespective of travel direction, are treated as barriers.
#' * The resulting travel rates are then multiplied by the vegetation type-based relative conductance factors, defined by [reclass_veg()].
#'
#' @param slope SpatRaster. The LANDFIRE slope dataset derived from successful execution of [download_lf()].
#' @param elev SpatRaster. The LANDFIRE elevation dataset derived from successful execution of [download_lf()].
#' @param veg_type SpatRaster. The OpenStreetMap transportation data-augmented LANDFIRE vegetation type dataset derived from successful execution of [download_lf()], [download_trans()], and [reclass_veg()].
#' @return A TransitionLayer representing the pedestrian conductance values (in meters per second) associated with every 8-directional adjacent cell pairing throughout the study area.
#' @examples
#' # read in aoi
#' aoi <- vect("C:/temp/study_area.shp")
#'
#' # download landfire data using a 15 minute evacuation time
#' lf <- download_lf(aoi, 15, "C:/temp")
#'
#' # access outputs
#' slope <- lf$slope
#' elev <- lf$elev
#' veg_type <- lf$veg_type
#'
#' # download transportation data using a 15 minute evacuation time
#' trans <- download_trans(aoi, 15, "C:/temp")
#'
#' # reclassify vegetation type
#' veg_type_rc <- reclass_veg(veg_type, trans, "C:/temp")
#'
#' # generate transition matrix
#' tm <- gen_tm(slope, elev, veg_type_rc)
#' @export
gen_tm <- function(slope, elev, veg_type){
  slope <- terra::ifel(slope >= 45, 0, 1)
  slope <- raster::raster(slope)
  tm_slope <- gdistance::transition(slope, min, 8, symm = T)
  veg_type <- raster::raster(veg_type)
  tm_veg <- gdistance::transition(veg_type, mean, 8, symm = T)
  elev <- raster::raster(elev)
  tm_elev <- gdistance::transition(elev, function(x) {x[2] - x[1]}, 8, symm = F)
  tm_elev <- gdistance::geoCorrection(tm_elev)
  tm_elev <- atan(tm_elev) * 180/pi
  adj <- raster::adjacent(elev, cells = 1:raster::ncell(elev), pairs = T, directions = 8)
  tm_elev[adj] <- slope_rate(tm_elev[adj])
  tm <- tm_slope * tm_veg * tm_elev
  tm <- gdistance::geoCorrection(tm)
  return(tm)
}
