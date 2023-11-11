#' Slope-Travel Rate Calculation
#'
#' @description
#' Estimate pedestrian travel rate, in meters per second, as a function of the terrain slope, in degrees.
#'
#' @details
#' The equation that forms the basis of this calculation emerged from an analysis of crowdsourced *AllTrails* GPS data. It represents the median slope-travel rate model among a broad-ranging population of individuals on a similarly diverse range of trails throughout the US.
#'
#' @param slope Numeric. Terrain slope, in degrees. You can provide a single slopes or a vector of multiple slopes.
#' @return Returns a numeric vector the same length as that which is provided to the `slope` argument. The estimated pedestrian travel rate, in meters per second, resulting from traversing terrain with the defined slope(s).
#' @examples
#' travel_rate <- slope_rate(0)
#' travel_rate <- slope_rate(c(-30,-15,0,15,30))
#' @export
slope_rate <- function(slope){
  r <- (0.006512 * slope ^ 2 + 80.88869) /
    (0.140214 * slope ^ 2 + 70.389209)
  return(r)
}
