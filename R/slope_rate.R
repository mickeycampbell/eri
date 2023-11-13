#' Slope-Travel Rate Calculation
#'
#' @description
#' Estimate pedestrian travel rate, in meters per second, as a function of the terrain slope, in degrees.
#'
#' @details
#' The equation that forms the basis of this calculation emerged from [an analysis of crowdsourced *AllTrails* GPS data](https://content.csbs.utah.edu/~pdennison/reprints/denn/2022_Campbell_etal_CEUS.pdf). It represents the median slope-travel rate model among a broad-ranging population of individuals on a similarly diverse range of trails throughout the US.
#'
#' @param slope numeric. Terrain slope, in degrees. You can provide a single slopes or a vector of multiple slopes.
#' @return A numeric vector the same length as that which is provided to the `slope` argument. The estimated pedestrian travel rate, in meters per second, resulting from traversing terrain with the defined slope(s).
#' @examples
#' # estimate travel rate on a flat slope
#' travel_rate <- slope_rate(0)
#'
#' # estimate travel rate on a vector of several downhill and uphill slopes
#' travel_rate <- slope_rate(c(-30,-15,0,15,30))
#' @export
slope_rate <- function(slope){
  a <- -1.4579
  b <- 22.0787
  c <- 76.3271
  d <- 0.0525
  e <- -3.2002e-4
  r <- c*(1/(pi*b*(1+((slope-a)/b)^2)))+d+e*slope
  return(r)
}
