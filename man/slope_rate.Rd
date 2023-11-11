% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope_rate.R
\name{slope_rate}
\alias{slope_rate}
\title{Slope-Travel Rate Calculation}
\usage{
slope_rate(slope)
}
\arguments{
\item{slope}{Numeric. Terrain slope, in degrees. You can provide a single slopes or a vector of multiple slopes.}
}
\value{
Returns a numeric vector the same length as that which is provided to the `slope` argument. The estimated pedestrian travel rate, in meters per second, resulting from traversing terrain with the defined slope(s).
}
\description{
Estimate pedestrian travel rate, in meters per second, as a function of the terrain slope, in degrees.
}
\details{
The equation that forms the basis of this calculation emerged from an analysis of crowdsourced *AllTrails* GPS data. It represents the median slope-travel rate model among a broad-ranging population of individuals on a similarly diverse range of trails throughout the US.
}
\examples{
travel_rate <- slope_rate(0)
travel_rate <- slope_rate(c(-30,-15,0,15,30))
}