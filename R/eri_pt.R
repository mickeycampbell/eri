#' Calculate ERI for a Single Point
#'
#' @description
#' Runs the Escape Route Index (ERI) algorithm for a single starting point. ERI is the ratio between the distance one can travel within a given time frame, factoring in the surrounding landscape impediments, and the distance one could travel in that same time frame absent any impediments. It ranges from 0 (cannot evacuate from a particular location) to 1 (can evacuate with ideal speed).
#'
#' @details
#' * In two dimensions, ERI is simple: the ratio between the distance one can travel through an existing landscape and the distance one could travel through an imaginary landscape without any impedance. If the former was 100 m and the latter was 200 m, ERI would be equal to 0.5.
#' * However, ERI is calculated in two dimensions, as follows:
#'   1. Given a starting point...
#'   2. An accumulative cost surface (representing travel time from that starting point) is generated using slope, elevation, and vegetation type-based transition matrix from [gen_tm()].
#'   3. A threshold is applied to that cost surface, defined by the evacuation time frame of interest, parameterized by `time`.
#'   4. The edge of the resulting area represents a travel time isochrone (a line of equal travel time) from the starting point.
#' * Accordingly, ERI can, and is, calculated in a number of ways:
#'   - ERImean is the mean distance to all points along the isochrone divided by the unimpeded distance one can travel in the same time frame, representing average evacuation ability in all directions
#'   - ERImin is the distance to the closest point along the isochrone divided by the unimpeded distance one can travel in the same time frame, representing evacuation ability in the worst-case scenario
#'   - ERImax is the distance to the furthest point along the isochrone divided by the unimpeded distance one can travel in the same time frame, representing evacuation ability in the best-case scenario
#'   - ERIaz is the azimuth between the starting point and the furthest point along the isochrone, representing the best evacuation direction
#' * The primary purpose of this function is to drive the [eri_area()] function, which maps ERImean, ERImin, ERImax, and ERIaz across your entire AOI.
#' * However, you can also run this on a single SpatVector starting point. If you intend to do so for a point of interest, you may consider setting `spatial = TRUE`. In addition to returning a list of numeric ERI metrics, this will return a series of SpatVector and SpatRaster objects that are useful for visualization and can be fed into [plot_eri()] for a simple visualization.
#' * Note that, in order to avoid erroneous calculations at the edges of your AOI, `time` should be equal to or less than the same argument provided to both [download_lf()] and [download_trans()].
#'
#' @param veg_type SpatRaster. The OpenStreetMap transportation data-augmented LANDFIRE vegetation type dataset derived from successful execution of [download_lf()], [download_trans()], and [reclass_veg()].
#' @param tm TransitionLayer. The transition matrix generated using [gen_tm()].
#' @param start_pt SpatVector. The starting point from which evacuation will be simulated and ERI will be calculated.
#' @param time numeric. The total simulated evacuation time, in minutes, within which you aim to evaluate ERI.
#' @param spatial Boolean. Defines whether you only want quantitative ERI metrics (`spatial = FALSE`) or if you also want SpatVector/SpatRaster representations of the spatial features used to calculate ERI (`spatial = T`).
#' @return
#' If `spatial == FALSE`, a list with four numeric objects:
#' * `eri_mean`
#' * `eri_min`
#' * `eri_max`
#' * `eri_az`
#'
#' If `spatial == TRUE`, all of the above, plus:
#' * `start_pt`: SpatVector. The same starting point supplied by the user in the function call
#' * `max_poly`: SpatVector. A polygon representing the distance one can travel in the absence of travel impediments. Takes the form of an octagon due to the 8-directional nature of cell adjacency in the accumulative cost simulations.
#' * `act_poly`: SpatVector. A polygon representing the distance one can travel in the presence of travel impediments.
#' * `tt_rast`: SpatRaster. The accumulative cost raster, clipped to the extent of `act_poly`, where each cell represents travel time in seconds from the starting point.
#' * `pt_max`: SpatVector. A polygon representation of the travel time raster cell furthest from the starting point used as the basis of ERImax and ERIaz calculations.
#' * `pt_min`: SpatVector. A polygon representation of the travel time raster cell closest from the starting point used as the basis of ERImin calculation.
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
#'
#' # generate a random starting point within the aoi
#' start_pt <- terra::project(aoi, veg_type_rc) |>
#'   terra::spatSample(1)
#'
#' # calculate eri from that point using a 15 minute evacuation time
#' eri_pt_result <- eri_pt(veg_type_rc, tm, start_pt, 15, T)
#' @export
eri_pt <- function(veg_type, tm, start_pt, time, spatial = F){
  dist <- slope_rate(0) * time * 60
  buff <- terra::buffer(start_pt, dist + 30)
  clipper <- terra::ext(buff)
  clip <- terra::crop(veg_type, clipper)
  cells <- terra::cells(veg_type, clipper)
  tm_sub <- new("TransitionLayer",
                nrows = as.integer(terra::nrow(clip)),
                ncols = as.integer(terra::ncol(clip)),
                extent = raster::extent(terra::ext(clip)[1:4]),
                crs = terra::crs(clip),
                transitionMatrix = Matrix::Matrix(0, terra::ncell(clip), terra::ncell(clip)),
                transitionCells = 1:terra::ncell(clip))
  gdistance::transitionMatrix(tm_sub) <- gdistance::transitionMatrix(tm)[cells,cells]
  xy <- terra::geom(start_pt)[,c("x","y")]
  acc <- gdistance::accCost(tm_sub, xy) |> terra::rast()
  terra::crs(acc) <- terra::crs(veg_type)
  tt_poly <- terra::ifel(acc <= (time * 60), 1, NA)
  tt_edge <- terra::boundaries(tt_poly, falseval = NA)
  tt_dist <- terra::distance(tt_edge, start_pt) |>
    terra::mask(tt_edge)
  tt_max <- terra::where.max(tt_dist)
  tt_max_xy <- terra::xyFromCell(tt_dist, tt_max[1,2][[1]]) -
    terra::crds(start_pt)
  eri_az <- atan2(tt_max_xy[,2], tt_max_xy[,1]) * (180 / pi)
  eri_az <- (360 - eri_az + 90) %% 360
  eri_denom <- 0.9303937 * dist
  eri_mean <- terra::global(tt_dist, mean, na.rm = T) / eri_denom
  eri_min <- terra::global(tt_dist, min, na.rm = T) / eri_denom
  eri_max <- tt_max[1,3][[1]] / eri_denom
  if (spatial == F){
    return(list(eri_mean = eri_mean[[1]],
                eri_min = eri_min[[1]],
                eri_max = eri_max[[1]],
                eri_az = eri_az[[1]]))
  } else {
    xy0 <- terra::crds(start_pt)
    x0 <- xy0[,1]
    y0 <- xy0[,2]
    for (i in c(0,45,90,135,180,225,270,315)){
      rad <- i * pi / 180
      xi <- x0 + dist * sin(rad)
      yi <- y0 + dist * cos(rad)
      if (i == 0){
        max_poly <- rbind(c(1, 1, xi, yi, 0))
      } else {
        max_poly <- rbind(max_poly, c(1, 1, xi, yi, 0))
      }
    }
    colnames(max_poly) <- c("object", "part", "x", "y", "hole")
    max_poly <- terra::vect(max_poly, "polygons", crs = terra::crs(start_pt))
    act_poly <- terra::as.polygons(tt_poly)
    tt_rast <- terra::mask(acc, act_poly)
    pt_max <- terra::xyFromCell(tt_dist, tt_max[,2][[1]]) |>
      terra::vect(crs = terra::crs(tt_dist))
    pt_min <- terra::xyFromCell(tt_dist,
                                terra::where.min(tt_dist)[,2][[1]]) |>
      terra::vect(crs = terra::crs(tt_dist))
    return(list(start_pt = start_pt,
                max_poly = max_poly,
                act_poly = act_poly,
                tt_rast = tt_rast,
                pt_max = pt_max,
                pt_min = pt_min,
                eri_mean = eri_mean[[1]],
                eri_min = eri_min[[1]],
                eri_max = eri_max[[1]],
                eri_az = eri_az[[1]]))
  }
}
