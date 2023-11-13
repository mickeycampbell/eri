#' Calculate ERI for an Area
#'
#' @description
#' Runs the Escape Route Index (ERI) algorithm for every 30 x 30m cell within an area of interest. ERI is the ratio between the distance one can travel within a given time frame, factoring in the surrounding landscape impediments, and the distance one could travel in that same time frame absent any impediments. It ranges from 0 (cannot evacuate from a particular location) to 1 (can evacuate with ideal speed).
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
#' * Note that, in order to avoid erroneous calculations at the edges of your AOI, `time` should be equal to or less than the same argument provided to both [download_lf()] and [download_trans()].
#' * Successful execution of this function will result in the saving of a file named "eri.tif" located within the user-defined `out_dir`.
#'
#' @param veg_type SpatRaster. The OpenStreetMap transportation data-augmented LANDFIRE vegetation type dataset derived from successful execution of [download_lf()], [download_trans()], and [reclass_veg()].
#' @param tm TransitionLayer. The transition matrix generated using [gen_tm()].
#' @param aoi SpatVector. The area of interest within which you aim to map ERI.
#' @param time numeric. The total simulated evacuation time, in minutes, within which you aim to evaluate ERI.
#' @param ncores numeric. The number of cores you want to use to run the algorithm. Defaults to half of the available cores.
#' @param out_dir character. The output directory to which the resulting ERI raster will be saved.
#' @return A SpatRaster with four layers: (1) `eri_mean`; (2) `eri_min`; (3) `eri_max`; (4) `eri_az`.
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
#' # calculate eri within the aoi using a 15 minute evacuation time
#' eri_area_result <- eri_area(veg_type_rc, tm, aoi, 15, 8, "C:/temp")
#' @export
eri_area <- function(veg_type, tm, aoi, time,
                     ncores = floor(parallel::detectCores()/2), out_dir){
  veg_type_file <- terra::sources(veg_type)
  aoi <- terra::project(aoi, veg_type)
  cells_all <- terra::cells(veg_type, aoi)[,2]
  cell_chunks <- split(cells_all, sort(cells_all %% ncores))
  clust <- parallel::makeCluster(ncores)
  chunk_dir <- file.path(out_dir, "eri_chunks")
  dir.create(chunk_dir)
  doParallel::registerDoParallel(clust)
  `%dopar%` <- foreach::`%dopar%`
  eri_chunk_files <- foreach::foreach(i = 1:ncores,
                                      .export = c("eri_pt", "slope_rate"),
                                      .packages = c("terra", "raster", "gdistance", "Matrix")) %dopar%
    {
      veg_type <- terra::rast(veg_type_file)
      eri_mean_rast <- terra::rast(veg_type)
      eri_min_rast <- terra::rast(veg_type)
      eri_max_rast <- terra::rast(veg_type)
      eri_az_rast <- terra::rast(veg_type)
      terra::values(eri_mean_rast) <- NA
      terra::values(eri_min_rast) <- NA
      terra::values(eri_max_rast) <- NA
      terra::values(eri_az_rast) <- NA
      cells <- cell_chunks[[i]]
      for (cell in cells){
        start_pt <- terra::xyFromCell(eri_mean_rast, cell) |>
          terra::vect(type = "points", crs = terra::crs(eri_mean_rast))
        eri_vals <- eri_pt(veg_type, tm, start_pt, time)
        eri_mean_rast[cell] <- eri_vals$eri_mean
        eri_min_rast[cell] <- eri_vals$eri_min
        eri_max_rast[cell] <- eri_vals$eri_max
        eri_az_rast[cell] <- eri_vals$eri_az
      }
      eri_chunk <- c(eri_mean_rast, eri_min_rast, eri_max_rast, eri_az_rast)
      names(eri_chunk) <- c("eri_mean", "eri_min", "eri_max", "eri_az")
      chunk_file <- file.path(chunk_dir, paste0("eri_chunk_", i, ".tif"))
      terra::writeRaster(eri_chunk, chunk_file, overwrite = T)
      return(chunk_file)
    }
  parallel::stopCluster(clust)
  eri_chunk_rasts <- lapply(eri_chunk_files, terra::rast)
  eri_sprc <- terra::sprc(eri_chunk_rasts)
  eri_rast <- terra::merge(eri_sprc)
  out_file <- file.path(out_dir, "eri.tif")
  terra::writeRaster(eri_rast, out_file, overwrite = T)
  return(eri_rast)
}
