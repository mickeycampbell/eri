#' Download LANDFIRE Data
#'
#' @description
#' ERI relies on three raster input datasets, each of which are provided by LANDFIRE: (1) slope in degrees; (2) elevation in meters; and (3) existing vegetation type. This function queries LANDFIRE's API with the user-provided area of interest (AOI) and tries to download the three datasets.
#'
#' @details
#' * The ERI algorithm is very processing-intensive. Accordingly, users should avoid using very large AOIs (e.g., >>100 sqkm).
#' * In addition to AOI area, the `time` argument also partially dictates study area size and processing time. It is used in conjunction with the [slope_rate()] function to generate a buffer distance around the AOI within which data will be downloaded. This is to ensure that there are no edge effects in the calculation of ERI within a particular evacuation time frame.
#' * Users should specify the same AOI and travel time for both this function and [download_trans()].
#' * LANDFIRE's API has a small delay betweeen when requests are submitted and when data are returned. This function tries several times to download the returned data, with a 5-second delay between each try. If it reaches 10 tries, it stops, under the assumption that there is a server-side error or delay that is preventing successful data downloading.
#' * In the current version of the `eri` package, LF2020 data are provided. Note that vegetation may have been disturbed since the release of these data and thus ERI calculations may not reflect updated landscape conditions.
#' * Note that this function will save three TIF files to `out_dir`: (1) "slope.tif"; (2) "elevation.tif"; and (3) "veg_type.tif".
#' * By default, LANDFIRE's API returns data layers in a local Albers projection centered on the centroid of the AOI that maintains both pixel area and N-S-E-W orientation.
#'
#' @param aoi SpatVector. The area of interest (AOI) within which you aim to evaluate ERI. Can be defined using `terra::vect()`.
#' @param time numeric. The total simulated evacuation time, in minutes, within which you aim to evaluate ERI.
#' @param out_dir character. The output directory to which LANDFIRE data will be downloaded.
#' @return A list with three SpatRasters: `slope` (terrain slope in degrees); `elev` (terrain elevation in meters); and `veg_type` (LANDFIRE Existing Vegetation Type)
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
#' @export
download_lf <- function(aoi, time, out_dir){
  message(paste0(print_time(), " Trying to download LANDFIRE data..."))
  aoi <- terra::project(aoi, "epsg:4326")
  dist <- slope_rate(0) * time * 60 + 30
  buff <- terra::buffer(aoi, dist)
  e <- terra::ext(buff)
  xmin <- e[1]
  xmax <- e[2]
  ymin <- e[3]
  ymax <- e[4]
  base_url <- "https://lfps.usgs.gov/arcgis/rest/services/LandfireProductService/GPServer/LandfireProductService/submitJob?"
  prod_url <- "Layer_List=SLPD2020;ELEV2020;220EVT"
  aoi_url <- paste0("&Area_Of_Interest=",
                    xmin, "%20", ymin, "%20", xmax, "%20", ymax)
  full_url <- paste0(base_url, prod_url, aoi_url)
  resp <- httr::GET(full_url)
  if (resp$status_code != 200){
    status <- httr::http_status(resp)$message
    msg <- paste0(print_time(),
                  " The LANDFIRE API returned a status code of: ", status,
                  ". Please try again later.")
    stop(msg)
  }
  job_url <- resp[[1]]
  output_url <- paste0(job_url, "/results/Output_File?f=pjson")
  output_url
  success <- F
  n_tries <- 0
  while (success == F & n_tries < 10){
    Sys.sleep(5)
    n_tries <- n_tries + 1
    message(paste0(print_time(), " Try #", n_tries))
    resp <- httr::GET(output_url)
    j <- rawToChar(resp$content) |> jsonlite::fromJSON()
    if ("error" %in% names(j)) next
    message(paste0(print_time(), "   Success! Downloading file...\n"))
    zip_file <- file.path(out_dir, paste0("lf_",
                                          format(Sys.time(), "%Y%m%d%H%M%S"), ".zip"))
    dl_url <- j$value$url
    download.file(dl_url, zip_file, method = "libcurl")
    success <- T
  }
  if (!file.exists(zip_file)){
    stop(paste0(print_time(), " Maximum number of tries reached. ",
                "Try a smaller study area or try again later."))
  }
  temp_dir <- file.path(out_dir, "temp")
  unzip(zip_file, exdir = temp_dir)
  tif_file <- list.files(temp_dir, "*.tif$", full.names = T)
  r <- terra::rast(tif_file)
  out_files <- file.path(out_dir, c("slope.tif", "elevation.tif", "veg_type.tif"))
  terra::writeRaster(r, out_files, overwrite = T)
  unlink(temp_dir, recursive = T)
  unlink(zip_file)
  tif_slp <- terra::rast(file.path(out_dir, "slope.tif"))
  tif_elv <- terra::rast(file.path(out_dir, "elevation.tif"))
  tif_evt <- terra::rast(file.path(out_dir, "veg_type.tif"))
  return(list(slope = tif_slp, elev = tif_elv, veg_type = tif_evt))
}
