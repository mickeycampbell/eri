#' Download Transportation Data
#'
#' @description
#' Transportation features (roads, trails, etc.) represent easily walkable portions of the landscape, some of which may be missed in the LANDFIRE Existing Vegetation Type classification. This function downloads transportation linear features (roads, trails, etc.) from OpenStreetMap for use in the ERI modeling process.
#'
#' @details
#' * The ERI algorithm is very processing-intensive. Accordingly, users should avoid using very large AOIs (e.g., >>100 sqkm).
#' * In addition to AOI area, the `time` argument also partially dictates study area size and processing time. It is used in conjunction with the [slope_rate()] function to generate a buffer distance around the AOI within which data will be downloaded. This is to ensure that there are no edge effects in the calculation of ERI within a particular evacuation time frame.
#' * Users should specify the same AOI and travel time for both this function and [download_lf()].
#' * Note that this function will save the following shapefile to `out_dir`: "transportation.shp".
#' * Although LANDFIRE Existing Vegetation Type captures the dominant land cover conditions, smaller roads, trails, and other minor transportation features may still be classified as restrictive vegetation types (i.e., a trail through the woods being classified as forest). Given the ease with which pedestrians can traverse these constructed features, this function provides the opportunity to download OpenStreetMap transportation linear features. In a subsequent function ([reclass_veg()]), the LANDFIRE Existing Vegetation Type and the OpenStreetMap transportation features are merged to create a combined land cover travel impedance raster layer.
#'
#' @param aoi SpatVector. The area of interest (AOI) within which you aim to evaluate ERI. Can be defined using `terra::vect()`.
#' @param time numeric. The total simulated evacuation time, in minutes, within which you aim to evaluate ERI.
#' @param out_dir character. The output directory to which transportation data will be downloaded.
#' @return A SpatVector of all transportation features clipped to the extent of the travel time-buffered AOI.
#' @examples
#' # read in aoi
#' aoi <- vect("C:/temp/study_area.shp")
#'
#' # download transportation data using a 15 minute evacuation time
#' trans <- download_trans(aoi, 15, "C:/temp")
#' @export
download_trans <- function(aoi, time, out_dir){
  aoi <- terra::project(aoi, "epsg:4326")
  dist <- slope_rate(0) * time * 60 + 30
  buff <- terra::buffer(aoi, dist)
  e <- terra::ext(buff)
  xmin <- e[1][[1]]
  xmax <- e[2][[1]]
  ymin <- e[3][[1]]
  ymax <- e[4][[1]]
  bb <- c(xmin, ymin, xmax, ymax)
  feats <- list("highway" = "motorway",
                "highway" = "trunk",
                "highway" = "primary",
                "highway" = "secondary",
                "highway" = "tertiary",
                "highway" = "unclassified",
                "highway" = "residential",
                "highway" = "motorway_link",
                "highway" = "trunk_link",
                "highway" = "primary_link",
                "highway" = "secondary_link",
                "highway" = "tertiary_link",
                "highway" = "living_street",
                "highway" = "service",
                "highway" = "pedestrian",
                "highway" = "track",
                "highway" = "bus_guideway",
                "highway" = "escape",
                "highway" = "raceway",
                "highway" = "road",
                "highway" = "busway",
                "highway" = "footway",
                "highway" = "bridleway",
                "highway" = "steps",
                "highway" = "corridor",
                "highway" = "path",
                "highway" = "cycleway",
                "footway" = "sidewalk",
                "footway" = "crossing",
                "footway" = "access_aisle",
                "footway" = "link",
                "footway" = "traffic_island",
                "footway" = "alley",
                "footway" = "lane",
                "cycleway" = "lane",
                "cycleway" = "opposite",
                "cycleway" = "opposite_lane",
                "cycleway" = "track",
                "cycleway" = "opposite_track",
                "cycleway" = "share_busway",
                "cycleway" = "opposite_share_busway",
                "cycleway" = "shared_lane",
                "busway" = "lane",
                "route" = "bicycle",
                "route" = "bus",
                "route" = "foot",
                "route" = "hiking",
                "route" = "horse",
                "route" = "mtb",
                "route" = "railway",
                "route" = "road",
                "route" = "running",
                "route" = "train",
                "route" = "tracks")
  q <- osmdata::opq(bb) |>
    osmdata::add_osm_features(feats) |>
    osmdata::osmdata_sf()
  lns <- q$osm_lines |>
    terra::vect() |>
    terra::crop(e)
  out_file <- file.path(out_dir, "transportation.shp")
  terra::writeVector(lns, out_file, overwrite = T)
  return(lns)
}
