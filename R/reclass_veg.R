#' Reclassify Vegetation Data
#'
#' @description
#' Combines the LANDFIRE Existing Vegetation Type and OpenStreetMap transportation data to create a single raster dataset representing the land cover-based relative conductance for pedestrian travel.
#'
#' @details
#' * The conductance values modeled in this function were gleaned from a generalization of the pedestrian impedance factors used in the [WFDSS Estimated Ground Evacuation Time](https://wfdss.usgs.gov/wfdss_help/WFDSSHelp_Est_Grd_Medevac_Time.html) layer. They are as follows:
#'   - Grass-dominated/non-burnable/transportation: 1x
#'   - Shrub-dominated: 0.5x
#'   - Tree-dominated: 0.25x
#'   - Water: 0x (acts as a barrier)
#' * Note that this function will save the following TIF to `out_dir`: "veg_type_rc.tif".
#'
#' @param veg_type SpatRaster. The LANDFIRE Existing Vegetation Type dataset derived from successful execution of [download_lf()].
#' @param trans SpatVector. The OpenStreetMap transportation dataset derived from successful execution of [download_trans()].
#' @param out_dir character. The output directory to which the resulting reclassified raster will be saved.
#' @return A SpatRaster with cell values representing relative pedestrian conductance (see Details).
#' @examples
#' # read in aoi
#' aoi <- vect("C:/temp/study_area.shp")
#'
#' # download landfire data using a 15 minute evacuation time
#' lf <- download_lf(aoi, 15, "C:/temp")
#'
#' # access output
#' veg_type <- lf$veg_type
#'
#' # download transportation data using a 15 minute evacuation time
#' trans <- download_trans(aoi, 15, "C:/temp")
#'
#' # reclassify vegetation type
#' veg_type_rc <- reclass_veg(veg_type, trans, "C:/temp")
#' @export
reclass_veg <- function(veg_type, trans, out_dir){
  lut_file <- list.files(system.file("extdata", package = "eri"), full.names = TRUE)
  lf_df <- read.csv(lut_file) |>
    na.omit()
  vals_tree <- lf_df$VALUE[lf_df$EVT_ORDER == "Tree-dominated"]
  vals_shrb <- lf_df$VALUE[lf_df$EVT_ORDER == "Shrub-dominated"]
  vals_grnb <- lf_df$VALUE[lf_df$EVT_ORDER %in% c("No dominant lifeform",
                                                  "Herbaceous / Nonvascular-dominated",
                                                  "No Dominant Life Form",
                                                  "Herbaceous/Non-vascular",
                                                  "No Dominant Lifeform",
                                                  "Non-vegetated")]
  vals_watr <- lf_df$VALUE[lf_df$EVT_NAME == "Open Water"]
  vals_grnb <- vals_grnb[!vals_grnb %in% vals_watr]
  rc_mat <- matrix(c(vals_tree, vals_shrb, vals_grnb, vals_watr,
                     rep(0.25, length(vals_tree)),
                     rep(0.50, length(vals_shrb)),
                     rep(1.00, length(vals_grnb)),
                     rep(0, length(vals_watr))),
                   ncol = 2, byrow = F)
  evt_rc <- terra::classify(veg_type, rc_mat)
  trans$rastval <- 1
  trans <- terra::project(trans, evt_rc)
  trans_rc <- terra::rasterize(trans, evt_rc, "rastval")
  evt_rc <- terra::ifel(trans_rc == 1, 1, evt_rc)
  out_file <- file.path(out_dir, "veg_type_rc.tif")
  terra::writeRaster(evt_rc, out_file, overwrite = T)
  r <- terra::rast(out_file)
  return(r)
}
