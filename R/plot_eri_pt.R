#' Plot ERI for a Single Point
#'
#' @description
#' Uses the outputs of [eri_pt()] to create a map that displays all of the spatial features that were used to calculate ERI.
#'
#' @param elev SpatRaster. The LANDFIRE elevation dataset derived from successful execution of [download_lf()].
#' @param veg_type SpatRaster. The OpenStreetMap transportation data-augmented LANDFIRE vegetation type dataset derived from successful execution of [download_lf()], [download_trans()], and [reclass_veg()].
#' @param eri_pt_result list. The list returned from the successful execution of [eri_pt()].
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
#'
#' # plot out the results
#' plot_eri_pt(elev, veg_type_rc, eri_pt_result)
#' @export
plot_eri_pt <- function(elev, veg_type, eri_pt_result){
  dev.off()
  clipper <- terra::buffer(eri_pt$max_poly, 100)
  elev <- terra::crop(elev, clipper)
  veg_type <- terra::crop(veg_type, clipper)
  ll_ext <- terra::project(veg_type, "epsg:4326") |>
    terra::ext()
  lat_vals <- pretty(c(ll_ext[3][[1]], ll_ext[4][[1]]))
  lon_vals <- pretty(c(ll_ext[1][[1]], ll_ext[2][[1]]))
  grat <- terra::graticule(lon_vals, lat_vals, terra::crs(clipper)) |>
    terra::crop(terra::ext(veg_type))
  pt_max <- eri_pt$pt_max
  pt_max$rastval <- 1
  pix_max <- terra::rasterize(pt_max, veg_type, "rastval") |>
    terra::as.polygons()
  pt_min <- eri_pt$pt_min
  pt_min$rastval <- 1
  pix_min <- terra::rasterize(pt_min, veg_type, "rastval") |>
    terra::as.polygons()
  veg_type <- terra::classify(veg_type, matrix(c(0,1,0.25,2,0.5,3,1,4),
                                               byrow = T, ncol = 2))
  terra::coltab(veg_type) <- data.frame(values = c(1,2,3,4),
                                        cols = c(rgb(23,165,194,126, maxColorValue = 255),
                                                 rgb(29,133,40,126, maxColorValue = 255),
                                                 rgb(219,134,59,126, maxColorValue = 255),
                                                 rgb(224,217,72,126, maxColorValue = 255)))
  terra::plot(grat, mar = c(1.25,2.25,1.5,10))
  terra::plot(veg_type, legend = F, add = T)
  terra::lines(grat)
  terra::contour(elev, add = T, col = "white")
  terra::lines(x = c(terra::crds(start_pt)[,1],
                     terra::crds(pt_max)[,1]),
               y = c(terra::crds(start_pt)[,2],
                     terra::crds(pt_max)[,2]),
               lwd = 3, col = 4)
  terra::lines(x = c(terra::crds(start_pt)[,1],
                     terra::crds(pt_min)[,1]),
               y = c(terra::crds(start_pt)[,2],
                     terra::crds(pt_min)[,2]),
               lwd = 3, col = 2)
  terra::plot(pix_max, col = 4, add = T)
  terra::plot(pix_min, col = 2, add = T)
  terra::plot(eri_pt$max_poly, lwd = 3, bg = "black", lty = 2, add = T)
  terra::plot(eri_pt$act_poly, lwd = 3, bg = "black", add = T)
  terra::plot(eri_pt$start_pt, pch = 21, bg = "white", cex = 3, add = T)
  legend(x = terra::ext(veg_type)[1],
         y = terra::ext(veg_type)[4],
         legend = c(bquote(ERI[mean]==.(round(eri_pt$eri_mean, 2))),
                    bquote(ERI[max]==.(round(eri_pt$eri_max, 2))),
                    bquote(ERI[min]==.(round(eri_pt$eri_min, 2))),
                    bquote(ERI[azim]==.(round(eri_pt$eri_az))*degree)),
         x.intersp = 0,
         bg = "lightgray")
  terra::sbar(xy = c(terra::ext(veg_type)[1] +
                       0.05 * (terra::ext(veg_type)[2] - terra::ext(veg_type)[1]),
                     terra::ext(veg_type)[3] +
                       0.05 * (terra::ext(veg_type)[4] - terra::ext(veg_type)[3])),
              below = "meters", lwd = 4, halo = F)
  par(xpd = NA)
  legend(x = terra::ext(veg_type)[2] +
           0.02 * (terra::ext(veg_type)[2] - terra::ext(veg_type)[1]),
         y = mean(terra::ext(veg_type)[3:4]),
         legend = c("Tree", "Shrub", "Grass/NB", "Water",
                    "Start Pt", expression(ERI[max]~Pt),
                    expression(ERI[min]~Pt), "Actual Distance",
                    "Optimal Distance"),
         pch = c(15,15,15,15,21,22,22,NA,NA),
         lty = c(NA,NA,NA,NA,NA,NA,NA,1,2),
         lwd = c(NA,NA,NA,NA,NA,NA,NA,3,3),
         col = c(rgb(29,133,40,126, maxColorValue = 255),
                 rgb(219,134,59,126, maxColorValue = 255),
                 rgb(224,217,72,126, maxColorValue = 255),
                 rgb(23,165,194,126, maxColorValue = 255),
                 1,1,1,1,1),
         pt.bg = c(NA,NA,NA,NA,"white",4,2,NA,NA),
         pt.cex = 2,
         yjust = 0.5, bty = "n")
  ll_start <- terra::project(start_pt, "epsg:4326") |>
    terra::crds()
  title <- paste0("Start Point: ", round(ll_start[1,2], 6), "?, ", round(ll_start[1,1], 6), "?")
  text(x = mean(terra::ext(veg_type)[1:2]),
       y = terra::ext(veg_type)[4],
       labels = title, pos = 3, offset = 0.5, font = 2)
  box <- terra::ext(veg_type) |>
    terra::as.polygons() |>
    terra::plot(lwd = 2, col = NA, add = T)
  par(xpd = T)
}
