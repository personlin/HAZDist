#' Function for get fault distances.
#'
#' \code{getfaultdist} returns distances of the fault.
#'
#' @param fault fault coordinates with list name of fZ, fLat, fLong, nDD, pt.line.
#' @param site data.frame with x and y.
#' @param seismoDepth seismogenic depth (fault minDepth in fault file).
#'
#' @return A data.frame will be return, including column of HWFlag, hypoDepth, ZTOR, RupWidth, RupLen,
#' seismoDepth, distJB, distRup, distSeism, distepi, disthypo, Rx, Ry, Ry0.
#'
#' @examples
#' site <- data.frame(x=121.001, y=24.7794)
#' fault <- getfaultcoord(FT, 1, "ID", "WIDTH_KM", "FAULT_DIP1")
#' getfaultdist(fault, site, seismoDepth=0)
#'
#'
#' @export
getfaultdist <- function(fault, site, seismoDepth=0){
  # Convert Coordinates
  xzy <- ConvertCoordinates2(site$x, site$y, fault$fZ, fault$fLat, fault$fLong, fault$nDD, fault$pt.line)
  # calcFltGrid
  grid <- calcfltgrid(xzy$xFlt, xzy$yFlt, xzy$zFlt, xzy$x0, xzy$y0, xzy$z0)
  # CalcDist
  dist <- CalcDist(grid$faultW, grid$faultLen, seismoDepth=seismoDepth, grid$fltGrid_x, grid$fltGrid_y, grid$fltGrid_z,
                   grid$fltGrid_x1, grid$fltGrid_y1, grid$fltGrid_z1,
                   grid$fltGrid_x2, grid$fltGrid_y2, grid$fltGrid_z2,
                   grid$fltGrid_x3, grid$fltGrid_y3, grid$fltGrid_z3,
                   grid$fltGrid_x4, grid$fltGrid_y4, grid$fltGrid_z4,
                   grid$fltGrid_Rrup, grid$fltGrid_Rjb)
  return(as.data.frame(dist))
}
