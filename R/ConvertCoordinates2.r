#' Function for converting latitude/longitude to x/y coordinates on the basis of a spherical earth
#' from Norman Abrahamson's HAZ program.
#'
#' \code{ConvertCoordinates2} returns fault coordinates in grid(km).
#'
#' @param siteX Site coordinate X (longitude).
#' @param siteY Site coordinate Y (latitude).
#' @param ft.Z Fault line coordinate (top of fault).
#' @param ft.Lat Fault line coordinates (latitude).
#' @param ft.Long Fault line coordinate (longitude).
#' @param nDD number of points on down dip direction.
#' @param pt.line number of points on fault line.
#'
#' @return A list will be return, including siteX, siteY, xFlt(,), yFlt(,), zFlt(,), x0, y0, z0.
#'
#' @examples
#' ConvertCoordinates2(121.5, 25, array(c(0,0),c(1,1,2)),
#'                     array(c(24.74554,24.7901),c(1,1,2)),
#'                     array(c(120.9609, 121.0747),c(1,1,2)), 2, 2)
#' ConvertCoordinates2(121.0, 24.5, array(c(0,0),c(1,1,2)),
#'                     array(c(24.74554,24.7901),c(1,1,2)),
#'                     array(c(120.9609, 121.0747),c(1,1,2)), 2, 2)
#'
#' @export
ConvertCoordinates2 <- function(siteX, siteY, ft.Z, ft.Lat, ft.Long, nDD, pt.line){
  # subroutine ConvertCoordinates2 (nfp, iFlt, iCoor, grid_n, sourceType, nDD,
  #                                 siteX, siteY, fLat, fLong, fZ, grid_lat,
  #                                 grid_long, grid_dlat, grid_dlong, nPts1, xFlt,
  #                                 yFlt, zFlt, grid_x, grid_y, grid_dx, grid_dy,
  #                                 x0, y0, z0)
  # nfp = number of points for fault line
  #
  iflt = 1
  icoor = 1       # for lon, lat coordinates
  grid_n = 0      # for sourcetype = 3 or 4
  sourcetype = 1  # 1 for fault sourcetype

  MAX_FLT = 1
  MAX_DD = 12
  MAX_SEG = 300
  MAX_GRID = 32000

  grid_lat = array(0, c(MAX_FLT, MAX_GRID))
  grid_long = array(0, c(MAX_FLT, MAX_GRID))
  grid_dlat = 0.1
  grid_dlong = 0.1
  nfp = c(pt.line,0)
  nPts1 = pt.line

  fLong <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fLat <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fZ <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  xFlt = array(0, c(MAX_DD, MAX_SEG))
  yFlt = array(0, c(MAX_DD, MAX_SEG))
  zFlt = array(0, c(MAX_DD, MAX_SEG))

  fLong[1,1:nDD,1:pt.line] <- ft.Long
  fLat[1,1:nDD,1:pt.line] <- ft.Lat
  fZ[1,1:nDD,1:pt.line] <- ft.Z

  retvals <- .Fortran("ConvertCoordinates2_f", nfp=as.integer(nfp), iFlt=as.integer(iflt), iCoor=as.integer(icoor),
                      grid_n=as.integer(grid_n), sourceType=as.integer(sourcetype), nDD=as.integer(nDD),
                      siteX=as.single(siteX), siteY=as.single(siteY),
                      fLat=as.single(fLat), fLong=as.single(fLong), fZ=as.single(fZ),
                      grid_lat=as.single(grid_lat), grid_long=as.single(grid_long),
                      grid_dlat=as.single(grid_dlat), grid_dlong=as.single(grid_dlong), nPts1=as.integer(nPts1),
                      xFlt=as.single(xFlt), yFlt=as.single(yFlt), zFlt=as.single(zFlt),
                      grid_x=as.single(0), grid_y=as.single(0), grid_dx=as.single(0), grid_dy=as.single(0),
                      x0=as.single(0), y0=as.single(0), z0=as.single(0))
  names(retvals) <- c("nfp", "iFlt", "iCoor", "grid_n", "sourceType", "nDD",
                      "siteX", "siteY", "fLat", "fLong", "fZ", "grid_lat", "grid_long",
                      "grid_dlat", "grid_dlong", "nPts1", "xFlt", "yFlt", "zFlt", "grid_x", "grid_y",
                      "grid_dx", "grid_dy", "x0", "y0", "z0")
  retvals$iFlt <- NULL
  retvals$iCoor <- NULL
  retvals$grid_n <- NULL
  retvals$sourceType <- NULL
  retvals$grid_lat <- NULL
  retvals$grid_long <- NULL
  retvals$grid_dlat <- NULL
  retvals$grid_dlong <- NULL
  retvals$grid_x <- NULL
  retvals$grid_y <- NULL
  retvals$grid_dx <- NULL
  retvals$grid_dy <- NULL
  retvals$xFlt <- array(retvals$xFlt, c(MAX_DD, MAX_SEG))
  retvals$yFlt <- array(retvals$yFlt, c(MAX_DD, MAX_SEG))
  retvals$zFlt <- array(retvals$zFlt, c(MAX_DD, MAX_SEG))
  out <- list(siteX = retvals$siteX, siteY = retvals$siteY,
              xFlt = retvals$xFlt[1:2,1:pt.line], yFlt = retvals$yFlt[1:2,1:pt.line],
              zFlt = retvals$zFlt[1:2,1:pt.line],
              x0 = retvals$x0, y0 = retvals$y0, z0 = retvals$z0)
  return(out)
}
