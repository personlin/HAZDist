#' Function for setting bottom of fault for standard faults (source type 1)
#' from Norman Abrahamson's HAZ program.
#'
#' \code{SetFltBottom} returns coordinates on bottom of fault.
#'
#' @param dip Dip angle of the fault plane.
#' @param faultwidth Faultwidth(km) = faultThick1 = seismogenic thickness (fault widths in fault file).
#' @param ft.Z Fault line coordinate (top of fault).
#' @param ft.Lat Fault line coordinates (latitude).
#' @param ft.Long Fault line coordinate (longitude).
#' @param pt.line Number of points of the fault line.
#'
#' @return A list will be return, including dip, faultwidth, fZ(,), fLat(,), fLong(,), nDD.
#'
#' @examples
#' SetFltBottom(45, 15, c(0,0), c(24.74554,24.7901), c(120.9609, 121.0747), 2)
#' SetFltBottom(30, 15, c(0,0), c(24.74554,24.7901), c(120.9609, 121.0747), 2)
#'
#'
#' @export
SetFltBottom <- function(dip, faultwidth, ft.Z, ft.Lat, ft.Long, pt.line){
  # subroutine SetFltBottom (iCoor, iFlt, nfp, dip2, faultThick, fZ,
  #                          flat, flong, nDD)
  #

  icoor = 1 # for lon, lat coordinates
  iflt = 1 # only deal with one fault at one time
  nDD = 1
  nfp = c(pt.line,0)
  nDD = c(1,0)

  MAX_FLT = 1
  MAX_DD = 12
  MAX_SEG = 300

  fLong <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fLat <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fZ <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))

  fLong[1,1,1:pt.line] <- ft.Long
  fLat[1,1,1:pt.line] <- ft.Lat
  fZ[1,1,1:pt.line] <- ft.Z

  retvals <- .Fortran("SetFltBottom_f", iCoor=as.integer(icoor), iFlt=as.integer(iflt), nfp=as.integer(nfp),
                      dip2=as.single(dip), faultThick=as.single(faultwidth),
                      fZ=as.single(fZ), flat=as.single(fLat), flong=as.single(fLong),
                      nDD=as.integer(nDD))
  names(retvals) <- c("iCoor", "iFlt", "nfp", "dip", "faultwidth", "fZ", "fLat", "fLong", "nDD")
  retvals$iCoor <- NULL
  retvals$iFlt <- NULL
  retvals$nfp <- NULL
  retvals$fZ <- array(retvals$fZ, c(MAX_FLT,MAX_DD,MAX_SEG))
  retvals$fLat <- array(retvals$fLat, c(MAX_FLT,MAX_DD,MAX_SEG))
  retvals$fLong <- array(retvals$fLong, c(MAX_FLT,MAX_DD,MAX_SEG))
  out <- list(dip = retvals$dip, faultwidth = retvals$faultwidth, fLong = retvals$fLong[1,1:2,1:pt.line],
              fLat = retvals$fLat[1,1:2,1:pt.line], fZ = retvals$fZ[1,1:2,1:pt.line],
              nDD = retvals$nDD[1])
  return(out)
}
