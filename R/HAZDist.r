#' HAZDist: Distance calculation function using Norman Abrahamson's HAZ fortran code.
#'
#' The HAZDist package provides distance calculation functions:
#' SetFltBottom, ConvertCoordinates2, calcfltgrid, CalcDist, getfaultcoord, getfaultdist
#'
#'
#' @docType package
#' @name HAZDist
#' @useDynLib HAZDist, .registration = TRUE
#' @param libpath library path

.onUnload <- function(libpath){
  library.dynam.unload("HAZDist", libpath)
}
