#' HAZDist: Distance calculation function using Norman Abrahamson's HAZ fortran code.
#'
#' The HAZDist package provides distance calculation functions:
#' SetFltBottom, ConvertCoordinates2, calcfltgrid, CalcDist, getfaultcoord, getfaultdist
#'
#' @docType package
#' @name HAZDist
#' @aliases HAZDist-package
#' @useDynLib HAZDist, .registration = TRUE
"_PACKAGE"

.onUnload <- function(libpath) {
  library.dynam.unload("HAZDist", libpath)
}
