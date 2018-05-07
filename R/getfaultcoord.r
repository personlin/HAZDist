#' Function for get fault coordinates from spatial data.frame.
#'
#' \code{getfaultcoord} returns coordinates of fault from spatial data frame.
#'
#' @param faultlines a spatial data frame of fault.
#' @param id ID from fault.
#' @param c.id Column name of the fault id.
#' @param c.fault.top Column name of depth of top of fault(km).
#' @param c.fault.thick Column name of the fault thick(km) (seismogenic thickness).
#' @param c.dip Column name of the fault dip(degree).
#'
#' @return A list will be return, including dip, faultwidth, fZ(,), fLat(,), fLong(,), nDD, pt.line.
#'
#' @examples
#' getfaultcoord(FT, 1, "ID", "WIDTH_KM", "FAULT_DIP1")
#'
#' @export
getfaultcoord <- function(faultlines, id, c.id, c.fault.top, c.fault.thick, c.dip){
  # fixed the R CMD Check
  ft <- NULL
  dip <- NULL
  faultwidth <- NULL

  # select the fault with id
  eval(parse(text=paste0("ft <- subset(faultlines, ",c.id, " == ",id,")")))

  # set parameters
  MAX_FLT = 1
  MAX_DD = 12
  MAX_SEG = 300
  eval(parse(text=paste0("ft.top <- ft$",c.fault.top)))
  # set array for data
  fLong <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fLat <- array(0, c(MAX_FLT,MAX_DD,MAX_SEG))
  fZ <- array(ft.top, c(MAX_FLT,MAX_DD,MAX_SEG))

  # get coordinate of fault line
  pt.line <- length(ft@lines[[1]]@Lines[[1]]@coords[,1])
  fLong[1,1,1:pt.line] <- ft@lines[[1]]@Lines[[1]]@coords[,1]
  fLat[1,1,1:pt.line] <- ft@lines[[1]]@Lines[[1]]@coords[,2]

  # get fault parameters
  ft.data <- ft@data
  eval(parse(text=paste0("faultwidth <- ft.data$",c.fault.thick)))
  #faultwidth <- ft@data$WIDTH_KM[1]
  eval(parse(text=paste0("dip <- ft.data$",c.dip)))
  #dip <- ft@data$FAULT_DIP1[1]

  nfp = c(pt.line,0)
  nDD = c(1,0)

  retvals <- .Fortran("SetFltBottom", iCoor=as.integer(1), iFlt=as.integer(1), nfp=as.integer(nfp),
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
              nDD = retvals$nDD[1],
              pt.line=pt.line)
  return(out)

}
