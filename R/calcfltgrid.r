#' Function for converting fault into grid from Norman Abrahamson's HAZ program.
#'
#' \code{calcfltgrid} returns xyz.
#'
#' @param ft.xFlt Fault grid x(nDD, npt), Numeric.
#' @param ft.yFlt Fault grid y(nDD, npt), Numeric.
#' @param ft.zFlt Fault grid z(nDD, npt), Numeric.
#' @param x0 grid origin x.
#' @param y0 grid origin y.
#' @param z0 grid origin z.
#' @param xstep xstep in fault input file.
#'
#' @return A list will be return, including fltGrid_x, fltGrid_y, fltGrid_z, fltGrid_x1, fltGrid_y1, fltGrid_z1,
#' fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4,
#' fltGrid_Rrup, fltGrid_Rjb
#'
#' @examples
#' calcfltgrid(ft.xFlt=t(array(c(-8.214971, -7.930996, -2.216696, 15.31437, 22.89828, 24.87581,
#'                       -1.370474, -1.087678,  4.623962, 22.15565, 29.74227, 31.71796), c(6,2))),
#'             ft.yFlt=t(array(c(-0.8335323,  4.558042, 15.116348, 27.63519, 32.26784, 34.74363,
#'                       -7.2026749, -1.810954,  8.750206, 21.27779, 25.91442, 28.39125), c(6,2))),
#'             ft.zFlt=(array(c( 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
#'                       16.17, 16.17, 16.17, 16.17, 16.17, 16.17), c(6,2))), x0=0, y0=0, z0=0, xstep=1.0)
#'
#' @export
calcfltgrid <- function(ft.xFlt, ft.yFlt, ft.zFlt, x0=0, y0=0, z0=0, xstep=1.0){
  # subroutine calcFltGrid ( xFlt, yFlt, zFlt, npts, nDD, fltGrid_x, fltGrid_y,
  #                          fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
  #                          fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, faultW,
  #                          step, fltGrid_fLen, fltGrid_x1, fltGrid_y1,
  #                          fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3,
  #                          fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )

  npts = dim(ft.xFlt)[2]
  nDD = dim(ft.xFlt)[1]

  MAX_DD = 12
  MAX_SEG = 300
  MAXFLT_DD=2000
  MAXFLT_AS=2000

  xFlt = array(0, c(MAX_DD, MAX_SEG))
  yFlt = array(0, c(MAX_DD, MAX_SEG))
  zFlt = array(0, c(MAX_DD, MAX_SEG))

  xFlt[1:nDD,1:npts] <- ft.xFlt
  yFlt[1:nDD,1:npts] <- ft.yFlt
  zFlt[1:nDD,1:npts] <- ft.zFlt

  fltGrid_x=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z=array(0, c(MAXFLT_DD,MAXFLT_AS))
  # nfltGrid=c(0,0)
  fltGrid_a=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_w=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_Rrup=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_Rjb=array(0, c(MAXFLT_DD,MAXFLT_AS))
  step=xstep  # xstep in fault input file
  fltGrid_fLen = array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_x1=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y1=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z1=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_x2=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y2=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z2=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_x3=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y3=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z3=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_x4=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y4=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z4=array(0, c(MAXFLT_DD,MAXFLT_AS))

  retvals <- .Fortran("calcFltGrid_f", xFlt=as.single(xFlt), yFlt=as.single(yFlt), zFlt=as.single(zFlt),
                      npts=as.integer(npts), nDD=as.integer(nDD),
                      fltGrid_x=as.single(fltGrid_x), fltGrid_y=as.single(fltGrid_y), fltGrid_z=as.single(fltGrid_z),
                      nfltGrid=as.integer(c(0,0)), fltGrid_a=as.single(fltGrid_a),
                      fltGrid_w=as.single(fltGrid_w),
                      x0=as.single(0), y0=as.single(0), z0=as.single(0),
                      fltGrid_Rrup=as.single(fltGrid_Rrup), fltGrid_Rjb=as.single(fltGrid_Rrup),
                      faultArea=as.single(0), faultLen=as.single(0), faultW=as.single(0),
                      step=as.single(step), fltGrid_fLen=as.single(fltGrid_fLen),
                      fltGrid_x1=as.single(fltGrid_x1), fltGrid_y1=as.single(fltGrid_y1), fltGrid_z1=as.single(fltGrid_z1),
                      fltGrid_x2=as.single(fltGrid_x2), fltGrid_y2=as.single(fltGrid_y2), fltGrid_z2=as.single(fltGrid_z2),
                      fltGrid_x3=as.single(fltGrid_x3), fltGrid_y3=as.single(fltGrid_y3), fltGrid_z3=as.single(fltGrid_z3),
                      fltGrid_x4=as.single(fltGrid_x4), fltGrid_y4=as.single(fltGrid_y4), fltGrid_z4=as.single(fltGrid_z4))
  names(retvals) <- c("xFlt", "yFlt", "zFlt", "npts", "nDD", "fltGrid_x", "fltGrid_y",
                      "fltGrid_z", "nfltGrid", "fltGrid_a", "fltGrid_w", "x0", "y0", "z0",
                      "fltGrid_Rrup", "fltGrid_Rjb", "faultArea", "faultLen", "faultW",
                      "step", "fltGrid_fLen", "fltGrid_x1", "fltGrid_y1",
                      "fltGrid_z1", "fltGrid_x2", "fltGrid_y2", "fltGrid_z2", "fltGrid_x3",
                      "fltGrid_y3", "fltGrid_z3", "fltGrid_x4", "fltGrid_y4", "fltGrid_z4" )
  nfltGrid <- retvals$nfltGrid
  retvals$fltGrid_x = array(retvals$fltGrid_x, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_y = array(retvals$fltGrid_y, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_z = array(retvals$fltGrid_z, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_x1 =array(retvals$fltGrid_x1, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_y1 =array(retvals$fltGrid_y1, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_z1 =array(retvals$fltGrid_z1, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_x2 =array(retvals$fltGrid_x2, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_y2 =array(retvals$fltGrid_y2, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_z2 =array(retvals$fltGrid_z2, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_x3 =array(retvals$fltGrid_x3, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_y3 =array(retvals$fltGrid_y3, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_z3 =array(retvals$fltGrid_z3, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_x4 =array(retvals$fltGrid_x4, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_y4 =array(retvals$fltGrid_y4, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_z4 =array(retvals$fltGrid_z4, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_a = array(retvals$fltGrid_a, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_w = array(retvals$fltGrid_w, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_Rrup =array(retvals$fltGrid_Rrup, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_Rjb = array(retvals$fltGrid_Rjb, c(MAXFLT_DD,MAXFLT_AS))
  retvals$fltGrid_fLen =array(retvals$fltGrid_fLen, c(MAXFLT_DD,MAXFLT_AS))

  out <- list(fltGrid_x = retvals$fltGrid_x[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_y = retvals$fltGrid_y[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_z = retvals$fltGrid_z[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_x1 = retvals$fltGrid_x1[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_y1 = retvals$fltGrid_y1[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_z1 = retvals$fltGrid_z1[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_x2 = retvals$fltGrid_x2[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_y2 = retvals$fltGrid_y2[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_z2 = retvals$fltGrid_z2[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_x3 = retvals$fltGrid_x3[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_y3 = retvals$fltGrid_y3[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_z3 = retvals$fltGrid_z3[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_x4 = retvals$fltGrid_x4[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_y4 = retvals$fltGrid_y4[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_z4 = retvals$fltGrid_z4[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_Rrup = retvals$fltGrid_Rrup[1:nfltGrid[1],1:nfltGrid[2]],
              fltGrid_Rjb = retvals$fltGrid_Rjb[1:nfltGrid[1],1:nfltGrid[2]],
              faultArea=retvals$faultArea, faultLen=retvals$faultLen, faultW=retvals$faultW,
              nfltGrid=retvals$nfltGrid)
  return(out)
}
