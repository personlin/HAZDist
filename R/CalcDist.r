#' Distance calculation function from Norman Abrahamson's HAZ program.
#'
#' \code{CalcDist} computes the distance parameters for the GMPEs.
#'                 (e.g. Rjb, Rrup, ZTOR, Rx, Ry, Ry0)
#'
#' @param faultW Fault width
#' @param faultLen Fault length.
#' @param seismoDepth seismogenic depth (fault minDepth in fault file).
#' @param ft.fltGrid_x Fault grid coordinate x.
#' @param ft.fltGrid_y Fault grid coordinate y.
#' @param ft.fltGrid_z Fault grid coordinate z.
#' @param ft.fltGrid_x1 Fault grid coordinate x1.
#' @param ft.fltGrid_y1 Fault grid coordinate y1.
#' @param ft.fltGrid_z1 Fault grid coordinate z1.
#' @param ft.fltGrid_x2 Fault grid coordinate x2.
#' @param ft.fltGrid_y2 Fault grid coordinate y2.
#' @param ft.fltGrid_z2 Fault grid coordinate z2.
#' @param ft.fltGrid_x3 Fault grid coordinate x3.
#' @param ft.fltGrid_y3 Fault grid coordinate y3.
#' @param ft.fltGrid_z3 Fault grid coordinate z3.
#' @param ft.fltGrid_x4 Fault grid coordinate x4.
#' @param ft.fltGrid_y4 Fault grid coordinate y4.
#' @param ft.fltGrid_z4 Fault grid coordinate z4.
#' @param ft.fltGrid_Rrup Fault grid rupture distance.
#' @param ft.fltGrid_Rjb Fault grid JB distance.
#' @param ystep ystep in fault input file.
#'
#' @return A list will be returned, includingHWFlag, hypoDepth, ZTOR, RupWidth, RupLen, seismoDepth,
#' distJB, distRup, distSeism, distepi, disthypo, Rx, Ry, Ry0.
#'
#' @examples
#' CalcDist(14.14, 8.72, seismoDepth=0, array(0, c(4,4)), array(0, c(4,4)), array(0, c(4,4)),
#' array(1, c(4,4)), array(1, c(4,4)), array(1, c(4,4)),
#' array(2, c(4,4)), array(2, c(4,4)), array(2, c(4,4)),
#' array(3, c(4,4)), array(3, c(4,4)), array(3, c(4,4)),
#' array(4, c(4,4)), array(4, c(4,4)), array(4, c(4,4)),
#' array(2, c(4,4)), array(2, c(4,4)), ystep=1.0)
#'
#'
#' @export
CalcDist <- function(faultW, faultLen, seismoDepth=0, ft.fltGrid_x, ft.fltGrid_y, ft.fltGrid_z,
                     ft.fltGrid_x1, ft.fltGrid_y1, ft.fltGrid_z1,
                     ft.fltGrid_x2, ft.fltGrid_y2, ft.fltGrid_z2,
                     ft.fltGrid_x3, ft.fltGrid_y3, ft.fltGrid_z3,
                     ft.fltGrid_x4, ft.fltGrid_y4, ft.fltGrid_z4,
                     ft.fltGrid_Rrup, ft.fltGrid_Rjb, ystep=1.0){
  # subroutine CalcDist (sourceType, pscorflag, nFltGrid, n1AS, iLocX, iLocY, n2AS,
  #                      iFltWidth, iFlt, iMag, ystep, grid_top, RupWidth, RupLen, r_horiz, seismoDepth,
  #                      fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1, fltgrid_y1, fltgrid_z1,
  #                      fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3, fltgrid_z3,
  #                      fltgrid_x4, fltgrid_y4, fltgrid_z4, fltgrid_Rrup, fltgrid_Rjb, dip, dipS7,
  #                      distS7, HWFlag, n1, n2, icellRupstrike, icellRupdip, hypoDepth, distJB,
  #                      distRup, ZTOR, distSeismo, distepi, disthypo, dipavgd, Rx, Ry, Ry0)

  MAX_FLT = 1
  MAX_SEG=300
  MAXFLT_DD=2000
  MAXFLT_AS=2000
  MAX_WIDTH=15
  MAX_S7=70000
  MAX_GRID = 32000

  sourcetype = 1  # 1 for fault sourcetype
  nFltGrid = dim(ft.fltGrid_x)
  ystep = ystep # ystep in fault input file

  grid_top=array(0, c(MAX_FLT,MAX_GRID))

  fltGrid_x=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_y=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_z=array(0, c(MAXFLT_DD,MAXFLT_AS))
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

  fltGrid_a=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_w=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_Rrup=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_Rjb=array(0, c(MAXFLT_DD,MAXFLT_AS))
  fltGrid_fLen = array(0, c(MAXFLT_DD,MAXFLT_AS))

  dip = array(0, c(MAX_FLT,MAX_WIDTH,MAX_SEG))
  distS7 = array(0, c(MAX_FLT,MAX_S7))
  dipS7 = array(0, c(MAX_FLT,MAX_S7))

  RupWidth = faultW
  RupLen = faultLen

  fltGrid_x[1:nFltGrid[1],1:nFltGrid[2]]  = ft.fltGrid_x
  fltGrid_y[1:nFltGrid[1],1:nFltGrid[2]]  = ft.fltGrid_y
  fltGrid_z[1:nFltGrid[1],1:nFltGrid[2]]  = ft.fltGrid_z
  fltGrid_x1[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_x1
  fltGrid_y1[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_y1
  fltGrid_z1[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_z1
  fltGrid_x2[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_x2
  fltGrid_y2[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_y2
  fltGrid_z2[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_z2
  fltGrid_x3[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_x3
  fltGrid_y3[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_y3
  fltGrid_z3[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_z3
  fltGrid_x4[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_x4
  fltGrid_y4[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_y4
  fltGrid_z4[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_z4
  fltGrid_Rrup[1:nFltGrid[1],1:nFltGrid[2]]=ft.fltGrid_Rrup
  fltGrid_Rjb[1:nFltGrid[1],1:nFltGrid[2]] =ft.fltGrid_Rjb

  retvals <- .Fortran("CalcDist_f", sourceType=as.integer(sourcetype), pscorflag=as.integer(1), nFltGrid=as.integer(nFltGrid),
                      n1AS=as.integer(1), iLocX=as.integer(1), iLocY=as.integer(1), n2AS=as.integer(nFltGrid[2]),
                      iFltWidth=as.integer(1), iFlt=as.integer(1), iMag=as.integer(1),
                      ystep=as.single(ystep), grid_top=as.single(grid_top), RupWidth=as.single(RupWidth), RupLen=as.single(RupLen),
                      r_horiz=as.single(0), seismoDepth=as.single(seismoDepth),
                      fltgrid_x=as.single(fltGrid_x), fltgrid_y=as.single(fltGrid_y), fltgrid_z=as.single(fltGrid_z),
                      fltgrid_x1=as.single(fltGrid_x1), fltgrid_y1=as.single(fltGrid_y1), fltgrid_z1=as.single(fltGrid_z1),
                      fltgrid_x2=as.single(fltGrid_x2), fltgrid_y2=as.single(fltGrid_y2), fltgrid_z2=as.single(fltGrid_z2),
                      fltgrid_x3=as.single(fltGrid_x3), fltgrid_y3=as.single(fltGrid_y3), fltgrid_z3=as.single(fltGrid_z3),
                      fltgrid_x4=as.single(fltGrid_x4), fltgrid_y4=as.single(fltGrid_y4), fltgrid_z4=as.single(fltGrid_z4),
                      fltgrid_Rrup=as.single(fltGrid_Rrup), fltgrid_Rjb=as.single(fltGrid_Rjb),
                      dip=as.single(dip), dipS7=as.single(dipS7), distS7=as.single(distS7),
                      HWFlag=as.integer(1), n1=as.integer(1), n2=as.integer(1), icellRupstrike=as.integer(1),
                      icellRupdip=as.integer(1),
                      hypoDepth=as.single(0), distJB=as.single(0), distRup=as.single(0), ZTOR=as.single(0),
                      distSeismo=as.single(0), distepi=as.single(0), disthypo=as.single(0), dipavgd=as.single(0),
                      Rx=as.single(0), Ry=as.single(0), Ry0=as.single(0))
  names(retvals) <- c("sourceType", "pscorflag", "nFltGrid", "n1AS", "iLocX", "iLocY", "n2AS", "iFltWidth",
                      "iFlt", "iMag", "ystep", "grid_top", "RupWidth", "RupLen", "r_horiz", "seismoDepth",
                      "fltgrid_x", "fltgrid_y", "fltgrid_z", "fltgrid_x1", "fltgrid_y1", "fltgrid_z1",
                      "fltgrid_x2", "fltgrid_y2", "fltgrid_z2", "fltgrid_x3", "fltgrid_y3", "fltgrid_z3",
                      "fltgrid_x4", "fltgrid_y4", "fltgrid_z4", "fltgrid_Rrup", "fltgrid_Rjb", "dip", "dipS7",
                      "distS7", "HWFlag", "n1", "n2", "icellRupstrike", "icellRupdip", "hypoDepth", "distJB",
                      "distRup", "ZTOR", "distSeismo", "distepi", "disthypo", "dipavgd", "Rx", "Ry", "Ry0")
  out <- list(HWFlag=retvals$HWFlag, hypoDepth=retvals$hypoDepth, ZTOR=retvals$ZTOR,
              RupWidth=retvals$RupWidth, RupLen=retvals$RupLen, seismoDepth=retvals$seismoDepth,
              distJB=retvals$distJB, distRup=retvals$distRup,
              distSeismo=retvals$distSeismo, distepi=retvals$distepi,  disthypo=retvals$disthypo,
              Rx=retvals$Rx, Ry=retvals$Ry, Ry0=retvals$Ry0)
  return(out)
}
