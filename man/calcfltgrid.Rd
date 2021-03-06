% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcfltgrid.r
\name{calcfltgrid}
\alias{calcfltgrid}
\title{Function for converting fault into grid from Norman Abrahamson's HAZ program.}
\usage{
calcfltgrid(ft.xFlt, ft.yFlt, ft.zFlt, x0 = 0, y0 = 0, z0 = 0,
  xstep = 1)
}
\arguments{
\item{ft.xFlt}{Fault grid x(nDD, npt), Numeric.}

\item{ft.yFlt}{Fault grid y(nDD, npt), Numeric.}

\item{ft.zFlt}{Fault grid z(nDD, npt), Numeric.}

\item{x0}{grid origin x.}

\item{y0}{grid origin y.}

\item{z0}{grid origin z.}

\item{xstep}{xstep in fault input file.}
}
\value{
A list will be return, including fltGrid_x, fltGrid_y, fltGrid_z, fltGrid_x1, fltGrid_y1, fltGrid_z1,
fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4,
fltGrid_Rrup, fltGrid_Rjb
}
\description{
\code{calcfltgrid} returns xyz.
}
\examples{
calcfltgrid(ft.xFlt=t(array(c(-8.214971, -7.930996, -2.216696, 15.31437, 22.89828, 24.87581,
                      -1.370474, -1.087678,  4.623962, 22.15565, 29.74227, 31.71796), c(6,2))),
            ft.yFlt=t(array(c(-0.8335323,  4.558042, 15.116348, 27.63519, 32.26784, 34.74363,
                      -7.2026749, -1.810954,  8.750206, 21.27779, 25.91442, 28.39125), c(6,2))),
            ft.zFlt=(array(c( 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                      16.17, 16.17, 16.17, 16.17, 16.17, 16.17), c(6,2))), x0=0, y0=0, z0=0, xstep=1.0)

}
