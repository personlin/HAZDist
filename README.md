[![Build Status](https://travis-ci.org/personlin/HAZDist.svg?branch=master)](https://travis-ci.org/personlin/HAZDist)

# HAZDist

`HAZDist` provides a set of R wrapper functions around the distance
calculation routines from Norman Abrahamson's
[HAZ](https://github.com/abrahamson/HAZ) Fortran code.  The package is
useful when computing fault–to–site distances required by ground motion
prediction equations (GMPEs).

## Features

* **SetFltBottom** – compute the bottom coordinates of a fault plane.
* **ConvertCoordinates2** – convert geographic coordinates to a local
  Cartesian system.
* **calcfltgrid** – create a fault grid used in distance calculations.
* **CalcDist** – return common GMPE distance measures such as `Rjb`,
  `Rrup`, `Rx`, `Ry` and `Ry0`.
* **getfaultcoord** – extract fault geometry information from a spatial
  data frame.
* **getfaultdist** – high level helper to compute distances from a fault
  to one or more sites.

## Installation

```r
install.packages("devtools")
# install from GitHub
devtools::install_github("personlin/HAZDist")
```

## Example

The package ships with a sample spatial data set `FT` describing the
Shanchiao fault.  The following code demonstrates how to compute the
fault–to–site distances for a single location:

```r
library(HAZDist)

# load example fault data
data(FT)

# extract geometry for the first fault
fault <- getfaultcoord(
  faultlines   = FT,
  id           = 1,
  c.id         = "ID",
  c.fault.top  = "DEPTH1_KM",
  c.fault.thick= "WIDTH_KM",
  c.dip        = "FAULT_DIP1"
)

# coordinates of the site (longitude, latitude)
site <- data.frame(x = 121.001, y = 24.7794)

# compute distance metrics
getfaultdist(fault, site, seismoDepth = 0)
```

## Data

`FT` is a spatial data frame containing the geometry and attributes of
part of the Shanchiao fault system.  The data include columns such as
`FAULT_NAME`, `FAULT_TYPE`, `LENGTH_KM`, `FAULT_DIP1`, `DEPTH1_KM` and
`WIDTH_KM`.

## License

HAZDist is released under the GPL (>= 3) license.
