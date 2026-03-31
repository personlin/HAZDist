[![R-CMD-check](https://github.com/personlin/HAZDist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/personlin/HAZDist/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# HAZDist

**HAZDist** is an R wrapper around the fault-distance routines from Norman Abrahamson's [HAZ](https://github.com/abrahamson/HAZ) Fortran program. It computes the source-to-site distance metrics required by Ground Motion Prediction Equations (GMPEs) for finite-fault seismic hazard analysis.

> 繁體中文說明請見 [README.zh-TW.md](README.zh-TW.md)

## Features

| Function | Description |
|---|---|
| `getfaultcoord()` | Extract fault geometry (top/bottom coordinates) from an `sf` object |
| `getfaultdist()` | High-level wrapper — given fault and site coordinates, returns all distance metrics |
| `SetFltBottom()` | Compute fault bottom coordinates from dip angle and down-dip width |
| `ConvertCoordinates2()` | Convert geographic coordinates (lon/lat) to a fault-strike local Cartesian system |
| `calcfltgrid()` | Build a finite-element fault grid for distance calculations |
| `CalcDist()` | Call the HAZ Fortran core to compute Rjb, Rrup, ZTOR, Rx, Ry, Ry0, etc. |

## System Requirements

- **R** ≥ 4.1.0
- **Fortran compiler** — gfortran must be available when installing from source
  - macOS: `brew install gcc` via [Homebrew](https://brew.sh/)
  - Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  - Linux: `sudo apt install gfortran` (Debian/Ubuntu)
- **GNU make**
- R package dependency: [`sf`](https://r-spatial.github.io/sf/)

## Installation

```r
# Install the sf dependency first
install.packages("sf")

# Install HAZDist from GitHub
remotes::install_github("personlin/HAZDist")
# or using pak
pak::pak("personlin/HAZDist")
```

## Quick Start

The package ships with `FT`, an example `sf` dataset of the Shanchiao Fault (山腳斷層) in Taiwan.

```r
library(HAZDist)

# Load the built-in fault dataset (sf object)
data(FT)

# Step 1: Extract geometry for fault segment 1
fault <- getfaultcoord(
  faultlines    = FT,
  id            = 1,
  c.id          = "ID",
  c.fault.top   = "DEPTH1_KM",   # depth to fault top (km)
  c.fault.thick = "WIDTH_KM",    # down-dip width (km)
  c.dip         = "FAULT_DIP1"   # dip angle (degrees)
)

# Step 2: Define the site location (WGS84 lon/lat)
site <- data.frame(x = 121.001, y = 24.7794)

# Step 3: Compute all distance metrics
result <- getfaultdist(fault, site, seismoDepth = 0)
print(result)
```

### Output fields from `getfaultdist()`

| Field | Description |
|---|---|
| `distJB` | Joyner–Boore distance Rjb (km) |
| `distRup` | Closest distance to the rupture surface Rrup (km) |
| `Rx` | Horizontal distance perpendicular to strike (km); positive on hanging-wall side |
| `Ry` / `Ry0` | Along-strike distance (km) |
| `ZTOR` | Depth to top of rupture (km) |
| `HWFlag` | Hanging-wall flag (1 = hanging-wall side) |
| `hypoDepth` | Assumed hypocentral depth (km) |
| `RupWidth` / `RupLen` | Rupture width / length (km) |
| `distSeism` | Seismogenic depth distance (km) |
| `distepi` / `disthypo` | Epicentral / hypocentral distance (km) |

## Dataset

`FT` is an `sf` LINESTRING object representing the Shanchiao Fault system. Key attribute columns:

| Column | Description |
|---|---|
| `ID` | Fault segment ID |
| `FAULT_NAME` | Fault name |
| `FAULT_TYPE` | Fault mechanism (reverse, normal, etc.) |
| `LENGTH_KM` | Segment length (km) |
| `FAULT_DIP1` | Dip angle (degrees) |
| `DEPTH1_KM` | Depth to fault top (km) |
| `WIDTH_KM` | Down-dip width (km) |

## Reference

Abrahamson, N. A. (2000). *HAZ: A computer program for the probabilistic seismic hazard analysis*. Pacific Gas and Electric Company.

## License

HAZDist is released under the [GPL (≥ 3)](LICENSE) license.
