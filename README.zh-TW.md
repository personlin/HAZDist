[![R-CMD-check](https://github.com/personlin/HAZDist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/personlin/HAZDist/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# HAZDist

**HAZDist** 是 Norman Abrahamson 的 [HAZ](https://github.com/abrahamson/HAZ) Fortran 距離計算程式的 R 包裝套件，提供計算地震動衰減方程式（GMPE）所需的有限斷層至測站距離參數，適用於機率性震災危害度分析（PSHA）。

> For English documentation, see [README.md](README.md)

## 功能

| 函式 | 說明 |
|---|---|
| `getfaultcoord()` | 從 `sf` 物件萃取斷層幾何資訊（頂部、底部坐標） |
| `getfaultdist()` | 高階函式：輸入斷層與測站坐標，一次返回所有距離參數 |
| `SetFltBottom()` | 根據傾角與斷層寬度計算斷層底部坐標 |
| `ConvertCoordinates2()` | 將地理坐標（經緯度）轉換成以斷層走向為基準的局部直角坐標系 |
| `calcfltgrid()` | 建立斷層有限元素格網，供距離計算使用 |
| `CalcDist()` | 呼叫 HAZ Fortran 核心，計算 Rjb、Rrup、ZTOR、Rx、Ry、Ry0 等距離量 |

## 系統需求

- **R** ≥ 4.1.0
- **Fortran 編譯器** — 從原始碼安裝時需有 gfortran
  - macOS：`brew install gcc`（透過 [Homebrew](https://brew.sh/)）
  - Windows：安裝 [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  - Linux（Debian/Ubuntu）：`sudo apt install gfortran`
- **GNU make**
- 相依 R 套件：[`sf`](https://r-spatial.github.io/sf/)

## 安裝

```r
# 先安裝相依套件
install.packages("sf")

# 從 GitHub 安裝 HAZDist
remotes::install_github("personlin/HAZDist")
# 或使用 pak
pak::pak("personlin/HAZDist")
```

## 快速上手

套件內附範例斷層資料集 `FT`，包含山腳斷層系統（Shanchiao Fault）的空間幾何與屬性資料，以 `sf` 物件格式儲存。

```r
library(HAZDist)

# 載入內建斷層資料（sf 物件）
data(FT)

# 步驟 1：萃取第一條斷層分段的幾何資訊
fault <- getfaultcoord(
  faultlines    = FT,
  id            = 1,
  c.id          = "ID",
  c.fault.top   = "DEPTH1_KM",   # 斷層頂部深度（公里）
  c.fault.thick = "WIDTH_KM",    # 順傾斜方向寬度（公里）
  c.dip         = "FAULT_DIP1"   # 斷層傾角（度）
)

# 步驟 2：定義測站位置（WGS84 經緯度）
site <- data.frame(x = 121.001, y = 24.7794)

# 步驟 3：計算所有距離參數
result <- getfaultdist(fault, site, seismoDepth = 0)
print(result)
```

### `getfaultdist()` 返回值說明

| 欄位 | 說明 |
|---|---|
| `distJB` | Joyner–Boore 距離 Rjb（公里） |
| `distRup` | 最短斷層面距離 Rrup（公里） |
| `Rx` | 垂直走向的水平距離（公里），正值為懸壁側 |
| `Ry` / `Ry0` | 沿走向距離（公里） |
| `ZTOR` | 斷層頂部深度（公里） |
| `HWFlag` | 懸壁效應旗標（1 = 懸壁側） |
| `hypoDepth` | 假設震源深度（公里） |
| `RupWidth` / `RupLen` | 破裂寬度／長度（公里） |
| `distSeism` | 地震深度距離（公里） |
| `distepi` / `disthypo` | 震央距離／震源距離（公里） |

## 資料集

`FT` 是山腳斷層系統的空間資料，以 `sf` LINESTRING 物件格式儲存。主要屬性欄位如下：

| 欄位 | 說明 |
|---|---|
| `ID` | 斷層分段編號 |
| `FAULT_NAME` | 斷層名稱 |
| `FAULT_TYPE` | 斷層類型（逆斷層、正斷層等） |
| `LENGTH_KM` | 斷層長度（公里） |
| `FAULT_DIP1` | 斷層傾角（度） |
| `DEPTH1_KM` | 斷層頂部深度（公里） |
| `WIDTH_KM` | 順傾斜方向寬度（公里） |

## 參考文獻

Abrahamson, N. A. (2000). *HAZ: A computer program for the probabilistic seismic hazard analysis*. Pacific Gas and Electric Company.

## 授權

HAZDist 依 [GPL (≥ 3)](LICENSE) 授權釋出。
