# ROKA - ROck slope Kinematic Analysis

[![MATLAB](https://img.shields.io/badge/MATLAB-R2018b%2B-e16737.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fapp11041698-orange.svg)](https://doi.org/10.3390/app11041698)

**ROKA** (ROck slope Kinematic Analysis) is a MATLAB-based algorithm designed to perform ad-hoc 3D kinematic stability analyses for rock slopes and cliffs using 3D point clouds and 3D fracture/discontinuity data.

By taking into account the local orientation of the slope alongside the position, dimension, and attitude of 3D mapped discontinuity planes, ROKA evaluates stability locally across all portions of a rock cliff. Additionally, the algorithm calculates real 3D discontinuity plane intersections to test for kinematic feasibility across multiple failure mechanisms.

---

## Key Features

- **Local Slope Orientation Analysis**: Computes kinematic susceptibility per point using local normal vectors from 3D point clouds.
- **Full Discontinuity Evaluation**: Evaluates individual 3D mapped discontinuity discs (Beacher's disc model).
- **Real Intersections Feasibility**: Tests real 3D intersections of mapped fracture planes for kinematic failure modes.
- **Multiple Failure Modes Tested**:
  - Planar Sliding (with or without lateral limits)
  - Flexural Toppling
  - Wedge Sliding
  - Direct Toppling
- **Vectorized & Parallelized Execution**: Leverages parallel computing for fast spatial calculation.
- **3D DXF Export**: Automatically generates 3D `.dxf` vector files for critical discontinuity planes and intersections.

---

## Requirements & Dependencies

To run ROKA, you need:

1. **MATLAB**
2. **Computer Vision Toolbox**: Required for plotting and spatial processing routines.
3. **Parallel Computing Toolbox**: Required for multi-core computation speedup.

---

## Input Data Formatting

ROKA requires two primary input files placed in the working directory:

### 1. Point Cloud File (`.txt`)
- **Format**: Space-delimited text file without header rows.
- **Column Order**: `x y z Nx Ny Nz`
  - `x, y, z`: Point coordinates.
  - `Nx, Ny, Nz`: Components of the normal vector for each point.
- **Coordinate System Recommendation**: A **local coordinate reference system** is strongly recommended over geographic coordinates (e.g., WGS84/UTM) to reduce digit length and simplify distance calculations.

### 2. Mapped Discontinuities File (`.xlsx`)
- **Format**: Excel spreadsheet (`.xlsx`).
- **Required Columns** (at least 9 headed columns, order-independent):
  - `Dip`: Discontinuity dip angle (degrees).
  - `DipDirection`: Discontinuity dip direction angle (degrees).
  - `Radius`: Discontinuity radius / Beacher's disc radius (meters).
  - `Xcenter`, `Ycenter`, `Zcenter`: 3D coordinates of the discontinuity center.
  - `Nx`, `Ny`, `Nz`: Components of the normal vector to the discontinuity plane.

#### Options to Generate the Discontinuity `.xlsx` File:

- **Method A (Included MATLAB Script)**:
  Run `dxf_plane_fit.m` (available in the main ROKA folder). This script fits 3D discs from a `.dxf` file containing discontinuity traces mapped on the 3D model and automatically exports the required Excel file.
  
- **Method B (CloudCompare Workflow)**:
  Map discontinuities in CloudCompare (e.g., via *Trace Polyline* or the *Compass* plugin by Thiele et al., 2017), fit planes (`Fit -> Plane`), and export plane info (`Batch Export -> Export plane info`).
  - *Note on Radius Calculation*: CloudCompare fits rectangular planes rather than discs. To convert rectangular height and length to the required Beacher's disc radius ($r$), apply Pythagoras' theorem:
    $$\text{radius} = \frac{\sqrt{\text{length}^2 + \text{height}^2}}{2}$$

---

## Calculation Parameters Definition

Before running the analysis, configure the calculation parameters inside the main code (`ROKAmain.m`):

| Parameter | Type / Unit | Description |
| :--- | :--- | :--- |
| `userSPr` | Meters (`m`) | Scan radius size (e.g., `0.15` for 15 cm). |
| `disccutoff` | Meters (`m`) | Cutoff value for intersection between discontinuity and scan-plane. Default is 90% of scan-plane diameter. |
| `latlimits` | Degrees (`°`) | Lateral limits angle for kinematic tests. |
| `frictionangle` | Degrees (`°`) | Friction angle of discontinuity surfaces. |
| `uselatlimits` | Binary (`1` / `0`) | Enforce lateral limits (`1` = active, `0` = inactive). |
| `export_AllIntersection` | Binary (`1` / `0`) | Export all real discontinuity intersection lines (`1` = yes, `0` = no). |

---

## How to Run

1. Open **MATLAB**.
2. Set the ROKA root directory as the MATLAB `Current Folder`.
3. Open `ROKAmain.m` in the Editor.
4. Set input file paths and adjust calculation parameters (`userSPr`, `disccutoff`, `frictionangle`, `latlimits`, etc.).
5. Execute the analysis by typing `ROKAmain` in the Command Window or clicking **Run** in the Editor toolbar.

---

## Output Data Description

All results are automatically saved in a directory named `KinAnalysis/` created within the input data path.

### 1. Evaluated Point Cloud (`.txt`)
An exported space-delimited, headed text file reporting the following columns per point:
```
X Y Z DipDir Dip AllDisc_Inters AllInt_inters PlanarSliding FlexuralToppling WedgeSliding DirectToppling
```
- `X Y Z`: Point coordinates.
- `DipDir Dip`: Local slope dip direction and dip angle.
- `AllDisc_Inters`: Number of discontinuity planes intersecting the scan neighborhood.
- `AllInt_inters`: Number of discontinuity intersections intersecting the scan neighborhood.
- `PlanarSliding`, `FlexuralToppling`, `WedgeSliding`, `DirectToppling`: Maximum local critical kinematic value among intersecting discontinuities for each failure mode.

### 2. Critical 3D Files (`.dxf`)
3D vector `.dxf` files representing critical discontinuity planes and intersections are stored in dedicated subfolders by failure mode:

- **Discontinuity Planes Path Example**: `KinAnalysis/Disc/PlanarSlidingLatLimits/`
  - Naming scheme: `{CriticalValue}_CriticDisc_{DiscIndex}.dxf`
  - *Example*: `100_CriticDisc_1.dxf` (Critical value: 100, Discontinuity plane index: 1).

- **Discontinuity Intersections Path Example**: `KinAnalysis/Inters/DirectToppling/`
  - Naming scheme: `{CriticalValue}_Inters_{IntersIndex}_Discs_{DiscIndex1}_{DiscIndex2}.dxf`
  - *Example*: `75_Inters_1_Discs_1_5.dxf` (Critical value: 75, Intersection index: 1 formed by discontinuity planes 1 and 5).

---

## Citation & Reference

If you use ROKA in your research, please cite the following paper:

> **Menegoni, N., Giordan, D., & Perotti, C. (2021).** ROKA: ROck slope Kinematic Analysis algorithm. *Applied Sciences*, 11(4), 1698.  
> DOI: [10.3390/app11041698](https://doi.org/10.3390/app11041698)

---

## Contact & Authors

- **Niccolò Menegoni** (Department of Earth and Environmental Sciences, University of Pavia)
  - Emails: `niccolo.menegoni01@universitadipavia.it`, `niccolo.menegoni@unipv.it`, `niccolo.menegoni@live.it`
- **Daniele Giordan** (CNR-IRPI)
  - Email: `daniele.giordan@irpi.cnr.it`
- **Cesare Perotti** (University of Pavia)
  - Email: `cesare.perotti@unipv.it`
