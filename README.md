# Oregon LiDAR Tree Age Analysis

This project calculates tree age from LiDAR-derived canopy height models using the Means & Sabin (1989) formula. The analysis tessellates the Siuslaw National Forest into 50-acre square grids and calculates age statistics for each grid cell.

## Features

- Downloads Canopy Height Model (CHM) data from Oregon DOGAMI's REST API
- Creates a square grid (50 acres per cell) covering the Siuslaw National Forest
- Detects individual trees using ForestTools' variable window filter
- Calculates tree age from height using the Means & Sabin (1989) formula
- Generates maps of tree age and height distribution

## Getting Started

### Prerequisites

- R 4.0.0 or higher
- Required R packages: sf, terra, dplyr, ForestTools

### Installation

1. Clone this repository
2. Install required R packages:

```R
install.packages(c("sf", "terra", "dplyr", "ForestTools"))
```

### Usage

To run the analysis:

```bash
Rscript oregon_tree_age_grid.R
```

## Output

The script produces the following outputs in the `output` directory:

- `siuslaw_boundary.gpkg`: The boundary of the Siuslaw National Forest
- `siuslaw_grid_sample.gpkg`: The sample of grid cells used for analysis
- `siuslaw_grid_visualization.png`: Visualization of the grid cells
- `siuslaw_tree_age.gpkg`: Grid cells with tree age and height statistics
- `siuslaw_age_map.png`: Map of tree age distribution

## How It Works

1. The script loads the Siuslaw National Forest boundary from the Forest Service shapefile
2. It creates a square grid (50 acres per cell) covering the forest
3. For each grid cell, it:
   - Downloads LiDAR data from Oregon DOGAMI
   - Detects individual trees using ForestTools
   - Calculates mean tree height and dominant height
   - Estimates tree age using the Means & Sabin formula
4. Results are saved as GeoPackage files and visualized as maps

## References

Means, J. E., & Sabin, T. E. (1989). Height growth and site index curves for Douglas-fir in the Siuslaw National Forest, Oregon. Western Journal of Applied Forestry, 4(4), 136-142. 