# Oregon LiDAR Tree Age Analysis
This 
This project calculates tree age from LiDAR-derived canopy height models using the Means & Sabin (1989) formula. The analysis tessellates the Siuslaw National Forest into 50-acre square grids and calculates age statistics for each grid cell.

## Features

- Downloads Canopy Height Model (CHM) data from Oregon DOGAMI's REST API
- Creates a square grid (2 acres per cell) covering the Siuslaw National Forest
  - With the later goal of updating to run for the entire state of Oregon
- Detects individual trees using ForestTools' variable window filter
- Calculates tree age from height using the Means & Sabin (1989) formula
- Generates maps of tree age and height distribution

## Getting Started

### Installation

1. Clone this repository
`git clone https://github.com/sailbot24/oregon_lidar`
2. Install required R packages:



### Usage

To run the analysis:


## Output

## How It Works

1. The script loads the Siuslaw National Forest boundary from the Forest Service shapefile
2. It creates a square grid (2 acres per cell) covering the forest
3. For each grid cell, it:
   - Downloads LiDAR data from Oregon DOGAMI
   - Detects individual trees using ForestTools
   - Calculates mean tree height
   - Estimates tree age using the Means & Sabin formula
4. Results are saved as GeoPackage files and visualized as maps

## References

Means, J. E., & Sabin, T. E. (1989). Height growth and site index curves for Douglas-fir in the Siuslaw National Forest, Oregon. Western Journal of Applied Forestry, 4(4), 136-142. 