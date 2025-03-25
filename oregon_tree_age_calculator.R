#' Oregon Tree Age Calculator
#' 
#' This script calculates tree age from LiDAR-derived canopy height models
#' using the Means & Sabin (1989) formula. The area is tessellated into
#' 50-acre hexagons, and age is calculated for each hexagon.
#' 
#' Author: Adams Geospatial
#' Date: Created on `r Sys.Date()`

# Load required libraries
library(sf)          # For spatial data handling
library(terra)       # For raster operations
library(dplyr)       # For data manipulation
library(ForestTools) # For tree detection
library(ggplot2)     # For visualization
library(tmap)        # For mapping
library(units)       # For unit conversion

# Source the Oregon CHM downloader function
source("oregon_chm_downloader.R")

#' Function to calculate Site Index (SI) from tree height and age
#' Based on Means & Sabin (1989)
calculate_site_index <- function(height, age) {
  # Parameters from Means & Sabin (1989)
  b0 <- 1.42956
  b1 <- 0.37900
  b2 <- 65.11
  b3 <- -0.022733
  b4 <- 1.27157

  # Compute Site Index (SI)
  SI <- 39.77 + (b0 / exp(age / 50) + b1 * exp(age / 300)) * 
        (height - b2 * (1 - exp(b3 * age))^b4)
  
  return(SI)
}

#' Function to estimate age from height using the Means & Sabin formula
#' This is an inverse of the site index formula, solved for age
#' We use an iterative approach to find the age that corresponds to a given height
estimate_age_from_height <- function(height, site_index = 40, max_age = 200, tolerance = 0.1) {
  # Convert height from feet to meters if needed
  height_m <- ifelse(height > 100, height * 0.3048, height)
  
  # Function to minimize: difference between observed height and predicted height for a given age
  height_diff <- function(age, height, si) {
    # Parameters from Means & Sabin (1989)
    b0 <- 1.42956
    b1 <- 0.37900
    b2 <- 65.11
    b3 <- -0.022733
    b4 <- 1.27157
    
    # Calculate predicted height for this age and SI
    predicted_si <- 39.77 + (b0 / exp(age / 50) + b1 * exp(age / 300)) * 
                   (height - b2 * (1 - exp(b3 * age))^b4)
    
    return(abs(predicted_si - si))
  }
  
  # Try a range of ages and find the one that minimizes the difference
  ages <- seq(10, max_age, by = 1)
  diffs <- sapply(ages, height_diff, height = height_m, si = site_index)
  best_age <- ages[which.min(diffs)]
  
  # If the best fit is still poor, return NA
  if (min(diffs) > tolerance) {
    warning("Could not find a good age estimate for height ", height_m, 
            " m. Best difference was ", min(diffs))
  }
  
  return(best_age)
}

#' Function to create hexagon grid covering an area of interest
#' Each hexagon will be approximately 50 acres
create_hexagon_grid <- function(aoi, area_acres = 50) {
  # Convert acres to square meters
  area_m2 <- area_acres * 4046.86
  
  # Calculate hexagon size (distance from center to vertex)
  # Area of a hexagon = 2.598 * size^2
  hex_size <- sqrt(area_m2 / 2.598)
  
  # Create hexagon grid
  hex_grid <- st_make_grid(aoi, cellsize = hex_size, square = FALSE)
  
  # Convert to sf object
  hex_grid_sf <- st_sf(id = 1:length(hex_grid), geometry = hex_grid)
  
  # Assign same CRS as AOI
  st_crs(hex_grid_sf) <- st_crs(aoi)
  
  # Clip to AOI
  hex_grid_clipped <- st_intersection(hex_grid_sf, st_union(aoi))
  
  return(hex_grid_clipped)
}

#' Function to process a single hexagon
#' Downloads CHM data, detects trees, and calculates age statistics
process_hexagon <- function(hexagon, id, save_dir = "output") {
  # Create output directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Download CHM for this hexagon
  chm_file <- file.path(save_dir, paste0("hex_", id, "_chm.tif"))
  
  tryCatch({
    # Download CHM data
    chm <- download_oregon_chm(hexagon, output_path = chm_file, mask = TRUE, quiet = TRUE)
    
    # Detect trees using variable window filter
    # Function for defining dynamic window size
    lin <- function(x){x * 0.05 + 0.6}
    
    # Detect treetops (minimum height 10 feet)
    ttops <- vwf(chm, winFun = lin, minHeight = 10)
    
    # If no trees detected, return empty result
    if (nrow(ttops) == 0) {
      return(data.frame(
        hex_id = id,
        mean_height = NA,
        mean_age = NA,
        tree_count = 0,
        dominant_height = NA,
        dominant_age = NA
      ))
    }
    
    # Calculate statistics
    mean_height <- mean(ttops$height, na.rm = TRUE)
    
    # Get dominant height (95th percentile)
    dominant_height <- quantile(ttops$height, 0.95, na.rm = TRUE)
    
    # Estimate ages
    mean_age <- estimate_age_from_height(mean_height)
    dominant_age <- estimate_age_from_height(dominant_height)
    
    # Return results
    return(data.frame(
      hex_id = id,
      mean_height = mean_height,
      mean_age = mean_age,
      tree_count = nrow(ttops),
      dominant_height = dominant_height,
      dominant_age = dominant_age
    ))
  }, error = function(e) {
    warning("Error processing hexagon ", id, ": ", e$message)
    return(data.frame(
      hex_id = id,
      mean_height = NA,
      mean_age = NA,
      tree_count = NA,
      dominant_height = NA,
      dominant_age = NA
    ))
  })
}

# Main execution function
run_tree_age_analysis <- function() {
  # Load National Forest boundaries
  nf_path <- "Veg_ShapeFiles/FS_National_Forests_Dataset_(US_Forest_Service_Proclaimed_Forests)/FS_National_Forests_Dataset_(US_Forest_Service_Proclaimed_Forests).shp"
  national_forests <- st_read(nf_path, quiet = TRUE)
  
  # Filter for Siuslaw National Forest
  siuslaw <- national_forests %>% 
    filter(FORESTNAME == "Siuslaw National Forest")
  
  # Create hexagon grid (50 acres each)
  message("Creating hexagon grid...")
  hex_grid <- create_hexagon_grid(siuslaw, area_acres = 50)
  
  # Save hexagon grid
  st_write(hex_grid, "output/siuslaw_hexagons.gpkg", delete_dsn = TRUE, quiet = TRUE)
  
  # Process each hexagon
  message("Processing hexagons...")
  results <- list()
  
  for (i in 1:nrow(hex_grid)) {
    message("Processing hexagon ", i, " of ", nrow(hex_grid))
    hexagon <- hex_grid[i, ]
    results[[i]] <- process_hexagon(hexagon, i)
  }
  
  # Combine results
  all_results <- bind_rows(results)
  
  # Join results to hexagon grid
  hex_grid_with_age <- left_join(hex_grid, all_results, by = c("id" = "hex_id"))
  
  # Save results
  st_write(hex_grid_with_age, "output/siuslaw_tree_age.gpkg", delete_dsn = TRUE, quiet = TRUE)
  
  # Create maps
  message("Creating maps...")
  
  # Age map
  age_map <- tm_shape(hex_grid_with_age) +
    tm_fill("mean_age", 
            title = "Mean Tree Age (years)",
            palette = "viridis",
            style = "cont") +
    tm_borders(col = "white", lwd = 0.1) +
    tm_layout(title = "Estimated Tree Age in Siuslaw National Forest",
              legend.position = c("right", "bottom"))
  
  # Height map
  height_map <- tm_shape(hex_grid_with_age) +
    tm_fill("mean_height", 
            title = "Mean Tree Height (ft)",
            palette = "viridis",
            style = "cont") +
    tm_borders(col = "white", lwd = 0.1) +
    tm_layout(title = "Tree Height in Siuslaw National Forest",
              legend.position = c("right", "bottom"))
  
  # Save maps
  tmap_save(age_map, "output/siuslaw_age_map.png")
  tmap_save(height_map, "output/siuslaw_height_map.png")
  
  message("Analysis complete. Results saved to output/siuslaw_tree_age.gpkg")
  
  return(hex_grid_with_age)
}

# Run the analysis if this script is executed directly
if (!interactive()) {
  run_tree_age_analysis()
} 