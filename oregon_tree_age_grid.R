#!/usr/bin/env Rscript

# Oregon Tree Age Calculator using Square Grids
# This script calculates tree age from LiDAR-derived canopy height models
# using the Means & Sabin (1989) formula. The area is tessellated into
# 50-acre square grids, and age is calculated for each grid cell.

# Load required libraries
library(sf)
library(dplyr)
library(terra)
library(ForestTools)

# Set options to avoid geometry issues
options(sf_use_s2 = FALSE)

# Function to calculate Site Index (SI) from tree height and age
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

# Function to estimate age from height using the Means & Sabin formula
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
  
  return(best_age)
}

# Create output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Load the National Forest boundaries
message("Loading National Forest boundaries...")
nf_path <- "Veg_ShapeFiles/FS_National_Forests_Dataset_(US_Forest_Service_Proclaimed_Forests)/FS_National_Forests_Dataset_(US_Forest_Service_Proclaimed_Forests).shp"
national_forests <- st_read(nf_path, quiet = TRUE)

# Filter for Siuslaw National Forest
message("Filtering for Siuslaw National Forest...")
siuslaw <- national_forests %>% 
  filter(FORESTNAME == "Siuslaw National Forest")

# Make sure the geometry is valid
siuslaw <- st_make_valid(siuslaw)

# Save the Siuslaw boundary for reference
st_write(siuslaw, "output/siuslaw_boundary.gpkg", delete_dsn = TRUE, quiet = TRUE)

# Calculate the side length for a 50-acre square (in meters)
# 50 acres = 202,343 square meters
# side length = sqrt(202,343) = 450 meters approximately
side_length <- sqrt(50 * 4046.86)

# Create a grid with a simplified approach
message("Creating square grid...")

# Get the bounding box of Siuslaw
bbox <- st_bbox(siuslaw)

# Create a simple rectangular polygon from the bounding box
bbox_poly <- st_as_sfc(bbox)

# Create a square grid over the bounding box
message("Creating grid over bounding box...")
grid <- st_make_grid(bbox_poly, cellsize = side_length, square = TRUE)

# Convert to sf object
grid_sf <- st_sf(id = 1:length(grid), geometry = grid)

# Identify which grid cells intersect with Siuslaw
message("Identifying grid cells that intersect with Siuslaw...")
intersects <- st_intersects(grid_sf, siuslaw)
grid_in_siuslaw <- grid_sf[lengths(intersects) > 0, ]

message("Created ", nrow(grid_in_siuslaw), " grid cells that intersect with Siuslaw")

# Take a sample of grid cells for processing
if(nrow(grid_in_siuslaw) > 5) {
  message("Taking a sample of 5 grid cells for processing...")
  sample_indices <- sample(1:nrow(grid_in_siuslaw), 5)
  grid_sample <- grid_in_siuslaw[sample_indices, ]
} else {
  grid_sample <- grid_in_siuslaw
}

# Save the grid sample for reference
st_write(grid_sample, "output/siuslaw_grid_sample.gpkg", delete_dsn = TRUE, quiet = TRUE)

# Create a visualization of the grid
message("Creating grid visualization...")
png("output/siuslaw_grid_visualization.png", width = 1200, height = 800)
plot(siuslaw$geometry, main = "Square Grid Cells (50 acres) in Siuslaw National Forest", 
     border = "darkgreen", lwd = 2)
plot(grid_sample$geometry, add = TRUE, border = "red", lwd = 1)
dev.off()

# Source the Oregon CHM downloader
source("oregon_chm_downloader.R")

# Process each grid cell
message("Processing grid cells...")
results <- data.frame()

for (i in 1:nrow(grid_sample)) {
  message("Processing grid cell ", i, " of ", nrow(grid_sample))
  
  # Convert sf to sp for compatibility with the downloader
  grid_sp <- as(grid_sample[i,], "Spatial")
  
  # Try to download CHM data
  tryCatch({
    # Download CHM data
    chm_file <- paste0("output/grid_", i, "_chm.tif")
    chm <- download_oregon_chm(grid_sp, output_path = chm_file, mask = TRUE, quiet = FALSE)
    
    # Check if CHM has valid data
    if (is.null(chm) || terra::global(chm, "max", na.rm = TRUE)[1] <= 0) {
      message("No valid CHM data for grid cell ", i)
      next
    }
    
    # Function for defining dynamic window size
    lin <- function(x){x * 0.05 + 0.6}
    
    # Detect treetops (minimum height 10 feet)
    ttops <- vwf(chm, winFun = lin, minHeight = 10)
    
    # If no trees detected, skip
    if (nrow(ttops) == 0) {
      message("No trees detected in grid cell ", i)
      next
    }
    
    # Calculate statistics
    mean_height <- mean(ttops$height, na.rm = TRUE)
    
    # Get dominant height (95th percentile)
    dominant_height <- quantile(ttops$height, 0.95, na.rm = TRUE)
    
    # Estimate ages
    mean_age <- estimate_age_from_height(mean_height)
    dominant_age <- estimate_age_from_height(dominant_height)
    
    # Add to results
    grid_result <- data.frame(
      grid_id = grid_sample$id[i],
      mean_height = mean_height,
      mean_age = mean_age,
      tree_count = nrow(ttops),
      dominant_height = dominant_height,
      dominant_age = dominant_age
    )
    
    results <- rbind(results, grid_result)
    
    message("Successfully processed grid cell ", i)
    print(grid_result)
    
  }, error = function(e) {
    message("Error processing grid cell ", i, ": ", e$message)
  })
}

# Save results
if (nrow(results) > 0) {
  # Join results to grid
  grid_with_age <- merge(grid_sample, results, by.x = "id", by.y = "grid_id", all.x = TRUE)
  
  # Save to file
  st_write(grid_with_age, "output/siuslaw_tree_age.gpkg", delete_dsn = TRUE, quiet = TRUE)
  
  message("Results saved to output/siuslaw_tree_age.gpkg")
  print(results)
  
  # Create a visualization of the results
  message("Creating age map visualization...")
  png("output/siuslaw_age_map.png", width = 1200, height = 800)
  plot(siuslaw$geometry, main = "Tree Age in Siuslaw National Forest (50-acre Square Grids)", 
       border = "darkgreen", lwd = 2)
  
  # Create a color palette
  if (nrow(results) > 1) {
    # Plot with color scale if we have multiple results
    age_range <- range(results$mean_age, na.rm = TRUE)
    age_breaks <- seq(age_range[1], age_range[2], length.out = 10)
    age_colors <- heat.colors(9)
    
    # Plot each grid cell with color based on age
    for (i in 1:nrow(results)) {
      grid_id <- results$grid_id[i]
      age <- results$mean_age[i]
      color_index <- findInterval(age, age_breaks)
      if (color_index < 1) color_index <- 1
      if (color_index > 9) color_index <- 9
      
      plot(grid_sample[grid_sample$id == grid_id, ]$geometry, add = TRUE, 
           col = age_colors[color_index], border = "white")
    }
    
    # Add a legend
    legend("topright", 
           legend = paste(round(age_breaks[-10]), "-", round(age_breaks[-1])), 
           fill = age_colors, 
           title = "Tree Age (years)",
           cex = 0.8)
  } else {
    # Simple plot if we only have one result
    plot(grid_with_age$geometry, add = TRUE, col = "red", border = "white")
    text(st_coordinates(st_centroid(grid_with_age$geometry))[,1], 
         st_coordinates(st_centroid(grid_with_age$geometry))[,2], 
         labels = paste("Age:", round(grid_with_age$mean_age)), 
         cex = 0.8)
  }
  
  dev.off()
  
  message("Visualization saved to output/siuslaw_age_map.png")
} else {
  message("No grid cells were successfully processed.")
} 