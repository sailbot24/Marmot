#!/usr/bin/env Rscript

# Run Tree Age Analysis
# This script runs the tree age analysis for the Siuslaw National Forest

# Source the tree age calculator
source("oregon_tree_age_calculator.R")

# Run the analysis
message("Starting tree age analysis for Siuslaw National Forest...")
results <- run_tree_age_analysis()
message("Analysis complete!")

# Print summary statistics
message("\nSummary Statistics:")
message("Number of hexagons: ", nrow(results))
message("Mean tree age: ", round(mean(results$mean_age, na.rm = TRUE), 1), " years")
message("Mean tree height: ", round(mean(results$mean_height, na.rm = TRUE), 1), " feet")
message("Total tree count: ", sum(results$tree_count, na.rm = TRUE))

message("\nResults saved to output/siuslaw_tree_age.gpkg")
message("Maps saved to output/siuslaw_age_map.png and output/siuslaw_height_map.png") 