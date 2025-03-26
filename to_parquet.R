#' Convert GIS file to GeoParquet format
#'
#' This function reads various GIS file formats and converts them to GeoParquet format.
#' It supports all formats readable by the sf package (Shapefiles, GeoJSON, GPKG, etc.).
#'
#' @param input_path Path to the input GIS file
#' @param output_path Path where the GeoParquet file should be saved (including .parquet extension)
#' @param compression Compression algorithm to use (default: "snappy")
#' @param overwrite Logical, whether to overwrite existing file (default: FALSE)
#' @return Invisible NULL, writes file to disk
#' @export
#' @examples
#' \dontrun{
#' # Convert a shapefile to GeoParquet
#' gis_to_geoparquet("data/shapes.shp", "output/shapes.parquet")
#'
#' # Convert a GeoJSON file with ZSTD compression
#' gis_to_geoparquet("data/features.geojson", "output/features.parquet", 
#'                   compression = "zstd")
#' }
gis_to_geoparquet <- function(input_path, output_path, 
                              compression = "snappy", overwrite = FALSE) {
  
  # Check if required packages are installed
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required but not installed. Please install with install.packages('sf')")
  }
  
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required but not installed. Please install with install.packages('arrow')")
  }
  
  # Check if output file exists and handle overwrite
  if (file.exists(output_path) && !overwrite) {
    stop("Output file already exists and overwrite = FALSE")
  }
  
  # Read the input file using sf
  sf_data <- sf::st_read(input_path, quiet = TRUE)
  
  # Validate that we have geometry data
  if (!inherits(sf_data, "sf")) {
    stop("Input file does not contain spatial data or could not be read as such")
  }
  
  # Convert to GeoParquet
  tryCatch({
    arrow::write_parquet(
      sf_data,
      sink = output_path,
      compression = compression
    )
    
    message(sprintf("Successfully converted %s to GeoParquet at %s", 
                    basename(input_path), output_path))
  }, error = function(e) {
    stop(sprintf("Failed to write GeoParquet: %s", e$message))
  })
  
  invisible(NULL)
}

