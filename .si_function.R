library(oRLiDAR)
library(sf)
library(terra)
### 
#' This is a function to calculate SI from USFS productivity layer for the Siuslaw National Forest. 
#' @author Gavin Adams 
#' #TODO other function like things go here
#' 
#' @return the SI column for the ttops form Forest Tools
#' 

si_calculator <- function(aoi, ttops,...){
  si_lookup <- c("0"= NA, "1" = 145, "2" = 125, "3" = 105, "4" = 85, "5" = 65, "6" = 45, "7" = 25)# this is -20 from the last one 
  
  plot(aoi, main = "Productivity Classes")
  plot(sf::st_geometry(ttops), add = TRUE, pch = 16, cex = 0.5, col = "red")
  
  prod_subset <- arcpullr::get_image_layer(url = 'https://di-usfsdata.img.arcgis.com/arcgis/rest/services/CONUS_site_productivity_2018_masked_202106032103033/ImageServer',
                                           sf_object = aoi)
  
 
  ttops$prod <- terra::extract(prod_subset, ttops)[, 2]
  ttops$SI <- si_lookup[as.character(ttops$prod)] #TODO as number 
  return (ttops$SI)
}

calculate_tree_age <- function(SI, height, bo = 123.25, b1 = 0.71698, b2 = -0.0001677, b3 =  0.95516, b4 = 0.0022182){
  (log((1-height/(bo+b1*SI))^(b3+b4*SI)))/(b2*SI)
}