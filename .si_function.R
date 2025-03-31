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
  
 
  prod_subset <- arcpullr::get_image_layer(url = 'https://di-usfsdata.img.arcgis.com/arcgis/rest/services/CONUS_site_productivity_2018_masked_202106032103033/ImageServer',
                                           sf_object = aoi)
  # Plot the productity class for the cell. 
  
  plot(aoi$geometry, main = "Productivity Classes")
  plot(prod_subset)
  plot(sf::st_geometry(ttops), add = TRUE, pch = 16, cex = 0.5, col = "red")
  
  ttops$prod <- terra::extract(prod_subset, ttops)[, 2]
  ttops$SI <- si_lookup[as.character(ttops$prod)] #TODO as number 
  return (ttops$SI)
}

#' this is the function to calculate the tree age from tree height 
#' as of now it works quite well for coastal forests in Oregon. 
#' refinement will be needed to see how well it works for forests in WA and other stands in OR when the forest is less doug fir 
#' this is based on the PAPER CITATION 
#' @author: "Gavin Adams" 
#' 

calculate_tree_age <- function(SI, height, bo = 123.25, b1 = 0.71698, b2 = -0.0001677, b3 =  0.95516, b4 = 0.0022182){
   age <-(log((1-height/(bo+b1*SI))^(b3+b4*SI)))/(b2*SI)
   return(age)
}