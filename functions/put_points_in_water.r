#This function checks whether points are inside or outside of our water mask and then changes the values of apparent non-water cells that are intersected to water by changing a value of 0 to 1
#
#
#input_raster<-list.of.trimmed.raster.objects$CATS_water_raster_10m  
#input_points <- CATS_sites.sf.utm  # sf with centroids for buffer
#field <- "Pop"
 
put_points_in_water <- function(input_raster, input_points, field) {
  
#First do a check as to whether points are in (1) or out (0) of the water layer by extracting the raster value that the points intersect.  Water is 1, land is 0.  
  points.in.and.out.of.water <- terra::extract(
    input_raster,  # raster layer
    input_points,  # sf with centroids for buffer
    method = "simple"
  )   
  
#If you have stations that are outside of the water mask, you need to correct the raster to convert the cells with a sampling point to be coded as water.
  
#Create a raster of the locations of the sampling points.
#First create a blank raster based on the water mask raster
  
blank_raster <- input_raster
values(blank_raster) <- 0
  
#rasterize the water stations based on the blank raster, counting the stations within each cell.
  
stationrast <- terra::rasterize(input_points, blank_raster, field = field, fun=sum)
  
  #convert the NA values to 0
  stationrast[is.na(stationrast)] <- 0
  
#Add the watermask raster with the raster that contains the number of stations per cell.
  # 
  corrected_trimmed_raster <- input_raster + stationrast
  
  #Change cells with values of 2 and 3 (i.e., there was >= 1 station in a cell that was already water) to 1 (indicating water), 
  corrected_trimmed_raster[corrected_trimmed_raster >1] <- 1
  
  return(corrected_trimmed_raster)
}