#This function trims our watermask rasters in which land is coded as 0 and water is coded as 1.

trim_rasters <- function(rast.name, trimvalue) {
  rast.name[rast.name == trimvalue] <- NA # convert values to be trimmed to NA
  rasttrim<-terra::trim(rast.name, value = NA) #trim NAs
  rasttrim[is.na(rasttrim)] <- trimvalue #gdistance won't work with NAs, 
                                          #replace the with trim values - typically 0's
  return(rasttrim)
}

