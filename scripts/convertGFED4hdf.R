library(rgdal)
library(gdalUtils)
library(raster)
library(ncdf4)

#not working currently

setwd('D:/Driving_C/GFED4/monthly/')

files <- dir(pattern = ".hdf")

filename <- substr(files,12,17)
filename <- paste0("GFED4BA_", filename, ".tif")


#test <- raster(files[10])

#test <- nc_open(files[1])

#test <- gdal_translate(files[1], dst_dataset='test.tif')

#test <- get_subdatasets(files[10])

for (i in 1:259){
  sds <- get_subdatasets(files[i])
  gdal_translate(sds[1], dst_dataset = filename[i])
}


#convert to annual