library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(scales)
library(ncdf4)
library(gridExtra)
library(rasterVis)
library(tls)

setwd('D:/Driving_C')

#read DGVM data and calculate veg carbon means per year 
#calculate TRENDY means without SDGVM, LPJ-wsl, LPJ Bern and VISIT

#4 dimensional netcdf
ncin <- nc_open(paste0("./DGVM/trendyv8_S3_cVeg_1901-2018.nc"))

print(ncin)

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)
 
cVeg <- ncvar_get(ncin,'cVeg')
fillvalue <- ncatt_get(ncin,'cVeg',"_FillValue")

cVeg[cVeg==fillvalue$value] <- NA

TRENDYstack <- cVeg[,,111:118,c(-10,-11,-15,-16)]

#get mean of all TRENDY models
TRENDYmeanstack <- apply(TRENDYstack,c(1,2,3),FUN=mean,na.rm=T)

TRENDYmeanbrick <- t(raster::flip(brick(TRENDYmeanstack),1))

extent(TRENDYmeanbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmeanbrick) <- CRS("+init=epsg:4326")

#get mean of all years per model
TRENDYmodelsmeanstack <- apply(TRENDYstack,c(1,2,4),FUN=mean,na.rm=T)

TRENDYmodelsmeanbrick <- t(raster::flip(brick(TRENDYmodelsmeanstack),1))

extent(TRENDYmodelsmeanbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelsmeanbrick) <- CRS("+init=epsg:4326")

#get trend slope per model
fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
TRENDYmodelstrendstack <- apply(TRENDYstack,c(1,2,4),FUN=fun1)

TRENDYmodelstrendbrick <- t(raster::flip(brick(TRENDYmodelstrendstack),1))

extent(TRENDYmodelstrendbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelstrendbrick) <- CRS("+init=epsg:4326")

nc_close(ncin)


writeRaster(TRENDYmeanbrick,'./DGVM/TRENDYcVeg2011_2018v3.tif')

writeRaster(TRENDYmodelsmeanbrick,'./DGVM/TRENDYpermodelcVeg2011_2018v3.tif')

writeRaster(TRENDYmodelstrendbrick,'./DGVM/TRENDYpermodeltrendcVeg2011_2018v3.tif',overwrite=T)
