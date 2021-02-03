library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(scales)
library("ncdf4", lib.loc="~/R/win-library/3.4")
library(gridExtra)
library(rasterVis)
library(tls)
#clean imports


#read DGVM data and calculate veg carbon means per year 
#calculate TRENDY means without SDGVM, LPJ-wsl, LPJ Bern and VISIT

#4 dimensional netcdf
ncin <- nc_open(paste0("D:/Driving_C/DGVM/trendyv8_S3_cVeg_1901-2018.nc"))

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

nc_close(ncin)

writeRaster(TRENDYmeanbrick,'D:/Driving_C/DGVM/TRENDYcVeg2011_2018v3.tif')
