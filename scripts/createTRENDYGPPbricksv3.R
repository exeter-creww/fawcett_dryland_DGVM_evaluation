
#read DGVM GPP data  for 2003 - 2018 and generate raster GPP files
#calculate TRENDY means without SDGVM, LPJ-wsl, LPJ Bern and VISIT

library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(trend)
library(ncdf4)
library(raster)

setwd('D:/Driving_C')

#4 dimensional netcdf
ncin <- nc_open(paste0("./DGVM/trendyv8_S3_gpp_1901-2018.nc"))

print(ncin)

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)

#time <- ncvar_get(ncin,'time')

trendygpp <- ncvar_get(ncin,'gpp',start=c(1,1,(100*12)+1,1),count=c(360,180,216,16))

fillvalue <- ncatt_get(ncin,'gpp',"_FillValue")

trendygpp[trendygpp==fillvalue$value] <- NA

TRENDYstack <- trendygpp[,,,c(-10,-11,-15,-16)]

TRENDYmeanstack <- apply(TRENDYstack,c(1,2,3),FUN=mean,na.rm=T)

TRENDYbrick <- t(raster::flip(brick(TRENDYmeanstack),1))

extent(TRENDYbrick) <- c(-180, 180, -90, 90)
projection(TRENDYbrick) <- CRS("+init=epsg:4326")

#get mean of all years per model
TRENDYmodelsmeanstack <- apply(TRENDYstack,c(1,2,4),FUN=mean,na.rm=T)

TRENDYmodelsmeanbrick <- t(raster::flip(brick(TRENDYmodelsmeanstack),1))

extent(TRENDYmodelsmeanbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelsmeanbrick) <- CRS("+init=epsg:4326")

nc_close(ncin)

monthyearindex <- rep(1:18, each=12)

TRENDYannualgpp <- stackApply(TRENDYbrick,monthyearindex,fun=mean)
TRENDYannualgpp <- TRENDYannualgpp*31556952 #from mean kg/m2/s to kg/m2/year

TRENDYmodelsannualgpp <- TRENDYmodelsmeanbrick*31556952 


writeRaster(TRENDYannualgpp,'./DGVM/TRENDYGPP_2001_2018v3.tif',overwrite=T)

writeRaster(TRENDYmodelsannualgpp,'./DGVM/TRENDYpermodelGPP2001_2018v3.tif',overwrite=T)

#get annual GPP trend per model

TRENDYmodelstrendstack <- stack()

fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}

for(i in 1:12){
  
  TRENDYmodelstack <- TRENDYstack[,,,i]
  TRENDYmodelbrick <- t(raster::flip(brick(TRENDYmodelstack),1))
  extent(TRENDYmodelbrick) <- c(-180, 180, -90, 90)
  projection(TRENDYmodelbrick) <- CRS("+init=epsg:4326")
  
  TRENDYmodelannualgpp <- stackApply(TRENDYmodelbrick,monthyearindex,fun=mean)
  TRENDYmodelannualgpp <- TRENDYmodelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year

  #TRENDYmodelstrendstack <- apply(TRENDYstack,c(1,2,4),FUN=fun1)
  
  TRENDYmodeltrendgpp <- calc(TRENDYmodelannualgpp,fun1)
  
  TRENDYmodelstrendstack <- stack(TRENDYmodelstrendstack,TRENDYmodeltrendgpp)
  
}

writeRaster(TRENDYmodelstrendstack,'./DGVM/TRENDYpermodeltrendGPP2001_2018v3.tif',overwrite=T)

