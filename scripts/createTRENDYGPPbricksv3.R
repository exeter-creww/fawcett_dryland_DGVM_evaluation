
#read DGVM GPP data  for 2003 - 2018 and generate raster GPP files
#calculate TRENDY means without SDGVM, LPJ-wsl, LPJ Bern and VISIT

library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(trend)
library(ncdf4)

setwd('D:/Driving_C')

#4 dimensional netcdf
ncin <- nc_open(paste0("./DGVM/trendyv8_S3_gpp_1901-2018.nc"))

print(ncin)

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)

#time <- ncvar_get(ncin,'time')

trendygpp <- ncvar_get(ncin,'gpp',start=c(1,1,(102*12)+1,1),count=c(360,180,192,16))

fillvalue <- ncatt_get(ncin,'gpp',"_FillValue")

trendygpp[trendygpp==fillvalue$value] <- NA

TRENDYstack <- trendygpp[,,,c(-10,-11,-15,-16)]

TRENDYmeanstack <- apply(TRENDYstack,c(1,2,3),FUN=mean,na.rm=T)

TRENDYbrick <- t(raster::flip(brick(TRENDYmeanstack),1))




extent(TRENDYbrick) <- c(-180, 180, -90, 90)
projection(TRENDYbrick) <- CRS("+init=epsg:4326")

nc_close(ncin)



monthyearindex <- rep(1:16,each=12)
 
TRENDYannualgpp <- stackApply(TRENDYbrick,monthyearindex,fun=mean)
TRENDYannualgpp <- TRENDYannualgpp*31556952 #from mean kg/m2/s to kg/m2/year


writeRaster(TRENDYannualgpp,'./DGVM/TRENDYGPP_2003_2018v3.tif',overwrite=T)
