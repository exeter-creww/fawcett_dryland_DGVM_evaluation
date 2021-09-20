library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(trend)
library("ncdf4", lib.loc="~/R/win-library/3.4")



#read DGVM GPP data  for 2003 - 2018 and generate raster GPP files

#4 dimensional netcdf
ncin <- nc_open(paste0("D:/Driving_C/DGVM/trendyv8_S3_gpp_1901-2018.nc"))

print(ncin)

modelnames <- ncatt_get(ncin,0,"models")
modelnames <- unlist(strsplit(modelnames$value,' '))

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)

#time <- ncvar_get(ncin,'time')

trendygpp <- ncvar_get(ncin,'gpp',start=c(1,1,(100*12)+1,1),count=c(360,180,216,16))

fillvalue <- ncatt_get(ncin,'gpp',"_FillValue")

trendygpp[trendygpp==fillvalue$value] <- NA
 
#iterate through models
for(i in 1:16){
TRENDYmodelstack <- trendygpp[,,,i]

TRENDYmodelmeanstack <- apply(TRENDYmodelstack,c(1,2,3),FUN=mean,na.rm=T)

TRENDYmodelbrick <- t(raster::flip(brick(TRENDYmodelmeanstack),1))


extent(TRENDYmodelbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelbrick) <- CRS("+init=epsg:4326")

monthyearindex <- rep(1:18,each=12)

TRENDYmodelannualgpp <- stackApply(TRENDYmodelbrick,monthyearindex,fun=mean)
TRENDYmodelannualgpp <- TRENDYmodelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year

writeRaster(TRENDYmodelannualgpp,paste0('D:/Driving_C/DGVM/TRENDYmodelsGPP/',modelnames[i],'_GPP_2001_2018v2.tif'),overwrite=T)
}
nc_close(ncin)