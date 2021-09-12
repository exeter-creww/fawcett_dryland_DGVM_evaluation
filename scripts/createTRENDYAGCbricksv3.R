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
 
modelnames <- ncatt_get(ncin,0,"models")
modelnames <- unlist(strsplit(modelnames$value,' '))

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)
 
cVeg <- ncvar_get(ncin,'cVeg')
fillvalue <- ncatt_get(ncin,'cVeg',"_FillValue")

cVeg[cVeg==fillvalue$value] <- NA

TRENDYstack <- cVeg[,,111:118,c(-10,-11,-15,-16)]*0.4 #0.4 as constant parameter to convert cVeg to AGC

#for models where AGC could be calculated from cVeg-cRoot, use this 
for(i in c(1,2,3,5)){
  for(j in 1:8){
  TRENDYmodelAGCstack <-   stack(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/",modelnames[i],"_AGC2011_2018_1deg.tif"))
#  TRENDYmodelAGCmean <- calc(TRENDYmodelAGCstack,mean,na.rm=T)
  TRENDYmodelAGCstack <- t(raster::flip(TRENDYmodelAGCstack,2))
  TRENDYmodelAGCmat <- as.matrix(TRENDYmodelAGCstack[[j]])
  TRENDYstack[,,j,i] <- TRENDYmodelAGCmat
 # TRENDYmodelsmeanbrick[[i]] <- TRENDYmodelAGCmean 
  }
}

#get mean of all TRENDY models
TRENDYmeanstack <- apply(TRENDYstack,c(1,2,3),FUN=mean,na.rm=T)

TRENDYmeanbrick <- t(raster::flip(brick(TRENDYmeanstack),1))

extent(TRENDYmeanbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmeanbrick) <- CRS("+init=epsg:4326")



# #get mean of all years per model
TRENDYmodelsmeanstack <- apply(TRENDYstack,c(1,2,4),FUN=mean,na.rm=T)

TRENDYmodelsmeanbrick <- t(raster::flip(brick(TRENDYmodelsmeanstack),1))

TRENDYmodelsmeanbrick <- TRENDYmodelsmeanbrick*0.4

extent(TRENDYmodelsmeanbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelsmeanbrick) <- CRS("+init=epsg:4326")


# 
# #get trend slope per model
fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
TRENDYmodelstrendstack <- apply(TRENDYstack,c(1,2,4),FUN=fun1)

TRENDYmodelstrendbrick <- t(raster::flip(brick(TRENDYmodelstrendstack),1))

extent(TRENDYmodelstrendbrick) <- c(-180, 180, -90, 90)
projection(TRENDYmodelstrendbrick) <- CRS("+init=epsg:4326")

nc_close(ncin)


writeRaster(TRENDYmeanbrick,'./DGVM/TRENDYAGC2011_2018v3.tif',overwrite=T)

writeRaster(TRENDYmodelsmeanbrick,'./DGVM/TRENDYpermodelAGC2011_2018v3.tif')

writeRaster(TRENDYmodelstrendbrick,'./DGVM/TRENDYpermodeltrendAGC2011_2018v3.tif',overwrite=T)
