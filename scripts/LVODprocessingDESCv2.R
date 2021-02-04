#convert daily vOD observations in nc to yearly mean, median and valid observation count rasters

#library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(fields)
library(viridis)

setwd('D:/Driving_C')

#read ncdf file

yearfnames <- list.files("./LVOD_WGS84/DESC/")


yeardaystack <- stack()
for(j in 1:11){

  
ncfnamesDESC <- list.files(paste0("./LVOD_WGS84/DESC/",yearfnames[j],"/"))


vodstack1 <- stack() #first filtering for mean calculation
vodstack2 <- stack() #SD filtering applied
voddaystack <- stack()
for(i in 1:length(ncfnamesDESC)){

  
ncin <- nc_open(paste0("./LVOD_WGS84/DESC/",yearfnames[j],"/",ncfnamesDESC[i]))

#get DOY
#obsdate <- as.Date(substr(ncfnamesDESC[i],20,27),'%Y%m%d')
#obsDOY <- as.numeric(strftime(obsdate, format = "%j"))

vod.array <- ncvar_get(ncin, "Optical_Thickness_Nad1")
sf.array <- ncvar_get(ncin, "Scene_Flags3")
rmse.array <- ncvar_get(ncin, "RMSE4")

sfmask <- sf.array<=1
sfmask[sfmask==0] <- NA
rmsemask <- rmse.array<=8 #rmse between observed and modelled TB [Kelvin]
rmsemask[rmsemask==0] <- NA

vodfiltered <- vod.array*sfmask*rmsemask

#vodobsday <- is.finite(vodfiltered)*obsDOY 
#vodobsday[vodobsday==0] <- NA

#voddaystack <- stack(voddaystack,raster(vodobsday))
vodstack1 <- stack(vodstack1,raster(vodfiltered))

lon <- ncvar_get(ncin, "lon")
nlon <- dim(lon)

lat <- ncvar_get(ncin, "lat")
nlat <- dim(lat)

nc_close(ncin)
}





 vodmean <- calc(vodstack1,mean,na.rm=T)
 vodsd <- calc(vodstack1,sd,na.rm=T)
 
# vodmedian <- calc(vodstack1,median,na.rm=T)
# vodmean <- t(flip(vodmean,1))
# vodmean[vodmean<0] <- NA
# extent(vodmean) <- c(min(lon), max(lon), min(lat), max(lat))
# projection(vodmean) <- CRS("+init=epsg:4326")
# writeRaster(vodmean,paste0("D:/Driving_C/LVOD_WGS84/composites/vodmean",yearfnames[j],"_DESC.tif"),overwrite=T)
#
#
 
 for(i in 1:length(ncfnamesDESC)){
   
   
   #get DOY
   obsdate <- as.Date(substr(ncfnamesDESC[i],20,27),'%Y%m%d')
   obsDOY <- as.numeric(strftime(obsdate, format = "%j"))
   
   
   vodfiltered <- vodstack1[[i]]#vod.array*sfmask*rmsemask
   
   voddiff <- abs(vodfiltered-vodmean)
   
   vodSDfiltered <- vodfiltered
   vodSDfiltered[voddiff>(2*vodsd)] <- NA
   
   vodobsday <- is.finite(vodSDfiltered)*obsDOY 
   vodobsday[vodobsday==0] <- NA
   
   voddaystack <- stack(voddaystack,vodobsday)
   vodstack2 <- stack(vodstack2,vodSDfiltered)
   
 }
 
vodmedian <- calc(vodstack2,median,na.rm=T)
vodmedian <- t(flip(vodmedian,1))
vodmedian[vodmedian<0] <- 0#NA
extent(vodmedian) <- c(-180, 180, -90, 90)
projection(vodmedian) <- CRS("+init=epsg:4326")
writeRaster(vodmedian,paste0("./LVOD_WGS84/composites/vodmedian",yearfnames[j],"_DESCv2.tif"),overwrite=T)

# vodmean <- calc(vodstack2,mean,na.rm=T)
# vodmean <- t(flip(vodmean,1))
# vodmean[vodmean<0] <- NA
# extent(vodmean) <- c(min(lon), max(lon), min(lat), max(lat))
# projection(vodmean) <- CRS("+init=epsg:4326")
# writeRaster(vodmean,paste0("D:/Driving_C/LVOD_WGS84/composites/vodmean",yearfnames[j],"_DESC.tif"),overwrite=T)

# #

#function to calculate circular mean for doy
# circdaymean <- function(d, ...){
#   if(length(d)==0){
#     return(NA)
#   }else{
# int <- 365
#  rad.d <- d*(360/int)*(pi/180)
#  sinr <- sum(sin(rad.d),na.rm=T)
#  cosr <- sum(cos(rad.d),na.rm=T)
#  circmean <- atan2(sinr, cosr)
#  x.deg <- circmean*(180/pi)+360
#  daymean <- (x.deg/(360/int))%%365
#  return(daymean)
#   }
# }
#
#
# vodmeanday <- calc(voddaystack,fun=function(x,na.rm=T){circdaymean(x,na.rm=T)},na.rm=T)
# vodmeanday <- t(flip(vodmeanday,1))
# vodmeanday[vodmeanday<0] <- NA
# extent(vodmeanday) <- c(min(lon), max(lon), min(lat), max(lat))
# projection(vodmeanday) <- CRS("+init=epsg:4326")
# yeardaystack <- stack(yeardaystack,vodmeanday)
#writeRaster(vodmedianday,paste0("D:/Driving_C/LVOD_WGS84/composites/vodmedianday",yearfnames[j],"_DESC.tif"),overwrite=T)

#
# vodcount <- calc(vodstack2,fun=function(x,na.rm=T){sum(!is.na(x))},na.rm=T)
# vodcount <- t(flip(vodcount,1))
# vodcount[vodcount<0] <- NA
# extent(vodcount) <- c(min(lon), max(lon), min(lat), max(lat))
# projection(vodcount) <- CRS("+init=epsg:4326")
# writeRaster(vodcount,paste0("D:/Driving_C/LVOD_WGS84/composites/vodcount",yearfnames[j],"_DESC.tif"),overwrite=T)


#image.plot(lon, lat, as.matrix(vodmean), col = viridis(100))

}
# yeardaystacksub <- yeardaystack[[2:10]]
# yeardaymeans <- calc(yeardaystacksub,fun=function(x,na.rm=T){circdaymean(x,na.rm=T)},na.rm=T)
#writeRaster(yeardaymeans,paste0("D:/Driving_C/LVOD_WGS84/composites/vodmeanday2010_2018_DESC.tif"),overwrite=T)
# rasterVis::levelplot(yeardaymeans,main="Circular mean DOY of VOD observation (2011-2018)",at=seq(1, 365, len = 100),col.regions=rainbow(100),margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey")))
