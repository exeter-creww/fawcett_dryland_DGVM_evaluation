library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(trend)
library(ggplot2)
library(gridExtra)
library(scales)
library(raster)
library(reshape2)
library(sf)
library(ncdf4)
  
setwd('D:/Driving_C')
#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = 'D:/Driving_C', layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = 'D:/Driving_C', layer = "NorthAmericaNoGreenland")
contsfordisp <- aggregate(continentshapes,dissolve=T)

 
yearlistGPP <- seq(2003,2018,1)
yearlistC <- seq(2011,2018,1)
yearlistmod <- seq(1901,2018,1)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- readOGR(dsn = getwd(), layer = "drylandsglobal")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")
drylandclassSub <- readOGR(dsn = getwd(), layer = "drylands4contsub")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")
drylandclasssf <- st_as_sfc(drylandclass) #spatialpolygonsdf to sfc for exactextractr
drylandclassSubsf <- st_as_sfc(drylandclassSub) #spatialpolygonsdf to sfc for exactextractr

#mask of dryland regions with reliable VOD data, generated from VOD data mask in ArcMap (error with R vectorisation) 
VODdatamaskdrylands <- readOGR(getwd(),'VODdatamaskdrylands')
VODdatamaskdryalndssf <- st_as_sfc(VODdatamaskdrylands) #spatialpolygonsdf to sfc for exactextractr


#preprocessing of PMLv2 GPP in GEE
GPPstack <- stack("D:/Driving_C/PMLV2sampled/PMLv2GPPstack10knew.tif")

#contnrlist <- c(1,6,3,4)#number in continent shapefiles

    
  gppmodelpath <- 'ORCHIDEE-CNP_S1_gpp.nc'
  cVegmodelpath <- 'ORCHIDEE-CNP_S1_cVeg.nc'
  cSoilmodelpath <- 'ORCHIDEE-CNP_S1_cSoil.nc'
  #lcpath <- 'ORCHIDEE-CNP_S3_oceanCoverFrac.nc'
  #4 dimensional netcdf
  ncingpp <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsGPP/origGrids/",gppmodelpath))
  ncincVeg <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/origGrids/",cVegmodelpath))
  ncincSoil <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscSoil/origGrids/",cSoilmodelpath))
  #lcncin <-  nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsLandCover/",lcpath))

  print(ncingpp)
  print(ncincVeg)
  #print(lcncin)
  
  lonDGVM <- ncvar_get(ncingpp, "lon")
  nlonDGVM <- dim(lonDGVM)
  
  latDGVM <- ncvar_get(ncingpp, "lat")
  nlatDGVM <- dim(latDGVM)
  
  #time <- ncvar_get(ncingpp,'time')
  
  modelgpp <- ncvar_get(ncingpp,'gpp',start=c(1,1,(302*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcVeg <- ncvar_get(ncincVeg,'cVeg')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcVegVODcomp <- modelcVeg[,,311:318]
  modelcVeg <- modelcVeg[,,201:318]
  modelcSoil <- ncvar_get(ncincSoil,'cSoil')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcSoil <- modelcSoil[,,201:318]
  fillvalue <- ncatt_get(ncingpp,'gpp',"_FillValue")
  
  modelgpp[modelgpp==fillvalue$value] <- NA
  modelcVegVODcomp[modelcVegVODcomp==fillvalue$value] <- NA
  modelcVeg[modelcVeg==fillvalue$value] <- NA
  modelcSoil[modelcSoil==fillvalue$value] <- NA
  
  modelgppbrick <- t((brick(modelgpp)))#no flip  needed for DLEM
  modelcVegVODcompbrick <-  t((brick(modelcVegVODcomp)))
  modelcVegbrick <-  t((brick(modelcVeg)))
  modelcSoilbrick <-  t((brick(modelcSoil)))
  
  extent(modelgppbrick) <- c(-180, 180, -90, 90)
  projection(modelgppbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegVODcompbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegVODcompbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcSoilbrick) <- c(-180, 180, -90, 90)
  projection(modelcSoilbrick) <- CRS("+init=epsg:4326")
  
  
  nc_close(ncingpp)
  nc_close(ncincVeg)
  nc_close(ncincSoil)
  
  monthyearindex <- rep(1:16,each=12)
  
  modelannualgpp <- stackApply(modelgppbrick,monthyearindex,fun=mean)
  modelannualgpp <- modelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year
  
  modelannualcVegVODcomp <- modelcVegVODcompbrick*10
  modelannualcVeg <- modelcVegbrick
  modelannualcSoil <- modelcSoilbrick
  
  
  GPPstacksresamp <- raster::resample(GPPstack,modelannualgpp[[1]])
  GPPfinbrick <- GPPstacksresamp/100
  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask[PMLdatamask ==0] <- NA
  
   modelannualgppmasked <- mask(modelannualgpp ,PMLdatamask)*10
   modelannualcVegVODcompmasked <- mask(modelannualcVegVODcomp ,VODdatamask)*10
   modelannualcVegmasked <- modelannualcVeg*10
   modelannualcSoilmasked <- modelannualcSoil*10
   
   arearaster <- area(modelannualgppmasked[[1]])*100#*lcfraster
   
   
   #GPP calc
   totalpercell <- arearaster*modelannualgppmasked 
   
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:16] <- totalglobalextract[,1:16]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:16]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ORCHIDEE-CNP_S1_dryland_GPP_2003_2018.csv",sep=",",row.names = F)
   
   #GPP calc cells touching
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:16] <- totalglobalextract[,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:16]
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ORCHIDEE-CNP_S1_dryland_GPP_2003_2018_extended.csv",sep=",",row.names = F)
   
   
   #GPP calc only cells contained
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
   
   totalglobalextract <- totalglobalextract[containedpixels,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ORCHIDEE-CNP_S1_dryland_GPP_2003_2018_contained.csv",sep=",",row.names = F)
   
   
   #cVeg VOD comp calc
   arearaster <- area(modelannualcVegVODcomp[[1]])*100#*lcfraster
   
   totalpercell <- arearaster*modelannualcVegVODcomp
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,VODdatamaskdryalndssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:8] <- totalglobalextract[,1:8]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ORCHIDEE-CNP_S1_dryland_cVeg_2011_2018.csv",sep=",",row.names = F)
   
   #cVeg VOD comp calc all cells touched
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:8] <- totalglobalextract[,1:8]#*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ORCHIDEE-CNP_S1_dryland_cVeg_2011_2018_extended.csv",sep=",",row.names = F)
   
   
   #cVeg VOD comp calc only cells with centre within
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
   
   totalglobalextract <- totalglobalextract[containedpixels,1:8]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ORCHIDEE-CNP_S1_dryland_cVeg_2011_2018_contained.csv",sep=",",row.names = F)
   
   
   
   #cVeg global 1901 calc
   totalpercell <- arearaster*modelannualcVegmasked 
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
   
   totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistmod,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ORCHIDEE-CNP_S1_dryland_cVeg_1901_2018.csv",sep=",",row.names = F)
   
   #cSoil global 1901 calc
   totalpercell <- arearaster*modelannualcSoilmasked 
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
   
   
   totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
   
   dfdrylandcSoil <- data.frame(year=yearlistmod,cSoil=totalglobalPgC)
   
   write.table(dfdrylandcSoil,"D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/ORCHIDEE-CNP_S1_dryland_cSoil_1901_2018.csv",sep=",",row.names = F)