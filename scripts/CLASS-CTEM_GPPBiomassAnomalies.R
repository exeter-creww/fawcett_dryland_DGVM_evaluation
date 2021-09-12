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

#TRENDY mean data bricks

TRENDYcVegbrick <- brick('./DGVM/TRENDYcVeg2011_2018v3.tif')*10 #kg C per m2 to Mg C per ha


#contnrlist <- c(1,6,3,4)#number in continent shapefiles

    
  gppmodelpath <- 'CLASS-CTEM_S3_gpp.nc'
  cVegmodelpath <- 'CLASS-CTEM_S3_cVeg.nc'
  #lcpath <- 'OCN_S3_oceanCoverFrac.nc'
  cSoilmodelpath <- 'CLASS-CTEM_S3_cSoil.nc'
  cRootmodelpath <- 'CLASS-CTEM_S3_cRoot.nc'
  
  #4 dimensional netcdf
  ncingpp <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsGPP/origGrids/",gppmodelpath))
  ncincVeg <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/origGrids/",cVegmodelpath))
  ncincSoil <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscSoil/origGrids/",cSoilmodelpath))
  ncincRoot <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscRoot/origGrids/",cRootmodelpath))
  
  #lcncin <-  nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsLandCover/",lcpath))

  print(ncingpp)
  print(ncincSoil)
 # print(lcncin)
  
  lonDGVM <- ncvar_get(ncingpp, "longitude")
  nlonDGVM <- dim(lonDGVM)
  
  latDGVM <- ncvar_get(ncingpp, "latitude")
  nlatDGVM <- dim(latDGVM)
  
  #time <- ncvar_get(ncingpp,'time')
  
  modelgpp <- ncvar_get(ncingpp,'gpp',start=c(1,1,(303*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelgpptotTS <- ncvar_get(ncingpp,'gpp',start=c(1,1,(201*12)+1),count=c(nlonDGVM,nlatDGVM,1416))
  
  modelcVeg <- ncvar_get(ncincVeg,'cVeg')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcVegVODcomp <- modelcVeg[,,312:319]
  modelcRoot <- ncvar_get(ncincRoot,'cRoot')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcRootVODcomp <- modelcRoot[,,312:319]
  modelcVeg <- modelcVeg[,,202:319]
  modelcSoil <- ncvar_get(ncincSoil,'cSoil')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcSoil <- modelcSoil[,,202:319]
  fillvalue <- ncatt_get(ncingpp,'gpp',"_FillValue")
  
  modelgpp[modelgpp==fillvalue$value] <- NA
  modelgpptotTS[modelgpptotTS==fillvalue$value] <- NA
  modelcVegVODcomp[modelcVegVODcomp==fillvalue$value] <- NA
  modelcRootVODcomp[modelcRootVODcomp==fillvalue$value] <- NA
  
  modelcVeg[modelcVeg==fillvalue$value] <- NA
  modelcSoil[modelcSoil==fillvalue$value] <- NA
  #get landcover fraction
 # landcoverfrac <- abs(ncvar_get(lcncin,'oceanCoverFrac',start=c(1,1),count=c(nlonDGVM,nlatDGVM))-1)
  
  #for DLEM this is area per grid cell not fraction
  #lcfraster <- t(raster::flip(raster(landcoverfrac),1))
  
  modelgppbrick <- t(raster::flip(brick(modelgpp),1))#no flip  needed for DLEM
  modelgpptotTSbrick <- t(raster::flip(brick(modelgpptotTS),1))#no flip  needed for DLEM
  modelcRootVODcompbrick <-  t(raster::flip(brick(modelcRootVODcomp),1))
  
  modelcVegVODcompbrick <-  t(raster::flip(brick(modelcVegVODcomp),1))
  modelcVegbrick <-  t(raster::flip(brick(modelcVeg),1))
  modelcSoilbrick <-  t(raster::flip(brick(modelcSoil),1))
  
  
  extent(modelgppbrick) <- c(-180, 180, -90, 90)
  projection(modelgppbrick) <- CRS("+init=epsg:4326")
  
  extent(modelgpptotTSbrick) <- c(-180, 180, -90, 90)
  projection(modelgpptotTSbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegVODcompbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegVODcompbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcRootVODcompbrick) <- c(-180, 180, -90, 90)
  projection(modelcRootVODcompbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcSoilbrick) <- c(-180, 180, -90, 90)
  projection(modelcSoilbrick) <- CRS("+init=epsg:4326")
  
 # extent(lcfraster) <- c(-180, 180, -90, 90)
 # projection(lcfraster) <- CRS("+init=epsg:4326")
  
  #export cRoot brick resampled to 1 degree
  modelcRootVODcompbrickresamp <- raster::resample(modelcRootVODcompbrick,TRENDYcVegbrick[[1]]) 
  modelcVegVODcompbrickresamp <- raster::resample(modelcVegVODcompbrick,TRENDYcVegbrick[[1]]) 
  
  AGCmap <- (modelcVegVODcompbrickresamp-modelcRootVODcompbrickresamp)
  AGCfracmap <- (modelcVegVODcompbrick-modelcRootVODcompbrick)/modelcVegVODcompbrick
  AGCfracdrylandsmap <- mask(AGCfracmap,drylandclass)
  AGCfracdrylandsmap[AGCfracdrylandsmap<0] <- NA
  
  
  writeRaster(modelcRootVODcompbrickresamp,'./DGVM/TRENDYmodelscRoot/CLASS-CTEM_cRoot2011_2018_1deg.tif',overwrite=T)
  writeRaster(AGCmap,'./DGVM/TRENDYmodelscVeg/CLASS-CTEM_AGC2011_2018_1deg.tif',overwrite=T)
  
  
  
  nc_close(ncingpp)
  nc_close(ncincVeg)
 # nc_close(lcncin)
  
  
  monthyearindex <- rep(1:16,each=12)
  
  modelannualgpp <- stackApply(modelgppbrick,monthyearindex,fun=mean)
  modelannualgpp <- modelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year
  
  monthyearindex <- rep(1:118,each=12)
  
  modelannualgpptotTS <- stackApply(modelgpptotTSbrick,monthyearindex,fun=mean)
  modelannualgpptotTS <- modelannualgpptotTS*31556952 #from mean kg/m2/s to kg/m2/year
  
  modelannualcVegVODcomp <- modelcVegVODcompbrick*10
  modelannualcRootVODcomp <- modelcRootVODcompbrick*10
  
  modelannualcVeg <- modelcVegbrick
  modelannualcSoil <- modelcSoilbrick

  
  #resample dryland classes to 1 deg (modal aggregation then neares neighbour resampling)
  #drylandclassag <- aggregate(drylandclass,fact=10,fun=modal)
  #drylandclassresamp <- resample(drylandclassag,modelannualgpp[[1]],method='ngb')
  #drylandclassresamp <- drylandclassresamp*longmask

  #drylandmask <- drylandclassresamp>0
  #drylandmask[drylandmask==0] <- NA
  
  #regiondrylandmask <- mask(drylandmask,continentshapes[contnrlist,])
  
  GPPstacksresamp <- raster::resample(GPPstack,modelannualgpp[[1]])
  GPPfinbrick <- GPPstacksresamp/100
  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask[PMLdatamask ==0] <- NA
  
  
   modelannualgppmasked <- mask(modelannualgpp ,PMLdatamask)*10
   modelannualgpptotTSmasked <- mask(modelannualgpptotTS ,PMLdatamask)*10
   
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
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/CLASS-CTEM_dryland_GPP_2003_2018.csv",sep=",",row.names = F)
   
   #GPP calc cells touching
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:16] <- totalglobalextract[,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:16]
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/CLASS-CTEM_dryland_GPP_2003_2018_extended.csv",sep=",",row.names = F)
   
   
   #GPP calc only cells contained
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
   
   totalglobalextract <- totalglobalextract[containedpixels,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/CLASS-CTEM_dryland_GPP_2003_2018_contained.csv",sep=",",row.names = F)
   
   
   #cVeg VOD comp calc
   arearaster <- area(modelannualcVegVODcomp[[1]])*100#*lcfraster
   
   totalpercell <- arearaster*modelannualcVegVODcomp
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,VODdatamaskdryalndssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:8] <- totalglobalextract[,1:8]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/CLASS-CTEM_dryland_cVeg_2011_2018.csv",sep=",",row.names = F)
   
   #AGC VOD comp calc
   arearaster <- area(modelannualcVegVODcomp[[1]])*100#*lcfraster
   
   totalpercell <- arearaster*(modelannualcVegVODcomp-modelannualcRootVODcomp)
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,VODdatamaskdryalndssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:8] <- totalglobalextract[,1:8]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/CLASS-CTEM_dryland_AGC_2011_2018.csv",sep=",",row.names = F)
   
   
   #cVeg VOD comp calc all cells touched
   
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:8] <- totalglobalextract[,1:8]#*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/CLASS-CTEM_dryland_AGC_2011_2018_extended.csv",sep=",",row.names = F)
   
   
   #cVeg VOD comp calc only cells with centre within
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
   
   totalglobalextract <- totalglobalextract[containedpixels,1:8]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/CLASS-CTEM_dryland_AGC_2011_2018_contained.csv",sep=",",row.names = F)
   
   #GPP global 1901 calc
   totalpercell <- arearaster*modelannualgpptotTSmasked 
   
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistmod,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/CLASS-CTEM_dryland_GPP_1901_2018.csv",sep=",",row.names = F)
   
   #cVeg global 1901 calc
   totalpercell <- arearaster*modelannualcVegmasked 
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
   
   totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistmod,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/CLASS-CTEM_dryland_cVeg_1901_2018.csv",sep=",",row.names = F)
   
   #cSoil global 1901 calc
   totalpercell <- arearaster*modelannualcSoilmasked 
   
   totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
   
   totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
   
   
   totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
   
   dfdrylandcSoil <- data.frame(year=yearlistmod,cSoil=totalglobalPgC)
   
   write.table(dfdrylandcSoil,"D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/CLASS-CTEM_dryland_cSoil_1901_2018.csv",sep=",",row.names = F)
   