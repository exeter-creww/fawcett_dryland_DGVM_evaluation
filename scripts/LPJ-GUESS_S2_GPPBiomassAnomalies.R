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
drylandclass <- readOGR(dsn = 'D:/Driving_C', layer = "drylandsglobal")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")
drylandclassSub <- readOGR(dsn = 'D:/Driving_C', layer = "drylands4contsub")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")


#preprocessing of VOD data to median annual composites can be found in LVODprocessing.R and LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("D:/Driving_C/LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#preprocessing of PMLv2 GPP in GEE
GPPstack <- stack("D:/Driving_C/PMLV2sampled/PMLv2GPPstack10knew.tif")

#2011 to 2018 (2010 does not have reliable values, extend to 2019 once TRENDY runs available)
VODstack <- VODannualstacktot[[2:9]]


#contnrlist <- c(1,6,3,4)#number in continent shapefiles

  gppmodelpath <- 'LPJ-GUESS_S2_gpp.nc'
  cVegmodelpath <- 'LPJ-GUESS_S2_cVeg.nc'
  #lcpath <- 'OCN_S2_oceanCoverFrac.nc'
  cSoilmodelpath <- 'LPJ-GUESS_S2_cSoil.nc'
  #4 dimensional netcdf
  ncingpp <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsGPP/origGrids/",gppmodelpath))
  ncincVeg <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/origGrids/",cVegmodelpath))
  ncincSoil <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscSoil/origGrids/",cSoilmodelpath))
  
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
  modelcVeg <- ncvar_get(ncincVeg,'cVeg')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcVegVODcomp <- modelcVeg[,,312:319]
  modelcVeg <- modelcVeg[,,202:319]
  modelcSoil <- ncvar_get(ncincSoil,'cSoil')#,start=c(1,1,(102*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcSoil <- modelcSoil[,,202:319]
  fillvalue <- ncatt_get(ncingpp,'gpp',"_FillValue")
  
  modelgpp[modelgpp==fillvalue$value] <- NA
  modelcVegVODcomp[modelcVegVODcomp==fillvalue$value] <- NA
  modelcVeg[modelcVeg==fillvalue$value] <- NA
  modelcSoil[modelcSoil==fillvalue$value] <- NA
  #get landcover fraction
 # landcoverfrac <- abs(ncvar_get(lcncin,'oceanCoverFrac',start=c(1,1),count=c(nlonDGVM,nlatDGVM))-1)
  
  #for DLEM this is area per grid cell not fraction
  #lcfraster <- t(raster::flip(raster(landcoverfrac),1))
  
  modelgppbrick <- t(raster::flip(brick(modelgpp),1))#no flip  needed for DLEM
  modelcVegVODcompbrick <-  t(raster::flip(brick(modelcVegVODcomp),1))
  modelcVegbrick <-  t(raster::flip(brick(modelcVeg),1))
  modelcSoilbrick <-  t(raster::flip(brick(modelcSoil),1))
  
  
  extent(modelgppbrick) <- c(-180, 180, -90, 90)
  projection(modelgppbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegVODcompbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegVODcompbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegbrick) <- c(-180, 180, -90, 90)
  projection(modelcVegbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcSoilbrick) <- c(-180, 180, -90, 90)
  projection(modelcSoilbrick) <- CRS("+init=epsg:4326")
  
 # extent(lcfraster) <- c(-180, 180, -90, 90)
 # projection(lcfraster) <- CRS("+init=epsg:4326")
  
  nc_close(ncingpp)
  nc_close(ncincVeg)
 # nc_close(lcncin)
  
  
  monthyearindex <- rep(1:16,each=12)
  
  modelannualgpp <- stackApply(modelgppbrick,monthyearindex,fun=mean)
  modelannualgpp <- modelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year
  
  modelannualcVegVODcomp <- modelcVegVODcompbrick*10
  modelannualcVeg <- modelcVegbrick
  modelannualcSoil <- modelcSoilbrick
  #create longitude based mask and apply to filter out northernmost and southernmost regions
  longmask <- modelannualgpp[[1]]
  longmask[,] <- 1
  longmask[1:68,] <- NA
  longmask[292:360,] <- NA
  
  
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
  
  VODstackresamp <- raster::resample(VODstack,modelannualcVegVODcomp[[1]])
  VODfinbrick <- VODstackresamp*37.522 
  VODdatamask <- !is.na(sum(VODfinbrick))
  VODdatamask[VODdatamask ==0] <- NA
   # TRENDYmodelGPPbrick <- brick(paste0('D:/Driving_C/DGVM/TRENDYmodelsGPP/',TRENDYmodelGPPnames[[j]]))
   modelannualgppmasked <- mask(modelannualgpp ,PMLdatamask)*10
   modelannualcVegVODcompmasked <- mask(modelannualcVegVODcomp ,VODdatamask)*10
   modelannualcVegmasked <- modelannualcVeg*10
   modelannualcSoilmasked <- modelannualcSoil*10
   
   arearaster <- area(modelannualgppmasked[[1]])*100#*lcfraster
   VODdatamaskdrylands <- readOGR('D:/Driving_C','VODdatamaskdrylands')
   
   
   #GPP calc
   totalpercell <- arearaster*modelannualgppmasked 
   
   
   totalglobalextract <- extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
   
   totalglobalextract[,2:17] <- totalglobalextract[,2:17]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:17]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/LPJ-GUESS_S2_dryland_GPP_2003_2018.csv",sep=",",row.names = F)
   
   #GPP calc only cells contained
   totalpercell <- arearaster*modelannualgppmasked 
   
   
   totalglobalextract <- extract(totalpercell,drylandclass,weights=F,small=F,df=T)
   
   #totalglobalextract[,2:17] <- totalglobalextract[,2:17]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:17]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/LPJ-GUESS_S2_dryland_GPP_2003_2018_contained.csv",sep=",",row.names = F)
   
   #GPP calc only cells contained
   totalpercell <- arearaster*modelannualgppmasked 
   
   
   totalglobalextract <- extract(totalpercell,drylandclass,weights=F,small=T,df=T)
   
   #totalglobalextract[,2:17] <- totalglobalextract[,2:17]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:17]
   
   
   totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
   
   dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
   
   write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/LPJ-GUESS_S2_dryland_GPP_2003_2018_extended.csv",sep=",",row.names = F)
   
   #cVeg VOD comp calc
   arearaster <- area(modelannualcVegVODcomp[[1]])*100#*lcfraster
   
   
   totalpercell <- arearaster*modelannualcVegVODcomp
   
   totalglobalextract <- extract(totalpercell,VODdatamaskdrylands,weights=T,normalizeWeights=F,df=T)
   
   totalglobalextract[,2:9] <- totalglobalextract[,2:9]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:9]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/LPJ-GUESS_S2_dryland_cVeg_2011_2018.csv",sep=",",row.names = F)
   
   #cVeg VOD comp calc only cells with centre within
   totalpercell <- arearaster*modelannualcVegVODcomp 
   
   totalglobalextract <- extract(totalpercell,VODdatamaskdrylands,weights=F,small=F,df=T)
   
   #totalglobalextract[,2:9] <- totalglobalextract[,2:9]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:9]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/LPJ-GUESS_S2_dryland_cVeg_2011_2018_contained.csv",sep=",",row.names = F)
   
   #cVeg VOD comp calc all cells touched
   totalpercell <- arearaster*modelannualcVegVODcomp
   
   totalglobalextract <- extract(totalpercell,VODdatamaskdrylands,weights=F,small=T,df=T)
   
   #totalglobalextract[,2:9] <- totalglobalextract[,2:9]*totalglobalextract$weight
   
   totalglobal <- colSums(totalglobalextract,na.rm=T)[2:9]
   
   totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
   
   dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
   
   write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/LPJ-GUESS_S2_dryland_cVeg_2011_2018_extended.csv",sep=",",row.names = F)
   
   
    
    #cVeg global 1901 calc
    totalpercell <- arearaster*modelannualcVegmasked 
    
    totalglobalextract <- extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    
    totalglobalextract[,2:119] <- totalglobalextract[,2:119]*totalglobalextract$weight
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[2:119]
    
    totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
    
    dfdrylandcVeg <- data.frame(year=yearlistmod,cVeg=totalglobalPgC)
    
    write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/LPJ-GUESS_S2_dryland_cVeg_1901_2018.csv",sep=",",row.names = F)
    
    #cSoil global 1901 calc
    totalpercell <- arearaster*modelannualcSoilmasked 
    
    totalglobalextract <- extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    
    totalglobalextract[,2:119] <- totalglobalextract[,2:119]*totalglobalextract$weight
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[2:119]
    
    totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
    
    dfdrylandcSoil <- data.frame(year=yearlistmod,cSoil=totalglobalPgC)
    
    write.table(dfdrylandcSoil,"D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/LPJ-GUESS_S2_dryland_cSoil_1901_2018.csv",sep=",",row.names = F)
    
    
