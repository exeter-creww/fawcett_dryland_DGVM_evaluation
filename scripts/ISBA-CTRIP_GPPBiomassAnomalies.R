library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(trend)
library(ggplot2)
library(gridExtra)
library(pals)
library(scales)
library(raster)
library(reshape2)
library(sf)
library(ncdf4)
 

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
    
  gppmodelpath <- 'ISBA-CTRIP_S3_gpp.nc'
  cVegmodelpath <- 'ISBA-CTRIP_S3_cVeg.nc'
  lcpath <- 'ISBA-CTRIP_sftlf.nc'
  cSoilmodelpath <- 'ISBA-CTRIP_S3_cSoil.nc'
  
  #4 dimensional netcdf
  ncingpp <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsGPP/origGrids/",gppmodelpath))
  ncincVeg <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/origGrids/",cVegmodelpath))
  ncincSoil <- nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelscSoil/origGrids/",cSoilmodelpath))
  lcncin <-  nc_open(paste0("D:/Driving_C/DGVM/TRENDYmodelsLandCover/",lcpath))

  print(ncingpp)
  print(ncincVeg)
  print(ncincSoil)
  print(lcncin)
  
  lonDGVM <- ncvar_get(ncingpp, "lon_FULL")
  nlonDGVM <- dim(lonDGVM)
  
  latDGVM <- ncvar_get(ncingpp, "lat_FULL")
  nlatDGVM <- dim(latDGVM)
  
  #time <- ncvar_get(ncingpp,'time')
  
  modelgpp <- ncvar_get(ncingpp,'gpp',start=c(1,1,(303*12)+1),count=c(nlonDGVM,nlatDGVM,192))
  modelcVeg <- ncvar_get(ncincVeg,'cVeg',start=c(1,1,(201*12)+1),count=c(nlonDGVM,nlatDGVM,1416))
  modelcVegVODcomp <- modelcVeg[,,((110*12)+1):length(modelcVeg[1,1,])]
  #modelcVeg <- modelcVeg[,,202:319]
  modelcSoil <- ncvar_get(ncincSoil,'cSoil',start=c(1,1,(201*12)+1),count=c(nlonDGVM,nlatDGVM,1416))
  #modelcSoil <- modelcSoil[,,((201*12)+1):(((201*12)+1)+1415)]
  #modelcSoil <- modelcSoil[,,202:319]
  fillvalue <- ncatt_get(ncingpp,'gpp',"_FillValue")
  
  modelgpp[modelgpp==fillvalue$value] <- NA
  modelcVegVODcomp[modelcVegVODcomp==fillvalue$value] <- NA
  modelcVeg[modelcVeg==fillvalue$value] <- NA
  modelcSoil[modelcSoil==fillvalue$value] <- NA
  
  #get landcover fraction
  landcoverfrac <- ncvar_get(lcncin,'sftlf',start=c(1,1),count=c(nlonDGVM,nlatDGVM))
  
  lcfraster <- t(raster::flip(raster(landcoverfrac),1))
  
  modelgppbrick <- t(raster::flip(brick(modelgpp),1))
  modelcVegVODcompbrick <-  t(raster::flip(brick(modelcVegVODcomp),1))
  modelcVegbrick <-  t(raster::flip(brick(modelcVeg),1))
  modelcSoilbrick <-  t(raster::flip(brick(modelcSoil),1))
  
  
  extent(modelgppbrick) <- c(-180, 180, -60, 90)
  projection(modelgppbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegVODcompbrick) <- c(-180, 180, -60, 90)
  projection(modelcVegVODcompbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcVegbrick) <- c(-180, 180, -60, 90)
  projection(modelcVegbrick) <- CRS("+init=epsg:4326")
  
  extent(modelcSoilbrick) <- c(-180, 180, -60, 90)
  projection(modelcSoilbrick) <- CRS("+init=epsg:4326")
  
  extent(lcfraster) <- c(-180, 180, -60, 90)
  projection(lcfraster) <- CRS("+init=epsg:4326")
  
  nc_close(ncingpp)
  nc_close(ncincVeg)
  nc_close(lcncin)
  nc_close(ncincSoil)
  
  
  monthyearindexgpp <- rep(1:16,each=12)
  monthyearindexcVegVODcomp <- rep(1:8,each=12)
  monthyearindexcVeg <- rep(1:118,each=12)
  monthyearindexcSoil <- rep(1:118,each=12)
  
  modelannualgpp <- stackApply(modelgppbrick,monthyearindexgpp,fun=mean)
  modelannualgpp <- modelannualgpp*31556952 #from mean kg/m2/s to kg/m2/year
  
  modelannualcVegVODcomp  <- stackApply(modelcVegVODcompbrick,monthyearindexcVegVODcomp,fun=mean)*10
  modelannualcVeg <-  stackApply(modelcVegbrick,monthyearindexcVeg,fun=mean)
  modelannualcSoil <- stackApply(modelcSoilbrick,monthyearindexcVeg,fun=mean)#modelcSoilbrick

  #resample dryland classes to 1 deg (modal aggregation then nearest neighbour resampling)

  GPPstacksresamp <- raster::resample(GPPstack,modelannualgpp[[1]])
  GPPfinbrick <- GPPstacksresamp/100
  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask[PMLdatamask ==0] <- NA

   modelannualgppmasked <- mask(modelannualgpp ,PMLdatamask)*10
   modelannualcVegmasked <- modelannualcVeg*10
   modelannualcSoilmasked <- modelannualcSoil*10
   
    arearaster <- area(modelannualgppmasked[[1]])*100*lcfraster
    
    
    #GPP calc
    totalpercell <- arearaster*modelannualgppmasked 
    
    
    totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:16] <- totalglobalextract[,1:16]*totalglobalextract$coverage_fraction
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:16]
    
    
    totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
    
    dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
    
    write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ISBA-CTRIP_dryland_GPP_2003_2018.csv",sep=",",row.names = F)
    
    #GPP calc cells touching
    
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:16] <- totalglobalextract[,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:16]
    
    totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
    
    dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
    
    write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ISBA-CTRIP_dryland_GPP_2003_2018_extended.csv",sep=",",row.names = F)
    
    
    #GPP calc only cells contained
    
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
    
    totalglobalextract <- totalglobalextract[containedpixels,1:16]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)
    
    totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
    
    dfdrylandGPP <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
    
    write.table(dfdrylandGPP,"D:/Driving_C/DGVM/DGVMdrylandTS/GPP/ISBA-CTRIP_dryland_GPP_2003_2018_contained.csv",sep=",",row.names = F)
    
    
    #cVeg VOD comp calc
    arearaster <- area(modelannualcVegVODcomp[[1]])*100*lcfraster
    
    totalpercell <- arearaster*modelannualcVegVODcomp
    
    totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,VODdatamaskdryalndssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:8] <- totalglobalextract[,1:8]*totalglobalextract$coverage_fraction
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
    
    
    totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
    
    dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
    
    write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ISBA-CTRIP_dryland_cVeg_2011_2018.csv",sep=",",row.names = F)
    
    #cVeg VOD comp calc all cells touched
    
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:8] <- totalglobalextract[,1:8]#*totalglobalextract$coverage_fraction
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
    
    totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
    
    dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
    
    write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ISBA-CTRIP_dryland_cVeg_2011_2018_extended.csv",sep=",",row.names = F)
    
    
    #cVeg VOD comp calc only cells with centre within
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    containedpixels <- totalglobalextract$coverage_fraction>=1 #only completely covered pixels
    
    totalglobalextract <- totalglobalextract[containedpixels,1:8]#*totalglobalextract$coverage_fraction #no weighting for touching pixels
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
    
    totalglobalPgC <- totalglobal/(10^9)*0.4 #from Mg to Pg, from cVeg to AGC
    
    dfdrylandcVeg <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
    
    write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ISBA-CTRIP_dryland_cVeg_2011_2018_contained.csv",sep=",",row.names = F)
    
    
    
    #cVeg global 1901 calc
    totalpercell <- arearaster*modelannualcVegmasked 
    
    totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
    
    totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
    
    dfdrylandcVeg <- data.frame(year=yearlistmod,cVeg=totalglobalPgC)
    
    write.table(dfdrylandcVeg,"D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/ISBA-CTRIP_dryland_cVeg_1901_2018.csv",sep=",",row.names = F)
    
    #cSoil global 1901 calc
    totalpercell <- arearaster*modelannualcSoilmasked 
    
    totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    
    totalglobalextract[,1:118] <- totalglobalextract[,1:118]*totalglobalextract$coverage_fraction
    
    totalglobal <- colSums(totalglobalextract,na.rm=T)[1:118]
    
    
    totalglobalPgC <- totalglobal/(10^9)#from Mg to Pg, from cVeg to AGC
    
    dfdrylandcSoil <- data.frame(year=yearlistmod,cSoil=totalglobalPgC)
    
    write.table(dfdrylandcSoil,"D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/ISBA-CTRIP_dryland_cSoil_1901_2018.csv",sep=",",row.names = F)
    
    
    #test month to year aggregation

    # totalpercell <- arearaster*modelcSoilbrick*10
    # totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
    # totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
    # totalglobalextract[,1:1416] <- totalglobalextract[,1:1416]*totalglobalextract$coverage_fraction
    # 
    # totalglobal <- colSums(totalglobalextract,na.rm=T)[1:1416]
    # totalglobalPgC <- totalglobal/(10^9)
    # 
    # cSoilDiff <- modelcSoilbrick[[1400]]-modelcSoilbrick[[1]]
    # 
    # 
    # #functions to create diverging plots in levelplot
    # getMinMax <- function(inraster){
    #   if(abs(cellStats(inraster,min,na.rm=T))>abs(cellStats(inraster,max,na.rm=T))){
    #     colbarbounds <- c(cellStats(inraster,min,na.rm=T),abs(cellStats(inraster,min,na.rm=T)))
    #   }else{
    #     colbarbounds <- c(-cellStats(inraster,max,na.rm=T),cellStats(inraster,max,na.rm=T))
    #   }
    #   return(colbarbounds)
    # }
    # 
    # diverge0 <- function(p, ramp) {
    #   # p: a trellis object resulting from rasterVis::levelplot
    #   # ramp: the name of an RColorBrewer palette (as character), a character 
    #   #       vector of colour names to interpolate, or a colorRampPalette.
    #   require(RColorBrewer)
    #   require(rasterVis)
    #   if(length(ramp)==1 && is.character(ramp) && ramp %in% 
    #      row.names(brewer.pal.info)) {
    #     ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
    #   } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    #     ramp <- colorRampPalette(ramp)
    #   } else if(!is.function(ramp)) 
    #     stop('ramp should be either the name of a RColorBrewer palette, ', 
    #          'a vector of colours to be interpolated, or a colorRampPalette.')
    #   rng <- range(p$legend[[1]]$args$key$at)
    #   s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
    #   i <- findInterval(rng[which.min(abs(rng))], s)
    #   zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
    #   p$legend[[1]]$args$key$at <- s[zlim]
    #   p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
    #   p
    # } 
    # 
    # 
    # my.settings <- list(par.main.text = list(font = 2, just = "left",  x = grid::unit(5, "mm")),panel.background=list(col="lightgrey"))
    # 
    # #VOD trends
    # 
    # trendMinMax <- getMinMax(cSoilDiff)
    # 
    # 
    # diverge0(levelplot(cSoilDiff,par.settings=my.settings,main="ISBA-CTRIP Mg C/ ha 2018 - 1901",at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(drylandclass,col='black'))
    # 
    