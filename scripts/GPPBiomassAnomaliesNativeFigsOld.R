library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(sf)
library(trend)
library(ggplot2)
library(gridExtra)
library(pals)
library(scales)
library(rgeos)
library(reshape2)
library(ncdf4)
library(raster)

#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = 'D:/Driving_C', layer = "WorldContinents")
contnrlist <- c(1,6,3,4)#number in continent shapefiles

NorthAmericaShape <- readOGR(dsn = 'D:/Driving_C', layer = "NorthAmericaNoGreenland")
contsfordisp <- aggregate(continentshapes,dissolve=T)

yearlistGPP <- seq(2003,2018,1)
yearlistC <- seq(2011,2018,1)
yearlistmod <- seq(1901,2018,1)

GPPdf <- data.frame(year=yearlistGPP)
cVegVODcompdf <- data.frame(year=yearlistC)
cVegalldf <- data.frame(year=yearlistmod)
#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclassraster <- raster("D:/Driving_C/Plots/drylandclassclipfin.tif")

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- readOGR(dsn = 'D:/Driving_C', layer = "drylandsglobal")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")
drylandclassSub <- readOGR(dsn = 'D:/Driving_C', layer = "drylands4contsub")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")


#preprocessing of VOD data to median annual composites can be found in LVODprocessing.R and LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("D:/Driving_C/LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#preprocessing of PMLv2 GPP in GEE
GPPstack <- stack("D:/Driving_C/PMLV2sampled/PMLv2GPPstack10knew.tif")

#2011 to 2018 (2010 does not have reliable values, extend to 2019 once TRENDY runs available)
VODstack <- VODannualstacktot[[2:9]]




#RS data bricks in Mg C ha

VODCarbonfinbrick <- VODstack*37.522 #factor based on regression LVOD vs Avitabile

GPPfinbrick <- GPPstack/100

#TRENDY mean data bricks

TRENDYcVegbrick <- brick('D:/Driving_C/DGVM/TRENDYcVeg2011_2018v3.tif')*10
TRENDYGPPbrick <- brick('D:/Driving_C/DGVM/TRENDYGPP_2003_2018v3.tif')*10


#replace this dryland masking by exactextractr dryland shapefile extraction

drylandmask <- drylandclassresamp>0
drylandmask[drylandmask==0] <- NA

#GPP dryland raster mask

drylandclassresamp <- raster::resample(drylandclassraster,GPPstack[[1]],method='ngb')
drylandclassPMLresamp <- drylandclassresamp

drylandmaskPML <- drylandclassPMLresamp>0
drylandmaskPML[drylandmaskPML==0] <- NA


#VOD dryland raster mask
drylandclassag <- aggregate(drylandclassraster,fact=10,fun=modal,na.rm=F)
drylandclassresamp <- raster::resample(drylandclassag,VODstack[[1]],method='ngb')

drylandclassVODresamp <- drylandclassresamp

drylandmaskVOD <- drylandclassVODresamp>0
drylandmaskVOD[drylandmaskVOD==0] <- NA





#contnamelist <- c('Africa','South America','Australia','North America')


modelmatrix =  matrix(NA, nrow = (2018-2003+1), ncol = 16)


  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask [PMLdatamask ==0] <- NA

  arearaster <- raster::mask(raster::mask(area(GPPfinbrick)*100,drylandmaskPML),PMLdatamask)# 100 hectares in 1 km2
  
  totalpercell <- arearaster*GPPfinbrick
  
  totalglobal <- cellStats(totalpercell,sum,na.rm=T)
  
  totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
  
  drylandGPPPML <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
  

  
  #TRENDY mean GPP

  
  arearaster <- area(TRENDYGPPbrick[[1]])*100#*lcfraster
  totalpercell <- arearaster*TRENDYGPPbrick 
  
  totalglobalextract <- extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
  
  totalglobalextract[,2:17] <- totalglobalextract[,2:17]*totalglobalextract$weight
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[2:17]
  
  totalglobalPgC <- (totalglobal/(10^9)) #from Mg to Pg, from cVeg to AGC
  
  TRENDYGPPTS <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
  
  #GPP extracted from native resolution models
  
  GPPPMLcompAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018.csv"))
  
  
  df <- cbind(TRENDY=TRENDYGPPTS$GPP,MODIS=drylandGPPPML$GPP,GPPPMLcompAllModels)
  
  melted <- melt(df,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  GPPplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    #geom_line(data=data.frame(years=yearlistGPP,value=PML,variable='test'),aes(color="Main Forecast"),fill="black",lwd=1.5)+
    #geom_line(data=data.frame(years=yearlistGPP,value=TRENDYmean,variable='test2'),colour="grey",lwd=1.5)+
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values=ltytest)+
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    #plot.margin=margin(1,1,1,1,'cm')
    #guides(colour = guide_legend(override.aes = list(lty = 1, fill = colpalette)))+
    labs(title='a)',y=bquote("GPP " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')
  #geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  plot(GPPplot)

  
  #subtract mean of first year
  dfvar <- df
  dfvar[,-3] <-  df[,-3]-as.list(colMeans(df[,-3]))
  dfvar <- dfvar[,-1]
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.75)}) 
  
  #TRENDYsd <- apply(dfvar[,3:14],1,FUN=sd)
  
#GPP plot minus mean (variance)

  melted <- melt(dfvar,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#000000",colpalette)
  
  lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  
  GPPplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    #geom_line(data=data.frame(years=yearlistGPP,value=PML,variable='test'),aes(color="Main Forecast"),fill="black",lwd=1.5)+
    #geom_line(data=data.frame(years=yearlistGPP,value=TRENDYmean,variable='test2'),colour="grey",lwd=1.5)+
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    #plot.margin=margin(1,1,1,1,'cm')
    #guides(colour = guide_legend(override.aes = list(lty = 1, fill = colpalette)))+
    labs(title='b)',y=bquote("norm. GPP IAV " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')+
   geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  plot(GPPplotvariance)
    
  
  
#ribbon plot variant
  
  dfribbon <- df
  dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3]))
  
  GPPribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbon$TRENDY,MODIS=dfribbon$MODIS,TRENDY_Quartlower,TRENDY_Quartupper)
  GPPplotribbon <- ggplot(GPPribbonplotdf, aes(x = years)) + #vegetation crown plot
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    
    xlab('Time [yr]') +
    ylab(bquote("norm. GPP IAV " ~ "["~ Pg ~ C  ~ yr^{-1}~ "]")) +
    labs(title='c)')+
    #scale_x_continuous(breaks=seq(1,12,1))+
    # ylim(c(0.4,1))+
    geom_line(aes(y = MODIS),lwd=2) +
    geom_line(aes(y = TRENDYmean),col='darkgrey',lwd=2) +
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
    geom_ribbon(aes(ymin =TRENDY_Quartlower,
                    ymax = TRENDY_Quartupper), alpha = 0.2,fill='darkgrey')
  
  print(GPPplotribbon)
  
  GPPplotlist <- list(GPPplot,GPPplotvariance,GPPplotribbon)
  grid.arrange(GPPplot,arrangeGrob(GPPplotvariance,GPPplotribbon),ncol=2)#,layout_matrix = c(1,1,2,3))
  
  #GPP bias and variance statistics
  
  modelBiasGPP <- colMeans(df[2:14]-drylandGPPPML$GPP)
  modelVarianceGPP <- colMeans(abs(t(apply(df[2:14],1,'-',colMeans(df[2:14])))-(drylandGPPPML$GPP-mean(drylandGPPPML$GPP))))
  
  GPPtempstatsresdf <- data.frame(model=names(df[2:14]),bias=modelBiasGPP,variance=modelVarianceGPP)
   
  write.table(GPPtempstatsresdf,"D:/Driving_C/stats/TRENDYmodelsGPPstats_2003_2018old.csv",row.names=F)
  
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  
  
  TStrends <- apply(df,2,fun2)
  
  write.table(TStrends,"D:/Driving_C/stats/TRENDYmodelsGPPtrends_2003_2018.csv",row.names=F)
  
  ##############
  
  Cplotlist <- list()
  #par(mfcol=c(5,2))
  #par(mar=c(4,2.25,2.5,1))
  #par(oma=c(0,0,0,0))
  #matrix of regression/correlation values

  VODdatamask <- !is.na(sum(VODCarbonfinbrick)) #only use pixels with data valid GPP data for all years
  VODdatamask[VODdatamask ==0] <- NA
  
  VODdatamaskpoly <- rasterToPolygons(VODdatamask,na.rm=T,dissolve=T)
  
  #writeOGR(VODdatamaskpoly,'D:/Driving_C/shapefiles','VODdatamask',driver='ESRI Shapefile')
  
  #intersection of VODdatamaskpoly and drylandsubconts done in ArcMap (error in R)
  VODdatamaskdrylands <- readOGR('D:/Driving_C','VODdatamaskdrylands')
  
  
  drylandVODpixels <- raster::mask(drylandclassraster,VODdatamaskdrylands)
  drylandVODpixelsN <- cellStats(!is.na(drylandVODpixels),sum,na.rm=T)
  drylandtotpixelsN <- cellStats(!is.na(raster::mask(drylandclassraster,drylandclassSub)),sum,na.rm=T)
  
  #drylandVODpixelsN/drylandtotpixelsN
  
  arearaster <- area(VODCarbonfinbrick)*100# 100 hectares in 1 km2
  
  totalpercell <- arearaster*VODCarbonfinbrick
  
  totalglobalextract <- extract(totalpercell,VODdatamaskdrylands,weights=T,normalizeWeights=F,df=T)
  
  totalglobalextract[,2:9] <- totalglobalextract[,2:9]*totalglobalextract$weight
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[2:9]
  
  totalglobalPgC <- (totalglobal/(10^9)) #from Mg to Pg, from cVeg to AGC
  
  
  drylandcVegVOD <- data.frame(year=yearlistC,LVOD=totalglobalPgC)
 
  
  #TRENDY mean cVeg

  
  arearaster <- area(TRENDYcVegbrick[[1]])*100#*lcfraster
  totalpercell <- arearaster*TRENDYcVegbrick#*0.4
  
  totalglobalextract <- extract(totalpercell,VODdatamaskdrylands,weights=T,normalizeWeights=F,df=T)
  
  totalglobalextract[,2:9] <- totalglobalextract[,2:9]*totalglobalextract$weight
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[2:9]
  
  totalglobalPgC <- (totalglobal/(10^9))*0.4 #from Mg to Pg, from cVeg to AGC
  
  TRENDYcVegTS <- data.frame(year=yearlistC,cVeg=totalglobalPgC)
  

  


  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  
  
  cVegVODcompAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018.csv"))
  
  
  df <- cbind(TRENDY=TRENDYcVegTS$cVeg,LVOD=drylandcVegVOD$LVOD,cVegVODcompAllModels)
  
  
  melted <- melt(df,id.vars='year')
  
  Cplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +#ggplot(data=melted,aes(x=years)) + 
   #geom_line(data=data.frame(years=yearlistC,value=VOD,variable='test'),colour="black",lwd=1.5)+
    #geom_line(data=data.frame(years=yearlistC,value=TRENDYmean,variable='test'),colour="grey",lwd=1.5)+
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    #plot.margin=margin(1,1,1,1,'cm')
    labs(title='a)',y=bquote("AGC " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
  geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplot)
  
  
  
  #subtract mean of first year
  dfvar <- df
  dfvar[,-3] <-  df[,-3]-as.list(colMeans(df[,-3]))
  dfvar <- dfvar[,-1]
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.75)}) 
  
  #TRENDYsd <- apply(dfvar[,3:14],1,FUN=sd)

  #GPP plot minus mean (variance)
  colpalette <- hue_pal()(12)
  colpalette <- c("#000000",colpalette)
  
  lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  melted <- melt(dfvar,id.vars='year')
  
  Cplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +#ggplot(data=melted,aes(x=years)) + 
    #geom_line(data=data.frame(years=yearlistC,value=VOD,variable='test'),colour="black",lwd=1.5)+
    #geom_line(data=data.frame(years=yearlistC,value=TRENDYmean,variable='test'),colour="grey",lwd=1.5)+
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    #plot.margin=margin(1,1,1,1,'cm')
    labs(title='b)',y=bquote("norm. AGC IAV " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplotvariance)
  
  
  #ribbon plot variant
  
  dfribbon <- df
  dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3]))
  
  
  
  Cribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbon$TRENDY,LVOD=dfribbon$LVOD,TRENDY_Quartlower,TRENDY_Quartupper)
  Cplotribbon <- ggplot(Cribbonplotdf, aes(x = years)) + #vegetation crown plot
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='c)')+
    xlab('Time [yr]') +
    ylab(bquote("norm. AGC IAV " ~ "["~ Pg ~ C ~ "]")) +
    #scale_x_continuous(breaks=seq(1,12,1))+
    # ylim(c(0.4,1))+
    geom_line(aes(y = LVOD),lwd=2) +
    geom_line(aes(y = TRENDYmean),col='darkgrey',lwd=2) +
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
    geom_ribbon(aes(ymin =TRENDY_Quartlower,
                    ymax = TRENDY_Quartupper), alpha = 0.2,fill='darkgrey')
  
  print(Cplotribbon)
  
  
  grid.arrange(Cplot,arrangeGrob(Cplotvariance,Cplotribbon),ncol=2)#,layout_matrix = c(1,1,2,3))
  
  
  
  #AGC bias and variance statistics
  
  modelBiasAGC <- colMeans(df[2:14]- drylandcVegVOD$LVOD)
  modelVarianceAGC <- colMeans(abs(t(apply(df[2:14],1,'-',colMeans(df[2:14])))-(drylandcVegVOD$LVOD-mean(drylandcVegVOD$LVOD))))
  
  AGCtempstatsresdf <- data.frame(model=names(df[2:14]),bias=modelBiasAGC,variance=modelVarianceAGC)
  
  write.table(AGCtempstatsresdf,"D:/Driving_C/stats/TRENDYmodelscVegstats_2011_2018.csv",row.names=F)
  
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  
  
  TStrends <- apply(df,2,fun2)
  
  write.table(TStrends,"D:/Driving_C/stats/TRENDYmodelsAGCtrends_2011_2018.csv",row.names=F)
  
  
  
  
  #cVeg all models native resolution
  
  cVegAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
  
  cVegAllModels[,2:13] <-  cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13])
  
  
  #dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
  melted <- melt(cVegAllModels,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c(colpalette,"#909090")
  
  ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)
  
  lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)
#dfallyears[,12] <- NA
  
  #get mean of first year
  
p1 <- ggplot(data=melted,aes(x=year,y=value,group=variable)) + 
   geom_line(aes(colour=variable,lwd=variable,lty=variable))+
   scale_colour_manual(values = colpalette)+
   scale_size_manual(values = lwdtest)+
   scale_linetype_manual(values = ltytest)+
   ylim(c(-55,55))+
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
   geom_hline(aes(yintercept=0),lty=2)+
    #plot.margin=margin(1,1,1,1,'cm')
    labs(y=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='a)')
  

  

#cSoil from model native resolutions

cSoilAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))

cSoilAllModels[,2:13] <-  cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13])


#dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
melted <- melt(cSoilAllModels,id.vars='year')

colpalette <- hue_pal()(12)
colpalette <- c(colpalette,"#909090")

ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)

lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)
#dfallyears[,12] <- NA

#get mean of first year

#firstyearmean <- mean(unlist(dfallyears[1,2:16]),na.rm=T)
#subtract mean of first year

#melted <- melt(dfallyears,id.vars='years')
p2 <- ggplot(data=melted,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-55,55))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='b)')

plot(test)

#Total Dryland biomass

cVegAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
cSoilAllModels <- data.frame(read.csv("D:/Driving_C/DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))


cTotalAllModels <- cSoilAllModels+cVegAllModels

cTotalAllModels[,2:13] <-  cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13])

cTotalAllModels[,1] <- cSoilAllModels[,1]

#dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
melted <- melt(cTotalAllModels,id.vars='year')

colpalette <- hue_pal()(12)
colpalette <- c(colpalette,"#909090")

ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)

lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)
#dfallyears[,12] <- NA

#melted <- melt(dfallyears,id.vars='years')
p3 <- ggplot(data=melted,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-55,55))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='c)')

grid.arrange(p1,p2,p3,ncol=3)#,layout_matrix = c(1,1,2,3))


#model comparison

cVegdf <- data.frame(cVegAllModels[,2:13])

cSoildf <- data.frame(cSoilAllModels[,2:13])

cTotaldf <- data.frame(cTotalAllModels[,2:13])

deltacVegdf <- cVegdf[118,]-cVegdf[1,]
deltacSoildf <- cSoildf[118,]-cSoildf[1,]
deltaTotaldf <- cTotaldf[118,]-cTotaldf[1,]


deltacdf <- data.frame(deltacVeg=unlist(deltacVegdf),deltacSoil=unlist(deltacSoildf),names=names(deltacVegdf),deltaTotaldf=unlist(deltaTotaldf),modelnr=seq(1,12,1))

test <- ggplot(data=deltacdf,aes(x=deltacVeg,y=deltacSoil,color=deltaTotaldf,label=modelnr))+
  geom_hline(aes(yintercept=0),lty=2,col='grey')+
  geom_vline(aes(xintercept=0),lty=2,col='grey')+
  geom_point()+
  xlim(c(-55,55))+
  ylim(c(-55,55))+
  geom_text(aes(label=modelnr),hjust=0.5, vjust=-1,col='black',size=3)+
 # geom_point(aes(colour=color))+
  scale_color_gradientn(colours = viridis(5))+
  labs(colour=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"))

plot(test)
  

  
nc_close(ncin)
               