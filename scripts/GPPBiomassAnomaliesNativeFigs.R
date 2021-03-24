#Temporal analysis of GPP and AGC over 1) drylands in four continents with LVOD data coverage or 2) globally (for GPP and model time series)
#creates time series plots per model and observation dataset
#includes comparison of cVeg and cSoil from model time series since 1901


library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(sf)
library(trend)
library(ggplot2)
library(gridExtra)
library(scales)
library(rgeos)
library(reshape2)
library(ncdf4)
library(raster)
library(grid)

setwd('D:/Driving_C')

#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
contnrlist <- c(1,6,3,4)#number in continent shapefiles

NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")
contsfordisp <- aggregate(continentshapes,dissolve=T)

yearlistGPP <- seq(2003,2018,1)
yearlistC <- seq(2011,2018,1)
yearlistmod <- seq(1901,2018,1)

GPPdf <- data.frame(year=yearlistGPP)
cVegVODcompdf <- data.frame(year=yearlistC)
cVegalldf <- data.frame(year=yearlistmod)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclassraster <- raster("./Plots/drylandclassclipfin.tif")

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- readOGR(dsn = getwd(), layer = "drylandsglobal")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")
drylandclassSub <- readOGR(dsn = getwd(), layer = "drylands4contsub")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")


#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#preprocessing of PMLv2 GPP in GEE
GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew.tif")

#2011 to 2018 (2010 does not have reliable values, extend to 2019 once TRENDY runs available)
VODstack <- VODannualstacktot[[2:9]]

#RS data bricks in Mg C ha

VODCarbonfinbrick <- VODstack*37.522 #factor based on regression LVOD vs Avitabile

GPPfinbrick <- GPPstack/100

#TRENDY mean data bricks

TRENDYcVegbrick <- brick('./DGVM/TRENDYcVeg2011_2018v3.tif')*10 #kg C per m2 to Mg C per ha
TRENDYGPPbrick <- brick('./DGVM/TRENDYGPP_2003_2018v3.tif')*10 #kg C per m2 to Mg C per ha

##############
#GPP plots

#GPP dryland raster mask, to replace with exactextractr!

drylandclassresamp <- resample(drylandclassraster,GPPstack[[1]],method='ngb')
drylandclassPMLresamp <- drylandclassresamp

drylandmaskPML <- drylandclassPMLresamp>0
drylandmaskPML[drylandmaskPML==0] <- NA


modelmatrix =  matrix(NA, nrow = (2018-2003+1), ncol = 16)


  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask [PMLdatamask ==0] <- NA

  arearaster <- mask(mask(area(GPPfinbrick)*100,drylandmaskPML),PMLdatamask)# 100 hectares in 1 km2
  
  totalpercell <- arearaster*GPPfinbrick
  
  totalglobal <- cellStats(totalpercell,sum,na.rm=T)
  
  totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
  
  drylandGPPPML <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
  

  
  #TRENDY mean GPP

  
  arearaster <- area(TRENDYGPPbrick[[1]])*100#*lcfraster
  totalpercell <- arearaster*TRENDYGPPbrick 
  
  totalglobalextract <- extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)#replace with exactectractr as done for spatial
  
  totalglobalextract[,2:17] <- totalglobalextract[,2:17]*totalglobalextract$weight
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[2:17]
  
  totalglobalPgC <- (totalglobal/(10^9)) #from Mg to Pg, from cVeg to AGC
  
  TRENDYGPPTS <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
  
  #GPP extracted from native resolution models
  
  GPPPMLcompAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018.csv"))
  
  
  df <- cbind(TRENDY=TRENDYGPPTS$GPP,MODIS=drylandGPPPML$GPP,GPPPMLcompAllModels)
  
  melted <- melt(df,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  #plot GPP time series 2003-2018 per model 
  GPPplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values=ltytest)+
    ylim(0, 35)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='a)',y=bquote("GPP " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')

  
  plot(GPPplot)

  
  #subtract mean of first year
  dfvar <- df
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.75)}) 
  

  #normalize to mean of time series
  dfvar[,-3] <-  df[,-3]-as.list(colMeans(df[,-3]))
  dfvar <- dfvar[,-1]
  
#GPP plot minus mean (variance), detrended

  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfvar,2,fun1)
  TStrendpval <- apply(dfvar,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfvar,2,median)-TStrendslope*median(seq(1:16))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[2] <- F #without year vec
  
  #detrend time series with significant trends
  dfvardetrend <- dfvar
  dfvardetrend[,sigtrendsmask] <- dfvardetrend[,sigtrendsmask]-t(t(matrix(seq(1:16), nrow=16, ncol=14, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  
  melted <- melt(dfvardetrend,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#000000",colpalette)
  
  lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  
  GPPplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='b)',y=bquote("norm. GPP IAV " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')+
   geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  plot(GPPplotvariance)
    
  
  
#GPP ribbon plot variant
  
  dfribbon <- df
 # dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3])) #decided against sho
  
  GPPribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbon$TRENDY,MODIS=dfribbon$MODIS,TRENDY_Quartlower,TRENDY_Quartupper)
  GPPplotribbon <- ggplot(GPPribbonplotdf, aes(x = years)) + #vegetation crown plot
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    
    xlab('Time [yr]') +
    ylab(bquote("norm. GPP IAV " ~ "["~ Pg ~ C  ~ yr^{-1}~ "]")) +
    labs(title='c)')+
    geom_line(aes(y = MODIS),lwd=2) +
    geom_line(aes(y = TRENDYmean),col='darkgrey',lwd=2) +
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
    geom_ribbon(aes(ymin =TRENDY_Quartlower,
                    ymax = TRENDY_Quartupper), alpha = 0.2,fill='darkgrey')
  
  print(GPPplotribbon)
  
  GPPplotlist <- list(GPPplot,GPPplotvariance,GPPplotribbon)
  grid.arrange(GPPplot,arrangeGrob(GPPplotvariance,GPPplotribbon),ncol=2)
  
  #GPP bias and variance statistics
  
  dfstats <- df[c(-2,-3)] #remove observed and year columns
  
  modelBiasGPP <- colMeans(dfstats[1:13]-drylandGPPPML$GPP)
  
  modelVarianceGPP <- colMeans(abs(t(apply(dfstats[1:13],1,'-',colMeans(dfstats[1:13])))-(drylandGPPPML$GPP-mean(drylandGPPPML$GPP))))
  
  GPPtempstatsresdf <- data.frame(model=names(dfstats[1:13]),bias=modelBiasGPP,variance=modelVarianceGPP)
   
  write.table(GPPtempstatsresdf,"./stats/TRENDYmodelsGPPstats_2003_2018.csv",row.names=F)
  
  #calculate sen's slope and pvalue for each series
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(df,2,fun1)
  TStrendpval <- apply(df,2,fun2)
  
  write.table(data.frame(modelname=names(df),slope=TStrendslope,pval=TStrendpval),"./stats/TRENDYmodelsGPPtrends_2003_2018.csv",row.names=F)
  
  ##############
  #AGC plots
  
  Cplotlist <- list()


  VODdatamask <- !is.na(sum(VODCarbonfinbrick)) #only use pixels with data valid GPP data for all years
  VODdatamask[VODdatamask ==0] <- NA
  
  VODdatamaskpoly <- rasterToPolygons(VODdatamask,na.rm=T,dissolve=T)
  
  #writeOGR(VODdatamaskpoly,'D:/Driving_C/shapefiles','VODdatamask',driver='ESRI Shapefile')
  
  #intersection of VODdatamaskpoly and drylandsubconts done in ArcMap (error in R)
  VODdatamaskdrylands <- readOGR(getwd(),'VODdatamaskdrylands')
  
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
  
  
  cVegVODcompAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018.csv"))
  
  
  df <- cbind(TRENDY=TRENDYcVegTS$cVeg,LVOD=drylandcVegVOD$LVOD,cVegVODcompAllModels)
  
  
  melted <- melt(df,id.vars='year')
  
  Cplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='a)',y=bquote("AGC " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
  geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplot)
  
  #subtract mean of first year
  dfvar <- df
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfvar[,3:14],1,FUN=function(x){quantile(x,0.75)}) 
  
  
  #normalize to mean of time series
  dfvar[,-3] <-  df[,-3]-as.list(colMeans(df[,-3]))
  dfvar <- dfvar[,-1]
  
  #GPP plot minus mean (variance), detrended
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfvar,2,fun1)
  TStrendpval <- apply(dfvar,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfvar,2,median)-TStrendslope*median(seq(1:8))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[2] <- F #without year vec
  
  #detrend time series with significant trends
  dfvardetrend <- dfvar
  dfvardetrend[,sigtrendsmask] <- dfvardetrend[,sigtrendsmask]-t(t(matrix(seq(1:8), nrow=8, ncol=14, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  
  melted <- melt(dfvardetrend,id.vars='year')
  
  #GPP plot minus mean (variance)
  colpalette <- hue_pal()(12)
  colpalette <- c("#000000",colpalette)
  
  lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  Cplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='b)',y=bquote("norm. AGC IAV " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplotvariance)
  
  
  #ribbon plot variant
  
  dfribbon <- df
 # dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3]))
  
  
  
  Cribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbon$TRENDY,LVOD=dfribbon$LVOD,TRENDY_Quartlower,TRENDY_Quartupper)
  Cplotribbon <- ggplot(Cribbonplotdf, aes(x = years)) + #vegetation crown plot
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='c)')+
    xlab('Time [yr]') +
    ylab(bquote("norm. AGC IAV " ~ "["~ Pg ~ C ~ "]")) +
    geom_line(aes(y = LVOD),lwd=2) +
    geom_line(aes(y = TRENDYmean),col='darkgrey',lwd=2) +
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
    geom_ribbon(aes(ymin =TRENDY_Quartlower,
                    ymax = TRENDY_Quartupper), alpha = 0.2,fill='darkgrey')
  
  print(Cplotribbon)
  
  #arrange all AGC time series plots
  grid.arrange(Cplot,arrangeGrob(Cplotvariance,Cplotribbon),ncol=2)
  
  
  
  #AGC bias and variance statistics
  
  dfstats <- df[c(-2,-3)] #remove observed and year columns
  
  modelBiasAGC <- colMeans(dfstats[1:13]- drylandcVegVOD$LVOD)
  modelVarianceAGC <- colMeans(abs(t(apply(dfstats[1:13],1,'-',colMeans(dfstats[1:13])))-(drylandcVegVOD$LVOD-mean(drylandcVegVOD$LVOD))))
  
  AGCtempstatsresdf <- data.frame(model=names(dfstats[1:13]),bias=modelBiasAGC,variance=modelVarianceAGC)
  
  write.table(AGCtempstatsresdf,"./stats/TRENDYmodelscVegstats_2011_2018.csv",row.names=F)
  
  
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(df,2,fun1)
  TStrendpval <- apply(df,2,fun2)
  
  write.table(data.frame(modelname=names(df),slope=TStrendslope,pval=TStrendpval),"./stats/TRENDYmodelsAGCtrends_2011_2018.csv",row.names=F)
  
  #############
  #model time series of cVeg and cSoil since 1901
  
  #cVeg all models native resolution
  
  cVegAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
  
  cVegAllModels[,2:13] <-  cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13])
  
  melted <- melt(cVegAllModels,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c(colpalette,"#909090")
  
  ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)
  
  lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)

  
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

cSoilAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))

cSoilAllModels[,2:13] <-  cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13])


#dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
melted <- melt(cSoilAllModels,id.vars='year')

colpalette <- hue_pal()(12)
colpalette <- c(colpalette,"#909090")

ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)

lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)

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

#plot(test)

#Total Dryland biomass

cVegAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
cSoilAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))


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
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
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


colPal <- colorRampPalette(c("red", "grey", "blue"))

legend_entry <- paste("  ", seq(1,12,1), names(deltacVegdf))


deltacdf <- data.frame(deltacVeg=unlist(deltacVegdf),deltacSoil=unlist(deltacSoildf),legend_entry,names=names(deltacVegdf),deltaTotaldf=unlist(deltaTotaldf),modelnr=seq(1,12,1))

ymin <- min(deltacdf$deltacSoil)
ymax <- 55#max(deltacdf$deltacSoil)

y_values_legend <- ymax-(ymax-ymin)*(1:nrow(deltacdf))/nrow(deltacdf)

p1 <- ggplot(data=deltacdf,aes(x=deltacVeg,y=deltacSoil,label=modelnr))+
  geom_hline(aes(yintercept=0),lty=2,col='grey')+
  geom_vline(aes(xintercept=0),lty=2,col='grey')+
  geom_point(aes(colour = deltaTotaldf))+
  theme_bw()+
  xlim(c(-55,55))+
  ylim(c(-55,55))+
  scale_colour_gradient2(low='red',mid='grey',high='blue')+
  #(colours = terrain.colors(10))+
  geom_text(aes(label=modelnr),hjust=0.5, vjust=-1,col='black',size=3)+
  geom_text(aes(label=names, x=Inf, y=y_values_legend, hjust=0)) +
  #geom_point(aes(colour=color))+
  #scale_color_gradientn(colours = viridis(5))+
  labs(colour=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"))

plot(p1)

# add numbered legend, final plot is WIP, currently created externally

p2 <- ggplot(data=deltacdf,aes(x=deltacVeg,y=deltacSoil,label=modelnr))+
  geom_hline(aes(yintercept=0),lty=2,col='grey')+
  geom_vline(aes(xintercept=0),lty=2,col='grey')+
  geom_point(aes(colour = deltaTotaldf))+
  theme_bw()+
  xlim(c(-55,55))+
  ylim(c(-55,55))+
  scale_colour_gradient2(low='red',mid='grey',high='blue')+
  #(colours = terrain.colors(10))+
  geom_text(aes(label=modelnr),hjust=0.5, vjust=-1,col='black',size=3)+
  geom_text(aes(label=legend_entry, x=Inf, y=y_values_legend, hjust=0)) +
  theme(legend.position='false',plot.margin = unit(c(1,15,1,1), "lines"))+
  #geom_point(aes(colour=color))+
  #scale_color_gradientn(colours = viridis(5))+
  labs(colour=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"))


gt <- ggplot_gtable(ggplot_build(p2))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)



#TRENDY mean stats comparison

#quantify divergence of models from the TRENDY mean (1901-2018)


cVegTRENDYmean <- apply(cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13]),1,mean)
cVegTRENDYSD <- apply(cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13]),1,sd)

cVegTRENDYSDend <- cVegTRENDYSD[118]
cVegTRENDYmeanend <- cVegTRENDYmean[118]

cSoilTRENDYmean <- apply(cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13]),1,mean)
cSoilTRENDYSD <- apply(cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13]),1,sd)

cSoilTRENDYSDend <- cSoilTRENDYSD[118]
cSoilTRENDYmeanend <- cSoilTRENDYmean[118]


cTotalTRENDYmean <- apply(cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13]),1,mean)
cTotalTRENDYSD <- apply(cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13]),1,sd)

cTotalTRENDYSDend <- cTotalTRENDYSD[118]
cTotalTRENDYmeanend <- cTotalTRENDYmean[118]

cVegdifftoTRENDY <- (deltacVegdf-cVegTRENDYmeanend)/cVegTRENDYSDend
cSoildifftoTRENDY <- (deltacSoildf-cSoilTRENDYmeanend)/cSoilTRENDYSDend
totaldifftoTRENDY <- (deltaTotaldf-cTotalTRENDYmeanend)/cTotalTRENDYSDend

write.table(data.frame(modelname=names(cVegdifftoTRENDY),cVegchange=unlist(deltacVegdf),cVegdeviation=unlist(cVegdifftoTRENDY),cSoilchange=unlist(deltacSoildf),cSoildeviation=unlist(cSoildifftoTRENDY),unlist(deltaTotaldf),cTotaldeviation=unlist(totaldifftoTRENDY)),"./stats/TRENDYmodelsdeviation_1901_2018.csv",row.names=F,sep=',')

nc_close(ncin)
               