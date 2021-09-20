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

yearlistGPP <- seq(2001,2018,1)
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
drylandclasssf <- st_as_sfc(drylandclass) #spatialpolygonsdf to sfc for exactextractr
drylandclassSubsf <- st_as_sfc(drylandclassSub) #spatialpolygonsdf to sfc for exactextractr



#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#preprocessing of PMLv2 GPP in GEE
GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew_2001_2018_v016.tif")

#2011 to 2018 (2010 does not have reliable values, extend to 2019 once TRENDY runs available)
VODstack <- VODannualstacktot[[2:9]]

#RS data bricks in Mg C ha

VODCarbonfinbrick <- VODstack*52.48 #calibration with Globbiomass #37.522 #factor based on regression LVOD vs Avitabile

GPPfinbrick <- GPPstack/100

##############
#GPP plots

legendlabels <- c('TRENDY-mean', 'MODIS','CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','ISAM','ISBA-CTRIP','JSBACH','JULES','LPJ-GUESS','OCN','ORCHIDEE','ORCHIDEE-CNP')

#GPP dryland raster mask, to replace with exactextractr!

drylandclassresamp <- raster::resample(drylandclassraster,GPPstack[[1]],method='ngb')
drylandclassPMLresamp <- drylandclassresamp

drylandmaskPML <- drylandclassPMLresamp>0
drylandmaskPML[drylandmaskPML==0] <- NA


modelmatrix =  matrix(NA, nrow = (2018-2001+1), ncol = 18)


  PMLdatamask <- !is.na(sum(GPPfinbrick)) #only use pixels with data valid GPP data for all years
  PMLdatamask [PMLdatamask ==0] <- NA
  
  
  arearaster <- raster::area(PMLdatamask)*100#*lcfraster
  
  
  #GPP calc
  
  totalpercell <- arearaster*GPPfinbrick
  
  totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,drylandclasssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
  totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
  
  totalglobalextract[,1:18] <- totalglobalextract[,1:18]*totalglobalextract$coverage_fraction
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[1:18]
  
  
  totalglobalPgC <- totalglobal/(10^9) #from Mg to Pg
  drylandGPPPML <- data.frame(year=yearlistGPP,GPP=totalglobalPgC)
  
  
  #GPP extracted from native resolution models
  
  GPPPMLcompAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2001_2018.csv"))
  
  newTRENDYmean <- rowMeans(GPPPMLcompAllModels[,2:13])
  
  df <- cbind(TRENDY=newTRENDYmean,MODIS=drylandGPPPML$GPP,GPPPMLcompAllModels)
  
  melted <- melt(df,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  #plot GPP time series 2001-2018 per model 
  GPPplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(labels=legendlabels,values = colpalette)+
    scale_size_manual(labels=legendlabels,values = lwdtest)+
    scale_linetype_manual(labels=legendlabels,values=ltytest)+
    ylim(0, 35)+
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='a)',y=bquote("GPP " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')
  
  plot(GPPplot)

  
  dfvar <- df
  
  #normalize to mean of time series
  dfvar[,-3] <-  dfvar[,-3]-as.list(colMeans(dfvar[,-3]))

#GPP plot minus mean (variance), detrended

  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfvar,2,fun1)
  TStrendpval <- apply(dfvar,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfvar,2,median)-TStrendslope*median(seq(1:18))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[] <- T #decided to detrend all instead
  sigtrendsmask[3] <- F #without year vec
  
  #detrend all time series (previously: with significant trends)
  dfvardetrend <- dfvar
  dfvardetrend[,sigtrendsmask] <- dfvardetrend[,sigtrendsmask]-t(t(matrix(seq(1:18), nrow=18, ncol=15, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  melted <- melt(dfvardetrend,id.vars='year')
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  alphavals <- c(1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  
  GPPplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable,alpha=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    scale_alpha_manual(values=alphavals)+
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='b)',y=bquote("norm. GPP IAV " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='Time [yr]')+
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  plot(GPPplotvariance)
  
  
  
  #GPP ribbon plot variant
  
  
  dfribbon <-  cbind(TRENDY=newTRENDYmean,MODIS=drylandGPPPML$GPP,GPPPMLcompAllModels)#cbind(TRENDY=TRENDYcVegTS$cVeg,LVOD=drylandcVegVOD$LVOD,cVegVODcompAllModels)
  
  #normalize to mean of time series
  dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3]))

  #GPP plot minus mean (variance), detrended
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfribbon,2,fun1)
  TStrendpval <- apply(dfribbon,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfribbon,2,median)-TStrendslope*median(seq(1:18))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[] <- T #decided to detrend all instead
  sigtrendsmask[3] <- F #without year vec
  
  #detrend time series with significant trends
  
  #detrend time series (previously: only with significant trends)
  dfribbondetrend <- dfribbon
  dfribbondetrend[,sigtrendsmask] <- dfribbondetrend[,sigtrendsmask]-t(t(matrix(seq(1:18), nrow=18, ncol=15, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfribbondetrend[,4:15],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfribbondetrend[,4:15],1,FUN=function(x){quantile(x,0.75)}) 
  
  dfribbondetrend$TRENDY <-  apply(dfribbondetrend[,4:15],1,mean)
  
  GPPribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbondetrend$TRENDY,MODIS=dfribbondetrend$MODIS,TRENDY_Quartlower,TRENDY_Quartupper)
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
  
  #old three panel plot
  #GPPplotlist <- list(GPPplot,GPPplotvariance,GPPplotribbon)
  #grid.arrange(GPPplot,arrangeGrob(GPPplotvariance,GPPplotribbon),ncol=2)
  
  #new two panel plot
  blank <- grid.rect(gp=gpar(col="white"))
  
  GPPplotlist <- list(GPPplot,blank,GPPplotvariance)
  
  lay <- rbind(c(1,2),
               c(1,3))
  grid.arrange(GPPplot,GPPplotvariance,blank,ncol=2,layout_matrix=lay)
  
  #GPP bias and variance statistics
  
  dfstats <- df[c(-2,-3)] #remove observed and year columns
  
  modelBiasGPP <- colMeans(dfstats[1:13]-drylandGPPPML$GPP)
  

  #mean of MAE between model and 
  modelVarianceGPP <- colMeans(abs(dfribbondetrend[,c(-2,-3)]-dfribbondetrend$MODIS)) 
  
  GPPtempstatsresdf <- data.frame(model=names(dfstats[1:13]),bias=modelBiasGPP,variance=modelVarianceGPP,modelmeans=colMeans(dfstats[1:13]),PMLmean=mean(drylandGPPPML$GPP))
   
  write.table(GPPtempstatsresdf,"./stats/TRENDYmodelsGPPstats_2001_2018_fin_fixed_PMLV2v016_truemean.csv",sep=',',row.names=F)
  
  #calculate sen's slope and pvalue for each series
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  fun3=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$conf.int[1]) }}
  fun4=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$conf.int[2]) }}
  
  
  TStrendslope <- apply(df,2,fun1)
  TStrendpval <- apply(df,2,fun2)
  TStrendCIlow <- apply(df,2,fun3)
  TStrendCIhigh <- apply(df,2,fun4)
  
  write.table(data.frame(modelname=names(df),slope=TStrendslope,CIlow=TStrendCIlow,CIhigh=TStrendCIhigh,pval=TStrendpval),"./stats/TRENDYmodelsGPPtrends_2001_2018_fin_PMLV2v016_truemean.csv",sep=',',row.names=F)
  
   
  ##############
  #AGC plots
  
  Cplotlist <- list()


  #VODdatamask <- !is.na(sum(VODCarbonfinbrick)) #only use pixels with data valid GPP data for all years
  #VODdatamask[VODdatamask ==0] <- NA
  
  #VODdatamaskpoly <- rasterToPolygons(VODdatamask,na.rm=T,dissolve=T)
  
  #export the data mask as shapefile
  #writeOGR(VODdatamaskpoly,'D:/Driving_C/shapefiles','VODdatamask',driver='ESRI Shapefile')
  
  #intersection of VODdatamaskpoly and drylandsubconts done in ArcMap (error in R)
  VODdatamaskdrylands <- readOGR(getwd(),'VODdatamaskdrylands')
  VODdatamaskdrylandssf <- st_as_sfc(VODdatamaskdrylands)
  
  arearaster <- raster::area(VODCarbonfinbrick)*100# 100 hectares in 1 km2
  
  totalpercell <- arearaster*VODCarbonfinbrick
  
 
  totalglobalextractperpoly <- exactextractr::exact_extract(totalpercell,VODdatamaskdrylandssf,force_df=T)#extract(totalpercell,drylandclass,weights=T,normalizeWeights=F,df=T)
  totalglobalextract <- do.call('rbind',totalglobalextractperpoly)
  
  totalglobalextract[,1:8] <- totalglobalextract[,1:8]*totalglobalextract$coverage_fraction
  
  totalglobal <- colSums(totalglobalextract,na.rm=T)[1:8]
  
  totalglobalPgC <- (totalglobal/(10^9)) #from Mg to Pg, from cVeg to AGC
  
  drylandcVegVOD <- data.frame(year=yearlistC,LVOD=totalglobalPgC)
  
  
  #TRENDY mean AGC

  
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
   
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  
  legendlabels <- c('TRENDY-mean', 'L-VOD','CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','ISAM','ISBA-CTRIP','JSBACH','JULES','LPJ-GUESS','OCN','ORCHIDEE','ORCHIDEE-CNP')
  
  cVegVODcompAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/AGC_drylands_2011_2018.csv"))
  
  newTRENDYmean <- rowMeans(cVegVODcompAllModels[,2:13])
  
  df <- cbind(TRENDY=newTRENDYmean,LVOD=drylandcVegVOD$LVOD,cVegVODcompAllModels)
  
  melted <- melt(df,id.vars='year')
  
  Cplot <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    geom_line(aes(colour=variable,lwd=variable,lty=variable))+
    scale_colour_manual(labels=legendlabels,values = colpalette)+
    scale_size_manual(labels=legendlabels,values = lwdtest)+
    scale_linetype_manual(labels=legendlabels,values = ltytest)+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='a)',y=bquote("AGC " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
  geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplot)
  

  dfvar <-  df
  
  #normalize to mean of time series
  dfvar[,-3] <-  dfvar[,-3]-as.list(colMeans(dfvar[,-3]))
  #dfvar <- dfvar[,-1]
  
  #GPP plot minus mean (variance), detrended
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfvar,2,fun1)
  TStrendpval <- apply(dfvar,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfvar,2,median)-TStrendslope*median(seq(1:8))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[] <- T #decided to detrend all instead
  sigtrendsmask[3] <- F #without year vec
  
  #detrend time series with significant trends

  #detrend time series (previously: only with significant trends)
  dfvardetrend <- dfvar
  dfvardetrend[,sigtrendsmask] <- dfvardetrend[,sigtrendsmask]-t(t(matrix(seq(1:8), nrow=8, ncol=15, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  
  melted <- melt(dfvardetrend,id.vars='year')
  
  #GPP plot minus mean (variance)
  colpalette <- hue_pal()(12)
  colpalette <- c("#909090","#000000",colpalette)
  
  lwdtest <- c(2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  ltytest <- c(1,1,1,2,1,2,1,2,1,2,1,2,1,2)
  alphavals <- c(1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  
  
  Cplotvariance <- ggplot(data=melted,aes(x=year,y=value,group=variable)) +
    
    geom_line(aes(colour=variable,lwd=variable,lty=variable,alpha=variable))+
    scale_colour_manual(values = colpalette)+
    scale_size_manual(values = lwdtest)+
    scale_linetype_manual(values = ltytest)+
    scale_alpha_manual(values=alphavals)+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
    labs(title='b)',y=bquote("norm. AGC IAV " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]')+
    geom_hline(yintercept=0,linetype='dashed',alpha=0.5)
  
  
  #plot VOD and modelled carbon 2011-2018
  plot(Cplotvariance)
  
  
  #ribbon plot variant
   
  dfribbon <- df
  
  
  dfribbon <-  cbind(TRENDY=newTRENDYmean,LVOD=drylandcVegVOD$LVOD,cVegVODcompAllModels)
  
  #normalize to mean of time series
  dfribbon[,-3] <-  dfribbon[,-3]-as.list(colMeans(dfribbon[,-3]))
  
  #GPP plot minus mean (variance), detrended
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  
  
  TStrendslope <- apply(dfribbon,2,fun1)
  TStrendpval <- apply(dfribbon,2,fun2)
  
  #calculation of TS intercept varies, this method using medians is taken from https://pubs.usgs.gov/tm/2006/tm4a7/pdf/USGSTM4A7.pdf
  TStrendintercept <- apply(dfribbon,2,median)-TStrendslope*median(seq(1:8))
  
  sigtrendsmask <- TStrendpval<0.05
  sigtrendsmask[] <- T #decided to detrend all instead
  sigtrendsmask[3] <- F #without year vec
  
  #detrend time series with significant trends
  
  #detrend time series (previously: only with significant trends)
  dfribbondetrend <- dfribbon
  dfribbondetrend[,sigtrendsmask] <- dfribbondetrend[,sigtrendsmask]-t(t(matrix(seq(1:8), nrow=8, ncol=15, byrow=F)[,sigtrendsmask])*TStrendslope[sigtrendsmask]+TStrendintercept[sigtrendsmask])
  
  
  #interquartile range
  TRENDY_Quartlower<- apply(dfribbondetrend[,4:15],1,FUN=function(x){quantile(x,0.25)})
  TRENDY_Quartupper <-  apply(dfribbondetrend[,4:15],1,FUN=function(x){quantile(x,0.75)}) 
  
  #dfribbon$TRENDY <-  apply(dfribbondetrend[,4:15],1,mean)
  
  Cribbonplotdf <- data.frame(years=dfribbon$year,TRENDYmean=dfribbondetrend$TRENDY,LVOD=dfribbondetrend$LVOD,TRENDY_Quartlower,TRENDY_Quartupper)
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
  
  #arrange all AGC time series plots (old 3 panel plot)
  #grid.arrange(Cplot,arrangeGrob(Cplotvariance,Cplotribbon),ncol=2)
  
  #arrange all AGC time series plots
  #grid.arrange(Cplot,Cplotvariance,ncol=1,heights=c(2,1))
  
  #new two panel plot
  blank <- grid.rect(gp=gpar(col="white"))
  
  lay <- rbind(c(1,2),
               c(1,3))
  grid.arrange(Cplot,Cplotvariance,blank,ncol=2,layout_matrix=lay)
  
  #AGC bias and variance statistics
  
  dfstats <- df[c(-2,-3)] #remove observed and year columns
  
  modelBiasAGC <- colMeans(dfstats[1:13]- drylandcVegVOD$LVOD)
  

  #sensitivity calculation using detrended time series
  modelVarianceAGC <- colMeans(abs(dfribbondetrend[,c(-2,-3)]-dfribbondetrend$LVOD)) 
  
  AGCtempstatsresdf <- data.frame(model=names(dfstats[1:13]),bias=modelBiasAGC,variance=modelVarianceAGC,modelmeans=colMeans(dfstats[1:13]),LVODmean=mean(drylandcVegVOD$LVOD))
  
  write.table(AGCtempstatsresdf,"./stats/TRENDYmodelsAGCstats_2011_2018_fin_truemean_test.csv",sep=',',row.names=F)
  
  fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
  fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}
  fun3=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$conf.int[1]) }}
  fun4=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$conf.int[2]) }}
  
  
  TStrendslope <- apply(df,2,fun1)
  TStrendpval <- apply(df,2,fun2)
  TStrendCIlow <- apply(df,2,fun3)
  TStrendCIhigh <- apply(df,2,fun4)
  
  write.table(data.frame(modelname=names(df),slope=TStrendslope,CIlow=TStrendCIlow,CIhigh=TStrendCIhigh,pval=TStrendpval),"./stats/TRENDYmodelsAGCtrends_2011_2018_fin_truemean.csv",sep=',',row.names=F)
  
  #############
  #model time series of cVeg and cSoil since 1901
  
  #cVeg all models native resolution
  
  cVegAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
  
  cVegAllModels[,2:13] <-  cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13])
  
  
  cVegAllModelsMean <- apply(cVegAllModels[,2:13],1,mean,na.rm=T)

  cVegAllModelsinclmean <- data.frame(mean=cVegAllModelsMean,cVegAllModels[,(-7)]) 
  
  melted <- melt(cVegAllModelsinclmean,id.vars='year')
  
   
  legendlabels <- c('TRENDY-mean', 'CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','ISAM','JSBACH','JULES','LPJ-GUESS','OCN','ORCHIDEE','ORCHIDEE-CNP')
  
  
  colpalette <- hue_pal()(12)
  colpalette <- c("#464646",colpalette[c(-6)])#remove colour for ISBA-CTRIP
  
  ltytest <- c(1,1,2,1,2,1,1,2,1,2,1,2)
  
  lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1)

  
p1 <- ggplot(data=melted) + 

   geom_line(aes(x=year,y=value,group=variable,colour=variable,lwd=variable,lty=variable))+
   scale_colour_manual(labels=legendlabels,values = colpalette)+
   scale_size_manual(labels=legendlabels,values = lwdtest)+
   scale_linetype_manual(labels=legendlabels,values = ltytest)+
   ylim(c(-25,25))+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
   geom_hline(aes(yintercept=0),lty=2)+
   
    labs(y=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='a)')
  
plot(p1)  


#cSoil from model native resolutions

cSoilAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))

cLitterAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cLitter/cLitter_drylands_1901_2018.csv"))
cCwdAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cCwd/cCwd_drylands_1901_2018.csv"))

#replace NA with 0 (missing cLitter info)
cLitterAllModels[is.na(cLitterAllModels)] <- 0

cCwdAllModels[is.na(cCwdAllModels)] <- 0

cSoilAllModels[,2:13] <- cSoilAllModels[,2:13]+cLitterAllModels[,2:13]+cCwdAllModels[,2:13] #add cLitter to cSoil

cSoilAllModels[,2:13] <-  cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13])

cSoilAllModelsMean <- apply(cSoilAllModels[,2:13],1,mean,na.rm=T)

cSoilAllModelsinclmean <- data.frame(mean=cSoilAllModelsMean,cSoilAllModels[,(-7)])#remove ISBA-CTRIP 


#dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
melted <- melt(cSoilAllModelsinclmean,id.vars='year')

colpalette <- hue_pal()(12)
colpalette <- c("#464646",colpalette[c(-6)])#remove colour for ISBA-CTRIP

ltytest <- c(1,1,2,1,2,1,1,2,1,2,1,2)

lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1)



p2 <- ggplot(data=melted) + 

  geom_line(aes(x=year,y=value,group=variable,colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(labels=legendlabels,values = colpalette)+
  scale_size_manual(labels=legendlabels,values = lwdtest)+
  scale_linetype_manual(labels=legendlabels,values = ltytest)+
  ylim(c(-25,25))+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  labs(y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='b)')

plot(p2)  

#plot(test)

#Total Dryland biomass

cVegAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
cSoilAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))
cLitterAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cLitter/cLitter_drylands_1901_2018.csv"))
cCwdAllModels <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cCwd/cCwd_drylands_1901_2018.csv"))

#replace NA with 0 (missing cLitter info)
cLitterAllModels[is.na(cLitterAllModels)] <- 0

cCwdAllModels[is.na(cCwdAllModels)] <- 0

cSoilAllModels[,2:13] <- cSoilAllModels[,2:13]+cLitterAllModels[,2:13]+cCwdAllModels[,2:13] #add cLitter to cSoil

cTotalAllModels <- cSoilAllModels+cVegAllModels

cTotalAllModels[,2:13] <-  cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13])

cTotalAllModels[,1] <- cSoilAllModels[,1]


cTotalAllModelsMean <- apply(cTotalAllModels[,2:13],1,mean,na.rm=T)

cTotalAllModelsinclmean <- data.frame(mean=cTotalAllModelsMean,cTotalAllModels[,(-7)]) 


#dfallyears <- cbind(dfallyears,TRENDY=TRENDYmean)
melted <- melt(cTotalAllModelsinclmean,id.vars='year')

colpalette <- hue_pal()(12)
colpalette <- c("#464646",colpalette[c(-6)])#remove colour for ISBA-CTRIP

ltytest <- c(1,1,2,1,2,1,1,2,1,2,1,2)

lwdtest <- c(2,1,1,1,1,1,1,1,1,1,1,1)
#dfallyears[,12] <- NA


p3 <- ggplot(data=melted) + 
  
  geom_line(aes(x=year,y=value,group=variable,colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(labels=legendlabels,values = colpalette)+
  scale_size_manual(labels=legendlabels,values = lwdtest)+
  scale_linetype_manual(labels=legendlabels,values = ltytest)+
  ylim(c(-25,25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  labs(y=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='c)')

plot(p3) 

#legend.position='none',

grid.arrange(p1,p2,p3,ncol=3)#,layout_matrix = c(1,1,2,3))

#ribbon plots of cVeg, cSoil and total since 1901

#95% CI range
# 
# 
# cSoilAllModels[,2:12] <-  cSoilAllModels[,2:12]-as.list(cSoilAllModels[1,2:12])
# cVegAllModels[,2:12] <-  cVegAllModels[,2:12]-as.list(cVegAllModels[1,2:12])
# 
# cTotalAllModelsMean <- apply(cTotalAllModels[,2:12],1,mean,na.rm=T)
# cTotalAllModels_Quartlower<- apply(cTotalAllModels[,2:12],1,FUN=function(x){quantile(x,0.025)})
# cTotalAllModels_Quartupper <-  apply(cTotalAllModels[,2:12],1,FUN=function(x){quantile(x,0.975)}) 
# 
# cSoilAllModelsMean <- apply(cSoilAllModels[,2:12],1,mean,na.rm=T)
# cSoilAllModels_Quartlower<- apply(cSoilAllModels[,2:12],1,FUN=function(x){quantile(x,0.025)})
# cSoilAllModels_Quartupper <-  apply(cSoilAllModels[,2:12],1,FUN=function(x){quantile(x,0.975)}) 
# 
# cVegAllModelsMean <- apply(cVegAllModels[,2:12],1,mean,na.rm=T)
# cVegAllModels_Quartlower<- apply(cVegAllModels[,2:12],1,FUN=function(x){quantile(x,0.025)})
# cVegAllModels_Quartupper <-  apply(cVegAllModels[,2:12],1,FUN=function(x){quantile(x,0.975)}) 
# 
# cVegAllModelsribbonplotdf <- data.frame(years=cVegAllModels$year,TRENDYmean=cVegAllModelsMean,TRENDY_Quartlower=cVegAllModels_Quartlower,TRENDY_Quartupper=cVegAllModels_Quartupper)
# 
# p1 <- ggplot(cVegAllModelsribbonplotdf, aes(x = years)) + 
#   theme_classic() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
#   ylim(c(-25,25))+
#   xlab('Time [yr]') +
#   ylab(bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]")) +
#   labs(title=bquote('d)'))+
#   geom_line(aes(y = TRENDYmean),col='black',lwd=1) +
#   geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
#   geom_ribbon(aes(ymin =TRENDY_Quartlower,
#                   ymax = TRENDY_Quartupper), alpha = 0.5,fill='darkgrey')
# 
# 
# cSoilAllModelsribbonplotdf <- data.frame(years=cSoilAllModels$year,TRENDYmean=cSoilAllModelsMean,TRENDY_Quartlower=cSoilAllModels_Quartlower,TRENDY_Quartupper=cSoilAllModels_Quartupper)
# 
# p2 <- ggplot(cSoilAllModelsribbonplotdf, aes(x = years)) + 
#   theme_classic() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
#   ylim(c(-25,25))+
#   xlab('Time [yr]') +
#   ylab(bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]")) +
#   labs(title=bquote('e)'))+
#   geom_line(aes(y = TRENDYmean),col='black',lwd=1) +
#   geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
#   geom_ribbon(aes(ymin =TRENDY_Quartlower,
#                   ymax = TRENDY_Quartupper), alpha = 0.5,fill='darkgrey')
# 
# 
# cTotalAllModelsribbonplotdf <- data.frame(years=cTotalAllModels$year,TRENDYmean=cTotalAllModelsMean,TRENDY_Quartlower=cTotalAllModels_Quartlower,TRENDY_Quartupper=cTotalAllModels_Quartupper)
# p3 <- ggplot(cTotalAllModelsribbonplotdf, aes(x = years)) + 
#   theme_classic() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
#   ylim(c(-25,25))+
#   xlab('Time [yr]') +
#   ylab(bquote(Delta ~" cTotal " ~ "["~ Pg ~ C ~ "]")) +
#   labs(title=bquote('f)'))+
#   geom_line(aes(y = TRENDYmean),col='black',lwd=1) +
#   geom_hline(yintercept=0,linetype='dashed',alpha=0.5)+
#   geom_ribbon(aes(ymin =TRENDY_Quartlower,
#                   ymax = TRENDY_Quartupper), alpha = 0.5,fill='darkgrey')
# 
# 
# grid.arrange(p1,p2,p3,ncol=3)#,layout_matrix = c(1,1,2,3))


#model comparison

cVegdf <- data.frame(cVegAllModels[,c(-1,-7)])

cSoildf <- data.frame(cSoilAllModels[,c(-1,-7)])

cTotaldf <- data.frame(cTotalAllModels[,c(-1,-7)])

deltacVegdf <- cVegdf[118,]-cVegdf[1,]
deltacSoildf <- cSoildf[118,]-cSoildf[1,]
deltaTotaldf <- cTotaldf[118,]-cTotaldf[1,]


colPal <- colorRampPalette(c("red", "grey", "blue"))

legendlabels <- c('CABLE-POP','CLASS-CTEM','CLM5.0','DLEM','ISAM','JSBACH','JULES','LPJ-GUESS','OCN','ORCHIDEE','ORCHIDEE-CNP')

legend_entry <- paste("  ", seq(1,11,1), legendlabels)


deltacdf <- data.frame(deltacVeg=unlist(deltacVegdf),deltacSoil=unlist(deltacSoildf),legend_entry,names=names(deltacVegdf),deltaTotaldf=unlist(deltaTotaldf),modelnr=seq(1,11,1))

ymin <- min(deltacdf$deltacSoil)
ymax <- 25#max(deltacdf$deltacSoil)

y_values_legend <- ymax-(ymax-ymin)*(1:nrow(deltacdf))/nrow(deltacdf)

p1 <- ggplot(data=deltacdf,aes(x=deltacVeg,y=deltacSoil,label=modelnr))+
  geom_hline(aes(yintercept=0),lty=2,col='grey')+
  geom_vline(aes(xintercept=0),lty=2,col='grey')+
  geom_point(aes(colour = deltaTotaldf))+
  theme_bw()+
  xlim(c(-25,25))+
  ylim(c(-25,25))+
  scale_colour_gradient2(low='red',mid='grey',high='blue')+
  geom_text(aes(label=modelnr),hjust=0.5, vjust=-1,col='black',size=3)+
  geom_text(aes(label=names, x=Inf, y=y_values_legend, hjust=0)) +
  coord_fixed()+

  labs(colour=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"))

plot(p1)

# add numbered legend, final plot is WIP, currently created externally


p2 <- ggplot(data=deltacdf,aes(x=deltacVeg,y=deltacSoil,label=modelnr))+
  geom_hline(aes(yintercept=0),lty=2,col='grey')+
  geom_vline(aes(xintercept=0),lty=2,col='grey')+
  geom_point(aes(colour = deltaTotaldf))+
  theme_bw()+
  xlim(c(-25,25))+
  ylim(c(-25,25))+
  scale_colour_gradient2(low='red',mid='grey',high='blue')+
  geom_text(aes(label=modelnr),hjust=0.5, vjust=-1,col='black',size=3)+
  geom_text(aes(label=legend_entry, x=Inf, y=y_values_legend, hjust=0)) +
  theme(legend.position='false',plot.margin = unit(c(1,15,1,1), "lines"))+
  labs(colour=bquote(Delta ~" cEco " ~ "["~ Pg ~ C ~ "]"),y=bquote(Delta ~" cSoil " ~ "["~ Pg ~ C ~ "]"),x=bquote(Delta ~" cVeg " ~ "["~ Pg ~ C ~ "]"))

gt <- ggplot_gtable(ggplot_build(p2))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)



#TRENDY mean stats comparison

#quantify divergence of models from the TRENDY mean (1901-2018)


cVegTRENDYmean <- apply(cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13]),1,mean,na.rm=T)
cVegTRENDYSD <- apply(cVegAllModels[,2:13]-as.list(cVegAllModels[1,2:13]),1,sd,na.rm=T)

cVegTRENDYSDend <- cVegTRENDYSD[118]
cVegTRENDYmeanend <- cVegTRENDYmean[118]

cSoilTRENDYmean <- apply(cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13]),1,mean,na.rm=T)
cSoilTRENDYSD <- apply(cSoilAllModels[,2:13]-as.list(cSoilAllModels[1,2:13]),1,sd,na.rm=T)

cSoilTRENDYSDend <- cSoilTRENDYSD[118]
cSoilTRENDYmeanend <- cSoilTRENDYmean[118]


cTotalTRENDYmean <- apply(cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13]),1,mean,na.rm=T)
cTotalTRENDYSD <- apply(cTotalAllModels[,2:13]-as.list(cTotalAllModels[1,2:13]),1,sd,na.rm=T)

cTotalTRENDYSDend <- cTotalTRENDYSD[118]
cTotalTRENDYmeanend <- cTotalTRENDYmean[118]

cVegdifftoTRENDY <- (deltacVegdf-cVegTRENDYmeanend)/cVegTRENDYSDend
cSoildifftoTRENDY <- (deltacSoildf-cSoilTRENDYmeanend)/cSoilTRENDYSDend
totaldifftoTRENDY <- (deltaTotaldf-cTotalTRENDYmeanend)/cTotalTRENDYSDend

cVegModelsmean <- apply(cVegAllModels[,c(-1,-7)],2,mean,na.rm=T)
cSoilModelsmean <- apply(cSoilAllModels[,c(-1,-7)],2,mean,na.rm=T)
cTotalModelsmean <- apply(cTotalAllModels[,c(-1,-7)],2,mean,na.rm=T)

write.table(data.frame(modelname=names(cVegdifftoTRENDY),cVegmean=unlist(cVegModelsmean),cSoilmean=unlist(cSoilModelsmean),cTotalmean=unlist(cTotalModelsmean),cVegchange=unlist(deltacVegdf),cVegdeviation=unlist(cVegdifftoTRENDY),cSoilchange=unlist(deltacSoildf),cSoildeviation=unlist(cSoildifftoTRENDY),cTotaldeviation=unlist(deltaTotaldf),cTotaldeviation=unlist(totaldifftoTRENDY)),"./stats/TRENDYmodelsdeviation_1901_2018_finV2.csv",row.names=F,sep=',')

nc_close(ncin)
               