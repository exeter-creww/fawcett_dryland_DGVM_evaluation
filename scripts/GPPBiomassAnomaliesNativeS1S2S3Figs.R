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


  #############
  #model time series of cVeg and cSoil since 1901
  
  #cVeg all models native resolution
  
  cVegAllModelsS3 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_1901_2018.csv"))
  cVegAllModelsS2 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_S2_drylands_1901_2018.csv"))
  cVegAllModelsS1 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_S1_drylands_1901_2018.csv"))
  

  cVegAllModelsS3[,2:13] <-  cVegAllModelsS3[,2:13]-as.list(cVegAllModelsS3[1,2:13])
  cVegAllModelsS2[,2:13] <-  cVegAllModelsS2[,2:13]-as.list(cVegAllModelsS2[1,2:13])
  cVegAllModelsS1[,2:13] <-  cVegAllModelsS1[,2:13]-as.list(cVegAllModelsS1[1,2:13])
  
  meltedS3 <- melt(cVegAllModelsS3,id.vars='year')
  
  meltedS2 <- melt(cVegAllModelsS2,id.vars='year')
  
  meltedS1 <- melt(cVegAllModelsS1,id.vars='year')
  
  meltedS2minusS1<- meltedS2
  meltedS2minusS1$value <- meltedS2$value-meltedS1$value
  
  meltedS3minusS2<- meltedS3
  meltedS3minusS2$value <- meltedS3$value-meltedS2$value
  
  
  colpalette <- hue_pal()(12)
  colpalette <- c(colpalette,"#909090")
  
  ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)
  
  lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)

  
p1 <- ggplot(data=meltedS1,aes(x=year,y=value,group=variable)) + 
   geom_line(aes(colour=variable,lwd=variable,lty=variable))+
   scale_colour_manual(values = colpalette)+
   scale_size_manual(values = lwdtest)+
   scale_linetype_manual(values = ltytest)+
   ylim(c(-40,40))+
    #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
    #geom_line(aes(x = years,y=JULES),colour="orange")+
    theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
   geom_hline(aes(yintercept=0),lty=2)+
    #plot.margin=margin(1,1,1,1,'cm')
    labs(y=bquote(Delta ~" cVeg S1" ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='a)')
  
p2 <- ggplot(data=meltedS2minusS1,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-40,40))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cVeg S2 - S1" ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='b)')

p3 <- ggplot(data=meltedS3minusS2,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-40,40))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cVeg S3 - S2" ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='c)')


grid.arrange(p1,p2,p3,ncol=3)#,layout_matrix = c(1,1,2,3))

  
#cSoil from model native resolutions

cSoilAllModelsS3 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_drylands_1901_2018.csv"))
cSoilAllModelsS2 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_S2_drylands_1901_2018.csv"))
cSoilAllModelsS1 <- data.frame(read.csv("./DGVM/DGVMdrylandTS/cSoil/cSoil_S1_drylands_1901_2018.csv"))


cSoilAllModelsS3[,2:13] <-  cSoilAllModelsS3[,2:13]-as.list(cSoilAllModelsS3[1,2:13])
cSoilAllModelsS2[,2:13] <-  cSoilAllModelsS2[,2:13]-as.list(cSoilAllModelsS2[1,2:13])
cSoilAllModelsS1[,2:13] <-  cSoilAllModelsS1[,2:13]-as.list(cSoilAllModelsS1[1,2:13])

meltedS3 <- melt(cSoilAllModelsS3,id.vars='year')

meltedS2 <- melt(cSoilAllModelsS2,id.vars='year')

meltedS1 <- melt(cSoilAllModelsS1,id.vars='year')

meltedS2minusS1<- meltedS2
meltedS2minusS1$value <- meltedS2$value-meltedS1$value

meltedS3minusS2<- meltedS3
meltedS3minusS2$value <- meltedS3$value-meltedS2$value

colpalette <- hue_pal()(12)
colpalette <- c(colpalette,"#909090")

ltytest <- c(1,2,1,2,1,2,1,2,1,2,1,2)

lwdtest <- c(1,1,1,1,1,1,1,1,1,1,1,1)

#melted <- melt(dfallyears,id.vars='years')
p1 <- ggplot(data=meltedS1,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-40,40))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cSoil S1 " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='a)')

#melted <- melt(dfallyears,id.vars='years')
p2 <- ggplot(data=meltedS2minusS1,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-40,40))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cSoil S2 - S1 " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='b)')

#melted <- melt(dfallyears,id.vars='years')
p3 <- ggplot(data=meltedS3minusS2,aes(x=year,y=value,group=variable)) + 
  geom_line(aes(colour=variable,lwd=variable,lty=variable))+
  scale_colour_manual(values = colpalette)+
  scale_size_manual(values = lwdtest)+
  scale_linetype_manual(values = ltytest)+
  ylim(c(-40,40))+
  #geom_line(aes(x = years,y=TRENDY),colour="red",lwd=2)+
  #geom_line(aes(x = years,y=JULES),colour="orange")+
  theme(legend.position='none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(face="bold"))+
  geom_hline(aes(yintercept=0),lty=2)+
  #plot.margin=margin(1,1,1,1,'cm')
  labs(y=bquote(Delta ~" cSoil S3 - S2 " ~ "["~ Pg ~ C ~ "]"),x='Time [yr]',title='c)')



grid.arrange(p1,p2,p3,ncol=3)#,layout_matrix = c(1,1,2,3))

