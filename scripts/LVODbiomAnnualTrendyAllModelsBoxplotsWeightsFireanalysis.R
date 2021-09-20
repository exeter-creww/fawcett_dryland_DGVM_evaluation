

#burnt area analysis of differences between AGC derived from TRENDY ensemble models VOD derived values, 
#compares differences for different burnt area frequency bins

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
library(epiR)
library(ggplot2)
library(weights)
library(deming)
library(cccrm)
library(exactextractr)
library(sf)
library(RColorBrewer)
library(latticeExtra)

setwd('D:/Driving_C')
 
#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")
contsfordisp <- aggregate(continentshapes,dissolve=T)
 
#extract four studied continents: North and South America, Africa and Australia
contnrlist <- c(1,6,3,4)#number in continent shapefiles
studycontshapes <- aggregate(continentshapes[contnrlist,],dissolve=T)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE), for visualisation only
drylandclass <- raster("./Plots/drylandclassclipfin.tif")

#load burned area 2001-2018

MODburnedarea <- raster("D:/Driving_C/BurnedArea/MCD64A1_2018daterange_burnsum_10km_ag.tif")

#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#2011 to 2018 (2010 does not have reliable values)
VODstack <- VODannualstacktot[[2:9]]


TRENDYmeanbrick <- brick('./DGVM/TRENDYAGC2011_2018v3.tif')


#resample VOD and burned area to 1 deg

VODstack1deg <- raster::resample(VODstack,TRENDYmeanbrick[[1]])
MODburnedarea1deg <- raster::resample(MODburnedarea,TRENDYmeanbrick[[1]])
MODburnedarea1deg <- MODburnedarea1deg/18 # #divide by years in record -> burn frequency

drylandclassSub <- readOGR(dsn = getwd(), layer = "drylands4contsub")#raster("D:/Driving_C/Plots/drylandclassclipfin.tif")

drylandclassSubsf <- st_as_sfc(drylandclassSub)

#convert from kg/m2 to Mg/ha (*10)

TRENDYCarbonfinbrick <- TRENDYmeanbrick*10

#convert VOD to biomass
VODCarbonfinbrick <- VODstack1deg*52.48  #calibration to Globbiomass

#calculate mean per pixel over all years
VODCarbonfinmeans <- calc(VODCarbonfinbrick,mean,na.rm=T)
TRENDYCarbonfinmeans <- calc(TRENDYCarbonfinbrick,mean,na.rm=T)


#plot boxplots of pixel TRENDY DGVM means and LVOD carbon density estimates

########

#list of boxplots
meanboxplotlist <- list()

#list of regression line plots
meanlineplotlist <- list()

#plot letter titles (currently unused)
titlenr <- 1
titlelist <- c('a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)')

DGVMstr <- 'TRENDY'

#get LVOD AGC pixel values and weights
VODCarbonyearmeansperpoly <- exactextractr::exact_extract(VODCarbonfinmeans,drylandclassSubsf,force_df=T)
VODCarbonyearmeans <- do.call('rbind',VODCarbonyearmeansperpoly)

#get Burned area pixel values and weights
MODburnedarea1degperpoly <- exactextractr::exact_extract(MODburnedarea1deg,drylandclassSubsf,force_df=T)
MODburnedarea1degmeans <- do.call('rbind',MODburnedarea1degperpoly)

#load regridded (1 deg) TRENDY models netcdf file and info
ncin <- nc_open(paste0("./DGVM/trendyv8_S3_cVeg_1901-2018.nc"))


modelnames <- ncatt_get(ncin,0,"models")
modelnames <- unlist(strsplit(modelnames$value,' '))

lonDGVM <- ncvar_get(ncin, "lon")
nlonDGVM <- dim(lonDGVM)

latDGVM <- ncvar_get(ncin, "lat")
nlatDGVM <- dim(latDGVM)

cVeg <- ncvar_get(ncin,'cVeg')
fillvalue <- ncatt_get(ncin,'cVeg',"_FillValue")

cVeg[cVeg==fillvalue$value] <- NA


modelvec <- seq(1,16,1)
modelvec <- modelvec[c(-10,-11,-15,-16)]

#matrix of regression/correlation values (remove CCC)
coeffmat <- matrix(NA,ncol=6,nrow=14)
coeffmat[1,] <- c('model','pearsonsr','pval','slope','intercept','npix')
rcount=2

#iterate through models

for(j in 1:12){#add individual models
  
  #read cVeg from netcdf, aggregate to mean raster
  
  if(j %in% c(1,2,3,5)){ #if model AGC can be calculated, use this, otherwise use 0.4
    modelindex <- modelvec[j]
    TRENDYmodelstack <- stack(paste0("D:/Driving_C/DGVM/TRENDYmodelscVeg/",modelnames[j],"_AGC2011_2018_1deg.tif"))
    TRENDYmodelCbrick <- TRENDYmodelstack#calc(TRENDYmodelstack,mean,na.rm=T)
    
  }else{
  
  modelindex <- modelvec[j]
  
  TRENDYmodelstack <- cVeg[,,111:118,modelindex]
  
  TRENDYmodelmeanstack <- apply(TRENDYmodelstack,c(1,2,3),FUN=mean,na.rm=T)
  
  TRENDYmodelCbrick <- t(raster::flip(brick(TRENDYmodelmeanstack),1))
  
  extent(TRENDYmodelCbrick) <- c(-180, 180, -90, 90)
  projection(TRENDYmodelCbrick) <- CRS("+init=epsg:4326")
  
  TRENDYmodelCbrick <- TRENDYmodelCbrick*0.4 #convert to AGC
  
  }
  
  #extract pixel values and convert to Mg C ha
  DGVMCarbonyearmeansperpoly <- exactextractr::exact_extract(calc(TRENDYmodelCbrick,mean,na.rm=T)*10,drylandclassSubsf,force_df=T)#,normalizeWeights=F,df=T)#mask(calc(TRENDYmodelCbrick,mean,na.rm=T),regionDrylandsmask)*0.4*10
  DGVMCarbonyearmeans <- do.call('rbind',DGVMCarbonyearmeansperpoly)

  DGVMCarbonDryClass <- DGVMCarbonyearmeans$value
  VODCarbonDryClass <- VODCarbonyearmeans$value
  MODburnedarea1degDryClass <- MODburnedarea1degmeans$value
  
  DeltaDGVMvsVODcarbon <- DGVMCarbonDryClass-VODCarbonDryClass
  
  carbondryclassdf <- data.frame(DGVMVODdelta=DeltaDGVMvsVODcarbon,MODburned=MODburnedarea1degDryClass)

  
  
  df <- data.frame(x = carbondryclassdf$MODburned, y = carbondryclassdf$DGVMVODdelta)
  
  df$bin <- cut(df$x, c(0,0.1, 0.2, 0.3, 0.4, 0.5, 1),include.lowest=T)
  
  
  DGVMname <- modelnames[modelindex]
  
  # define the summary function
  f <- function(x) {
    r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  meanboxplotlist[[j]] <- ggplot(df,aes(bin,y)) +
    stat_summary(fun.data = f, geom="boxplot")+
   # geom_boxplot(aes(bin, y),outlier.shape = NA)+
    
    #geom_point(aes(x, y),size=0.5,alpha=0.1) +
    #scale_color_identity() +
    theme_bw() +
    theme_classic() +
    theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
    labs(title=paste(titlelist[titlenr],DGVMname),x=bquote("Burned frequency"),y="modelled - LVOD C density")+#+
    geom_abline(slope=0,intercept=0,linetype='dashed')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    #coord_fixed()+
    ylim(c(-150,150))
  #xlim(c(0,1))
  titlenr <- titlenr+1
  



VODweights <- VODCarbonyearmeans$coverage_fraction


demingreg <- deming(DeltaDGVMvsVODcarbon~MODburnedarea1degDryClass,weights=VODweights)

pearsonsr <- wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsVODcarbon,weight=VODweights)[1,1]#extract person's r from matrix

pval <-  wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsVODcarbon,weight=VODweights)[1,4]#extract person's r from matrix

coeffmat[rcount,1] <- modelnames[modelindex]
coeffmat[rcount,2] <- pearsonsr
coeffmat[rcount,3] <- pval#epi.ccc(VODvals,DGVMvals)$rho.c[[1]]#to be replaced with weighted CCC
coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
coeffmat[rcount,6] <- sum(is.finite(MODburnedarea1degDryClass*DeltaDGVMvsVODcarbon))#number of pixels
rcount <- rcount+1

}

#add TRENDY mean comparison and stats


VODCarbonyearmeansperpoly <- exactextractr::exact_extract(VODCarbonfinmeans,drylandclassSubsf,force_df=T)#,normalizeWeights=F,df=T)#mask(calc(TRENDYmodelCbrick,mean,na.rm=T),regionDrylandsmask)*0.4*10
VODCarbonyearmeans <- do.call('rbind',VODCarbonyearmeansperpoly)

DGVMCarbonyearmeansperpoly <- exactextractr::exact_extract(TRENDYCarbonfinmeans,drylandclassSubsf,force_df=T)#,normalizeWeights=F,df=T)#mask(calc(TRENDYmodelCbrick,mean,na.rm=T),regionDrylandsmask)*0.4*10
DGVMCarbonyearmeans <- do.call('rbind',DGVMCarbonyearmeansperpoly)

DeltaDGVMvsVODcarbon <- DGVMCarbonyearmeans$value-VODCarbonyearmeans$value

carbondryclassdf <- data.frame(DGVMVODdelta=DeltaDGVMvsVODcarbon,MODburned=MODburnedarea1degDryClass)

df <- data.frame(x = carbondryclassdf$MODburned, y = carbondryclassdf$DGVMVODdelta)



df$bin <- cut(df$x, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1),include.lowest=T)

test <- df[order(df$bin),]

DGVMname <- 'TRENDY'

meanboxplotlist[[j+1]] <- ggplot(df,aes(bin,y)) +
 # geom_boxplot(aes(bin, y),outlier.shape = NA)+
  stat_summary(fun.data = f, geom="boxplot")+
  #geom_point(aes(x, y),size=0.5,alpha=0.1) +
  # scale_color_identity() +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  geom_abline(slope=0,intercept=0,linetype='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title=paste(titlelist[titlenr],DGVMname),x=bquote("Burned frequency"),y="modelled - LVOD C density")+
  ylim(c(-150,150))#+


#calculate statistics and add regression line to plot


demingreg <- deming(DeltaDGVMvsVODcarbon~MODburnedarea1degDryClass,weights=VODweights)


pearsonsr <- wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsVODcarbon,weight=VODweights)[1,1]#extract person's r from matrix

pval <-  wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsVODcarbon,weight=VODweights)[1,4]#extract person's r from matrix


coeffmat[rcount,1] <- 'TRENDY'
coeffmat[rcount,2] <- pearsonsr
coeffmat[rcount,3] <- pval#epi.ccc(VODvals,DGVMvals)$rho.c[[1]]#replace with weighted CCC
coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
coeffmat[rcount,6] <- sum(is.finite(MODburnedarea1degDryClass*DeltaDGVMvsVODcarbon))#number of pixels
rcount <- rcount+1



#make histogram


dfexclsmall <- data.frame(x = MODburnedarea1degDryClass, y = DeltaDGVMvsVODcarbon)
dfexclsmall[dfexclsmall$x<0.1,] <- NA
dfexclsmall <- dfexclsmall[complete.cases(dfexclsmall),]

dfexclsmall$bin <- cut(dfexclsmall$x, c(0.1, 0.2, 0.3, 0.4, 0.5, 1),include.lowest=F)

dfonlysmall <- data.frame(x = MODburnedarea1degDryClass, y = DeltaDGVMvsVODcarbon)


dfonlysmall$bin <- cut(dfonlysmall$x, c(0,0.1, 0.2, 0.3, 0.4, 0.5),include.lowest=T)

f <- function(x) {
  r <- length(x)
  r
}

stat_summary(fun.data=f)



barplot1 <- ggplot(dfexclsmall,aes(bin)) +
  geom_bar()+
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  #geom_abline(slope=0,intercept=0,linetype='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title=paste(titlelist[titlenr+1]),x=bquote("Burned frequency"),y="count")


barplot2 <- ggplot(dfonlysmall,aes(bin)) +
  geom_bar()+
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  #geom_abline(slope=0,intercept=0,linetype='dashed')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title=paste(titlelist[titlenr+1]),x=bquote("Burned frequency"),y="count")


meanboxplotlist[[j+2]] <- barplot1
meanboxplotlist[[j+3]] <- barplot2

#arrange and plot all boxplots
grid.arrange(grobs=meanboxplotlist,nrow=3,ncol=5)                     

#save stats table
write.csv(coeffmat,'./stats/deltaDGVM_LVOD_Carbon_BurnedArea_TrendStatsWeightedV2_GLOBBIOMASS_AGC.csv')



