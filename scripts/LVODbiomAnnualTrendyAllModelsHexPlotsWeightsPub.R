

#spatial comparisons of AGC derived from TRENDY ensemble models to VOD derived 
#values, makes scatterplots and calculates statistics #using exactextract and 
#weights for regressions and statistics to deal with partial pixel coverage
#contains code to make spatial maps of AGC and GPP (see end)

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
library(lattice)
library(raster)

setwd('D:/Driving_C')

#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")
contsfordisp <- aggregate(continentshapes,dissolve=T)
 
#extract four studied continents: North and South America, Africa and Australia
contnrlist <- c(1,6,3,4)#number in continent shapefiles
studycontshapes <- aggregate(continentshapes[contnrlist,],dissolve=T)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE, converted to vector)

drylandclassshp <- readOGR(dsn = getwd(), layer = "drylandsglobal")
drylandclasssfc <- st_as_sfc(drylandclassshp) #spatialpolygonsdf to sfc for exactextractr

drylandclassSub <- readOGR(dsn = getwd(), layer = "drylands4contsub")
drylandclassSubsf <- st_as_sfc(drylandclassSub)


#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#2011 to 2018 (2010 does not have reliable values)
VODstack <- VODannualstacktot[[2:9]]


TRENDYmeanbrick <- brick('./DGVM/TRENDYAGC2011_2018v3.tif')

#create longitude based mask and apply to filter out northernmost and southernmost regions
# longmask <- TRENDYmeanbrick[[1]]
# longmask[,] <- 1
# longmask[1:34,] <- NA
# longmask[146:180,] <- NA

#resample VOD to 1 deg

VODstack1deg <- raster::resample(VODstack,TRENDYmeanbrick[[1]])


#drylands mask as vectorised from 0.25 deg LVOD data
VODdatamaskdrylands <- readOGR(getwd(),'VODdatamaskdrylands')

VODdatamaskdrylandssf <- st_as_sfc(VODdatamaskdrylands) #spatialpolygonsdf to sfc for exactextractr
 
#convert from kg/m2 to Mg/ha (*10)

TRENDYCarbonfinbrick <- TRENDYmeanbrick*10

#convert VOD to biomass
VODCarbonfinbrick <- VODstack1deg*52.48#calibration to Globbiomass 

#calculate mean per pixel over all years
VODCarbonfinmeans <- calc(VODCarbonfinbrick,mean,na.rm=T)
TRENDYCarbonfinmeans <- calc(TRENDYCarbonfinbrick,mean,na.rm=T)


#plot scatterplots of pixel TRENDY DGVM means and LVOD carbon density estimates

########

#list of scatterplots
meanscatterplotlist <- list()

#list of regression line plots
meanlineplotlist <- list()

#plot letter titles (currently unused)
titlenr <- 3
titlelist <- c('a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)')

DGVMstr <- 'TRENDY'

#get LVOD AGC pixel values and weights
VODCarbonyearmeansperpoly <- exactextractr::exact_extract(VODCarbonfinmeans,drylandclassSubsf,force_df=T)
VODCarbonyearmeans <- do.call('rbind',VODCarbonyearmeansperpoly)

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
modelvec <- modelvec[c(-10,-11,-15,-16)] #remove models where some TRENDY inputs are missing
maxcountvec <- rep(1,13)
#matrix of regression/correlation values (remove CCC)
coeffmat <- matrix(NA,ncol=6,nrow=14)
coeffmat[1,] <- c('model','pearsonsr','linCCC','slope','intercept','npix')
rcount=2

#iterate through models and compare to LVOD AGC, make hexbin plots and extract stats

for(j in 1:12){#add individual models
  
  #read cVeg from netcdf, aggregate to mean raster
  
  if(j %in% c(1,2,3,5)){ #if model AGC can be calculated, use this, otherwise use 0.4
    modelindex <- modelvec[j]
    TRENDYmodelstack <- stack(paste0("./DGVM/TRENDYmodelscVeg/AGC_2011_2018_1deggrids/",modelnames[modelindex],"_AGC2011_2018_1deg.tif"))
    TRENDYmodelCbrick <- TRENDYmodelstack#calc(TRENDYmodelstack,mean,na.rm=T)
    TRENDYmodelCbrick[TRENDYmodelCbrick<0] <- 0
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
  
  carbondryclassdf <- data.frame(VOD=VODCarbonDryClass,DGVM=DGVMCarbonDryClass)

  
  df <- data.frame(x = carbondryclassdf$VOD, y = carbondryclassdf$DGVM)
  
  DGVMname <- modelnames[modelindex]
  
  
  #make scatterplot
  meanscatterplotlist[[j]] <- ggplot(df,aes(x=x,y=y)) +
    geom_hex(bins = 30, show.legend=F)+

    theme_bw() +
    theme_classic() +
    theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
    labs(title=paste(titlelist[titlenr],DGVMname),x="LVOD C density",y=bquote("modelled C density"))+#+#,title=paste0('Carbon density: ',DGVMname)) +
    coord_fixed()+
    ylim(c(-10,180))+
    xlim(c(-10,180))
    titlenr <- titlenr+1

#calculate statistics and add regression line to plot

    #get value of highest bin for adjusting colour scale
maxcountvec[j] <- max(ggplot_build(meanscatterplotlist[[j]])$data[[1]]$count)   
      
VODvals <- VODCarbonyearmeans$value
DGVMvals <- DGVMCarbonyearmeans$value

VODweights <- VODCarbonyearmeans$coverage_fraction


demingreg <- deming(DGVMvals~VODvals,weights=VODweights)


#add diagonal and regression line
meanscatterplotlist[[j]] <- meanscatterplotlist[[j]]+geom_abline(slope=1,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])

pearsonsr <- wtd.cor(VODvals,DGVMvals,weight=VODweights)[1,1]#extract person's r from matrix

coeffmat[rcount,1] <- modelnames[modelindex]
coeffmat[rcount,2] <- pearsonsr
coeffmat[rcount,3] <- NA#epi.ccc(VODvals,DGVMvals)$rho.c[[1]]#to be replaced with weighted CCC
coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
coeffmat[rcount,6] <- sum(is.finite(VODvals*DGVMvals))#number of pixels
rcount <- rcount+1

}

#add TRENDY mean comparison and stats


VODCarbonyearmeansperpoly <- exactextractr::exact_extract(VODCarbonfinmeans,drylandclassSubsf,force_df=T)#,normalizeWeights=F,df=T)#mask(calc(TRENDYmodelCbrick,mean,na.rm=T),regionDrylandsmask)*0.4*10
VODCarbonyearmeans <- do.call('rbind',VODCarbonyearmeansperpoly)

DGVMCarbonyearmeansperpoly <- exactextractr::exact_extract(TRENDYCarbonfinmeans,drylandclassSubsf,force_df=T)#,normalizeWeights=F,df=T)#mask(calc(TRENDYmodelCbrick,mean,na.rm=T),regionDrylandsmask)*0.4*10
DGVMCarbonyearmeans <- do.call('rbind',DGVMCarbonyearmeansperpoly)

carbondryclassdf <- data.frame(VOD=VODCarbonyearmeans$value,DGVM=DGVMCarbonyearmeans$value)

df <- data.frame(x = carbondryclassdf$VOD, y = carbondryclassdf$DGVM)


DGVMname <- 'TRENDY'

meanscatterplotlist[[j+1]] <- ggplot(df,aes(x=x,y=y)) +
  geom_hex(bins = 30, show.legend=F)+
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  labs(title=paste(titlelist[titlenr],DGVMname),x="LVOD C density",y=bquote("modelled C density"))+#,title=paste0('Carbon density: ',DGVMname)) +
  coord_fixed()+
  ylim(c(-10,180))+
  xlim(c(-10,180))

#get value of highest bin for adjusting colour scale
maxcountvec[j+1] <- max(ggplot_build(meanscatterplotlist[[j+1]])$data[[1]]$count)   

#calculate statistics and add regression line to plot

VODvals <- VODCarbonyearmeans$value
DGVMvals <- DGVMCarbonyearmeans$value

demingreg <- deming(DGVMvals~VODvals,weights=VODweights)

pearsonsr <- wtd.cor(VODvals,DGVMvals,weight=VODweights)[1,1]#extract person's r from matrix

#add diagonal and regression line
meanscatterplotlist[[j+1]] <- meanscatterplotlist[[j+1]]+geom_abline(slope=1,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])

coeffmat[rcount,1] <- 'TRENDY'
coeffmat[rcount,2] <- pearsonsr
coeffmat[rcount,3] <- NA#epi.ccc(VODvals,DGVMvals)$rho.c[[1]]#replace with weighted CCC
coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
coeffmat[rcount,6] <- sum(is.finite(VODvals*DGVMvals))#number of pixels
rcount <- rcount+1

#for each plot, define colour scale using max count of all plots
for(n in 1:13){
  meanscatterplotlist[[n]] <-  meanscatterplotlist[[n]]+scale_fill_viridis(name = "count", trans = "log",breaks = 10^(0:6),limits=c(1,max(maxcountvec)))
}

#arrange and plot all scatterplots
grid.arrange(grobs=meanscatterplotlist,nrow=3,ncol=5)                     

#save stats table
write.csv(coeffmat,'./stats/DGVMvsLVOD_CarbonTrendStatsWeightedV2_GLOBBIOMASS_AGC_nonVODdatafiltered.csv')


############
#maps of trends in dryland carbon

#functions to create diverging plots in levelplot
getMinMax <- function(inraster){
  if(abs(cellStats(inraster,min,na.rm=T))>abs(cellStats(inraster,max,na.rm=T))){
    colbarbounds <- c(cellStats(inraster,min,na.rm=T),abs(cellStats(inraster,min,na.rm=T)))
  }else{
    colbarbounds <- c(-cellStats(inraster,max,na.rm=T),cellStats(inraster,max,na.rm=T))
  }
  return(colbarbounds)
}

diverge0 <- function(p, ramp) { #source: https://rdrr.io/github/stevenpawley/Diverge0/src/R/diverge0.R

  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
  p
} 



my.settings <- list(par.main.text = list(font = 2, just = "left",  x = grid::unit(5, "mm")),panel.background=list(col="lightgrey"))


drylandsubshapes <- readOGR(dsn = getwd(), layer = "drylands4contsub")
#intersection of VODdatamaskpoly and drylandsubconts done in ArcMap (error in R)
VODdatamaskdrylands <- readOGR(getwd(),'VODdatamaskdrylands')

#AGC trends over time period
fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
#PMLannualdrygpptrend <- calc(PMLannualdrygpp,fun1)
VODCarbonfintrend <- calc(VODCarbonfinbrick,fun1)
TRENDYCarbonfintrend <- calc(TRENDYCarbonfinbrick,fun1)

#VOD trends


VODtrenddisp <- raster::mask(raster::mask(VODCarbonfintrend,drylandclassshp),studycontshapes)


trendMinMax <- getMinMax(VODtrenddisp)

VODbiomtrendplot <- diverge0(levelplot(VODtrenddisp,par.settings=my.settings,main="a) L-VOD",at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+latticeExtra::layer(sp.polygons(drylandsubshapes,col='darkgrey'))
plot(VODbiomtrendplot)

#TRENDY trends

TRENDYtrenddisp <- raster::mask(raster::mask(TRENDYCarbonfintrend,drylandclassshp),studycontshapes)

trendMinMax <- getMinMax(TRENDYtrenddisp)

TRENDYbiomtrendplot <- diverge0(levelplot(TRENDYtrenddisp,par.settings=my.settings,main="b) TRENDY-mean",at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+latticeExtra::layer(sp.polygons(drylandsubshapes,col='darkgrey'))

plot(TRENDYbiomtrendplot)



#GPP trends over time period

GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew_2001_2018_v016.tif")

TRENDYannualgpp <- brick('./DGVM/TRENDYGPP_2001_2018v3.tif')*10#from kgC per m2 to MgC per ha

GPPstackresamp <- raster::resample(GPPstack,TRENDYannualgpp[[1]])
PMLannualgpp <- GPPstackresamp/100 #from gC per m2 to MgC per ha

PMLannualdrygpp <- raster::mask(PMLannualgpp,drylandclassshp)#raster::mask(PMLannualgpp,drylandmask)
TRENDYannualdrygpp <- raster::mask(TRENDYannualgpp,drylandclassshp)

fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
#PMLannualdrygpptrend <- calc(PMLannualdrygpp,fun1)
PMLGPPfintrend <- calc(PMLannualdrygpp,fun1)
TRENDYGPPfintrend <- calc(TRENDYannualdrygpp,fun1)

#PMLv2 trends


trendMinMax <- getMinMax(PMLGPPfintrend)

PMLbiomtrendplot <- diverge0(levelplot(PMLGPPfintrend,par.settings=my.settings,main="a) MODIS PML-v2",at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))
plot(PMLbiomtrendplot)

#TRENDY trends

trendMinMax <- getMinMax(TRENDYGPPfintrend)

TRENDYbiomtrendplot <- diverge0(levelplot(TRENDYGPPfintrend,par.settings=my.settings,main="b) TRENDY-mean",at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))#+layer(sp.polygons(drylandsubshapes,col='darkgrey'))

plot(TRENDYbiomtrendplot)


#PER MODEL map figures of DGVM and VOD mean AGC and GPP 
#########

#continents for display (no antarctica)
contsfordispGPP <- aggregate(continentshapes[-7,],dissolve=T)

GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew_2001_2018_v016.tif")

TRENDYannualgpp <- brick('./DGVM/TRENDYGPP_2001_2018v3.tif')*10#from kgC per m2 to MgC per ha

#mask to include any cells intersecting drylands mask

# drylandclass_ras <- rasterize(drylandclass, PMLannualgpp[[1]], getCover=TRUE)
# drylandclass_ras[drylandclass_ras==0] <- NA

#means over time period
PMLannualdrygppmean <- calc(PMLannualdrygpp,mean,na.rm=T)
TRENDYannualdrygppmean <- calc(TRENDYannualdrygpp,mean,na.rm=T)

#trends over time period
fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates) }}
PMLannualdrygpptrend <- calc(PMLannualdrygpp,fun1)
VODCarbonfintrend <- calc(VODCarbonfinbrick,fun1)

#difference between trendy mean GPP and PMLv2 GPP
PMLTRENDYDiff <- TRENDYannualdrygppmean-PMLannualdrygppmean

trendMinMax <- getMinMax(PMLTRENDYDiff)

my.settings <- list(par.main.text = list(font = 2, just = "left",  x = grid::unit(5, "mm")),panel.background=list(col="lightgrey"))

#plot difference between PMLv2 GPP and TRENDY mean GPP
differencePMLTRENDY <- diverge0(levelplot(PMLTRENDYDiff,par.settings=my.settings,main=bquote("b) Mean GPP difference: TRENDY - MODIS " ~ "["~ Mg ~ C ~ ha^{-1} ~ y^{-1}~"]"),at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+latticeExtra::layer(sp.polygons(contsfordispGPP,col='black',lwd=1.5))

plot(differencePMLTRENDY)

cols <- colorRampPalette(brewer.pal(9,"YlGn"))
cols2 <- colorRampPalette(brewer.pal(9,"YlGnBu"))

VODCarbonfinmeans <- calc(VODCarbonfinbrick,mean,na.rm=T)

VODCarbonfinmeanstotalobs <- raster::mask(VODCarbonfinmeans,VODdatamaskdrylands)

VODTRENDYDiff <- TRENDYCarbonfinmeans-VODCarbonfinmeans 

drylandmask4conts <- raster::mask(drylandmask,studycontshapes)

trendMinMax <- getMinMax(raster::mask(VODTRENDYDiff,drylandclassSub))

differenceVODTRENDY <- diverge0(levelplot(raster::mask(VODTRENDYDiff,drylandclassSub),par.settings=my.settings,main=bquote("b) Mean AGC difference: TRENDY - L-VOD " ~ "["~ Mg ~ C ~ ha^{-1} ~"]"),at=seq(trendMinMax[1], trendMinMax[2], len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+latticeExtra::layer(sp.polygons(studycontshapes,col='black',lwd=1.5))

plot(differenceVODTRENDY)

#L-VOD biomass global plot

VODbiomplot <- levelplot(raster::mask(VODCarbonfinmeans,drylandclassSub),main=bquote("a) Mean AGC: L-VOD " ~ "["~ Mg ~ C ~ ha^{-1} ~ "]"),
                         par.settings=my.settings,at=seq(0, cellStats(raster::mask(VODCarbonfinmeans,drylandclassSub),max), len = 100),margin=FALSE,col.regions=cols,maxpixels = 2e10,add=T)+
                          latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+
                          latticeExtra::layer(sp.polygons(studycontshapes,col='black',lwd=1.5))
plot(VODbiomplot)

#PMLv2 GPP global plot
PMLGPPplot <- levelplot(PMLannualdrygppmean,main=bquote("a) Mean GPP: MODIS PML-v2 " ~ "["~ Mg ~ C ~ ha^{-1} ~ yr^{-1}~"]"),
                        par.settings=my.settings,at=seq(0, cellStats(PMLannualdrygppmean,max), len = 100),margin=FALSE,col.regions=cols2,maxpixels = 2e10,add=T)+
                          latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))+
                          latticeExtra::layer(sp.polygons(contsfordispGPP,col='black',lwd=1.5))
plot(PMLGPPplot)


TRENDYmodelsmeanGPP <- brick('./DGVM/TRENDYpermodelGPP2001_2018v3.tif')*10#from kgC per m2 to MgC per ha
TRENDYmodelsmeancVeg <- brick('./DGVM/TRENDYpermodelAGC2011_2018v3.tif')*10

TRENDYmodelstrendGPP <- brick('./DGVM/TRENDYpermodeltrendGPP2001_2018v3.tif')*10
TRENDYmodelstrendcVeg <- brick('./DGVM/TRENDYpermodeltrendAGC2011_2018v3.tif')*10

#models vs obs differences means
TRENDYmodelsmeanGPPdrylandsdiff <- raster::mask(TRENDYmodelsmeanGPP,drylandclassshp)-PMLannualdrygppmean 
names(TRENDYmodelsmeanGPPdrylandsdiff) <- modelnames[c(-10,-11,-15,-16)]
TRENDYmodelsmeancVegdrylandsdiff <- raster::mask(TRENDYmodelsmeancVeg,drylandclassshp)-VODCarbonfinmeans 
names(TRENDYmodelsmeancVegdrylandsdiff) <- modelnames[c(-10,-11,-15,-16)]

#models vs obs differences trends
TRENDYmodelstrendGPPdrylandsdiff <- raster::mask(TRENDYmodelstrendGPP,drylandclassshp)-PMLannualdrygpptrend
names(TRENDYmodelstrendGPPdrylandsdiff) <- modelnames[c(-10,-11,-15,-16)]
TRENDYmodelstrendcVegdrylandsdiff <- raster::mask(TRENDYmodelstrendcVeg,drylandclassshp)-VODCarbonfintrend
names(TRENDYmodelstrendcVegdrylandsdiff) <- modelnames[c(-10,-11,-15,-16)]


###GPP difference global plots for all models

#mean
GPPdiffmaxval <- cellStats(calc(TRENDYmodelsmeanGPPdrylandsdiff,max),max)
GPPdiffminval <- cellStats(calc(TRENDYmodelsmeanGPPdrylandsdiff,min),min)
if(abs(GPPdiffmaxval)>abs(GPPdiffminval)){
  GPPdiffextremeval <- abs(GPPdiffmaxval)
}else{
  GPPdiffextremeval <- abs(GPPdiffminval)

}

diverge0(levelplot(TRENDYmodelsmeanGPPdrylandsdiff,par.settings=my.settings,main=bquote("Model - MODIS PML-v2 mean annual GPP " ~ "["~ Mg ~ C ~ ha^{-1} ~ yr^{-1}~"]"),at=seq(-GPPdiffextremeval, GPPdiffextremeval, len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))

#trend

GPPdiffmaxval <- cellStats(calc(TRENDYmodelstrendGPPdrylandsdiff,max),max)
GPPdiffminval <- cellStats(calc(TRENDYmodelstrendGPPdrylandsdiff,min),min)
if(abs(GPPdiffmaxval)>abs(GPPdiffminval)){
  GPPdiffextremeval <- abs(GPPdiffmaxval)
}else{
  GPPdiffextremeval <- abs(GPPdiffminval)
  
}

diverge0(levelplot(TRENDYmodelstrendGPPdrylandsdiff,par.settings=my.settings,main=bquote("Model - MODIS PML-v2 annual GPP trend " ~ "["~ Mg ~ C ~ ha^{-1} ~ yr^{-1}~"]"),at=seq(-GPPdiffextremeval, GPPdiffextremeval, len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))

###AGC difference global plots for all models

#mean

cVegdiffmaxval <- cellStats(calc(TRENDYmodelsmeancVegdrylandsdiff,max),max)
cVegdiffminval <- cellStats(calc(TRENDYmodelsmeancVegdrylandsdiff,min),min)
if(abs(cVegdiffmaxval)>abs(cVegdiffminval)){
  cVegdiffextremeval <- abs(cVegdiffmaxval)
}else{
  cVegdiffextremeval <- abs(cVegdiffminval)
  
}

diverge0(levelplot(raster::mask(TRENDYmodelsmeancVegdrylandsdiff,studycontshapes),par.settings=my.settings,main=bquote("Mean AGC difference: Model - L-VOD " ~ "["~ Mg ~ C ~ ha^{-1}~"]"),at=seq(-cVegdiffextremeval, cVegdiffextremeval, len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))

#trend

cVegdiffmaxval <- cellStats(calc(TRENDYmodelstrendcVegdrylandsdiff,max),max)
cVegdiffminval <- cellStats(calc(TRENDYmodelstrendcVegdrylandsdiff,min),min)
if(abs(cVegdiffmaxval)>abs(cVegdiffminval)){
  cVegdiffextremeval <- abs(cVegdiffmaxval)
}else{
  cVegdiffextremeval <- abs(cVegdiffminval)
}

diverge0(levelplot(raster::mask(TRENDYmodelstrendcVegdrylandsdiff,studycontshapes),par.settings=my.settings,main=bquote("Mean AGC trend difference: Model - L-VOD " ~ "["~ Mg ~ C ~ ha^{-1}~yr^{-1}~"]"),at=seq(-cVegdiffextremeval, cVegdiffextremeval, len = 100),margin=FALSE,maxpixels = 2e10),colorRampPalette(c('red','white','blue')))+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))
