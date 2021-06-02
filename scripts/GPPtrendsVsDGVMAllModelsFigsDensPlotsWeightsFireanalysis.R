#spatial comparisons of GPP derived from TRENDY ensemble models to MODIS PMLv2 derived values, makes scatterplots and calculates statistics
#using exactextract and weights for regressions and statistics to deal with partial pixel coverage
#NOTE: weighted CCC currently not implemented


library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(epiR)
library(trend)
library(ggplot2)
library(gridExtra)
library(ncdf4)
library(deming) 
library(weights)

setwd('D:/Driving_C')

#continent outlines for plotting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")

contsfordisp <- aggregate(continentshapes,dissolve=T)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- readOGR(dsn = getwd(), layer = "drylandsglobal")
drylandclasssfc <- st_as_sfc(drylandclass) #spatialpolygonsdf to sfc for exactextractr

#load burned area 2000-2021

MODburnedarea <- raster("D:/Driving_C/BurnedArea/MCD64A1_totaldaterange_burnsum_10km_ag.tif")

#GPP stack PMLv2

GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew.tif")

years <- seq(2003,2018,1)

TRENDYannualgpp <- brick('./DGVM/TRENDYGPP_2003_2018v3.tif')*10 #from kgC per m2 to MgC per ha

GPPstackresamp <- raster::resample(GPPstack,TRENDYannualgpp[[1]])
PMLannualgpp <- GPPstackresamp/100 #from gC per m2 to MgC per ha

MODburnedarea1deg <- raster::resample(MODburnedarea,TRENDYannualgpp[[1]])

#means over time period
PMLGPPfinmeans <- calc(PMLannualgpp,mean,na.rm=T)
TRENDYGPPfinmeans <- calc(TRENDYannualgpp,mean,na.rm=T)

#PML_TRENDYdiff <- TRENDYannualdrygppmean-PMLannualdrygppmean 
#PML_JULESdiff <- JULESannualdrygppmean-PMLannualdrygppmean

#plot scatterplots of pixel TRENDY GPP means and PML carbon density estimates
########

#list of scatterplots
meanscatterplotlist <- list()

#list of regression line plots
meanlineplotlist <- list()

#letter titles not used currently
titlenr <- 3
titlelist <- c('a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)')

#get PMLv2 GPP pixel values and weights
PMLGPPyearmeansperpoly <- exactextractr::exact_extract(PMLGPPfinmeans,drylandclasssfc,force_df=T)
PMLGPPyearmeans <- do.call('rbind',PMLGPPyearmeansperpoly)

#get Burned area pixel values and weights
MODburnedarea1degperpoly <- exactextractr::exact_extract(MODburnedarea1deg,drylandclasssfc,force_df=T)
MODburnedarea1degmeans <- do.call('rbind',MODburnedarea1degperpoly)

  
#get TRENDY mean GPP pixel values and weights
TRENDYGPPyearmeansperpoly <- exactextractr::exact_extract(TRENDYGPPfinmeans,drylandclasssfc,force_df=T)
TRENDYGPPyearmeans <- do.call('rbind',TRENDYGPPyearmeansperpoly)

  carbonyearmeansdf <- data.frame(PML=PMLGPPyearmeans$value,DGVM=TRENDYGPPyearmeans$value)
  
  #load regridded (1 deg) TRENDY models netcdf file and info
  ncin <- nc_open(paste0("./DGVM/trendyv8_S3_cVeg_1901-2018.nc"))
  
  modelnames <- ncatt_get(ncin,0,"models")
  modelnames <- unlist(strsplit(modelnames$value,' '))
  
  nc_close(ncin)
  
  TRENDYmodelGPPnames <- list.files('./DGVM/TRENDYmodelsGPP')
  
  #matrix of regression/correlation values (remove CCC)
  coeffmat <- matrix(NA,ncol=6,nrow=14)
  coeffmat[1,] <- c('model','pearsonsr','pval','slope','intercept','npix')
  rcount=2
  
  modelvec <- seq(1,16,1)
  modelvec <- modelvec[c(-10,-11,-15,-16)]
  
  for(j in 1:12){#add individual models
    
    #read precomputed annual GPP per model from raster
    
    modelindex <- modelvec[j]
    
    TRENDYmodelGPPbrick <- brick(paste0('./DGVM/TRENDYmodelsGPP/',TRENDYmodelGPPnames[[ modelindex ]]))
    TRENDYGPPfinbrick <- TRENDYmodelGPPbrick*10 # kg per m2 to Mg C per ha
    
    DGVMGPPyearmeansperpoly <- exactextractr::exact_extract(calc(TRENDYGPPfinbrick,mean,na.rm=T),drylandclasssfc,force_df=T)
    DGVMGPPyearmeans <- do.call('rbind',DGVMGPPyearmeansperpoly)
    
    DGVMGPPDryClass <- DGVMGPPyearmeans$value
    PMLGPPyearmeansDryClass <- PMLGPPyearmeans$value
    MODburnedarea1degDryClass <- MODburnedarea1degmeans$value
    
    DeltaDGVMvsPMLGPP <- DGVMGPPDryClass-PMLGPPyearmeansDryClass 
  

    GPPdryclassdf <- data.frame(DGVMPMLdelta=DeltaDGVMvsPMLGPP,MODburned=MODburnedarea1degDryClass)
    
    df <- data.frame(x = GPPdryclassdf$MODburned, y = GPPdryclassdf$DGVMPMLdelta)
    
    DGVMname <- modelnames[modelindex]
    
    meanscatterplotlist[[j]] <- ggplot(df) +
      geom_point(aes(x, y),size=0.5,alpha=0.1) +
      scale_color_identity() +
      theme_bw() +
      theme_classic() +
      theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
      labs(title=paste(titlelist[titlenr],DGVMname),x=bquote("Burned frac. sum"),y="modelled - PMLv2 GPP")+
      #coord_fixed()+
      ylim(c(-40,40))+
      xlim(c(0,20))
    titlenr <- titlenr+1

    PMLvals <- PMLGPPyearmeans$value
    DGVMvals <- DGVMGPPyearmeans$value
    
    
    
    #calculate statistics and add regression line to plot
    
   # PMLvals <- PMLGPPyearmeans$value
    #DGVMvals <- DGVMGPPyearmeans$value
    
    PMLweights <- PMLGPPyearmeans$coverage_fraction
    
    demingreg <- deming(DeltaDGVMvsPMLGPP~MODburnedarea1degDryClass,weights=PMLweights)
    
    #add diagonal and regression line
    meanscatterplotlist[[j]] <- meanscatterplotlist[[j]]+geom_abline(slope=0,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])
    
    pearsonsr <- wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsPMLGPP,weight=PMLweights)[1,1]#extract person's r from matrix
    
    pval <- wtd.cor(MODburnedarea1degDryClass,DeltaDGVMvsPMLGPP,weight=PMLweights)[1,4]#extract person's r from matrix
    
    
    coeffmat[rcount,1] <- modelnames[modelindex]
    coeffmat[rcount,2] <- pearsonsr
    coeffmat[rcount,3] <- pval#to be replaced with weighted CCC
    coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
    coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
    coeffmat[rcount,6] <- sum(is.finite(DeltaDGVMvsPMLGPP*MODburnedarea1degDryClass))#number of pixels
    rcount <- rcount+1
    }

  #add TRENDY mean comparison and stats
  TRENDYPMLdeltaGPP = carbonyearmeansdf$DGVM-carbonyearmeansdf$PML
  
  df <- data.frame(x = MODburnedarea1degDryClass, y = TRENDYPMLdeltaGPP)#,

  
  DGVMname <- 'TRENDY'
  
  meanscatterplotlist[[j+1]] <- ggplot(df) +
    geom_point(aes(x, y),size=0.5,alpha=0.1) +
    scale_color_identity() +
    theme_bw() +
    theme_classic() +
    theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
    labs(title=paste(titlelist[titlenr],DGVMname),x=bquote("Burned frac. sum"),y="modelled - PMLv2 GPP")+
   # coord_fixed()+
    ylim(c(-40,40))+
    xlim(c(0,20))


  demingreg <- deming(TRENDYPMLdeltaGPP~MODburnedarea1degDryClass,weights=PMLweights)
  
  
  #add diagonal and regression line
  meanscatterplotlist[[j+1]] <- meanscatterplotlist[[j+1]]+geom_abline(slope=0,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])
  
  pearsonsr <- wtd.cor(MODburnedarea1degDryClass,TRENDYPMLdeltaGPP,weight=PMLweights)[1,1]#extract person's r from matrix
  
  pval <-  wtd.cor(MODburnedarea1degDryClass,TRENDYPMLdeltaGPP,weight=PMLweights)[1,4]#extract person's r from matrix

  
  coeffmat[rcount,1] <- 'TRENDY'
  coeffmat[rcount,2] <- pearsonsr
  coeffmat[rcount,3] <- pval
  coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
  coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
  coeffmat[rcount,6] <- sum(is.finite(MODburnedarea1degDryClass*TRENDYPMLdeltaGPP))#number of pixels
  rcount <- rcount+1
 
  #arrange and plot all scatterplots
grid.arrange(grobs=meanscatterplotlist,nrow=3,ncol=5)

#save stats table
write.csv(coeffmat,'./stats/DGVMvsMODIS_GPP_BurnedArea_TrendStatsWeightedV1.csv')


#make burned area map

diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
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

MODburnedarea1degmasked <- mask(MODburnedarea1deg,drylandclass)
my.settings <- list(par.main.text = list(font = 2, just = "left",  x = grid::unit(5, "mm")),panel.background=list(col="lightgrey"))


levelplot(MODburnedarea1degmasked,main=bquote("a) Area burned 2000-2021 " ~ "["~ Mg ~ C ~ ha^{-1} ~ "]"),
                         par.settings=my.settings,at=seq(0, cellStats(TRENDYCarbonfinmeans,max), len = 100),margin=FALSE,maxpixels = 2e10,add=T)+
  latticeExtra::layer(sp.polygons(contsfordisp,col='grey'))#+
  #latticeExtra::layer(sp.polygons(studycontshapes,col='black',lwd=1.5))



