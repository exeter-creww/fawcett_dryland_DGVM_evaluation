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


#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- readOGR(dsn = getwd(), layer = "drylandsglobal")
drylandclasssfc <- st_as_sfc(drylandclass) #spatialpolygonsdf to sfc for exactextractr
#GPP stack PMLv2

GPPstack <- stack("./PMLV2sampled/PMLv2GPPstack10knew.tif")

years <- seq(2003,2018,1)

TRENDYannualgpp <- brick('./DGVM/TRENDYGPP_2003_2018v3.tif')*10 #from kgC per m2 to MgC per ha

GPPstackresamp <- raster::resample(GPPstack,TRENDYannualgpp[[1]])
PMLannualgpp <- GPPstackresamp/100 #from gC per m2 to MgC per ha


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
  coeffmat[1,] <- c('model','pearsonsr','linCCC','slope','intercept','npix')
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
    
    GPPdryclassdf <- data.frame(PML=PMLGPPyearmeans$value,DGVM=DGVMGPPyearmeans$value)
    
    df <- data.frame(x = GPPdryclassdf$PML, y = GPPdryclassdf$DGVM)
    
    DGVMname <- modelnames[modelindex]
    
    meanscatterplotlist[[j]] <- ggplot(df) +
      geom_point(aes(x, y),size=0.5,alpha=0.1) +
      scale_color_identity() +
      theme_bw() +
      theme_classic() +
      theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
      labs(title=paste(titlelist[titlenr],DGVMname),x="MODIS GPP",y=bquote("modelled GPP"))+
      coord_fixed()+
      ylim(c(0,40))+
      xlim(c(0,40))
    titlenr <- titlenr+1

    PMLvals <- PMLGPPyearmeans$value
    DGVMvals <- DGVMGPPyearmeans$value
    
    
    
    #calculate statistics and add regression line to plot
    
    PMLvals <- PMLGPPyearmeans$value
    DGVMvals <- DGVMGPPyearmeans$value
    
    PMLweights <- PMLGPPyearmeans$coverage_fraction
    
    demingreg <- deming(DGVMvals~PMLvals,weights=PMLweights)
    
    #add diagonal and regression line
    meanscatterplotlist[[j]] <- meanscatterplotlist[[j]]+geom_abline(slope=1,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])
    
    pearsonsr <- wtd.cor(PMLvals,DGVMvals,weight=PMLweights)[1,1]#extract person's r from matrix
    
    coeffmat[rcount,1] <- modelnames[j]
    coeffmat[rcount,2] <- pearsonsr
    coeffmat[rcount,3] <- NA#to be replaced with weighted CCC
    coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
    coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
    coeffmat[rcount,6] <- sum(is.finite(PMLvals*DGVMvals))#number of pixels
    rcount <- rcount+1
    }

  #add TRENDY mean comparison and stats
  
  
  df <- data.frame(x = carbonyearmeansdf$PML, y = carbonyearmeansdf$DGVM)#,

  
  DGVMname <- 'TRENDY'
  
  meanscatterplotlist[[j+1]] <- ggplot(df) +
    geom_point(aes(x, y),size=0.5,alpha=0.1) +
    scale_color_identity() +
    theme_bw() +
    theme_classic() +
    theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
    labs(title=paste(titlelist[titlenr],DGVMname),x="MODIS GPP",y=bquote("modelled GPP"))+
    coord_fixed()+
    ylim(c(0,40))+
    xlim(c(0,40))

  
  PMLvals <- carbonyearmeansdf$PML
  DGVMvals <- carbonyearmeansdf$DGVM
  
  PMLweights <- PMLGPPyearmeans$coverage_fraction
  
  demingreg <- deming(DGVMvals~PMLvals,weights=PMLweights)
  
  #add diagonal and regression line
  meanscatterplotlist[[j+1]] <- meanscatterplotlist[[j+1]]+geom_abline(slope=1,intercept=0,linetype='dashed')+geom_abline(slope=demingreg$coefficients[2],intercept=demingreg$coefficients[1])
  
  pearsonsr <- wtd.cor(PMLvals,DGVMvals,weight=PMLweights)[1,1]#extract person's r from matrix
  
  coeffmat[rcount,1] <- 'TRENDY'
  coeffmat[rcount,2] <- pearsonsr
  coeffmat[rcount,3] <- NA#to be replaced with weighted CCC
  coeffmat[rcount,4] <- demingreg$coefficients[2]#slope
  coeffmat[rcount,5] <- demingreg$coefficients[1]#intercept
  coeffmat[rcount,6] <- sum(is.finite(PMLvals*DGVMvals))#number of pixels
  rcount <- rcount+1
 
  #arrange and plot all scatterplots
grid.arrange(grobs=meanscatterplotlist,nrow=3,ncol=5)

#save stats table
write.csv(coeffmat,'./stats/DGVMvsMODIS_GPPTrendStatsWeightedV1.csv')





