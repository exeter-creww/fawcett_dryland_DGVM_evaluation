library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(scales)
library("ncdf4", lib.loc="~/R/win-library/3.4")
library(gridExtra)
library(rasterVis)
library(tls)

 
#Combine individual years into VOD stack, filtering out pixels with large differences in ASC-DESC and pixels with low observation counts

#composite type 
#comptype <- 'median'

ASCfiles <- list.files(paste0("D:/Driving_C/LVOD_WGS84/composites/median/ASC/"))
countfiles <- list.files(paste0("D:/Driving_C/LVOD_WGS84/composites/count/ASC/"))
DESCfiles <- list.files(paste0("D:/Driving_C/LVOD_WGS84/composites/median/DESC/"))


VODstack2 <- stack()

#use ASC values and mask pixels where ASC and DESC significantly different (indicating RFI) or too few observations
for(i in 1:(length(ASCfiles))){
  
  ASCrasterout <- raster(paste0("D:/Driving_C/LVOD_WGS84/composites/median/ASC/",ASCfiles[i]))
  ASCraster <- ASCrasterout
  DESCraster <- raster(paste0("D:/Driving_C/LVOD_WGS84/composites/median/DESC/",DESCfiles[i]))
  countraster <- raster(paste0("D:/Driving_C/LVOD_WGS84/composites/count/ASC/",countfiles[i]))
  
    #filter out pixels with less than 20 observations, replace pixels with 50 observations with long-term mean (Brandt et al., 2018)
  #ASCrasterout[countraster<50] <- VODASCmeans[countraster<50]#approach of Brandt et al. 2018, we decided not to use
  ASCrasterout[countraster<20] <- NA
  
  #calculate difference in mean ASC and DESC values, filter out differences greater than 0.05, indicating RFI (pers. comm. Amen Al-Yaari)
  ASCDESCdifmask <- ((ASCraster-DESCraster)<0.05)&((ASCraster-DESCraster)>-0.05)
  ASCDESCdifmask[ASCDESCdifmask==0] <- NA
  vodmasked <- mask(ASCrasterout,ASCDESCdifmask)
  
  VODstack2 <- stack(VODstack2,vodmasked)
}

writeRaster(VODstack2,"D:/Driving_C/LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif",overwrite=T)

