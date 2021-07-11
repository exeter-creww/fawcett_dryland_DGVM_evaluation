#derive LVOD to AGC calibration from Avitabile biomass map

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
library(MASS)
library(epiR)
library(ggplot2)


setwd('D:/Driving_C')

#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")
contsfordisp <- raster::aggregate(continentshapes,dissolve=T)

#extract four studied continents: North and South America, Africa and Australia
contnrlist <- c(1,6,3,4)#number in continent shapefiles
studycontshapes <- raster::aggregate(continentshapes[contnrlist,],dissolve=T)

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- raster("./Plots/drylandclassclipfin.tif")
 

#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#use 2011 as reliable calibration year (see Fan et al. 2019)
VOD2011 <- VODannualstacktot[[2]]

#Avitabile biomass dataset

AvitabileBiom <- raster("./Avitabile_AGB_Map/Avitabile_AGB_Map.tif")

AvitabileBiomResamp <- raster::resample(AvitabileBiom,VOD2011)
AvitabileCResamp <- AvitabileBiomResamp*0.47

#resample dryland classes to 1 deg (modal aggregation then nearest neighbour resampling)
drylandclassag <- raster::aggregate(drylandclass,fact=5,fun=modal)
drylandclassresamp <- raster::resample(drylandclassag,VOD2011,method='ngb')

drylandmask <- drylandclassresamp>0
drylandmask[drylandmask==0] <- NA

VOD2011masked <- VOD2011*drylandmask


df <- data.frame(x = getValues(VOD2011masked), y = getValues(AvitabileCResamp),
                 d = densCols(getValues(VOD2011masked), getValues(AvitabileCResamp),nbin=200, colramp = colorRampPalette(rev(c('yellow','orange','turquoise4','dodgerblue4')))))#colorRampPalette(rev(rainbow(10, end = 4/6)))))

linmodzerointerc <- lm(y~0+x,df)
p <- ggplot(df) +
  geom_point(aes(x, y,col=d),size=0.5) +
  scale_color_identity() +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=18))+
  geom_abline(intercept= 0, slope=linmodzerointerc$coefficients[1],lwd=1.5)+
  #geom_abline(intercept = 0, slope = 116.56,lty=2) + #compare old Brandt et al based model
  xlab(bquote("LVOD"))+
  ylab(bquote("Avitabile AGC" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))+
  ylim(c(0,120))+
  xlim(c(0,1))#+

print(p)

summary(linmodzerointerc)

