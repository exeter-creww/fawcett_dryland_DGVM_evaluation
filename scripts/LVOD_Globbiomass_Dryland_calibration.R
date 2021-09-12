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
library(raster)
library(rmatio)
 

setwd('D:/Driving_C')

#continent outlines for plotting and region subsetting
continentshapes <- readOGR(dsn = getwd(), layer = "WorldContinents")
NorthAmericaShape <- readOGR(dsn = getwd(), layer = "NorthAmericaNoGreenland")
contsfordisp <- raster::aggregate(continentshapes,dissolve=T)

#extract four studied continents: North and South America, Africa and Australia
contnrlist <- c(1,6,3,4)#number in continent shapefiles
studycontshapes <- raster::aggregate(continentshapes[contnrlist,],dissolve=T)

drylands4contsub= readOGR(dsn="./Drylandboundaries",layer="drylands4contsub")

#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclass <- raster("./Plots/drylandclassclipfin.tif")
 
drylandclass <- raster::mask(drylandclass,drylands4contsub)
#preprocessing of VOD data to median annual composites can be found in LVODprocessingAnnualStackFin.R

VODannualstacktot <- stack("./LVOD_WGS84/composites/median/VOD_ASC_annual_median_filtv2.tif")

#use 2011 as reliable calibration year (see Fan et al. 2019)
VOD2011 <- VODannualstacktot[[2]]

#Globbiomass dataset

GlobBiom <- read.mat("./GlobBiomass/GloBiomass.mat")
GlobBiomraster <- raster(GlobBiom$GloBiomass)
projection(GlobBiomraster) <- projection(VOD2011)


extent(GlobBiomraster) <- c(-180, 180, -90, 90)
projection(GlobBiomraster) <- CRS("+init=epsg:4326")


GlobBiomResamp <- raster::resample(GlobBiomraster,VOD2011)
GlobBiomCResamp <- GlobBiomResamp*0.47

#resample dryland classes to 1 deg (modal aggregation then nearest neighbour resampling)
drylandclassag <- raster::aggregate(drylandclass,fact=5,fun=modal)
drylandclassresamp <- raster::resample(drylandclassag,VOD2011,method='ngb')

drylandmask <- drylandclassresamp>0
drylandmask[drylandmask==0] <- NA

#drylandmasknodesert <- 

VOD2011masked <- VOD2011*drylandmask

#step currently not in calibration
#VOD2011masked <- raster::mask(VOD2011masked,studycontshapes)
#GlobBiomCResamp  <- raster::mask(GlobBiomCResamp,studycontshapes)
plot(VOD2011masked,col=viridis(100))
plot(drylands4contsub,lwd=1,add=T)

df <- data.frame(x = getValues(VOD2011masked), y = getValues(GlobBiomCResamp),
                 d = densCols(getValues(VOD2011masked), getValues(GlobBiomCResamp),nbin=200, colramp = colorRampPalette(rev(c('yellow','orange','turquoise4','dodgerblue4')))))#colorRampPalette(rev(rainbow(10, end = 4/6)))))
pearsonsr <- cor(df$x,df$y,'complete.obs')
linmodzerointerc <- lm(y~0+x,df)
p <- ggplot(df) +
  geom_point(aes(x, y,col='blue'),size=0.5) +
  scale_color_identity() +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=18))+
  geom_abline(intercept= 0, slope=linmodzerointerc$coefficients[1],lwd=1.5)+
  #geom_abline(intercept = 0, slope = 116.56,lty=2) + #compare old Brandt et al based model
  xlab(bquote("LVOD"))+
  ylab(bquote("GlobBiom AGC" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))+
  annotate(
    "text", label = paste0("R = ",round(pearsonsr,2)),
    x = 0, y = 140, size = 5, colour = "black",hjust='left'
  )+
  annotate(
    "text", label = paste0("AGC = ",round(linmodzerointerc$coefficients[1],2)," * LVOD"),
    x = 0, y = 130, size = 5, colour = "black",hjust='left'
  )+
  theme_bw()+
  ylim(c(0,150))+
  xlim(c(0,1))

print(p)

#density plot
ggplot(df, aes(x=x, y=y,fill = ..density..) ) +
  geom_bin2d(binwidth=c(0.02,3)) +
  scale_fill_continuous(type = "viridis",trans='log10',name='log10(density)') +
  geom_abline(intercept= 0, slope=linmodzerointerc$coefficients[1],lwd=1.5)+
  annotate(
    "text", label = paste0("R = ",round(pearsonsr,2)),
    x = 0.05, y = 120, size = 4, colour = "black",hjust='left'
  )+
  annotate(
    "text", label = paste0("AGC = ",round(linmodzerointerc$coefficients[1],2)," * LVOD"),
    x = 0.05, y = 110, size = 4, colour = "black",hjust='left'
  )+
  xlab(bquote("LVOD"))+
  ylab(bquote("GlobBiom AGC" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))+
  theme_bw()+
  #ylim(c(-10,150))+
  #xlim(c(-0.05,1))+
  scale_y_continuous(expand = c(0, 0))+#,limits=c(-0.5,150))+
  scale_x_continuous(expand = c(0, 0))#,limits=c(-0.01,1))

summary(linmodzerointerc)

