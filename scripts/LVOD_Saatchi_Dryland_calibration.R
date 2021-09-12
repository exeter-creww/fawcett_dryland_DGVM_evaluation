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
library("rhdf5")

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

#use 2015 as reliable calibration year (see Fan et al. 2019)
VOD2015 <- VODannualstacktot[[6]]

#Saatchi biomass dataset

h5ls("./Saatchi/AGB_Saatchi_2015.mat")
Saatchi <- h5read("./Saatchi/AGB_Saatchi_2015.mat",'AGB_Saatchi')
Saatchiraster <- raster(Saatchi)
projection(Saatchiraster) <- projection(VOD2015)


extent(Saatchiraster) <- c(-180, 180, -90, 90)
projection(Saatchiraster) <- CRS("+init=epsg:4326")


SaatchiResamp <- raster::resample(Saatchiraster,VOD2015)
SaatchiCResamp <- SaatchiResamp*0.47

#resample dryland classes to 1 deg (modal aggregation then nearest neighbour resampling)
drylandclassag <- raster::aggregate(drylandclass,fact=5,fun=modal)
drylandclassresamp <- raster::resample(drylandclassag,VOD2015,method='ngb')

drylandmask <- drylandclassresamp>1
drylandmask[drylandmask==0] <- NA

#not masking desert areas
drylandmask <- drylandclassresamp>0
drylandmask[drylandmask==0] <- NA

#drylandmasknodesert <- 

VOD2015masked <- VOD2015*drylandmask

#step currently not in calibration
#VOD2015masked <- raster::mask(VOD2015masked,studycontshapes)
#SaatchiCResamp  <- raster::mask(SaatchiCResamp,studycontshapes)
plot(VOD2015masked,col=viridis(100))
plot(drylands4contsub,lwd=1,add=T)

df <- data.frame(x = getValues(VOD2015masked), y = getValues(SaatchiCResamp),
                 d = densCols(getValues(VOD2015masked), getValues(SaatchiCResamp),nbin=200, colramp = colorRampPalette(rev(c('yellow','orange','turquoise4','dodgerblue4')))))#colorRampPalette(rev(rainbow(10, end = 4/6)))))
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
  ylab(bquote("Saatchi AGC" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))+
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
  geom_bin2d(binwidth=c(0.02,4)) +
  scale_fill_continuous(type = "viridis",trans='log10',name='log10(density)') +
  geom_abline(intercept= 0, slope=linmodzerointerc$coefficients[1],lwd=1.5)+
  annotate(
    "text", label = paste0("R = ",round(pearsonsr,2)),
    x = 0, y = 140, size = 5, colour = "black",hjust='left'
  )+
  annotate(
    "text", label = paste0("AGC = ",round(linmodzerointerc$coefficients[1],2)," * LVOD"),
    x = 0, y = 130, size = 5, colour = "black",hjust='left'
  )+
  xlab(bquote("LVOD"))+
  ylab(bquote("Saatchi AGC" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))+
  theme_bw()+
  ylim(c(-0.1,150))+
  xlim(c(0,1))



summary(linmodzerointerc)

