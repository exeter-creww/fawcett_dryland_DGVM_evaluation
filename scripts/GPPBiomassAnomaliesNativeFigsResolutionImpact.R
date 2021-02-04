library(mblm)
library(viridis)
library(rgdal)
library(Kendall) 
library(sf)
library(trend)
library(ggplot2)
library(gridExtra)
library(pals)
library(scales)
library(rgeos)
library(reshape2)
library(ncdf4)

setwd('D:/Driving_C')

#Compares the outputs from using different weighting for extracting cVeg and GPP from native model resolutions (done in scripts for individual models)

#standard: weighted by pixel area covered
#contained: only consider pixels whose centre is contained in mask
#extended: considers all pixels touched by mask

#read GPP model native resolutions, different exctration methods
modelGPPregular<- data.frame(GPPregular=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018.csv"))[2:13])
modelGPPcontained<- data.frame(GPPcontained=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018_contained.csv"))[2:13])
modelGPPextended<- data.frame(GPPextended=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018_extended.csv"))[2:13])

modelGPP <- data.frame(modelGPPregular,modelGPPcontained,modelGPPextended,modelnames=row.names(modelGPPregular))

melted <- melt(modelGPP,id.vars='modelnames')

p<-ggplot(data=melted,aes(x=variable, y=value))+#, aes(x=dose, y=len)) +
  geom_bar(aes(), color='black',stat="identity")+
labs(y=bquote("GPP " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='method')+
  scale_x_discrete(labels=c("GPPregular" = "exact", "GPPcontained" = "contained",
                            "GPPextended" = "extended"))+
 
  facet_wrap(~ modelnames)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90))
  

plot(p)


#read cVeg model native resolutions, different exctration methods
modelcVegregular<- data.frame(cVegregular=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018.csv"))[2:13])
modelcVegcontained<- data.frame(cVegcontained=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018_contained.csv"))[2:13])
modelcVegextended<- data.frame(cVegextended=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018_extended.csv"))[2:13])

modelcVeg <- data.frame(modelcVegregular,modelcVegcontained,modelcVegextended,modelnames=row.names(modelcVegregular))

melted <- melt(modelcVeg,id.vars='modelnames')

p<-ggplot(data=melted,aes(x=variable, y=value))+#, aes(x=dose, y=len)) +
  geom_bar(aes(), color='black',stat="identity")+
  labs(y=bquote("AGC " ~ "["~ Pg ~ C ~ "]"),x='method')+
  scale_x_discrete(labels=c("cVegregular" = "exact", "cVegcontained" = "contained",
                            "cVegextended" = "extended"))+
  
  facet_wrap(~ modelnames)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90))


plot(p)


resvec <- c(1,2.8125,1.25,0.5,0.5,1,2,1.875,0.5,1.2,0.5,2)

plot(resvec,modelcVeg$cVegregular/modelcVeg$cVegextended,xlab='resolution',ylab='exact/extended')
plot(resvec,modelcVeg$cVegcontained/modelcVeg$cVegregular,xlab='resolution',ylab='contained/exact')

