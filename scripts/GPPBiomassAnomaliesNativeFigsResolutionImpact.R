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

setwd('D:/Driving_C')

#Compares the outputs from using different weighting for extracting cVeg and GPP from native model resolutions (done in scripts for individual models)

#standard: weighted by pixel area covered
#contained: only consider pixels whose centre is contained in mask
#extended: considers all pixels touched by mask

#read GPP model native resolutions, different exctration methods
modelGPPregular<- data.frame(GPPregular=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018.csv"))[2:13])
modelGPPcontained<- data.frame(GPPcontained=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018_contained.csv"))[2:13])
modelGPPextended<- data.frame(GPPextended=colMeans(read.csv("./DGVM/DGVMdrylandTS/GPP/GPP_drylands_2003_2018_extended.csv"))[2:13])

modelGPP <- data.frame(modelGPPcontained,modelGPPregular,modelGPPextended,modelnames=row.names(modelGPPregular))

melted <- melt(modelGPP,id.vars='modelnames')

p<-ggplot(data=melted,aes(x=variable, y=value))+#, aes(x=dose, y=len)) +
  geom_bar(aes(), color='black',stat="identity")+
labs(y=bquote("GPP " ~ "["~ Pg ~ C ~ yr^{-1}~ "]"),x='method')+
  scale_x_discrete(labels=c( "GPPcontained" = "contained","GPPregular" = "exact",
                            "GPPextended" = "extended"))+
 
  facet_wrap(~ modelnames)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90))
  

plot(p)


#read cVeg model native resolutions, different exctration methods
modelcVegregular<- data.frame(cVegregular=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018.csv"))[2:13])
modelcVegcontained<- data.frame(cVegcontained=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018_contained.csv"))[2:13])
modelcVegextended<- data.frame(cVegextended=colMeans(read.csv("./DGVM/DGVMdrylandTS/cVeg/cVeg_drylands_2011_2018_extended.csv"))[2:13])

modelcVeg <- data.frame(modelcVegcontained,modelcVegregular,modelcVegextended,modelnames=row.names(modelcVegregular))

melted <- melt(modelcVeg,id.vars='modelnames')

p<-ggplot(data=melted,aes(x=variable, y=value))+#, aes(x=dose, y=len)) +
  geom_bar(aes(), color='black',stat="identity")+
  labs(y=bquote("AGC " ~ "["~ Pg ~ C ~ "]"),x='method')+
  scale_x_discrete(labels=c("cVegcontained" = "contained","cVegregular" = "exact", 
                            "cVegextended" = "extended"))+
  
  facet_wrap(~ modelnames)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90))


plot(p)


resvec <- c(1,2.8125,1.25,0.5,0.5,1,2,1.875,0.5,1.2,0.5,2)

#plot(resvec,modelcVeg$cVegregular/modelcVeg$cVegextended,xlab='resolution',ylab='exact/extended')
#plot(resvec,modelcVeg$cVegcontained/modelcVeg$cVegregular,xlab='spatial resolution [degrees]',ylab='contained vs exact ratio')

modelcVegresimpactdf <- data.frame(resvec,ratio=modelcVeg$cVegcontained/modelcVeg$cVegregular)

plot1 <- ggplot(modelcVegresimpactdf,aes(x=resvec, y=ratio)) +
  geom_point() +
  #scale_color_identity() +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  labs(title='a) cVeg',x='spatial resolution [degrees]',y=bquote('contained to exact ratio'))+#,title=paste0('Carbon density: ',DGVMname)) +
  ylim(c(0,1))+
  xlim(c(0,3))


resvec <- c(1,2.8125,1.25,0.5,0.5,1,2,1.875,0.5,1.2,0.5,2)

#plot(resvec,modelGPP$GPPregular/modelGPP$GPPextended,xlab='resolution',ylab='exact/extended')
#plot(resvec,modelGPP$GPPcontained/modelGPP$GPPregular,xlab='resolution',ylab='contained/exact')

modelGPPresimpactdf <- data.frame(resvec,ratio=modelGPP$GPPcontained/modelGPP$GPPregular)

plot2 <- ggplot(modelGPPresimpactdf,aes(x=resvec, y=ratio)) +
  geom_point() +
  #scale_color_identity() +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=12),plot.title = element_text(face="bold",size=12))+
  labs(title='b) GPP',x='spatial resolution [degrees]',y=bquote('contained to exact ratio'))+#,title=paste0('Carbon density: ',DGVMname)) +
  ylim(c(0,1))+
  xlim(c(0,3))

grid.arrange(grobs=list(plot1,plot2),nrow=1,ncol=2)      
