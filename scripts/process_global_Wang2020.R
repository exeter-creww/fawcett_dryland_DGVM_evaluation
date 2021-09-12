library(ncdf4)
library(fields)
library(raster)

#source("~/R/func.R")

setwd('D:/Driving_C')

lon = seq(-179.975,179.975,length.out = 7200)
lat = seq(-89.975,89.975,length.out = 3600) 

dimmat <- raster(matrix(1,nrow=3600,ncol=7200))

extent(dimmat) <- c(-180, 180, -90, 90)
projection(dimmat) <- CRS("+init=epsg:4326")

#areas = area.mask.ext(lon,lat)



#dryland classes according to EU JRC dryland definition, also Yao et al. 2020 (exported from GEE)
drylandclassraster <- raster("./Plots/drylandclassclipfin.tif")
drylandclassmask <- drylandclassraster>0
#drylandclassmask[drylandclassmask==0] <- NA

drylandclassmaskresamp <- raster::resample(drylandclassmask,dimmat,method='ngb')

drylandclassmaskresamp[drylandclassmaskresamp==0] <- NA  

dir = "./GPP_datasets/Wang_2020/"

years = seq(16)+2002

dat_full = array(NA,dim=c(16))

for(y in 1:16){
  files = list.files(paste0(dir,years[y]),pattern="nc",full.names=T)
  dat_ann = array(NA,dim=c(length(files)))
  for(f in 1:length(files)){
    print(paste0("year=",years[y],", file=",f,"/12"))
    ncin = nc_open(files[f])
    dat = ncvar_get(ncin,"GPP")
    dat[dat == -9999] = NA
    #dat = rotate_clockwise(dat)
    nc_close(ncin)

    datraster = raster(dat)
    extent(datraster) <- c(-180, 180, -90, 90)
    projection(datraster) <- CRS("+init=epsg:4326")
    
    datrastermasked = datraster*drylandclassmaskresamp
    
    dat_ann[f] = sum(dat*areas,na.rm=T) * 1e-15 * 0.001 # PgC/day (monthly mean)
  }
  dat_full[y] = mean(dat_ann)*365 # PgC/yr
}

GPP_global_Wang = dat_full

save(GPP_global_Wang,file="GPP_global_Wang_82-16.RData")
