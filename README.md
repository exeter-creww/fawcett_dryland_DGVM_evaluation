# fawcett_dryland_DGVM_evaluation

R scripts used for final result generation:

**LVODprocessingASCv2.R:** Filter and aggregated daily netCDF datasets of LVOD, generating annual median and count rasters.
**LVODprocessingDESCv2.R:** Same as above but for DESC datasets (used for masking ASC DESC differences)
**LVODprocessingAnnualStackFin.R:** Filters out unreliable pixels based on observation count and ASC-DESC differences, creates filtered annual stack

**createTRENDYGPPbricksv3.R:** Extracts GPP data for 2003-2018 from netCDF and writes raster bricks with converted GPP values (kg/m2/y) for further analysis.
**createTRENDYcVegbricksv3.R:** Extracts vegetation carbon data for 2011-2018 from netCDF and writes raster bricks with veg. carbon values (kg/m2, above and below ground!) for further analysis.

**LVODbiomAnnualTrendyMeanFigsDensPlotsWeights.R:** Generate scatterplot and linear model figures for TRENDY models and LVOD biomass carbon comparison

**GPPtrendsVsDGVMFigsDensPlotsWeights.R:** Generate scatterplot figures for TRENDY models and PMLv2 GPP comparison

**[MODELNAME]_GPPBiomassAnomalies.R:** Extracts GPP and Biomass anomalies across time series (from 1901 and 2003/2011) from model at native resolution

**GPPBiomassAnomaliesNativeFigs.R:** Generate plots of GPP and Biomass anomalies across time series (from 1901 and 2003/2011)

**GPPBiomassAnomaliesNativeFigsResolutionImpact.R:** Investigate differences between extraction methods (bar plots)
