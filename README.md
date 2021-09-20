

[![DOI](https://zenodo.org/badge/333717763.svg)](https://zenodo.org/badge/latestdoi/333717763)


# Code for 'Assessing the predictions of carbon dynamics in global drylands'

This repository contains R code for processing, analysing and visualizing the 
data presented in 'Assessing model predictions of carbon dynamics in global 
drylands' by Dominic Fawcett, Andrew Cunliffe, Karen Anderson, Richard Brazier, 
Tim Hill, Stephen Sitch, Mike O’Sullivan and Others.

Contact: df332@exeter.ac.uk or a.cunliffe@exeter.ac.uk

## The manuscript is In Preparation


INCLUDE DATA REFERENCE
This code analyses the following datasets:
1) VOD dataset supplied by Jean Pierre Wigneron (SMOS-IC V2)
2) Output from the Trendy land surface model intercomparison (v8; https://sites.exeter.ac.uk/trendy/).
3) PML-v2 MODIS based GPP data (Zhang et al. 2019 https://github.com/gee-hydro/gee_PML)

Google Earth Engine Repository for exporting dryland classes, burned area and PML-v2 data: https://earthengine.googlesource.com/users/dfawcett/DRIVING_C_RS_publication)


R scripts used for final result generation:

**LVODprocessingASCv2.R:** Filter and aggregated daily netCDF datasets of LVOD, generating annual median and count rasters.

**LVODprocessingDESCv2.R:** Same as above but for DESC datasets (used for masking ASC DESC differences)

**LVODprocessingAnnualStackFin.R:** Filters out unreliable pixels based on observation count and ASC-DESC differences, creates filtered annual stack

**LVODprocessingAnnualStackFin_noDESCfilt.R:** Filters out unreliable pixels based on observation count, NO ASC-DESC data filtering, creates filtered annual stack

**createTRENDYGPPbricksv3.R:** Extracts GPP data for 2001-2018 from netCDF and writes raster bricks with converted GPP values (kg/m2/y) for further analysis.

**createTRENDYcVegbricksv3.R:** Extracts vegetation carbon data for 2011-2018 from netCDF and writes raster bricks with veg. carbon values (kg/m2, above and below ground!) for further analysis.

**LVODbiomAnnualTrendyAllModelsHexPlotsWeightsPub.R:** Generate hexbin scatterplot and linear model figures for TRENDY models and LVOD biomass carbon comparison, export spatial statistics and create spatial maps of AGC and GPP distributions and trends

**GPPtrendsVsDGVMAllModelsFigsHexPlotsWeightsPub.R:** Generate hexbin scatterplot and linear model figures for TRENDY models and PML-v2 GPP comparison, export spatial statistics

**LVODbiomAnnualTrendyAllModelsBoxplotsWeightsFireanalysis.R:** Generate boxplots of LVOD and TRENDY AGC difference per burned area bin 

**[MODELNAME]_GPPBiomassAnomalies.R:** Extracts GPP (2001-2018), aboveground carbon, vegetation carbon and soil carbon across (1901 and 2003/2011) time series from model at native resolution

**[MODELNAME]_GPPBiomassAnomaliesSI.R:** Extracts GPP (1901-2018), litter and coarse woody debris (1901-2018) across time series from model at native resolution

**GPPBiomassAnomaliesNativeFigsFinalPub.R:** Generate plots of GPP and Biomass anomalies across time series (from 1901 and 2001/2011), export statistics

**GPPBiomassAnomaliesNativeFigsResolutionImpact.R:** Investigate differences between extraction methods (bar plots)
