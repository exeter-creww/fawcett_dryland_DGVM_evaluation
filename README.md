

[![DOI](https://zenodo.org/badge/333717763.svg)](https://zenodo.org/badge/latestdoi/333717763)


# Code for 'Assessing the predictions of carbon dynamics in global drylands'

This repository contains R code for processing, analysing and visualizing the 
data presented in 'Assessing the predictions of carbon dynamics in global 
drylands' by Dominic Fawcett, Andrew Cunliffe, Karen Anderson, Richard Brazier, 
Tim Hill, Stephen Sitch, Mike O’Sullivan, Jean Pierre Wigneron, and Others.

Contact: df332@exeter.ac.uk or a.cunliffe@exeter.ac.uk

## The manuscript is In Preparation


INCLUDE DATA REFERENCE
This code analyses the following datasets:
1) VOD dataset supplied by Jean Pierre Wigneron (REFERENCE)
2) Output from the Trendy land surface model intercomparison (v8; https://sites.exeter.ac.uk/trendy/).
3) Dryland classes imported from Google Earth Engine Repository: https://earthengine.googlesource.com/users/dfawcett/DRIVING_C/+/refs/heads/master


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
