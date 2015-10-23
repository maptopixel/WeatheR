# WeatheR package #
## Some R functions for wrangling Weather Underground, Met Office, Environment Agency rainfall data##

The scripts process data from:

* UK MIDAS datasets, as retrieved from NERC BADC Web Processing Service

* Weather Underground API

The functions are:

* loadAndJoinMidas.R - use for taking a shapefile of points and a table (e.g. CSV) of MIDAS data and creating a joined dataset of observations chosen based on MIDAS Quality Control parameters.

* loadAndJoinWU.R - use for taking a directory of CSVs (one CSV per station) of data produced by WU API and returning a data frame of the observations.

* processMidasQC.R - use for taking a table (e.g. CSV) of MIDAS and producing a table of 1 observation per interval (hourly or daily).

* getHistoryForLocationsWU.R - use for downloading and processing (aggregating) Weather Underground data. Data are downloaded (and cached) as JSON files and processed in the same step.

* getLocationsWU.R - creates a spatial dataset of personal weather station points from the WU api

* arealInterpolation.R - computes the weighted mean of polygon attribute values, where weighting is based on area contribution

The demo scripts are:

* demo_processMidasQC.R - use for testing and generating observations for one MIDAS station

* demo_rainfall_interpolate - use for demonstrating the processing of MIDAS and WU observations