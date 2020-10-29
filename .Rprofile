# WaveLightGLS .Rprofile
# loading project environment
# A. Roy
#----------------------------------------------------

# clean environnement
rm(list = ls())

# # loading libraries
library(GeoLight)
library(SGAT)
library(TwGeos)
library(tidyverse)
library(lubridate)
library(maptools)
library(rgdal)
library(MASS)
library(WaveletComp)
library(viridis)
library(ncdf4)
library(plot3D)
library(gridExtra)

# load home-made functions
a = lapply(list.files(path= "toolbox", pattern = "[.]R$", 
                      full.names = TRUE, recursive = TRUE), source)
rm(a)

# set seed
set.seed(1)

# Presentation
cat("\014")
cat("----------------------------------------------------\n")
cat("WaveLightGLS - tools for GLS data analysis\n")
cat("A. Roy\n")
cat("----------------------------------------------------\n")
