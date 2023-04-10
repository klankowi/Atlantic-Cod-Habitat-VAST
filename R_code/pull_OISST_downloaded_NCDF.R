# Adapted from code written by Kimberly Bastille and Sarah Gaichas
# Theoretically should be able to download NetCDF files directly into R from
# NOAA, but recently the files have downloaded as corrupt and not openable.
# Manual downloading of yearly NetCDF files from 
# (https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html) is slower
# but more reliable. This code requires those NetCDF files to already be 
# manually downloaded from NOAA's website and placed into the project directory.

# Clear workspace
rm(list=ls())

# Load libraries
library(ncdf4)
library(raster)
library(tidyverse)
library(here)

# Load functions
source(here("R_code/utilities/nc_to_raster_func.R"))
source(here("R_code/utilities/raster_to_sstdf_func.R"))

# Set variable name to sea surface temperature
varname <- 'sst'

# Set period of years to work with
years <- 1982:2022

# Loop through years
for(i in years) {
  # Set file name of .nc file
  name <- paste0("sst.day.mean.", i, ".nc")
  # Set file name of output .gri, .grd files
  filename <- here::here("Data","Density_Covariates", "SST", 
                         paste0("test_", i, ".grd"))
  # Convert specified NetCDF file to raster
  text <- knitr::knit_expand(text = "test_{{year}} <- nc_to_raster(nc = name, varname = varname)
                                     raster::writeRaster(test_{{year}}, filename = filename, overwrite=TRUE)",
                             year = i)
  try(eval(parse(text = text)))
  # Remove nc file to save space (will be deleted from directory)
  unlink(name)
  # Set file name of raster (.gri, .grd files)
  name <- get(paste0("test_",i))
  # Set file name of output .rds files
  filename <- here::here("Data","Density_Covariates", "SST", 
                         paste0("sst", i, ".rds"))
  # Save rasters as .rds files
  text <- knitr::knit_expand(text = "sst{{year}} <- raster_to_sstdf(brick = name)
                                     saveRDS(sst{{year}}, filename)",
                             year = i)
  try(eval(parse(text = text)))
# Terminate outer loop
}