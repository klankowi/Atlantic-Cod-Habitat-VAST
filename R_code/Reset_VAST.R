rm(list=ls())

remove.packages('VAST')

rstudioapi::restartSession()

working_dir <- 'C:/Users/klankowicz/AppData/Local/R/win-library/4.4/VAST'
# Create working directory if it doesn't exist
if(!dir.exists(working_dir)) {
  dir.create(working_dir)
}

setwd(working_dir)

library(devtools)

install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")

# Load package
library(VAST)

# load data set
# see `?load_example` for list of stocks with example data 
# that are installed automatically with `FishStatsUtils`. 
example = load_example( data_set="EBS_pollock" )

setwd('C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/VAST_runs/test')

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 100, 
                          Region = example$Region, 
                          purpose = "index2", 
                          bias.correct = FALSE )

# Run model
fit = fit_model( settings = settings, 
                 Lat_i = example$sampling_data[,'Lat'], 
                 Lon_i = example$sampling_data[,'Lon'], 
                 t_i = example$sampling_data[,'Year'], 
                 b_i = example$sampling_data[,'Catch_KG'], 
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 working_dir = 'C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/VAST_runs/test')

# Plot results
#plot( fit )