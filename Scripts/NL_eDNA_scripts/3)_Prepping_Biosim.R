###########################################################################
#------Script for prepping Biosim climate model data----------

###########################################################################


# ----Prepping Biosim data for historical comparisons-------
library(here)
biosim20 <- read.csv(here("BioSim","ArcticChar2000_2020.csv"))
biosim21_22 <- read.csv(here("BioSim","ArcticChar2021_2022.csv"))
biosim60 <- read.csv(here("BioSim","climate1960_1980.csv"))
biosim_full <- rbind(biosim20,biosim21_22,biosim60)
# Biosim data is currently in daily format, need to convert to yearly averages
# Calculate averages for each year in a 20 yr period for each time point (1960-1980, and 2000-2020)
biosim_full$KeyID <- as.factor(biosim_full$KeyID)
# Calculate 20 yr average ending at 2021:
library(dplyr)
biosim_20year_avgs <- biosim_full %>% filter(Year==2001 |
                                               Year==2002 | Year==2003 |
                                               Year==2004 | Year==2005 | Year==2006 | Year==2007 |
                                               Year==2008 | Year==2009 | Year==2010 | Year==2011 |
                                               Year==2012 | Year==2013 | Year==2014 | Year==2015 |
                                               Year==2016 | Year==2017 | Year==2018 | Year==2019 |
                                               Year==2020) %>% 
  group_by(KeyID, Latitude, Longitude, Elevation,Year) %>%
  summarise(mean_air_temp=mean(Air.Temperature), #mean overall temperature over each year
            max_air_temp=max(Maximum.Air.Temperature), # max air temp in each year
            min_air_temp=min(Minimum.Air.Temperature), # min air temp in each year
            mean_tot_precip=mean(Total.Precipitation), # Mean overall precip each year
            temp_seasonality=100*sd(Air.Temperature),# sd of temp over each year
            temp_range=(max(Maximum.Air.Temperature)-min(Minimum.Air.Temperature)), # temp range each year
            precip_seasonality=100*sd(Total.Precipitation), # sd of precip over each year
            mean_rel_humid=mean(Relative.Humidity), # mean rel humididty over each year
            mean_wind_speed_2m=mean(Wind.Speed.at.2.meters), # mean wind speed each year
            mean_solar_rad=mean(Solar.Radiation), #mean solar rad each year
            mean_atm_press=mean(Atmospheric.Pressure), #mean atm pressure each year
            mean_snow_precip=mean(Snow.Precipitation), # mean snow precip each year
            mean_snow_accum=mean(Snow.Depth.Accumulation)) %>% 
  group_by(KeyID, Latitude,Longitude,Elevation) %>%
  summarise(mean_air_temp=mean(mean_air_temp), #mean overall temps over 20 years
            mean_max_air_temp=mean(max_air_temp),#mean of yearly max air temps over 20 years
            mean_min_air_temp=mean(min_air_temp), #mean of yearly min air temps over 20 years
            mean_tot_precip=mean(mean_tot_precip), #mean overall precip over 20 years
            mean_temp_seasonality=mean(temp_seasonality), #mean of yearly temp seasonalities over 20 years
            mean_temp_range=mean(temp_range), #mean of yearly temp ranges over 20 years
            mean_precip_seasonality=mean(precip_seasonality), #mean of yearly precip seasonalities over 20 years
            mean_rel_humid=mean(mean_rel_humid),
            mean_wind_speed_2m=mean(mean_wind_speed_2m),
            mean_solar_rad=mean(mean_solar_rad),
            mean_atm_press=mean(mean_atm_press),
            mean_snow_precip=mean(mean_snow_precip),
            mean_snow_accum=mean(mean_snow_accum)) %>%
  mutate(year="2001-2020") %>%
  as.data.frame()

# Calculate 20 yr average ending at 1980:
biosim_80year_avgs <- biosim_full %>% filter(Year==1961 | Year==1962 |
                                               Year==1963 | Year==1964 | Year==1965 | Year==1966 |
                                               Year==1967 | Year==1968 | Year==1969 | Year==1970 |
                                               Year==1971 | Year==1972 | Year==1973 | Year==1974 |
                                               Year==1975 | Year==1976 | Year==1977 | Year==1978 |
                                               Year==1979 | Year==1980) %>% 
  group_by(KeyID, Latitude, Longitude, Elevation,Year) %>%
  summarise(mean_air_temp=mean(Air.Temperature), #mean overall temperature over each year
            max_air_temp=max(Maximum.Air.Temperature), # max air temp in each year
            min_air_temp=min(Minimum.Air.Temperature), # min air temp in each year
            mean_tot_precip=mean(Total.Precipitation), # Mean overall precip each year
            temp_seasonality=100*sd(Air.Temperature),# sd of temp over each year
            temp_range=(max(Maximum.Air.Temperature)-min(Minimum.Air.Temperature)), # temp range each year
            precip_seasonality=100*sd(Total.Precipitation), # sd of precip over each year
            mean_rel_humid=mean(Relative.Humidity), # mean rel humididty over each year
            mean_wind_speed_2m=mean(Wind.Speed.at.2.meters), # mean wind speed each year
            mean_solar_rad=mean(Solar.Radiation), #mean solar rad each year
            mean_atm_press=mean(Atmospheric.Pressure), #mean atm pressure each year
            mean_snow_precip=mean(Snow.Precipitation), # mean snow precip each year
            mean_snow_accum=mean(Snow.Depth.Accumulation)) %>% 
  group_by(KeyID, Latitude,Longitude,Elevation) %>%
  summarise(mean_air_temp=mean(mean_air_temp), #mean overall temps over 20 years
            mean_max_air_temp=mean(max_air_temp),#mean of yearly max air temps over 20 years
            mean_min_air_temp=mean(min_air_temp), #mean of yearly min air temps over 20 years
            mean_tot_precip=mean(mean_tot_precip), #mean overall precip over 20 years
            mean_temp_seasonality=mean(temp_seasonality), #mean of yearly temp seasonalities over 20 years
            mean_temp_range=mean(temp_range), #mean of yearly temp ranges over 20 years
            mean_precip_seasonality=mean(precip_seasonality), #mean of yearly precip seasonalities over 20 years
            mean_rel_humid=mean(mean_rel_humid),
            mean_wind_speed_2m=mean(mean_wind_speed_2m),
            mean_solar_rad=mean(mean_solar_rad),
            mean_atm_press=mean(mean_atm_press),
            mean_snow_precip=mean(mean_snow_precip),
            mean_snow_accum=mean(mean_snow_accum)) %>%
  mutate(year="1961-1980") %>%
  as.data.frame()

# ----b) Getting uncorrelated enviro variables----
# Code from Cam Nugent (https://github.com/CNuge/salmon-genomic-vulnerability/blob/main/scripts/03_pair_data_for_rda.r)
# load necessary functions
# load the necessary libraries
library(tidyverse)

# load necessary functions
remove_zero_variance_columns <- function(dat) {
  out = lapply(dat, function(x) length(unique(x)))
  z_var = names(which(!out > 1))
  return(dat[!names(dat) %in% z_var])
}

#' Take a data frame with ONLY a set of predictor variables
#' Create a new data frame with only a set of variables for which the 
#' correlation coefficients do not exceede the specified threshold (default is 0.7).
get_uncorrelated_predictors = function(env_vals, 
                                       threshold = 0.75, drop_zero_variance = TRUE){
  if(drop_zero_variance == TRUE){
    env_vals = remove_zero_variance_columns(env_vals)}
  tmp = cor(env_vals)
  tmp[upper.tri(tmp)] = 0
  diag(tmp) = 0
  env_uncorr = env_vals[, !apply(tmp, 2, function(x) any(abs(x) > 0.75, 
                                                         na.rm = TRUE))]
  return(env_uncorr)
}

all_climate <- rbind(biosim_20year_avgs,biosim_80year_avgs)
env_uncorr_all <- get_uncorrelated_predictors(all_climate[,5:17])
# This returns 7 uncorrelated variables- check to see which ones make sense biologically

# -----------Combine enviro variables back location data-------
#Bind back together:
all_climate2 <- cbind(all_climate[,c(1,18)],env_uncorr_all)
colnames(all_climate2)[1] <- "RiverCode"
colnames(all_climate2)[2] <- "source" # This refers to data source of corresponding ecological data (eDNA vs. historical)
all_climate2 <- all_climate2 %>% 
  mutate(source=case_when(source=="2001-2020"~"eDNA",
                          source=="1961-1980"~"historical"))


# RiverCode CRK is listed as CKR in climate data- switch first to match properly
all_climate2$RiverCode <- as.character(all_climate2$RiverCode)
all_climate2 <- all_climate2 %>%
  mutate(RiverCode=case_when(RiverCode=="CKR" ~ "CRK",
                             TRUE~ RiverCode))
write.csv(all_climate2,"climate_avgs.csv",row.names=FALSE)