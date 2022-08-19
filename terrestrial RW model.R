library(tidyverse)
library(tsibble)
library(fable)

#read in the targets data
targets <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_daily/terrestrial_daily-targets.csv.gz", guess_max = 1e6)

# Make the targets into a tsibble with explicit gaps
targets_ts <- targets %>%
  as_tsibble(key = c('variable', 'site_id'), index = 'time') %>%
  # add NA values up to today (index)
  fill_gaps(.end = Sys.Date())

# Read in the site_data information
# site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
#   dplyr::filter(terrestrial == 1)

# # For the RW model need to start the forecast at the last non-NA day and then run to 35 days in the future
# forecast_starts <- targets %>%
#   filter(!is.na(observed)) %>%
#   group_by(site_id, variable) %>%
#   # Start the day after the most recent non-NA value
#   summarise(start_date = max(time) + 1) %>% # Date 
#   mutate(h = (Sys.Date() - start_date) + 35) %>% # Horizon value 
#   ungroup()
# 
# ggplot(targets, aes(x=time, y=observed)) +
#   geom_point() +
#   facet_grid(variable~site_id, scales = 'free')

# Function carry out a random walk forecast
RW_forecast <- function(site, var, h,
                        bootstrap = FALSE, boot_number = 200, 
                        transformation = 'none', verbose = TRUE,...) {
  # Work out when the forecast should start
  forecast_starts <- targets %>%
    filter(!is.na(observed) & site_id == site & variable == var) %>%
    # Start the day after the most recent non-NA value
    summarise(start_date = max(time) + 1) %>% # Date
    mutate(h = (Sys.Date() - start_date) + h) %>% # Horizon value
    ungroup()
  
  if (verbose == T) {
    message(
      site,
      ' ',
      var,
      ' forecast with transformation = ',
      transformation,
      ' and bootstrap = ',
      bootstrap
    )
  }
  
  # filter the targets data set to the site_var pair
  targets_use <- targets %>%
    filter(site_id == site,
           variable == var) %>%
    as_tsibble(key = c('variable', 'site_id'), index = 'time') %>%
    # add NA values up to today (index)
    fill_gaps(.end = Sys.Date()) %>%
    # Remove the NA's put at the end, so that the forecast starts from the last day with an observation,
    # rather than today
    filter(time < forecast_starts$start_date)
  
  if (transformation == 'log') {
    RW_model <- targets_use %>%
      model(RW = RW(log(observed)))
  }
  if (transformation == 'log1p') {
    RW_model <- targets_use %>%
      model(RW = RW(log1p(observed)))
  }
  if (transformation == 'none') {
    RW_model <- targets_use %>%
      model(RW = RW(observed))
  }
  
  if (bootstrap == T) {
    forecast <-
      RW_model %>% generate(
        h = as.numeric(forecast_starts$h),
        bootstrap = T,
        times = boot_number
      )
  }  else
    forecast <-
    RW_model %>% forecast(h = as.numeric(forecast_starts$h))
  return(forecast)
  
}

# 1. Run through each via a loop
# Generate a table with each combination of the site and variable
# site_var_combinations <- expand.grid(site = unique(targets$site_id),
#                                      var = unique(targets$variable)) %>%
#   # assign the transformation depending on the variable
#   mutate(transformation = ifelse(var == 'le', 'log', 'none'))
# 
# 
# for (i in 1:nrow(site_var_combinations)) {
# 
#   forecast <- RW_forecast(site_var_combinations$site[i], 
#               site_var_combinations$var[i], 
#               h  = 35, 
#               bootstrap = T, 
#               boot_number = 200, 
#               verbose = T, 
#               transformation =  site_var_combinations$transformation[i])
#   
#   if (!exists("RW_forecasts")) {
#     RW_forecasts <- forecast
#   } else {
#     RW_forecasts <- bind_rows(forecast, RW_forecasts)
#   }
#     
# }
# 

# 2. Run through each via map
site_var_combinations <- expand.grid(site = unique(targets$site_id),
                                     var = unique(targets$variable)) %>%
  # assign the transformation depending on the variable
  mutate(transformation = ifelse(var == 'le', 'log', 'none')) %>%
  mutate(boot_number = 200,
         h = 35,
         bootstrap = T, 
         verbose = T)

# runs the RW forecast for each combination of variable and site_id
RW_forecasts <- purrr::pmap_dfr(site_var_combinations, RW_forecast) 

# convert the output into EFI standard
RW_forecasts_EFI <- RW_forecasts %>%
  rename(ensemble = .rep,
         predicted = .sim) %>%
  group_by(site_id, variable) %>%
  mutate(start_time = min(time) - 1) %>%
  select(time, start_time, site_id, ensemble, variable, predicted)
