
# functions

# 1-1. load and rotate raster (for daily data)
load_and_rotate_raster <- function(file_paths, level = NULL) {
  
  # Helper function to load and rotate a single raster file
  load_single <- function(file) {
    # Load the raster brick
    raster_data <- if (!is.null(level)) {
      brick(file, level = level)  # Load specific level if provided
    } else {
      brick(file)  # Load all levels
    }
    
    # Rotate the raster to fix coordinates
    rotated_data <- rotate(raster_data)
    return(rotated_data)
  }
  
  # Apply the helper function to all file paths
  raster_list <- lapply(file_paths, load_single)
  
  return(raster_list)
}

# 1-2. load data (for monthly data)
load_raster <- function(file_path, level = NULL) {
  # Load the raster brick
  if (!is.null(level)) {
    raster_data <- brick(file_path, level = level)
  } else {
    raster_data <- brick(file_path)
  }
  
  return(raster_data)
}

# 2. convert raster to df 

raster_to_df <- function(raster_data) {
  
 df_list <- lapply(raster_list, function(raster) as.data.frame(raster, xy = TRUE))
 return(df_list)
  
}

# 3. Extract climate data for country 

extract_clim_data <- function(clim_raster, rhum_raster, shape_path, country, start_date, end_date) {
  
  # load shape files
  world <- st_read(shape_path)
  country_sf <- st_as_sf(world[world$sovrgnt == country, ])
  
  # extract temp and rhum 
  temp_extract <- exact_extract(clim_raster, country_sf, include_xy = T)
  temp_df <- do.call(rbind, temp_extract) - 273.15
  
  rh_extract <- exact_extract(rhum_raster, country_sf, include_xy = T)
  rh_df <- do.call(rbind, rh_extract)
  
  # create dataframe
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  n_days <- length(dates)
  
  country_data <- data.frame(
    date = as.character(dates),
    T    = colMeans(temp_df[, 1:n_days], na.rm = T),
    H    = colMeans(rh_df[, 1:n_days], na.rm = T)
  )
  
  return(country_data)
  
}

# 3-1. Extract climate data from specific country using country shape file
extract_clim_data_local <- function(clim_raster, rhum_raster, shape_path, start_date, end_date) {
  
  # load shape files
  country <- st_read(shape_path)
  
  # extract unique adm values
  uniq_adm <- unique(country$id)
  
  # Create date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  n_days <- length(dates)
  
  # Initialize list to store data frames for each adm
  all_adm_data <- list()
  
  # loop
  for(adm in uniq_adm) {
    
    country_sf <- st_as_sf(country[country$id == adm, ]) 
    # extract temp and rhum 
    temp_extract <- exact_extract(clim_raster, country_sf, include_xy = T)
    temp_df <- do.call(rbind, temp_extract) - 273.15
    rh_extract <- exact_extract(rhum_raster, country_sf, include_xy = T)
    rh_df <- do.call(rbind, rh_extract)
    
    # Calculate daily averages for temperature and humidity
    temp_avg <- colMeans(temp_df[, 1:n_days], na.rm = TRUE)
    rh_avg <- colMeans(rh_df[, 1:n_days], na.rm = TRUE)
    
    # Combine into a data frame
    adm_data <- data.frame(
      date = as.character(dates),
      T = temp_avg,
      H = rh_avg,
      adm = adm  # Add the admin unit as a column
    )
    
    # Store the data frame in the list
    all_adm_data[[as.character(adm)]] <- adm_data
  }
  
  return(all_adm_data)
  
}



# 4. Cumulative FoI 

cum_foi <- function(cum_indexP_list, annual_foi) {
  
  foi_list <- list()
  
  for(i in seq_along(cum_indexP_list)){
    year_data <- cum_indexP_list[[i]]
    
    # Normalize and compute FOI using mid, lo, hi values
    cum_foi_mid <- annual_foi * (year_data$mid / max(year_data$mid, na.rm = TRUE))
    cum_foi_lo  <- annual_foi * (year_data$lo / max(year_data$lo, na.rm = TRUE))
    cum_foi_hi  <- annual_foi * (year_data$hi / max(year_data$hi, na.rm = TRUE))
    
    foi_list[[i]] <- data.frame(
      date = year_data$date,
      foi_mid = cum_foi_mid,
      foi_lo  = cum_foi_lo,
      foi_hi  = cum_foi_hi
    )
    
  }
  
  return(foi_list)
}

# 4-1. loop through admin
cum_foi_adm <- function(cum_indexP_list, annual_foi) {
  
  # Initialize an empty list to store results for all years
  foi_list <- lapply(seq_along(cum_indexP_list), function(year_index) {
    year_data_list <- cum_indexP_list[[year_index]]
    
    # Initialize a list to store FOI for each admin unit within the year
    lapply(seq_along(year_data_list), function(admin_index) {
      admin_data <- year_data_list[[admin_index]]
      
      # Normalize and compute FOI using mid, lo, hi values
      cum_foi_mid <- annual_foi * (admin_data$mid / max(admin_data$mid, na.rm = TRUE))
      cum_foi_lo  <- annual_foi * (admin_data$lo / max(admin_data$lo, na.rm = TRUE))
      cum_foi_hi  <- annual_foi * (admin_data$hi / max(admin_data$hi, na.rm = TRUE))
      
      # Create a data frame for the FOI values and return it
      data.frame(
        date = admin_data$date,
        foi_mid = cum_foi_mid,
        foi_lo  = cum_foi_lo,
        foi_hi  = cum_foi_hi,
        admin = admin_data$admin[1]  # Add admin name for reference
      )
    })
  })
  
  # Set names for the outer and inner lists based on the structure of cum_indexP_list
  names(foi_list) <- names(cum_indexP_list)
  lapply(seq_along(foi_list), function(i) {
    names(foi_list[[i]]) <- names(cum_indexP_list[[i]])
  })
  
  return(foi_list)
}

# 5. Daily downscale FoI

downscale_foi <- function(cum_foi_list) {
  
  daily_foi_list <- list()
  
  for(i in seq_along(cum_foi_list)){
    
    cum_foi <- cum_foi_list[[i]] # for each year
    
    # create a dataframe to store daily foi (how each year looks like)
    daily_foi <- data.frame(
      date = cum_foi$date,
      foi_mid = numeric(nrow(cum_foi)),
      foi_lo  = numeric(nrow(cum_foi)),
      foi_hi  = numeric(nrow(cum_foi))
    )
    
    for (j in 1:nrow(cum_foi)) {
      if (j == 1) {
        # First row: Daily FOI is the same as the cumulative FOI for that day
        daily_foi[j, "foi_mid"] <- cum_foi[j, "foi_mid"]
        daily_foi[j, "foi_lo"]  <- cum_foi[j, "foi_lo"]
        daily_foi[j, "foi_hi"]  <- cum_foi[j, "foi_hi"]
      } else {
        # Subsequent rows: Compute the difference with the previous row
        daily_foi[j, "foi_mid"] <- cum_foi[j, "foi_mid"] - cum_foi[j - 1, "foi_mid"]
        daily_foi[j, "foi_lo"] <- cum_foi[j, "foi_lo"] - cum_foi[j - 1, "foi_lo"]
        daily_foi[j, "foi_hi"] <- cum_foi[j, "foi_hi"] - cum_foi[j - 1, "foi_hi"]
      }
    }
    daily_foi_list[[i]] <- daily_foi
  }
  
  # Return the list of daily FOI dataframes
  return(daily_foi_list)
}

# 5-1. loop through admin
downscale_foi_adm <- function(cum_foi_list) {
  
  # Outer `lapply` to process each year
  daily_foi_list <- lapply(cum_foi_list, function(year_data_list) {
    
    # Inner `lapply` to process each admin unit within the year
    lapply(year_data_list, function(cum_foi) {
      
      # Create a dataframe to store daily FOI
      daily_foi <- data.frame(
        date = cum_foi$date,
        admin = cum_foi$admin, 
        foi_mid = numeric(nrow(cum_foi)),
        foi_lo  = numeric(nrow(cum_foi)),
        foi_hi  = numeric(nrow(cum_foi))
      )
      
      # Compute daily FOI based on cumulative FOI
      daily_foi <- within(daily_foi, {
        foi_mid[1] <- cum_foi$foi_mid[1]
        foi_lo[1]  <- cum_foi$foi_lo[1]
        foi_hi[1]  <- cum_foi$foi_hi[1]
        
        foi_mid[-1] <- diff(cum_foi$foi_mid)
        foi_lo[-1]  <- diff(cum_foi$foi_lo)
        foi_hi[-1]  <- diff(cum_foi$foi_hi)
      })
      
      return(daily_foi)  # Return the processed data frame for this admin unit
    })
  })
  
  # Set the names for outer and inner lists
  names(daily_foi_list) <- names(cum_foi_list)
  for (i in seq_along(daily_foi_list)) {
    names(daily_foi_list[[i]]) <- names(cum_foi_list[[i]])
  }
  
  # Return the list of daily FOI data frames
  return(daily_foi_list)
}

