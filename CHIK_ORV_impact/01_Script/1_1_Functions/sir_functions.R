
# Function to extract data for a given scenario
extract_scenario_data <- function(scenario, all_results) {
  scenario_data <- lapply(all_results, function(iteration) {
    if (scenario %in% names(iteration)) {
      iteration[[scenario]]  # Extract the data frame
    } else {
      NULL  # Handle missing scenario gracefully
    }
  })
  
  # Remove NULL values
  scenario_data <- scenario_data[!sapply(scenario_data, is.null)]
  
  # Combine all iterations for this scenario into one data frame
  combined_scenario_data <- do.call(rbind, lapply(scenario_data, function(sim_result) {
    as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
      mutate(
        Scenario = scenario,
        AgeGroup = as.numeric(Var1),
        Week = as.numeric(Var2)
      )
  }))
  
  return(combined_scenario_data)
}


compute_median <- function(scenario_data) {
  scenario_data %>%
    group_by(Week, AgeGroup) %>%  # Group by Week and AgeGroup
    summarise(
      MedianCases = median(Cases, na.rm = TRUE),  # Compute median
      .groups = "drop"
    ) %>%
    mutate(Scenario = unique(scenario_data$Scenario))  # Add Scenario column
}


## functions to extract
extract_cuminf_data <- function(scenario, all_results_v2) {
  scenario_data <- lapply(all_results_v2, function(iteration) {
    if (scenario %in% names(iteration)) {
      iteration[[scenario]]  # Extract the data frame
    } else {
      NULL  # Handle missing scenario gracefully
    }
  })
  
  # Remove NULL values
  scenario_data <- scenario_data[!sapply(scenario_data, is.null)]
  
  # Combine all iterations for this scenario into one data frame
  combined_scenario_data <- do.call(rbind, lapply(scenario_data, function(sim_result) {
    
    R_data <- as.data.frame.table(sim_result$R, responseName = "R")
    RV_data <- as.data.frame.table(sim_result$RV, responseName = "RV")
    
    combined_data <- merge(R_data, RV_data, by = c("Var1", "Var2"))
    
    combined_data <- combined_data %>%
      mutate(
        CumInf = (R + RV) * rho,  # Sum R and RV
        Scenario = scenario,
        AgeGroup = as.numeric(Var1),
        Week = as.numeric(Var2)
      )
    combined_data <- combined_data %>%
      select(AgeGroup, Week, CumInf, Scenario)
    
    return(combined_data)
    
  }))
  
  return(combined_scenario_data)
}

compute_median_cuminf <- function(scenario_data) {
  scenario_data %>%
    group_by(Week, AgeGroup) %>%  # Group by Week and AgeGroup
    summarise(
      MedianCumInf = median(CumInf, na.rm = TRUE),  # Compute median
      .groups = "drop"
    ) %>%
    mutate(Scenario = unique(scenario_data$Scenario))  # Add Scenario column
}

extract_raw_allocation <- function(scenario, all_results) {
  raw_allocation_data <- lapply(all_results, function(iteration) {
    if (scenario %in% names(iteration)) {
      iteration[[scenario]]$raw_allocation  # Extract raw_allocation directly
    } else {
      NULL  # Handle missing scenario gracefully
    }
  })
  
  # Remove NULL values
  raw_allocation_data <- raw_allocation_data[!sapply(raw_allocation_data, is.null)]
  
  # Combine all iterations for this scenario into one data frame
  combined_raw_allocation <- do.call(rbind, lapply(seq_along(raw_allocation_data), function(i) {
    raw_allocation_df <- as.data.frame(raw_allocation_data[[i]])  # Convert matrix to data frame
    colnames(raw_allocation_df) <- paste0("Week_", seq_len(ncol(raw_allocation_df)))  # Rename columns for weeks
    raw_allocation_df$AgeGroup <- seq_len(nrow(raw_allocation_df))  # Add AgeGroup as a column
    raw_allocation_df$Iteration <- i  # Add Iteration identifier
    raw_allocation_df$Scenario <- scenario  # Add Scenario identifier
    raw_allocation_df  # Return modified data frame
  }))
  
  return(combined_raw_allocation)
}


extract_fatal_data <- function(scenario, all_results) {
  scenario_data <- lapply(all_results, function(iteration) {
    if (scenario %in% names(iteration)) {
      iteration[[scenario]]  # Extract the data frame
    } else {
      NULL  # Handle missing scenario gracefully
    }
  })
  
  # Remove NULL values
  scenario_data <- scenario_data[!sapply(scenario_data, is.null)]
  
  # Combine all iterations for this scenario into one data frame
  combined_scenario_data <- do.call(rbind, lapply(scenario_data, function(sim_result) {
    
    IFatal_data  <- as.data.frame.table(sim_result$IFatal_detect, responseName = "IFatal")
    IVFatal_data <- as.data.frame.table(sim_result$IVFatal_detect, responseName = "IVFatal")
    
    combined_data <- merge(IFatal_data, IVFatal_data, by = c("Var1", "Var2"))
    
    combined_data <- combined_data %>%
      mutate(
        Fatal = (IFatal + IVFatal),  
        Scenario = scenario,
        AgeGroup = as.numeric(Var1),
        Week = as.numeric(Var2)
      )
    combined_data <- combined_data %>%
      select(AgeGroup, Week, Fatal, Scenario)
    
    return(combined_data)
    
  }))
  
  return(combined_scenario_data)
}


compute_median_fatal <- function(scenario_data) {
  scenario_data %>%
    group_by(Week, AgeGroup) %>%  # Group by Week and AgeGroup
    summarise(
      MedianFatal = median(Fatal, na.rm = TRUE),  # Compute median
      .groups = "drop"
    ) %>%
    mutate(Scenario = unique(scenario_data$Scenario))  # Add Scenario column
}

# for a single scenario 
extract_single_scenario_data <- function(i, scenario_result) {
  # Access the specific scenario data
  scenario_data <- scenario_result[[i]]
  
  # Extract age_stratified_cases and convert it to a data frame
  sim_result <- as.data.frame.table(scenario_data$age_stratified_cases, responseName = "Cases") %>%
    mutate(
      Scenario = paste0("Scenario_", i),  # Add Scenario label
      AgeGroup = as.numeric(Var1),       # Convert Var1 to AgeGroup
      Week = as.numeric(Var2)            # Convert Var2 to Week
    )
  
  return(sim_result)
}

extract_single_fatal_data <- function(i, scenario_result) {
  # Access the specific scenario data
  scenario_data <- scenario_result[[i]]
  
  # Extract age_stratified_cases and convert it to a data frame
  sim_result <- as.data.frame.table(scenario_data$IFatal_detect, responseName = "IFatal_Cases") %>%
    left_join(
      as.data.frame.table(scenario_data$IVFatal_detect, responseName = "IVFatal_Cases"),
      by = c("Var1", "Var2")
    ) %>%
    mutate(
      Cases = IFatal_Cases + IVFatal_Cases,
      Scenario = paste0("Scenario_", i),  # Add Scenario label
      AgeGroup = as.numeric(Var1),       # Convert Var1 to AgeGroup
      Week = as.numeric(Var2)            # Convert Var2 to Week
    )
  
  return(sim_result)
}


compute_median <- function(scenario_data) {
  scenario_data %>%
    group_by(Week, AgeGroup) %>%  # Group by Week and AgeGroup
    summarise(
      MedianCases = median(Cases, na.rm = TRUE),  # Compute median
      .groups = "drop"
    ) %>%
    mutate(Scenario = unique(scenario_data$Scenario))  # Add Scenario column
}


