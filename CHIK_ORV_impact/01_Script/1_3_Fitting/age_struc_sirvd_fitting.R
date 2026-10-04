sirv_sim_coverageSwitch_Fatal <- function(
    # Model dimensions
  T, A,
  N,                       # population by age, length A
  r,                       # aging rates, length A
  
  # Epidemic parameters
  base_beta,               # length T
  I0_draw,                 # initial infected, length A
  R0,                      # initial prportion of immunity (underlying immunity)
  rho,                     # detection probability
  gamma,                   # recovery rate
  
  # Vaccination parameters
  delay,                   # day (index) to start Phase 1
  vaccine_end,             # day (index) to stop all vaccination
  coverage_threshold,      # e.g. 0.7 (70%)
  
  # target_age can be:
  #   1) a 0/1 vector of length A (e.g. c(0,0,1,1,...) )
  #   2) or a vector of actual indices (c(3,4,5)) -- but here we assume 0/1
  target_age,              # vector of length A, with 1 if that age is prioritized
  
  weekly_supply_phase1,    # total "weekly" supply in Phase 1
  weekly_supply_phase2,    # total "weekly" supply in Phase 2
  
  VE_block,                  # vaccine efficacy (infection blocking leaky)
  VE_p,                      # vaccine efficacy (disease blocking)
  FR_infection
){
  # -------------------------
  # 1) Initialize compartments
  # -------------------------
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  IV <- matrix(0, nrow = A, ncol = T)
  RV <- matrix(0, nrow = A, ncol = T)
  IFatal <- matrix(0, nrow = A, ncol = T)
  IVFatal <- matrix(0, nrow = A, ncol = T)
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  IFatal_detect <- matrix(0, nrow = A, ncol = T)
  IVFatal_detect <- matrix(0, nrow = A, ncol = T)
  
  # Additional outputs
  # raw_allocation[a,t] = how many doses allocated to age a at time t
  raw_allocation <- matrix(0, nrow = A, ncol = T)
  
  # coverage_frac[a,t] = fraction of age a that is in V+RV at the end of time t
  coverage_frac <- matrix(0, nrow = A, ncol = T)
  
  # coverage_target[t] = overall coverage in the target group at the end of time t
  coverage_target <- rep(0, T)
  
  # Initial conditions at t=1
  for(a in seq_len(A)) {
    S[a,1] <- N[a] - I0_draw[a]
    I[a,1] <- I0_draw[a]
    R_[a,1] <- 0
    V[a,1]  <- 0
    IV[a,1] <- 0
    RV[a,1] <- 0
    IFatal[a,1] <- 0
    IVFatal[a,1] <- 0 
    age_stratified_cases[a,1] <- rho * (I[a,1] + IV[a,1])
    IFatal_detect[a,1] <- rho * IFatal[a,1]
    IVFatal_detect[a,1] <- rho * IVFatal[a,1]
  }
  
  # Identify which ages are target
  target_indices <- which(target_age == 1)  
  target_pop <- sum(N[target_indices])
  total_pop  <- sum(N)
  
  # If no target group, skip coverage-based logic
  no_target <- (length(target_indices) == 0 || target_pop == 0)
  
  # coverage_switch: has the coverage threshold been met?
  coverage_switch <- FALSE
  switch_day <- NA  # if we never switch
  
  # Record initial coverage at t=1
  for(a in seq_len(A)) {
    coverage_frac[a, 1] <- (V[a,1] + RV[a,1]) / N[a]
  }
  if(!no_target) {
    coverage_target[1] <- sum(V[target_indices, 1] + RV[target_indices, 1]) / target_pop
  } else {
    coverage_target[1] <- 0  # or NA
  }
  
  # ---------------------------
  # 2) Forward simulation t=2..T
  # ---------------------------
  for(t_ in 2:T) {
    
    # 2.A) Check coverage in the target group at end of t_-1
    #      skip if no_target = TRUE
    if(!no_target) {
      target_vaccinated_prev <- sum(V[target_indices, t_-1] + RV[target_indices, t_-1])
      coverage_now <- target_vaccinated_prev / target_pop
      
      # If we haven't switched yet AND coverage >= coverage_threshold, switch
      if(!coverage_switch && coverage_now >= coverage_threshold) {
        coverage_switch <- TRUE
        switch_day <- t_
      }
    } else {
      coverage_now <- 0
    }
    
    # 2.B) Force of infection
    infected_prev <- sum(I[, t_-1] + IV[, t_-1])
    region_pop    <- sum(N)  # or sum(S+I+R+V+IV+RV) if dynamic
    phi <- base_beta[t_-1] * (infected_prev / region_pop)
    
    # 2.C) Update compartments by age
    for(a in seq_len(A)) {
      
      # Decide vacc_eff
      # If no_target and weekly_supply_phase1=0 & weekly_supply_phase2=0 => no vaccine
      if(t_ < delay || t_ > vaccine_end) {
        # Outside vaccination window
        raw_allocation[a, t_] <- 0
        vacc_eff <- 0
        
      } else {
        # We are within [delay, vaccine_end]
        if(!coverage_switch) {
          # PHASE 1: target only
          if(!no_target && target_age[a] == 1) {
            # Proportionally allocate weekly_supply_phase1 to this age
            pop_frac_a <- N[a] / target_pop
            raw_allocation_a <- weekly_supply_phase1 * pop_frac_a
            raw_allocation[a, t_] <- raw_allocation_a
            
            # Convert to per-capita rate
            vacc_eff <- raw_allocation_a / N[a]
            # Optionally cap
            # vacc_eff <- min(vacc_eff, 1.0)
          } else {
            raw_allocation[a, t_] <- 0
            vacc_eff <- 0
          }
        } else {
          # PHASE 2: all ages
          pop_frac_a <- N[a] / total_pop
          raw_allocation_a <- weekly_supply_phase2 * pop_frac_a
          raw_allocation[a, t_] <- raw_allocation_a
          
          vacc_eff <- raw_allocation_a / N[a]
          # Optionally cap
          # vacc_eff <- min(vacc_eff, 1.0)
        }
      }
      
      # -----------------------
      # Standard SIRV transitions
      # -----------------------
      if(a == 1) {
        # age 1
        S[a,t_]  <- S[a,t_-1] - phi * S[a,t_-1] -
          vacc_eff * S[a,t_-1] -
          r[a] * S[a,t_-1]
        I[a,t_]  <- I[a,t_-1] + phi * S[a,t_-1] -
          gamma * I[a,t_-1] -
          r[a] * I[a,t_-1]
        R_[a,t_] <- R_[a,t_-1] + gamma * I[a,t_-1] * (1 - FR_infection[a]) -
          r[a] * R_[a,t_-1]
        IFatal[a,t_] <- IFatal[a,t_-1] + I[a, t_-1] * FR_infection[a]
        V[a,t_]  <- V[a,t_-1] + vacc_eff * S[a,t_-1] -
          (1 - VE_block) * phi * V[a,t_-1] -
          r[a] * V[a,t_-1]
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block) * phi * V[a,t_-1] -
          gamma * IV[a,t_-1] -
          r[a] * IV[a,t_-1]
        RV[a,t_] <- RV[a,t_-1] + gamma * IV[a,t_-1] * (1 - FR_infection[a]) -
          r[a] * RV[a,t_-1]
        IVFatal[a,t_] <- IVFatal[a,t_-1] + IV[a, t_-1] * FR_infection[a] * (1 - VE_p)
      } else {
        # age 2..A
        S[a,t_]  <- S[a,t_-1] - phi * S[a,t_-1] -
          vacc_eff * S[a,t_-1] -
          r[a] * S[a,t_-1] +
          r[a-1] * S[a-1, t_-1]
        I[a,t_]  <- I[a,t_-1] + phi * S[a,t_-1] -
          gamma * I[a,t_-1] -
          r[a] * I[a,t_-1] +
          r[a-1] * I[a-1, t_-1]
        R_[a,t_] <- R_[a,t_-1] + gamma * I[a,t_-1] * (1 - FR_infection[a]) -
          r[a] * R_[a,t_-1] +
          r[a-1] * R_[a-1, t_-1]
        V[a,t_]  <- V[a,t_-1] + vacc_eff * S[a,t_-1] -
          (1 - VE_block) * phi * V[a,t_-1] -
          r[a] * V[a,t_-1] +
          r[a-1] * V[a-1, t_-1]
        IFatal[a,t_] <- IFatal[a,t_-1] + I[a, t_-1] * FR_infection[a]
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block) * phi * V[a,t_-1] -
          gamma * IV[a,t_-1] -
          r[a] * IV[a,t_-1] +
          r[a-1] * IV[a-1, t_-1]
        RV[a,t_] <- RV[a,t_-1] + gamma * IV[a,t_-1] * (1 - FR_infection[a]) -
          r[a] * RV[a,t_-1] +
          r[a-1] * RV[a-1, t_-1]
        IVFatal[a,t_] <- IVFatal[a,t_-1] + IV[a, t_-1] * FR_infection[a] * (1 - VE_p)
      }
      
      # clamp at 0
      S[a,t_]  <- pmax(0, S[a,t_])
      I[a,t_]  <- pmax(0, I[a,t_])
      R_[a,t_] <- pmax(0, R_[a,t_])
      V[a,t_]  <- pmax(0, V[a,t_])
      IFatal[a,t_] <- pmax(0, IFatal[a,t_])
      IV[a,t_] <- pmax(0, IV[a,t_])
      RV[a,t_] <- pmax(0, RV[a,t_])
      IVFatal[a,t_] <- pmax(0, IVFatal[a,t_])
    }
    
    # 2.D) After updating compartments, record coverage fractions & reported cases
    for(a in seq_len(A)) {
      coverage_frac[a, t_] <- (V[a,t_] + RV[a,t_]) / N[a]
      age_stratified_cases[a,t_] <- rho * (I[a,t_] + IV[a,t_])
      IFatal_detect[a,t_] <- rho * IFatal[a,t_]
      IVFatal_detect[a,t_] <- rho * IVFatal[a,t_]
    }
    if(!no_target) {
      coverage_target[t_] <- sum(V[target_indices, t_] + RV[target_indices, t_]) / target_pop
    } else {
      coverage_target[t_] <- 0  # or NA
    }
  }
  
  # Return compartments and additional tracking
  list(
    S = S, I = I, R = R_,
    V = V, IV = IV, RV = RV,
    IFatal = IFatal, IVFatal = IVFatal,
    
    # epidemiological metrics
    age_stratified_cases = age_stratified_cases,
    IFatal_detect = IFatal_detect, IVFatal_detect = IVFatal_detect,
    # vaccination logic
    raw_allocation = raw_allocation,
    coverage_frac  = coverage_frac,     # coverage fraction per age over time
    coverage_target = coverage_target,  # coverage in the target group over time
    
    # phase-switch logic
    coverage_switch = coverage_switch,
    switch_day = switch_day
  )
}


fatal <- hosp_fatal_all$fatal_rate
fatal <- c(rep(fatal[1],2),
           rep(fatal[2],2),
           rep(fatal[3],2),
           rep(fatal[4],2),
           rep(fatal[5],2),
           rep(fatal[6],2),
           rep(fatal[7],2),
           rep(fatal[8],2),
           rep(fatal[9],2))

sim_result_3 <- sirv_sim_coverageSwitch_Fatal(
  T = 52, 
  A = 18,
  N = N_region1,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  rho = rho,
  gamma = gamma,
  
  delay = 1,                   # Vaccination starts from Week x
  VE_block = 0.75,             # Vaccine efficacy
  target_age = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),  # Targeting specific age groups
  coverage_threshold = 1,
  total_coverage = 0.6,
  weekly_delivery_speed = 0.05,
  VE_p = 0.8, 
  FR_infection = fatal
)

sim_result_novacc <- sirv_sim_coverageSwitch_Fatal(
  T = 52, 
  A = 18,
  N = N_region1,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  rho = rho,
  gamma = gamma,
  
  delay = 52,                   # Vaccination starts from Week x
  VE_block = 0,             # Vaccine efficacy
  target_age = c(rep(0,18)),  # Targeting specific age groups
  vaccine_end = 0,
  coverage_threshold = 0,
  weekly_supply_phase1 = 0, # x of pop in y weeks
  weekly_supply_phase2 = 0,
  VE_p = 0, 
  FR_infection = fatal
)

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Median")
summary_Ifatal_pre <- as.data.frame.table(sim_result_novacc$IFatal_detect, responseName = "IFatal")
summary_IVfatal_pre <- as.data.frame.table(sim_result_novacc$IVFatal_detect, responseName = "IVFatal")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_Ifatal_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_IVfatal_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_cases_pre$Scenario <- "Pre-Vaccination"

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_Ifatal_pre <- summary_Ifatal_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_Ifatal_pre$Scenario <- "Pre-Vaccination"

summary_Ifatal_pre_all <- summary_Ifatal_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

#-------------------------------------------------------------------------------
# grid search
#-------------------------------------------------------------------------------
total_supply_rates <- seq(0.1, 1, by = 0.1)         # Total supply rate (0% to 100%)
weekly_delivery_speeds <- seq(0.01, 0.05, by = 0.01) # Weekly delivery speed (1% to 5%)
VE_blocks <- seq(0.1, 1, by = 0.1)               # Transmission blocking efficacy
VE_ps <- seq(0.1, 1, by = 0.1)

param_grid <- expand.grid(
  TotalSupplyRate = total_supply_rates,
  WeeklyDeliverySpeed = weekly_delivery_speeds,
  VE_block = VE_blocks,
  VE_p = VE_ps
)
param_grid <- param_grid %>%
  mutate(
    TotalSupplyRate = as.numeric(TotalSupplyRate),
    WeeklyDeliverySpeed = as.numeric(WeeklyDeliverySpeed),
    VE_block = as.numeric(VE_block),
    VE_p = as.numeric(VE_p)
  )
results <- data.frame(
  CombinationID = seq_len(nrow(param_grid)),
  TotalSupplyRate = param_grid$TotalSupplyRate,
  WeeklyDeliverySpeed = param_grid$WeeklyDeliverySpeed,
  VE_block = param_grid$VE_block,
  VE_p = param_grid$VE_p,
  Impact = numeric(nrow(param_grid)),          # To store impact or relevant metric
  TotalIFatal = numeric(nrow(param_grid)),     # To store total detected fatalities
  TotalIVFatal = numeric(nrow(param_grid))     # To store total vaccinated fatalities
)

all_results_v4 <- list()  # To store detailed simulation summaries for each iteration

for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %.0f of %.0f", i_iter, num_iters)) 
  
  iteration_result <- list()
  iteration_sim_dfs <- list()
  
  # Iterate through all combinations in param_grid
  for (param_index in seq_len(nrow(param_grid))) {
    params <- param_grid[param_index, ]
    
    total_supply <- sum(N_region1) * params$TotalSupplyRate
    weekly_supply_phase1 <- total_supply * params$WeeklyDeliverySpeed
    total_vaccine_duration_weeks <- ceiling(total_supply / weekly_supply_phase1)
    unused_supply <- 0
    
    if (total_vaccine_duration_weeks > 52) {
      unused_supply <- total_supply - (weekly_supply_phase1 * 52)
      total_vaccine_duration_weeks <- 52  # Cap duration at 52 weeks
    }
    
    vaccine_end <- 1 + total_vaccine_duration_weeks

    for (scenario_index in seq_along(target_age_list)) {
      target <- target_age_list[[scenario_index]]
      
      print(sprintf(
        "Iteration %d, Scenario %d: Total Supply %.1f%%, Weekly Speed %.1f%%, VE_block %.2f, VE_p %.2f, Duration %d weeks",
        i_iter, scenario_index, params$TotalSupplyRate * 100,
        params$WeeklyDeliverySpeed * 100, params$VE_block, params$VE_p,
        total_vaccine_duration_weeks
      ))
      
      sim_result <- sirv_sim_coverageSwitch_Fatal(
        T = 52,
        A = 18,
        N = N_region1,
        r = rep(0, 18),
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        delay = 1,
        vaccine_end = vaccine_end,
        VE_block = params$VE_block,
        VE_p = params$VE_p,
        coverage_threshold = 0.9,
        FR_infection = fatal,
        target_age = target,
        weekly_supply_phase1 = weekly_supply_phase1,
        weekly_supply_phase2 = weekly_supply_phase1
      )
      
      # Summarize detected fatalities
      total_IFatal_detect <- sum(sim_result$IFatal_detect)
      total_IVFatal_detect <- sum(sim_result$IVFatal_detect)
      
      # Update cumulative results
      results$TotalIFatal[param_index] <- results$TotalIFatal[param_index] + 
        total_IFatal_detect
      results$TotalIVFatal[param_index] <- results$TotalIVFatal[param_index] + 
        total_IVFatal_detect
      
      # Summarize `age_stratified_cases` into `sim_df`
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
        mutate(
          Iteration = i_iter,
          Scenario = scenario_index,
          TotalSupplyRate = params$TotalSupplyRate,
          WeeklyDeliverySpeed = params$WeeklyDeliverySpeed,
          VE_block = params$VE_block,
          VE_p = params$VE_p,
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store summarized results in `iteration_sim_dfs`
      iteration_sim_dfs[[length(iteration_sim_dfs) + 1]] <- sim_df
      
      # Store simulation results in `iteration_result`
      iteration_result[[paste0(
        "Scenario_", scenario_index, 
        "_Supply_", params$TotalSupplyRate,
        "_Speed_", params$WeeklyDeliverySpeed,
        "_VEblock_", params$VE_block,
        "_VEp_", params$VE_p
      )]] <- sim_result
    }
  }
  
  # Store detailed results for this iteration
  all_results_v4[[i_iter]] <- iteration_result
}

unique_scenarios <- unique(unlist(lapply(all_results_v4, function(iteration) names(iteration))))

summary_list_v4 <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_scenario_data(scenario, all_results_v4)
  
  # Compute median across iterations
  median_summary <- compute_median(scenario_data)
  
  # Store the summarized data
  summary_list_v4[[scenario]] <- median_summary
}

final_summ_v4 <- lapply(summary_list_v4, function(df) {
  df %>%
    group_by(Week, Scenario) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(MedianCases, na.rm = TRUE),  # Summarize MedianCases
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_cases_pre_all$Median,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100
    ) %>% mutate(
      Supply = as.numeric(sub("Supply_", "", stringr::str_extract(Scenario, "Supply_\\d+(\\.\\d+)?"))),
      Speed = as.numeric(sub("Speed_", "", stringr::str_extract(Scenario, "Speed_\\d+(\\.\\d+)?"))),
      VE_block = as.numeric(sub("VEblock_", "", stringr::str_extract(Scenario, "VEblock_\\d+(\\.\\d+)?"))),
      VE_p = as.numeric(sub("VEp_", "", stringr::str_extract(Scenario, "VEp_\\d+(\\.\\d+)?")))
    )
})

final_summ_s1_ve_v4 <- final_summ_v4[startsWith(names(final_summ_v4), "Scenario_4_")]

total_impact_s1_ve_v4 <- lapply(final_summ_s1_ve_v4, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s1_ve_v4 <- do.call(rbind, total_impact_s1_ve_v4)


summary_table_s1_ve_v4 <- final_total_impact_s1_ve_v4 %>%
  # Extract Delay and Duration using regex
  mutate(
    VEblock   = as.numeric(sub("VEblock_", "", stringr::str_extract(Scenario, "VEblock_\\d+(\\.\\d+)?"))),
    VEp       = as.numeric(sub("VEp_", "", stringr::str_extract(Scenario, "VEp_\\d+(\\.\\d+)?"))),
    Supply    = as.numeric(sub("Supply_", "", str_extract(Scenario, "Supply_\\d+(\\.\\d+)?"))),
    Speed     = as.numeric(sub("Speed_", "", str_extract(Scenario, "Speed_\\d+(\\.\\d+)?")))
  )

ggplot(data = summary_table_s1_ve_v4, aes(x = VEblock, y = Supply, z = tot_impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_grid(Speed ~ VEp, labeller = label_both)+
  labs(caption = "Note: VEp = disease blocking vaccine efficacy")+
  theme(
    plot.caption = element_text(hjust = 0, size = 7)  # Left-align caption
  )

#-------------------------------------------------------------------------------
# fatal cases
#-------------------------------------------------------------------------------
fatal_list_v4 <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_fatal_data(scenario, all_results_v4)
  
  # Compute median across iterations
  median_summary <- compute_median_fatal(scenario_data)
  
  # Store the summarized data
  fatal_list_v4[[scenario]] <- median_summary
}

fatal_summ_v4 <- lapply(fatal_list_v4, function(df) {
  df %>%
    group_by(Week, Scenario) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(MedianFatal, na.rm = TRUE),  # Summarize MedianCases
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_Ifatal_pre_all$Median,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100
    ) %>% mutate(
      Supply = as.numeric(sub("Supply_", "", stringr::str_extract(Scenario, "Supply_\\d+(\\.\\d+)?"))),
      Speed = as.numeric(sub("Speed_", "", stringr::str_extract(Scenario, "Speed_\\d+(\\.\\d+)?"))),
      VE_block = as.numeric(sub("VEblock_", "", stringr::str_extract(Scenario, "VEblock_\\d+(\\.\\d+)?"))),
      VE_p = as.numeric(sub("VEp_", "", stringr::str_extract(Scenario, "VEp_\\d+(\\.\\d+)?")))
    )
})

fatal_summ_s1_ve_v4 <- fatal_summ_v4[startsWith(names(fatal_summ_v4), "Scenario_4_")]

fatal_impact_s1_ve_v4 <- lapply(fatal_summ_s1_ve_v4, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
              pre_vacc = sum(pre_vacc),
              tot_impact = tot_diff / pre_vacc * 100
              )
})

fatal_total_impact_s1_ve_v4 <- do.call(rbind, fatal_impact_s1_ve_v4)

fatal_table_s1_ve_v4 <- fatal_total_impact_s1_ve_v4 %>%
  # Extract Delay and Duration using regex
  mutate(
    VEblock   = as.numeric(sub("VEblock_", "", stringr::str_extract(Scenario, "VEblock_\\d+(\\.\\d+)?"))),
    VEp       = as.numeric(sub("VEp_", "", stringr::str_extract(Scenario, "VEp_\\d+(\\.\\d+)?"))),
    Supply    = as.numeric(sub("Supply_", "", str_extract(Scenario, "Supply_\\d+(\\.\\d+)?"))),
    Speed     = as.numeric(sub("Speed_", "", str_extract(Scenario, "Speed_\\d+(\\.\\d+)?")))
  )

ggplot(data = fatal_table_s1_ve_v4, aes(x = VEblock, y = Supply, z = tot_impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_grid(Speed ~ VEp, labeller = label_both)+
  labs(caption = "Note: VEp = disease blocking vaccine efficacy")+
  theme(
    plot.caption = element_text(hjust = 0, size = 7)  # Left-align caption
  )

#-------------------------------------------------------------------------------
# comparison within same condition: speed 0.05, VEp = 0.8
#-------------------------------------------------------------------------------
ve80_speed5 <- Filter(function(df) any(df$Speed == 0.05 & df$VE_p == 0.8), final_summ_v4)
ve80_speed5 <- do.call(rbind, ve80_speed5) %>% mutate(
  Scenario = as.numeric(sub("Scenario_", "", str_extract(Scenario, "Scenario_\\d+(\\.\\d+)?")))
)

ve80_speed5_summ <- ve80_speed5 %>%
  group_by(Scenario, Supply, Speed, VE_block, VE_p) %>%  # Group by these variables
  summarise(
    Median = sum(Median),
    pre_vacc = sum(pre_vacc),
    diff = pre_vacc - Median,
    impact = diff / pre_vacc * 100,
    .groups = "drop"  # Avoid unnecessary grouping in the output
  )

ggplot(data = ve80_speed5_summ, aes(x = VE_block, y = Supply, z = impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~Scenario, labeller = label_both)+
  labs(caption = "Note: VEp = disease blocking vaccine efficacy")+
  theme(
    plot.caption = element_text(hjust = 0, size = 7)  # Left-align caption
  )
