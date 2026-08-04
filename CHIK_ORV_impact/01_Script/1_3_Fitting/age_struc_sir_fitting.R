sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,
  N,                       # population by age, length A
  r,                       # aging rates, length A
  
  # Epidemic parameters
  base_beta,               # length T
  I0_draw,                 # initial infected, length A
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
  
  VE_block                  # vaccine efficacy (block rate)
){
  # -------------------------
  # 1) Initialize compartments
  # -------------------------
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  IV <- matrix(0, nrow = A, ncol = T)
  RV <- matrix(0, nrow = A, ncol = T)
  
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  
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
    R[a,1] <- 0
    V[a,1]  <- 0
    IV[a,1] <- 0
    RV[a,1] <- 0
    
    age_stratified_cases[a,1] <- rho * (I[a,1] + IV[a,1])
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
        R[a,t_] <- R[a,t_-1] + gamma * I[a,t_-1] -
          r[a] * R[a,t_-1]
        V[a,t_]  <- V[a,t_-1] + vacc_eff * S[a,t_-1] -
          (1 - VE_block) * phi * V[a,t_-1] -
          r[a] * V[a,t_-1]
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block) * phi * V[a,t_-1] -
          gamma * IV[a,t_-1] -
          r[a] * IV[a,t_-1]
        RV[a,t_] <- RV[a,t_-1] + gamma * IV[a,t_-1] -
          r[a] * RV[a,t_-1]
        
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
        R[a,t_] <- R[a,t_-1] + gamma * I[a,t_-1] -
          r[a] * R[a,t_-1] +
          r[a-1] * R[a-1, t_-1]
        V[a,t_]  <- V[a,t_-1] + vacc_eff * S[a,t_-1] -
          (1 - VE_block) * phi * V[a,t_-1] -
          r[a] * V[a,t_-1] +
          r[a-1] * V[a-1, t_-1]
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block) * phi * V[a,t_-1] -
          gamma * IV[a,t_-1] -
          r[a] * IV[a,t_-1] +
          r[a-1] * IV[a-1, t_-1]
        RV[a,t_] <- RV[a,t_-1] + gamma * IV[a,t_-1] -
          r[a] * RV[a,t_-1] +
          r[a-1] * RV[a-1, t_-1]
      }
      
      # clamp at 0
      S[a,t_]  <- pmax(0, S[a,t_])
      I[a,t_]  <- pmax(0, I[a,t_])
      R[a,t_] <- pmax(0, R[a,t_])
      V[a,t_]  <- pmax(0, V[a,t_])
      IV[a,t_] <- pmax(0, IV[a,t_])
      RV[a,t_] <- pmax(0, RV[a,t_])
    }
    
    # 2.D) After updating compartments, record coverage fractions & reported cases
    for(a in seq_len(A)) {
      coverage_frac[a, t_] <- (V[a,t_] + RV[a,t_]) / N[a]
      age_stratified_cases[a,t_] <- rho * (I[a,t_] + IV[a,t_])
    }
    if(!no_target) {
      coverage_target[t_] <- sum(V[target_indices, t_] + RV[target_indices, t_]) / target_pop
    } else {
      coverage_target[t_] <- 0  # or NA
    }
  }
  
  # Return compartments and additional tracking
  list(
    S = S, I = I, R = R,
    V = V, IV = IV, RV = RV,
    
    # epidemiological metrics
    age_stratified_cases = age_stratified_cases,
    
    # vaccination logic
    raw_allocation = raw_allocation,
    coverage_frac  = coverage_frac,     # coverage fraction per age over time
    coverage_target = coverage_target,  # coverage in the target group over time
    
    # phase-switch logic
    coverage_switch = coverage_switch,
    switch_day = switch_day
  )
}

# should i put fixed weekly throughput rate (e.g. 7%) or specify total coverage target (e.g. 70% of total population over x weeks )
# or make the model halt vaccination once it reaches certain coverage target? 
# varying coverage of varying per-week rate
# varying coverage approach during fixed durations : sum(N)*0.7 /5 (coverage / duration)
# varying per-week rate: sum(N)*0.07 (total pop * per week rate)

sim_result_2 <- sirv_sim_coverageSwitch(
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
  vaccine_end = 25,
  coverage_threshold = 1,
  weekly_supply_phase1 = sum(N_region1)*0.07, # x of pop in y weeks
  weekly_supply_phase2 = sum(N_region1)*0.07
)

sim_result_novacc_2 <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_region1,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  rho = rho,
  gamma = gamma,
  
  delay = 53,                   # Vaccination starts from Week 1
  VE_block = 0,             # Vaccine efficacy
  target_age = rep(0, 18),    # Targeting specific age groups
  vaccine_end = 0,
  coverage_threshold = 1,
  weekly_supply_phase1 = 0, # x of pop in y weeks
  weekly_supply_phase2 = 0
)

# post-vacc
summary_cases <- as.data.frame.table(sim_result_2$age_stratified_cases, responseName = "Median")

# Rename columns for clarity
colnames(summary_cases) <- c("AgeGroup", "Week", "Median")

# Convert columns to numeric for correct ordering
summary_cases <- summary_cases %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )
summary_cases$Scenario <- "Post-Vaccination"

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc_2$age_stratified_cases, responseName = "Median")
summary_R_pre <- as.data.frame.table(sim_result_novacc_2$R, responseName = "R")
summary_RV_pre <- as.data.frame.table(sim_result_novacc_2$RV, responseName = "RV")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_R_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_RV_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_cases_pre$Scenario <- "Pre-Vaccination"

summary_R_pre <- summary_R_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_R_pre$Scenario <- "Pre-Vaccination"

summary_RV_pre <- summary_RV_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_RV_pre$Scenario <- "Pre-Vaccination"

RRV <- summary_R_pre$Median + summary_RV_pre$Median

summary_R_pre$RRV_all <- RRV

ggplot()+
  geom_line(data = summary_cases, aes(x = Week, y = Median, colour = Scenario))+
  geom_line(data = summary_cases_pre, aes(x = Week, y = Median, colour = Scenario))+
  facet_wrap(~AgeGroup)+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_R_pre_all <- summary_R_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(RRV_all)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  ) %>% mutate(
    tot_pop = sum(N_region1),
    Median_rho = Median * rho,
    prop_inf = Median_rho / tot_pop * 100
  )


#-------------------------------------------------------------------------------
# cumulative infections (R+RV) by different prioritised agegroups

target_age_list <- list(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1))

all_results_v2 <- list()  # Initialize storage for all iterations
summary_dfs_v2 <- list()  # Initialize storage for summarized data frames

coverage_delay_combinations <- expand.grid(
  total_coverage = seq(0.1, 1, by = 0.1),  
  delay = seq(1, 52, by = 1)  
)

# Fixed vaccine duration (26 weeks)
vaccine_duration_weeks <- 26

# Run for 1 iteration (use median values)
num_iters <- 1

for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  # Iterate through all combinations of coverage and delay
  for (comb in seq_len(nrow(coverage_delay_combinations))) {
    
    # Extract coverage and delay values
    coverage <- coverage_delay_combinations[comb, "total_coverage"]
    delay <- coverage_delay_combinations[comb, "delay"]
    
    print(sprintf("Iteration: %d, Coverage: %.1f%%, Delay: %d weeks", 
                  i_iter, coverage * 100, delay))
    
    # Inner loop for vaccination scenarios
    for (i in seq_along(target_age_list)) {
      
      target <- target_age_list[[i]]
      
      print(sprintf("Iteration %d, Scenario %d: Target Age = %s", 
                    i_iter, i, toString(target)))
      
      # Simulate the SIRV model
      sim_result <- sirv_sim_coverageSwitch(
        T = 52, 
        A = 18,
        N = N,
        r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        delay = delay,               # Vaccination starts from Week x
        VE_block = 0.75,             # Vaccine efficacy
        target_age = target,         # Targeting specific age groups by scenario
        vaccine_end = delay + vaccine_duration_weeks,
        coverage_threshold = 1,
        weekly_supply_phase1 = sum(N) * coverage / vaccine_duration_weeks, 
        weekly_supply_phase2 = sum(N) * coverage / vaccine_duration_weeks
      )
    
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = paste0("Scenario_", i, "_Delay_", delay, 
                            "_Coverage_", coverage),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs_v2[[length(summary_dfs_v2) + 1]] <- sim_df
      
      # Store results for this combination
      iteration_result[[paste0("Scenario_", i, "_Coverage_", coverage, "_Delay_", delay)]] <- sim_result
    }
  }
  
  # Store all results from this iteration
  all_results_v2[[i_iter]] <- iteration_result
}

# Extract unique scenarios from the nested structure
unique_scenarios <- unique(unlist(lapply(all_results_v2, function(iteration) names(iteration))))

# Initialize a list to store results for all scenarios
scenario_summary_list <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_scenario_data(scenario, all_results_v2)
  
  # Compute median across iterations
  median_summary <- compute_median(scenario_data)
  
  # Store the summarized data
  scenario_summary_list[[scenario]] <- median_summary
}

final_summ <- lapply(scenario_summary_list, function(df) {
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
      Delay    = as.numeric(sub("Delay_", "", stringr::str_extract(Scenario, "Delay_\\d+"))),
      Coverage = as.numeric(sub("Coverage_", "", stringr::str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
    )
})

final_summ_s1 <- final_summ[startsWith(names(final_summ), "Scenario_1_")]

total_impact_s1 <- lapply(final_summ_s1, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s1 <- do.call(rbind, total_impact_s1)

summary_table_s1 <- final_total_impact_s1 %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s1$tot_impact <- as.numeric(summary_table_s1$tot_impact)

ggplot(summary_table_s1, aes(x = Coverage, y = Delay, fill = tot_impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% cases averted") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))


# scenario2
final_summ_s2 <- final_summ[startsWith(names(final_summ), "Scenario_2_")]

total_impact_s2 <- lapply(final_summ_s2, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s2 <- do.call(rbind, total_impact_s2)

summary_table_s2 <- final_total_impact_s2 %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s2$tot_impact <- as.numeric(summary_table_s2$tot_impact)

ggplot(summary_table_s2, aes(x = Coverage, y = Delay, fill = tot_impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% cases averted") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

# scenario3
final_summ_s3 <- final_summ[startsWith(names(final_summ), "Scenario_3_")]

total_impact_s3 <- lapply(final_summ_s3, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s3 <- do.call(rbind, total_impact_s3)

summary_table_s3 <- final_total_impact_s3 %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s3$tot_impact <- as.numeric(summary_table_s3$tot_impact)

ggplot(summary_table_s3, aes(x = Coverage, y = Delay, fill = tot_impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% cases averted") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

# scenario4
final_summ_s4 <- final_summ[startsWith(names(final_summ), "Scenario_4_")]

total_impact_s4 <- lapply(final_summ_s4, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s4 <- do.call(rbind, total_impact_s4)

summary_table_s4 <- final_total_impact_s4 %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s4$tot_impact <- as.numeric(summary_table_s4$tot_impact)

ggplot(summary_table_s4, aes(x = Coverage, y = Delay, fill = tot_impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% cases averted") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

#-------------------------------------------------------------------------------
# estimate cumulative new infections (r+rv) by each scenario 

cum_inf_list <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_cuminf_data(scenario, all_results_v2)
  
  # Compute median across iterations
  median_summary <- compute_median_cuminf(scenario_data)
  
  # Store the summarized data
  cum_inf_list[[scenario]] <- median_summary
}

cuminf_summ <- lapply(cum_inf_list, function(df) {
  df %>%
    group_by(Week, Scenario) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(MedianCumInf, na.rm = TRUE),  # Summarize MedianCases
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_R_pre_all$Median_rho,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100,
      # -- extract Delay & Coverage from the Scenario name --
      Delay    = as.numeric(sub("Delay_", "", stringr::str_extract(Scenario, "Delay_\\d+"))),
      Coverage = as.numeric(sub("Coverage_", "", stringr::str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
    )
})

final_summ_s1_cuminf <- cuminf_summ[startsWith(names(cuminf_summ), "Scenario_1_")]

total_impact_s1_cuminf <- lapply(final_summ_s1_cuminf, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
                pre_vacc = sum(pre_vacc),
                tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s1_cuminf <- do.call(rbind, total_impact_s1_cuminf)

summary_table_s1 <- final_total_impact_s1 %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s1$tot_impact <- as.numeric(summary_table_s1$tot_impact)

#-------------------------------------------------------------------------------
# filter 30% supply for all scenarios
#-------------------------------------------------------------------------------

cum_inf_30 <- Filter(function(df) any(df$Coverage == 0.3 & df$Delay == 1), cuminf_summ)
cum_inf_30_summ <- do.call(rbind, cum_inf_30) %>% mutate(
  Scenario = as.numeric(sub("Scenario_", "", str_extract(Scenario, "Scenario_\\d+(\\.\\d+)?")))
)

cum_inf_30_summ <- cum_inf_30_summ %>% mutate(
  tot_pop   = sum(N_region1),
  inf_prop  = Median / tot_pop,
  reduction = (pre_vacc - Median)/pre_vacc * 100 
)

ggplot()+
  geom_line(data = cum_inf_30_summ, aes(x = Week, y = Median, color = factor(Scenario)))+
  #geom_line(data = cum_inf_30_summ, aes(x = Week, y = inf_prop, color = factor(Scenario)))+
  geom_line(data = summary_R_pre_all, aes(x = Week, y = Median_rho, color = Scenario), linetype = "dashed")+
  #scale_y_continuous(labels = scales::percent_format(scale = 0.0001))+
  theme_bw()

ggplot()+
  geom_line(data = cum_inf_30_summ, aes(x = Week, y = reduction, color = factor(Scenario)))+
  #geom_line(data = cum_inf_30_summ, aes(x = Week, y = inf_prop, color = factor(Scenario)))+
  #scale_y_continuous(labels = scales::percent_format(scale = 0.0001))+
  theme_bw()

#-------------------------------------------------------------------------------
# epi curve by scenario
#-------------------------------------------------------------------------------

inf_30 <- Filter(function(df) any(df$Coverage == 0.3 & df$Delay == 1), final_summ)
inf_30_summ <- do.call(rbind, inf_30) %>% mutate(
  Scenario = as.numeric(sub("Scenario_", "", str_extract(Scenario, "Scenario_\\d+(\\.\\d+)?")))
)

inf_30_summ <- inf_30_summ %>% mutate(
  tot_pop   = sum(N_region1),
  inf_prop  = Median / tot_pop,
  reduction = (pre_vacc - Median)/pre_vacc * 100 
)

ggplot()+
  geom_line(data = inf_30_summ, aes(x = Week, y = Median, color = factor(Scenario)))+
  geom_line(data = summary_cases_pre_all, aes(x = Week, y = Median, color = Scenario), linetype = "dashed")+
  scale_y_continuous(labels = scales::comma_format(), breaks = seq(0, 30000, by = 5000))+
  theme_bw()


## extract delay = 1 
inf_delay1 <- Filter(function(df) any(df$Delay == 1), final_summ)
inf_delay1_summ <- do.call(rbind, inf_delay1) %>% mutate(
  Scenario = as.numeric(sub("Scenario_", "", str_extract(Scenario, "Scenario_\\d+(\\.\\d+)?")))
)

inf_delay1_summ <- inf_delay1_summ %>% group_by(Coverage, Scenario) %>%
                       summarise(Median   = sum(Median),
                                 pre_vacc = sum(pre_vacc)) %>% mutate(
                                 reduction = (pre_vacc - Median) / pre_vacc * 100
                                 )

ggplot()+
  geom_line(data = inf_delay1_summ, aes(x = Coverage, y = reduction, 
                                        color = factor(Scenario)))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))+
  theme_bw()
  
## extract coverage 30
alloc_summary_list <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_raw_allocation(scenario, all_results_v2)
  
  # Store the summarized data
  alloc_summary_list[[scenario]] <- scenario_data
}

alloc_summary <- do.call(rbind, alloc_summary_list)%>% mutate(
  Delay    = as.numeric(sub("Delay_", "", stringr::str_extract(Scenario, "Delay_\\d+"))),
  Coverage = as.numeric(sub("Coverage_", "", stringr::str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
)

alloc_summary <- alloc_summary %>% filter(Allocation!=0 & Week == 2)
alloc_summary <- alloc_summary[!duplicated(alloc_summary), ]

alloc_summary_cov0.3 <- alloc_summary[alloc_summary$Coverage==0.3 & alloc_summary$Delay ==1, ]
alloc <- alloc_summary_cov0.3$Allocation
alloc_summary_cov0.3$tot_pop <- c(rep(sum(alloc[1:4]),4),
                                rep(sum(alloc[5:12]),8),
                                rep(sum(alloc[5:18]),14),
                                rep(sum(alloc[15:18]),4))

alloc_summary_cov0.3 <- alloc_summary_cov0.3 %>% mutate(
  alloc_prop = Allocation / tot_pop
)

ggplot() +
  geom_bar(data = alloc_summary_cov0.3, aes(x = AgeGroup, y = alloc_prop, fill = Scenario), stat = "identity", position = "dodge") +
  facet_wrap(~Scenario, ncol = 1, scales = "free_y") +
  theme_bw()+
  scale_x_continuous(breaks = 1:18)+
  scale_y_continuous(labels = scales::percent)+
  ylab("Vaccine allocation (%)")+
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(size = 10), # Adjust facet label text size
    axis.text.x = element_text(angle = 90, hjust = 1) 
  )

#-------------------------------------------------------------------------------
# evaluate by different VE
#-------------------------------------------------------------------------------

target_age_list <- list(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1))

all_results_v3 <- list()  # Initialize storage for all iterations
summary_dfs_v3 <- list()  # Initialize storage for summarized data frames

coverage_ve_combinations <- expand.grid(
  total_ve = seq(0.1, 1, by = 0.1),  
  coverage = seq(0.1, 1, by = 0.1)  
)

# Fixed vaccine duration (26 weeks)
vaccine_duration_weeks <- 26

# Run for 1 iteration (use median values)
num_iters <- 1

for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  # Iterate through all combinations of coverage and delay
  for (comb in seq_len(nrow(coverage_ve_combinations))) {
    
    # Extract coverage and delay values
    ve       <- coverage_ve_combinations[comb, "total_ve"]
    coverage <- coverage_ve_combinations[comb, "coverage"]
    
    print(sprintf("Iteration: %d, VE: %.1f%%, Coverage: %.1f%%", 
                  i_iter, ve * 100, coverage * 100))
    
    # Inner loop for vaccination scenarios
    for (i in seq_along(target_age_list)) {
      
      target <- target_age_list[[i]]
      
      print(sprintf("Iteration %d, Scenario %d: Target Age = %s", 
                    i_iter, i, toString(target)))
      
      # Simulate the SIRV model
      sim_result <- sirv_sim_coverageSwitch(
        T = 52, 
        A = 18,
        N = N_region1,
        r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        delay = 1,                   # Vaccination starts from Week x
        VE_block = ve,               # Vaccine efficacy
        target_age = target,         # Targeting specific age groups by scenario
        vaccine_end = delay + vaccine_duration_weeks,
        coverage_threshold = 1,
        weekly_supply_phase1 = sum(N) * coverage / vaccine_duration_weeks, 
        weekly_supply_phase2 = sum(N) * coverage / vaccine_duration_weeks
      )
      
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = paste0("Scenario_", i, "_VE_", ve, 
                            "_Coverage_", coverage),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs_v3[[length(summary_dfs_v3) + 1]] <- sim_df
      
      # Store results for this combination
      iteration_result[[paste0("Scenario_", i, "_VE_", ve, "_Coverage_", coverage)]] <- sim_result
    }
  }
  
  # Store all results from this iteration
  all_results_v3[[i_iter]] <- iteration_result
}

# Extract unique scenarios from the nested structure
unique_scenarios <- unique(unlist(lapply(all_results_v3, function(iteration) names(iteration))))

# Initialize a list to store results for all scenarios
ve_summary_list <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_scenario_data(scenario, all_results_v3)
  
  # Compute median across iterations
  median_summary <- compute_median(scenario_data)
  
  # Store the summarized data
  ve_summary_list[[scenario]] <- median_summary
}

final_summ_ve <- lapply(ve_summary_list, function(df) {
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
      VE       = as.numeric(sub("VE_", "", stringr::str_extract(Scenario, "VE_\\d+(\\.\\d+)?"))),
      Coverage = as.numeric(sub("Coverage_", "", stringr::str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
    )
})

final_summ_s1_ve <- final_summ_ve[startsWith(names(final_summ_ve), "Scenario_1_")]

total_impact_s1_ve <- lapply(final_summ_s1_ve, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
              pre_vacc = sum(pre_vacc),
              tot_impact = tot_diff / pre_vacc * 100
              )
})

final_total_impact_s1_ve <- do.call(rbind, total_impact_s1_ve)

summary_table_s1_ve <- final_total_impact_s1_ve %>%
  # Extract Delay and Duration using regex
  mutate(
    VE       = as.numeric(sub("VE_", "", stringr::str_extract(Scenario, "VE_\\d+(\\.\\d+)?"))),
    Coverage = as.numeric(sub("Coverage_", "", str_extract(Scenario, "Coverage_\\d+(\\.\\d+)?")))
  )

summary_table_s1_ve$tot_impact <- as.numeric(summary_table_s1_ve$tot_impact)
summary_table_s1_ve$tot_impact <- summary_table_s1_ve$tot_impact * 100

ggplot(data = summary_table_s1_ve, aes(x = VE, y = Coverage, z = tot_impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5)
