sirv_sim_coverageSwitch_Fatal <- function(
    # Model dimensions
  T, A,
  N,                       # population by age, length A
  r,                       # aging rates, length A
  
  # Epidemic parameters
  base_beta,               # length T
  I0_draw,                 # initial infected, length A
  R0,                      # initial underlying immunity 
  rho,                     # detection probability
  gamma,                   # recovery rate
  
  # Vaccination parameters
  delay,                   # day (index) to start Phase 1
  coverage_threshold,      # e.g. 0.7 (70%)
  
  # target_age can be:
  #   1) a 0/1 vector of length A (e.g. c(0,0,1,1,...) )
  #   2) or a vector of actual indices (c(3,4,5)) -- but here we assume 0/1
  target_age,              # vector of length A, with 1 if that age is prioritized
  total_coverage,
  weekly_delivery_speed,
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
  raw_allocation <- matrix(0, nrow = A, ncol = T)
  coverage_frac <- matrix(0, nrow = A, ncol = T)
  coverage_target <- rep(0, T)
  
  # vaccine params 
  total_supply <- sum(N)*total_coverage
  
  # target pop 
  target_pop <- sum(N[target_age == 1])
  
  # total vaccine supply by targeted age group
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # weekly supply by targeted age group 
  weekly_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(total_available_supply_age[a] > 0) {
      weekly_supply_age[a] <- total_available_supply_age[a] * weekly_delivery_speed
    } else {
      weekly_supply_age[a] <- 0
    }
  }
  
  # total duration of vaccines (this can vary by delivery speed so we are tracking here)
  total_vaccine_duration_age <- numeric(A)
  
  for(a in seq_len(A)) {
    if(weekly_supply_age[a] > 0) {
      tmp_weeks <- ceiling(total_available_supply_age[a] / weekly_supply_age[a])
    } else {
      tmp_weeks <- 0 
    }
    if (tmp_weeks > 52) {
      tmp_weeks <- 52
    }
    total_vaccine_duration_age[a] <- tmp_weeks
  }

  vaccine_end_age <- total_vaccine_duration_age + 1
  total_supply_age <- numeric(A)
  unused_supply_age <- total_available_supply_age  # Start with all supply unused

  # Initial conditions at t=1
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I0_draw[a] - R_[a, 1]
    V[a, 1]  <- 0
    IV[a, 1] <- 0
    RV[a, 1] <- 0
    IFatal[a, 1] <- 0
    IVFatal[a, 1] <- 0 
    age_stratified_cases[a, 1] <- rho * (I[a,1] + IV[a,1])
    IFatal_detect[a, 1] <- rho * IFatal[a,1]
    IVFatal_detect[a, 1] <- rho * IVFatal[a,1]
    coverage_frac[a, 1] <- (V[a, 1] + RV[a, 1]) / N[a]
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
      if(t_ < delay || t_ > vaccine_end_age[a]) {
        # Outside vaccination window
        raw_allocation[a, t_] <- 0
        vacc_eff <- 0
        
      } else {
        # We are within [delay, vaccine_end]
        if(!coverage_switch) {
          # PHASE 1: target only
          if(!no_target && target_age[a] == 1) {
            doses_planned <- weekly_supply_age[a]
          } else {
            doses_planned <- 0
          }
        } else {
          # PHASE 2: vaccinate all ages
          doses_planned <- weekly_supply_age[a]
        }
        
        # --------------------------------------
        # CLAMP the planned doses so we don't exceed:
        #   1) unused supply in that age group
        #   2) the # of susceptible individuals
        # --------------------------------------
        if(doses_planned > 0) {
          # The max we can allocate this step:
          # - If S[a,t_-1] is 0, there's nobody left to vaccinate
          # - If unused_supply_age[a] is 0, supply is depleted
          # 
          # Since each step is 1 week, "doses_planned" is the actual # of doses,
          # not a fraction.
          # 
          # If your code is interpreting "weekly_supply_age[a]" as
          # "fraction of the population," then adjust accordingly.
          
          # Actual feasible doses (integer or real):
          doses_actual <- min(
            doses_planned,
            unused_supply_age[a],
            S[a, t_-1]  # you cannot vaccinate more than the susceptible count
          )
          
          raw_allocation[a, t_] <- doses_actual
          
          # Vaccination rate per total population in this age
          if(N[a] > 0) {
            vacc_eff <- doses_actual / N[a]
          } else {
            vacc_eff <- 0
          }
          
          # Update supply usage
          total_supply_age[a]  <- total_supply_age[a]  + doses_actual
          unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
          
        } else {
          raw_allocation[a, t_] <- 0
          vacc_eff <- 0
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
    # main compartments
    S = S, I = I, R = R_,
    V = V, IV = IV, RV = RV,
    IFatal = IFatal, IVFatal = IVFatal,
    
    # epidemiological metrics
    age_stratified_cases = age_stratified_cases,
    IFatal_detect = IFatal_detect,
    IVFatal_detect = IVFatal_detect,
    
    # vaccination logic
    raw_allocation = raw_allocation,
    coverage_frac  = coverage_frac,      # coverage fraction per age over time
    coverage_target = coverage_target,   # coverage in the target group over time
    
    # phase-switch logic
    coverage_switch = coverage_switch,
    switch_day = switch_day,
    
    # final supply/tracking
    total_supply = total_supply,
    total_available_supply_age = total_available_supply_age,
    weekly_supply_age = weekly_supply_age,
    total_vaccine_duration_age = total_vaccine_duration_age,
    vaccine_end_age = vaccine_end_age,
    
    total_supply_age = total_supply_age,
    unused_supply_age = unused_supply_age
  )
}

## version 2. leaky vaccine scenario
sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,
  N,                       # population by age, length A
  r,                       # aging rates, length A
  
  # Epidemic parameters
  base_beta,               # length T
  I0_draw,                 # initial infected, length A
  R0,                      # initial underlying immunity 
  rho,                     # detection probability
  gamma,                   # recovery rate
  
  # Vaccination parameters
  delay,                   # day (index) to start Phase 1
  coverage_threshold,      # e.g. 0.7 (70%)
  
  # target_age can be:
  #   1) a 0/1 vector of length A (e.g. c(0,0,1,1,...) )
  #   2) or a vector of actual indices (c(3,4,5)) -- but here we assume 0/1
  target_age,              # vector of length A, with 1 if that age is prioritized
  total_coverage,
  weekly_delivery_speed,
  VE_block                 # vaccine efficacy (infection blocking, "leaky")
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
  
  # For reporting, coverage, etc.
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  raw_allocation       <- matrix(0, nrow = A, ncol = T)
  coverage_frac        <- matrix(0, nrow = A, ncol = T)
  coverage_target      <- rep(0, T)
  
  # Vaccine supply parameters
  total_supply <- sum(N) * total_coverage
  
  # Identify target population
  target_pop <- sum(N[target_age == 1])
  
  # total vaccine supply allocated to each age group (among target ages)
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # weekly supply by targeted age group
  # weekly delivery speed as a fraction of target population per week
  weekly_dose_total <- total_supply * weekly_delivery_speed
  weekly_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(total_available_supply_age[a] > 0) {
      weekly_supply_age[a] <- weekly_dose_total * (N[a] / target_pop)
    } else {
      weekly_supply_age[a] <- 0
    }
  }
  
  # total duration (in weeks) we have supply, capped at 52 for convenience
  total_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(weekly_supply_age[a] > 0) {
      tmp_weeks <- ceiling(total_available_supply_age[a] / weekly_supply_age[a])
      tmp_weeks <- min(tmp_weeks, 52)
      total_vaccine_duration_age[a] <- tmp_weeks
    } else {
      total_vaccine_duration_age[a] <- 0
    }
  }
  
  vaccine_end_age   <- total_vaccine_duration_age + 1
  total_supply_age  <- numeric(A)
  unused_supply_age <- total_available_supply_age  # Start with all supply unused
  
  # Initial conditions at t=1
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I0_draw[a] - R_[a, 1]
    
    V[a, 1]  <- 0
    IV[a, 1] <- 0
    RV[a, 1] <- 0
    
    age_stratified_cases[a, 1] <- rho * (I[a,1] + IV[a,1])
    coverage_frac[a, 1]        <- (V[a, 1] + RV[a, 1]) / N[a]
  }
  
  # Identify which ages are target
  target_indices <- which(target_age == 1)  
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # coverage_switch: has coverage threshold been met?
  coverage_switch <- FALSE
  switch_day      <- NA  # if we never switch
  
  if(!no_target) {
    coverage_target[1] <- sum(V[target_indices, 1] + RV[target_indices, 1]) / target_pop
  } else {
    coverage_target[1] <- 0
  }
  
  # ---------------------------
  # 2) Forward simulation t=2..T
  # ---------------------------
  
  phi_vec <- numeric(T)
  
  for(t_ in 2:T) {
    
    # 2.A) Check coverage in target group (end of t_-1) to see if we switch to Phase 2
    if(!no_target) {
      target_vaccinated_prev <- sum(V[target_indices, t_-1] + RV[target_indices, t_-1])
      coverage_now <- target_vaccinated_prev / target_pop
      
      if(!coverage_switch && coverage_now >= coverage_threshold) {
        coverage_switch <- TRUE
        switch_day      <- t_
      }
    } else {
      coverage_now <- 0
    }
    
    # 2.B) Force of infection
    infected_prev <- sum(I[, t_-1] + IV[, t_-1])
    region_pop    <- sum(N)  
    phi           <- base_beta[t_-1] * (infected_prev / region_pop)
    
    phi_vec[t_] <- phi
    
    # 2.C) Update compartments by age
    for(a in seq_len(A)) {
      # Decide how many doses to allocate in this time step
      if(t_ < (delay + 2) || t_ > vaccine_end_age[a]) {
        # Outside vaccination window
        raw_allocation[a, t_] <- 0
        vacc_eff <- 0
        
      } else {
        # Within [delay, vaccine_end]
        if(!coverage_switch) {
          # PHASE 1: vaccinate only target ages
          if(!no_target && target_age[a] == 1) {
            doses_planned <- weekly_supply_age[a]
          } else {
            doses_planned <- 0
          }
        } else {
          # PHASE 2: vaccinate all ages
          doses_planned <- weekly_supply_age[a]
        }
        
        # Clamp planned doses so we don't exceed supply or susceptible
        if(doses_planned > 0) {
          
          susceptibles_now <- S[a, t_-1]
          remaining_supply <- unused_supply_age[a]
          
          feasible <- min(
            doses_planned, remaining_supply, susceptibles_now
          )
          
          feasible <- max(0, feasible)
          
          raw_allocation[a, t_] <- feasible
          vacc_eff              <- if(N[a] > 0) feasible / N[a] else 0
          
          # Update supply usage
          total_supply_age[a]  <- total_supply_age[a] + feasible
          unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
          
        } else {
          raw_allocation[a, t_] <- 0
          vacc_eff <- 0
        }
      }
      
      # SIRV transitions
      if(a == 1) {
        # age 1
        S[a,t_] <- S[a,t_-1] - phi*S[a,t_-1] -
          vacc_eff*S[a,t_-1] - r[a]*S[a,t_-1]
        
        I[a,t_] <- I[a,t_-1] + phi*S[a,t_-1] -
          gamma*I[a,t_-1] - r[a]*I[a,t_-1]
        
        R_[a,t_] <- R_[a,t_-1] + gamma*I[a,t_-1] - r[a]*R_[a,t_-1]
        
        V[a,t_] <- V[a,t_-1] + vacc_eff*S[a,t_-1] -
          (1 - VE_block)*phi*V[a,t_-1] - r[a]*V[a,t_-1]
        
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block)*phi*V[a,t_-1] -
          gamma*IV[a,t_-1] - r[a]*IV[a,t_-1]
        
        RV[a,t_] <- RV[a,t_-1] + gamma*IV[a,t_-1] - r[a]*RV[a,t_-1]
        
      } else {
        # age 2..A
        S[a,t_] <- S[a,t_-1] - phi*S[a,t_-1] -
          vacc_eff*S[a,t_-1] - r[a]*S[a,t_-1] +
          r[a-1]*S[a-1,t_-1]
        
        I[a,t_] <- I[a,t_-1] + phi*S[a,t_-1] -
          gamma*I[a,t_-1] - r[a]*I[a,t_-1] +
          r[a-1]*I[a-1,t_-1]
        
        R_[a,t_] <- R_[a,t_-1] + gamma*I[a,t_-1] -
          r[a]*R_[a,t_-1] + r[a-1]*R_[a-1,t_-1]
        
        V[a,t_] <- V[a,t_-1] + vacc_eff*S[a,t_-1] -
          (1 - VE_block)*phi*V[a,t_-1] - r[a]*V[a,t_-1] +
          r[a-1]*V[a-1,t_-1]
        
        IV[a,t_] <- IV[a,t_-1] + (1 - VE_block)*phi*V[a,t_-1] -
          gamma*IV[a,t_-1] - r[a]*IV[a,t_-1] +
          r[a-1]*IV[a-1,t_-1]
        
        RV[a,t_] <- RV[a,t_-1] + gamma*IV[a,t_-1] -
          r[a]*RV[a,t_-1] + r[a-1]*RV[a-1,t_-1]
      }
      
      # Clamp at 0
      S[a,t_]  <- pmax(0, S[a,t_])
      I[a,t_]  <- pmax(0, I[a,t_])
      R_[a,t_] <- pmax(0, R_[a,t_])
      V[a,t_]  <- pmax(0, V[a,t_])
      IV[a,t_] <- pmax(0, IV[a,t_])
      RV[a,t_] <- pmax(0, RV[a,t_])
    }
    
    # 2.D) After updating compartments, record coverage & reported cases
    for(a in seq_len(A)) {
      coverage_frac[a, t_]        <- (V[a,t_] + RV[a,t_]) / N[a]
      age_stratified_cases[a,t_]  <- rho * (I[a,t_] + IV[a,t_])
    }
    if(!no_target) {
      coverage_target[t_] <- sum(V[target_indices, t_] + RV[target_indices, t_]) / target_pop
    } else {
      coverage_target[t_] <- 0
    }
  }
  
  # Return compartments and additional tracking
  list(
    # main compartments
    S = S, I = I, R = R_,
    V = V, IV = IV, RV = RV,
    
    # epidemiological metrics
    age_stratified_cases = age_stratified_cases,
    phi = phi_vec,
    
    # vaccination logic
    raw_allocation = raw_allocation,
    coverage_frac  = coverage_frac,      # coverage fraction per age over time
    coverage_target = coverage_target,   # coverage in the target group over time
    
    # phase-switch logic
    coverage_switch = coverage_switch,
    switch_day = switch_day,
    
    # final supply/tracking
    total_supply = total_supply,
    total_available_supply_age = total_available_supply_age,
    weekly_supply_age = weekly_supply_age,
    total_vaccine_duration_age = total_vaccine_duration_age,
    vaccine_end_age = vaccine_end_age,
    
    total_supply_age = total_supply_age,
    unused_supply_age = unused_supply_age
  )
}

## version 2. all-or-nothing vaccine scenario
sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,
  N,                       # population by age, length A
  r,                       # aging rates, length A
  
  # Epidemic parameters
  base_beta,               # length T
  I0_draw,                 # initial infected, length A
  R0,                      # initial underlying immunity 
  rho,                     # detection probability
  gamma,                   # recovery rate
  
  # Vaccination parameters
  delay,                   # day (index) to start Phase 1
  coverage_threshold,      # e.g. 0.7 (70%)
  
  # target_age can be:
  #   1) a 0/1 vector of length A (e.g. c(0,0,1,1,...) )
  #   2) or a vector of actual indices (c(3,4,5)) -- but here we assume 0/1
  target_age,              # vector of length A, with 1 if that age is prioritized
  total_coverage,
  weekly_delivery_speed,
  VE_block                 # vaccine efficacy (infection blocking, "leaky")
){
  # -------------------------
  # 1) Initialize compartments
  # -------------------------
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)

  # For reporting, coverage, etc.
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  raw_allocation       <- matrix(0, nrow = A, ncol = T)
  coverage_frac        <- matrix(0, nrow = A, ncol = T)
  coverage_target      <- rep(0, T)
  
  # Vaccine supply parameters
  total_supply <- sum(N) * total_coverage
  
  # Identify target population
  target_pop <- sum(N[target_age == 1])
  
  # total vaccine supply allocated to each age group (among target ages)
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # weekly supply by targeted age group
  # weekly delivery speed as a fraction of target population per week
  weekly_dose_total <- total_supply * weekly_delivery_speed
  weekly_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(total_available_supply_age[a] > 0) {
      weekly_supply_age[a] <- weekly_dose_total * (N[a] / target_pop)
    } else {
      weekly_supply_age[a] <- 0
    }
  }
  
  # total duration (in weeks) we have supply, capped at 52 for convenience
  total_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(weekly_supply_age[a] > 0) {
      tmp_weeks <- ceiling(total_available_supply_age[a] / weekly_supply_age[a])
      tmp_weeks <- min(tmp_weeks, 52)
      total_vaccine_duration_age[a] <- tmp_weeks
    } else {
      total_vaccine_duration_age[a] <- 0
    }
  }
  
  vaccine_end_age   <- total_vaccine_duration_age + 1
  total_supply_age  <- numeric(A)
  unused_supply_age <- total_available_supply_age  # Start with all supply unused
  
  # Initial conditions at t=1
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I0_draw[a] - R_[a, 1]
    
    V[a, 1]  <- 0

    age_stratified_cases[a, 1] <- rho * (I[a,1])
    coverage_frac[a, 1]        <- (V[a, 1]) / N[a]
  }
  
  # Identify which ages are target
  target_indices <- which(target_age == 1)  
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # coverage_switch: has coverage threshold been met?
  coverage_switch <- FALSE
  switch_day      <- NA  # if we never switch
  
  if(!no_target) {
    coverage_target[1] <- sum(V[target_indices, 1]) / target_pop
  } else {
    coverage_target[1] <- 0
  }
  
  # ---------------------------
  # 2) Forward simulation t=2..T
  # ---------------------------
  
  phi_vec <- numeric(T)
  
  for(t_ in 2:T) {
    
    # 2.A) Check coverage in target group (end of t_-1) to see if we switch to Phase 2
    if(!no_target) {
      target_vaccinated_prev <- sum(V[target_indices, t_-1])
      coverage_now <- target_vaccinated_prev / target_pop
      
      if(!coverage_switch && coverage_now >= coverage_threshold) {
        coverage_switch <- TRUE
        switch_day      <- t_
      }
    } else {
      coverage_now <- 0
    }
    
    # 2.B) Force of infection
    infected_prev <- sum(I[, t_-1])
    region_pop    <- sum(N)  
    phi           <- base_beta[t_-1] * (infected_prev / region_pop)
    
    phi_vec[t_] <- phi
    
    # 2.C) Update compartments by age
    for(a in seq_len(A)) {
      # Decide how many doses to allocate in this time step
      if(t_ < (delay + 2) || t_ > vaccine_end_age[a]) {
        # Outside vaccination window
        raw_allocation[a, t_] <- 0
        vacc_eff <- 0
        
      } else {
        # Within [delay, vaccine_end]
        if(!coverage_switch) {
          # PHASE 1: vaccinate only target ages
          if(!no_target && target_age[a] == 1) {
            doses_planned <- weekly_supply_age[a]
          } else {
            doses_planned <- 0
          }
        } else {
          # PHASE 2: vaccinate all ages
          doses_planned <- weekly_supply_age[a]
        }
        
        # Clamp planned doses so we don't exceed supply or susceptible
          if(doses_planned > 0) {
          #feasible <- min(
          #  doses_planned
          #)
          
          susceptibles_now <- S[a, t_-1]
          remaining_supply <- unused_supply_age[a]
          
          feasible <- min(doses_planned, remaining_supply, susceptibles_now)
          feasible <- max(0, feasible)
          
          raw_allocation[a, t_] <- feasible
          vacc_eff              <- if(N[a] > 0) feasible / N[a] else 0
          
          # Update supply usage
          total_supply_age[a]  <- total_supply_age[a] + feasible
          unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
          
        } else {
          raw_allocation[a, t_] <- 0
          vacc_eff <- 0
        }
      }
      
      # SIRV transitions
      if(a == 1) {
        # age 1
        S[a,t_] <- S[a,t_-1] - phi*S[a,t_-1] -
          VE_block * vacc_eff * S[a,t_-1] - r[a]*S[a,t_-1]
        
        V[a,t_] <- V[a,t_-1] + VE_block * vacc_eff * S[a, t_-1] -
          r[a] * V[a,t_-1]
        
        I[a,t_] <- I[a,t_-1] + phi*S[a,t_-1] -
          gamma*I[a,t_-1] - r[a]*I[a,t_-1]
        
        R_[a,t_] <- R_[a,t_-1] + gamma*I[a,t_-1] - r[a]*R_[a,t_-1]
        
      } else {
        # age 2..A
        S[a,t_] <- S[a,t_-1] - phi*S[a,t_-1] -
          VE_block * vacc_eff * S[a,t_-1] - r[a]*S[a,t_-1] +
          r[a-1]*S[a-1,t_-1]
        
        V[a,t_] <- V[a,t_-1] + VE_block * vacc_eff * S[a, t_-1] - r[a]*V[a,t_-1] +
          r[a-1]*V[a-1,t_-1]
        
        I[a,t_] <- I[a,t_-1] + phi*S[a,t_-1] -
          gamma*I[a,t_-1] - r[a]*I[a,t_-1] +
          r[a-1]*I[a-1,t_-1]
        
        R_[a,t_] <- R_[a,t_-1] + gamma*I[a,t_-1] -
          r[a]*R_[a,t_-1] + r[a-1]*R_[a-1,t_-1]
        
      }
      
      # Clamp at 0
      S[a,t_]  <- pmax(0, S[a,t_])
      I[a,t_]  <- pmax(0, I[a,t_])
      R_[a,t_] <- pmax(0, R_[a,t_])
      V[a,t_]  <- pmax(0, V[a,t_])
    }
    
    # 2.D) After updating compartments, record coverage & reported cases
    for(a in seq_len(A)) {
      coverage_frac[a, t_]        <- (V[a,t_]) / N[a]
      age_stratified_cases[a,t_]  <- rho * (I[a,t_])
    }
    if(!no_target) {
      coverage_target[t_] <- sum(V[target_indices, t_]) / target_pop
    } else {
      coverage_target[t_] <- 0
    }
  }
  
  # Return compartments and additional tracking
  list(
    # main compartments
    S = S, I = I, R = R_,
    V = V,
    
    # epidemiological metrics
    age_stratified_cases = age_stratified_cases,
    phi = phi_vec,
    
    # vaccination logic
    raw_allocation = raw_allocation,
    coverage_frac  = coverage_frac,      # coverage fraction per age over time
    coverage_target = coverage_target,   # coverage in the target group over time
    
    # phase-switch logic
    coverage_switch = coverage_switch,
    switch_day = switch_day,
    
    # final supply/tracking
    total_supply = total_supply,
    total_available_supply_age = total_available_supply_age,
    weekly_supply_age = weekly_supply_age,
    total_vaccine_duration_age = total_vaccine_duration_age,
    vaccine_end_age = vaccine_end_age,
    
    total_supply_age = total_supply_age,
    unused_supply_age = unused_supply_age
  )
}




### updated
sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,          # T = number of weeks, A = number of age groups
  N,             # population by age
  r,             # aging rates
  
  # Epidemic parameters
  base_beta,     # length T
  I0_draw,       # initial infected by age
  R0,            # initial underlying immunity fraction by age
  rho,           # detection probability
  gamma,         # recovery rate
  
  # Vaccination parameters
  delay,                 # the week index to START physical vaccination
  coverage_threshold,    # e.g., 0.7 => 70%
  target_age,            # 0/1 vector of length A
  total_coverage,        # fraction of total population eventually vaccinated
  weekly_delivery_speed, # fraction of total supply used per week
  VE_block,              # vaccine efficacy
  
  # NEW: toggle coverage switch
  use_coverage_switch = FALSE  # set FALSE to disable Phase 2 reallocation
){
  # ------------------------------------------------------------------------
  # 1) Initialize compartments and data structures
  # ------------------------------------------------------------------------
  
  # S, I, R, V each dimension [A x T]
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  
  # Delayed vaccination storage: how many got vaccinated at t_, but immune 2 weeks later
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # For reporting
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)  # weekly detected cases by age
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)  # how many doses allocated each week
  coverage_frac        <- matrix(0, nrow = A, ncol = T)  # fraction of each age group in V
  coverage_target      <- numeric(T)                     # coverage fraction in target group
  phi_vec              <- numeric(T)                     # force of infection each week
  
  # Identify target group
  target_indices <- which(target_age == 1)
  target_pop     <- sum(N[target_indices])
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # Total supply (when we define coverage as a fraction of total pop, or fraction of targeted pop the definition changes)
  total_supply <- target_pop  * total_coverage
  
  # ------------------------------------------------------------------------
  # 1A) Phase 1 supply: target only
  # ------------------------------------------------------------------------
  
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  # In Phase 1, only target ages get a share of weekly_dose_total
  weekly_supply_age <- numeric(A)
  if(!no_target) {
    for(a in seq_len(A)) {
      if(target_age[a] == 1) {
        weekly_supply_age[a] <- weekly_dose_total * (N[a] / target_pop)
      } else {
        weekly_supply_age[a] <- 0
      }
    }
  } else {
    weekly_supply_age <- rep(0, A)
  }
  
  # total_available_supply_age: maximum doses that can be allocated to each age a
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1 && !no_target) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # Track usage vs leftover supply
  total_supply_age  <- numeric(A) # how many doses have actually been used for age [a]
  unused_supply_age <- total_available_supply_age
  
  # coverage_switch = TRUE once coverage_threshold is reached in target group
  coverage_switch <- FALSE
  switch_day      <- NA
  
  # ------------------------------------------------------------------------
  # 1B) Vaccination duration tracking
  # ------------------------------------------------------------------------
  vacc_start_week <- rep(NA, A)
  vacc_end_week   <- rep(NA, A)
  
  # ------------------------------------------------------------------------
  # 1C) Set initial conditions at t=1
  # ------------------------------------------------------------------------
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a,1] - R_[a,1]
    V[a, 1]  <- 0
    
    age_stratified_cases[a, 1] <- rho * I[a,1]
    coverage_frac[a, 1]       <- V[a, 1] / N[a]
  }
  
  coverage_target[1] <- if(!no_target) 0 else 0
  
  # ------------------------------------------------------------------------
  # 2) Main loop, t_ = 2..T
  # ------------------------------------------------------------------------
  for(t_ in 2:T) {
    
    # 2A) Coverage in target group (cumulative # of doses physically given)
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:(t_-1)])
      coverage_now           <- target_vacc_cumulative / target_pop
      
      if(use_coverage_switch) {
        if(!coverage_switch && coverage_now >= coverage_threshold) {
          coverage_switch <- TRUE
          switch_day      <- t_
          
          # Reallocate leftover supply to ALL ages if coverage_switch
          # vaccinate everyone (including both targeted + untargeted) who is still susceptible
          leftover_total <- sum(unused_supply_age)
          if(leftover_total > 0) {
            total_pop <- sum(N)
            for(a in seq_len(A)) {
              share_a <- N[a] / total_pop
              total_available_supply_age[a] <- total_supply_age[a] + leftover_total * share_a
              unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
            }
          }
          
          # Now define weekly_supply_age for ALL ages
          for(a in seq_len(A)) {
            weekly_supply_age[a] <- weekly_dose_total * (N[a] / sum(N))
          }
        }
      }
      
    } else {
      coverage_now <- 0
    }
    
    # 2B) Force of infection
    infected_prev <- sum(I[, t_-1])
    region_pop    <- sum(N)
    phi           <- base_beta[t_-1] * (infected_prev / region_pop)
    phi_vec[t_]   <- phi
    
    # ------------------------------------------------
    # 2C) Vaccination logic: feasible allocations
    # ------------------------------------------------
    for(a in seq_len(A)) {
      # Start vaccinating at exactly t_=delay (not delay+1).
      # => if t_ < delay, no vaccination
      if(t_ < delay) {
        doses_planned <- 0
      } else {
        doses_planned <- weekly_supply_age[a]
      }
      
      susceptibles_now <- S[a, t_-1]
      remaining_supply <- unused_supply_age[a]
      
      # determine how many vaccines can actually be administered this week 
      feasible <- min(doses_planned, remaining_supply, susceptibles_now)
      feasible <- max(feasible, 0)
      
      raw_allocation_age[a, t_] <- feasible
      
      if(feasible > 0) {
        # Mark start if first time
        if(is.na(vacc_start_week[a])) {
          vacc_start_week[a] <- t_
        }
        # Update last time
        vacc_end_week[a] <- t_
        
        # Use supply
        total_supply_age[a]  <- total_supply_age[a] + feasible
        unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
        
        # Store them in vacc_delayed => immune at (t_ + 2)
        vacc_delayed[a, t_] <- feasible
      }
    }
    
    # ------------------------------------------------
    # 2D) S-I-R transitions (no new V yet)
    # ------------------------------------------------
    for(a in seq_len(A)) {
      if(a == 1) {
        S[a, t_] <- S[a, t_-1] - phi*S[a, t_-1] - r[a]*S[a, t_-1]
        I[a, t_] <- I[a, t_-1] + phi*S[a, t_-1] - gamma*I[a, t_-1] - r[a]*I[a, t_-1]
        R_[a, t_] <- R_[a, t_-1] + gamma*I[a, t_-1] - r[a]*R_[a, t_-1]
        V[a, t_] <- V[a, t_-1] - r[a]*V[a, t_-1]
        
      } else {
        S[a, t_] <- S[a, t_-1] - phi*S[a, t_-1] - r[a]*S[a, t_-1] + r[a-1]*S[a-1,t_-1]
        I[a, t_] <- I[a, t_-1] + phi*S[a, t_-1] - gamma*I[a, t_-1] - r[a]*I[a, t_-1] + r[a-1]*I[a-1,t_-1]
        R_[a, t_] <- R_[a, t_-1] + gamma*I[a, t_-1] - r[a]*R_[a, t_-1] + r[a-1]*R_[a-1,t_-1]
        V[a, t_] <- V[a, t_-1] - r[a]*V[a, t_-1] + r[a-1]*V[a-1,t_-1]
      }
      
      # clamp at 0
      S[a, t_]  <- pmax(0, S[a, t_])
      I[a, t_]  <- pmax(0, I[a, t_])
      R_[a, t_] <- pmax(0, R_[a, t_])
      V[a, t_]  <- pmax(0, V[a, t_])
    }
    
    # ------------------------------------------------
    # 2E) Apply 2-week immunity delay
    # ------------------------------------------------
    # Everyone who received vaccine in week (t_-2) => become immune at t_
    if(t_ - 2 >= 1) {
      t_delay <- t_ - 2
      for(a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t_delay]
        # clamp if S[a, t_] < newly_effective
        actually_moving <- min(S[a, t_], newly_effective)
        
        S[a, t_] <- S[a, t_] - actually_moving
        V[a, t_] <- V[a, t_] + actually_moving
      }
    }
    
    # ------------------------------------------------
    # 2F) Post-update coverage & cases
    # ------------------------------------------------
    for(a in seq_len(A)) {
      coverage_frac[a, t_]       <- V[a, t_] / N[a]
      age_stratified_cases[a, t_] <- rho * I[a, t_]
    }
    
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:t_])
      coverage_target[t_]    <- target_vacc_cumulative / target_pop
    } else {
      coverage_target[t_] <- 0
    }
  }
  
  # ------------------------------------------------------------------------
  # 3) Compute total vaccination duration by age
  # ------------------------------------------------------------------------
  real_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(!is.na(vacc_start_week[a]) && !is.na(vacc_end_week[a])) {
      real_vaccine_duration_age[a] <- 
        vacc_end_week[a] - vacc_start_week[a] + 1
    } else {
      real_vaccine_duration_age[a] <- 0
    }
  }
  
  # ------------------------------------------------------------------------
  # 4) Return results
  # ------------------------------------------------------------------------
  return(list(
    # compartments
    S = S,
    I = I,
    R = R_,
    V = V,
    
    # 2-week delayed vaccination
    vacc_delayed = vacc_delayed,
    
    # supply usage
    total_supply               = total_supply,
    total_available_supply_age = total_available_supply_age,
    total_supply_age           = total_supply_age,
    unused_supply_age          = unused_supply_age,
    
    # coverage switch
    coverage_switch = coverage_switch,
    switch_day      = switch_day,
    
    # epidemiological reporting
    age_stratified_cases = age_stratified_cases,
    coverage_frac        = coverage_frac,
    coverage_target      = coverage_target,
    phi = phi_vec,
    
    # vaccination window
    vacc_start_week            = vacc_start_week,
    vacc_end_week              = vacc_end_week,
    real_vaccine_duration_age  = real_vaccine_duration_age,
    
    # raw allocation by age each week
    raw_allocation_age = raw_allocation_age
  ))
}


### v2.2
sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,          # T = number of weeks, A = number of age groups
  N,             # population by age
  r,             # aging rates
  
  # Epidemic parameters
  base_beta,     # length T
  I0_draw,       # initial infected by age
  R0,            # initial underlying immunity fraction by age
  rho,           # detection probability
  gamma,         # recovery rate
  
  # Vaccination parameters
  delay,                 # the week index to START physical vaccination
  coverage_threshold,    # e.g., 0.7 => 70%
  target_age,            # 0/1 vector of length A
  total_coverage,        # fraction of total population eventually vaccinated
  weekly_delivery_speed, # fraction of total supply used per week
  VE_block,              # vaccine efficacy
  
  # NEW: toggle coverage switch
  use_coverage_switch = FALSE  # set FALSE to disable Phase 2 reallocation
){
  # ------------------------------------------------------------------------
  # 1) Initialize compartments and data structures
  # ------------------------------------------------------------------------
  
  # S, I, R, V each dimension [A x T]
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  
  # Delayed vaccination storage: doses administered at time t that become effective 2 weeks later.
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # For reporting
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)  # weekly detected cases by age
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)  # how many doses allocated each week
  coverage_frac        <- matrix(0, nrow = A, ncol = T)  # fraction of each age group in V
  coverage_target      <- numeric(T)                     # coverage fraction in target group
  phi_vec              <- numeric(T)                     # force of infection each week
  
  # Identify target group
  target_indices <- which(target_age == 1)
  target_pop     <- sum(N[target_indices])
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # Total supply (when we define coverage as a fraction of total pop, or fraction of targeted pop the definition changes)
  total_supply <- target_pop  * total_coverage # target pop or sum(N)
  
  # ------------------------------------------------------------------------
  # 1A) Phase 1 supply: target only
  # ------------------------------------------------------------------------
  
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  # In Phase 1, only target ages get a share of weekly_dose_total
  weekly_supply_age <- numeric(A)
  if(!no_target) {
    for(a in seq_len(A)) {
      if(target_age[a] == 1) {
        weekly_supply_age[a] <- weekly_dose_total * (N[a] / target_pop)
      } else {
        weekly_supply_age[a] <- 0
      }
    }
  } else {
    weekly_supply_age <- rep(0, A)
  }
  
  # total_available_supply_age: maximum doses that can be allocated to each age a
  total_available_supply_age <- numeric(A)
  for(a in seq_len(A)) {
    if(target_age[a] == 1 && !no_target) {
      total_available_supply_age[a] <- total_supply * (N[a] / target_pop)
    } else {
      total_available_supply_age[a] <- 0
    }
  }
  
  # Track usage vs leftover supply
  total_supply_age  <- numeric(A) # how many doses have actually been used for age [a]
  unused_supply_age <- total_available_supply_age
  
  # coverage_switch becomes TRUE once coverage_threshold is reached in the target group
  coverage_switch <- FALSE
  switch_day      <- NA
  
  # ------------------------------------------------------------------------
  # 1B) Vaccination duration tracking
  # ------------------------------------------------------------------------
  vacc_start_week <- rep(NA, A)
  vacc_end_week   <- rep(NA, A)
  
  # ------------------------------------------------------------------------
  # 1C) Set initial conditions at t = 1
  # ------------------------------------------------------------------------
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a,1] - R_[a,1]
    V[a, 1]  <- 0
    
    age_stratified_cases[a, 1] <- rho * I[a,1]
    coverage_frac[a, 1]       <- V[a, 1] / N[a]
  }
  
  coverage_target[1] <- 0
  
  # ------------------------------------------------------------------------
  # 2) Main loop, t = 2..T
  # ------------------------------------------------------------------------
  for(t_ in 2:T) {
    
    # ============================================================
    # New Step 1: Apply 2-week vaccination effect BEFORE computing FOI
    # ============================================================
    # Vaccines administered at time (t_ - 2) become effective now.
    if(t_ - 2 >= 1) {
      for(a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t_ - 2]
        # Move at most what is available in susceptibles from last time step.
        actual_move <- min(S[a, t_ - 1], newly_effective)
        S[a, t_ - 1] <- S[a, t_ - 1] - actual_move
        V[a, t_ - 1] <- V[a, t_ - 1] + actual_move
      }
    }
    
    # ============================================================
    # 2A) Update coverage in target group (cumulative doses administered)
    # ============================================================
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:(t_ - 1)])
      coverage_now <- target_vacc_cumulative / target_pop
      
      if(use_coverage_switch) {
        if(!coverage_switch && coverage_now >= coverage_threshold) {
          coverage_switch <- TRUE
          switch_day <- t_
          
          # Compute the leftover supply from phase 1:
          leftover_total <- sum(unused_supply_age)
          
          # Reallocate leftover supply: now, for the non-target groups only.
          non_target_indices <- which(target_age == 0)
          total_non_target_pop <- sum(N[non_target_indices])
          
          if(total_non_target_pop > 0 && leftover_total > 0) {
            # For each non-target age group, set its maximum available supply
            for(a in non_target_indices) {
              total_available_supply_age[a] <- leftover_total * (N[a] / total_non_target_pop)
              unused_supply_age[a] <- total_available_supply_age[a]  # reset usage for these groups
            }
            
            # Now, recalculate the weekly dose total for Phase 2 based on leftover supply:
            new_weekly_dose_total <- leftover_total * weekly_delivery_speed
            
            # Set the weekly supply for non-target groups only:
            for(a in seq_len(A)) {
              if(a %in% non_target_indices) {
                weekly_supply_age[a] <- new_weekly_dose_total * (N[a] / total_non_target_pop)
              } else {
                weekly_supply_age[a] <- 0  # no further allocation to already-targeted groups
              }
            }
          }
        }
      }
      
    } else {
      coverage_now <- 0
    }
    
    # ============================================================
    # 2B) Compute Force of Infection (using updated state from t_-1)
    # ============================================================
    infected_prev <- sum(I[, t_ - 1])
    region_pop <- sum(N)
    phi <- base_beta[t_ - 1] * (infected_prev / region_pop)
    phi_vec[t_] <- phi
    
    # ============================================================
    # 2C) Vaccination allocation for week t_
    # ============================================================
    for(a in seq_len(A)) {
      # Vaccination starts exactly at t = delay.
      if(t_ < delay) {
        doses_planned <- 0
      } else {
        doses_planned <- weekly_supply_age[a]
      }
      
      susceptibles_now <- S[a, t_ - 1]
      remaining_supply <- unused_supply_age[a]
      
      feasible <- min(doses_planned, remaining_supply, susceptibles_now)
      feasible <- max(feasible, 0)
      
      raw_allocation_age[a, t_] <- feasible
      
      if(feasible > 0) {
        if(is.na(vacc_start_week[a])) {
          vacc_start_week[a] <- t_
        }
        vacc_end_week[a] <- t_
        
        total_supply_age[a] <- total_supply_age[a] + feasible
        unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
        
        # Record the doses so that they become effective 2 weeks later.
        vacc_delayed[a, t_] <- feasible
      }
    }
    
    # ============================================================
    # 2D) S-I-R transitions for week t_ (using state from t_-1)
    # ============================================================
    for(a in seq_len(A)) {
      if(a == 1) {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1]
      } else {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1] + r[a - 1] * S[a - 1, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1] + r[a - 1] * I[a - 1, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1] + r[a - 1] * R_[a - 1, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1] + r[a - 1] * V[a - 1, t_ - 1]
      }
      
      # Ensure no compartment goes negative.
      S[a, t_] <- pmax(0, S[a, t_])
      I[a, t_] <- pmax(0, I[a, t_])
      R_[a, t_] <- pmax(0, R_[a, t_])
      V[a, t_] <- pmax(0, V[a, t_])
    }
    
    # ============================================================
    # 2F) Post-update: update coverage and recorded cases for week t_
    # ============================================================
    for(a in seq_len(A)) {
      coverage_frac[a, t_] <- V[a, t_] / N[a]
      age_stratified_cases[a, t_] <- rho * I[a, t_]
    }
    
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:t_])
      coverage_target[t_] <- target_vacc_cumulative / target_pop
    } else {
      coverage_target[t_] <- 0
    }
    
  }  # end of main loop
  
  # ------------------------------------------------------------------------
  # 3) Compute total vaccination duration by age
  # ------------------------------------------------------------------------
  real_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(!is.na(vacc_start_week[a]) && !is.na(vacc_end_week[a])) {
      real_vaccine_duration_age[a] <- vacc_end_week[a] - vacc_start_week[a] + 1
    } else {
      real_vaccine_duration_age[a] <- 0
    }
  }
  
  # ------------------------------------------------------------------------
  # 4) Return results
  # ------------------------------------------------------------------------
  return(list(
    # compartments
    S = S,
    I = I,
    R = R_,
    V = V,
    
    # delayed vaccination (to become effective 2 weeks later)
    vacc_delayed = vacc_delayed,
    
    # supply usage
    total_supply               = total_supply,
    total_available_supply_age = total_available_supply_age,
    total_supply_age           = total_supply_age,
    unused_supply_age          = unused_supply_age,
    
    # coverage switch info
    coverage_switch = coverage_switch,
    switch_day      = switch_day,
    
    # epidemiological reporting
    age_stratified_cases = age_stratified_cases,
    coverage_frac        = coverage_frac,
    coverage_target      = coverage_target,
    phi                  = phi_vec,
    
    # vaccination window info
    vacc_start_week           = vacc_start_week,
    vacc_end_week             = vacc_end_week,
    real_vaccine_duration_age = real_vaccine_duration_age,
    
    # raw allocation by age each week
    raw_allocation_age = raw_allocation_age
  ))
}

# tracking global stockpile
sirv_sim_coverageSwitch <- function(
    # Model dimensions
  T, A,          # T = number of weeks, A = number of age groups
  N,             # population by age (vector of length A)
  r,             # aging rates (vector of length A)
  
  # Epidemic parameters
  base_beta,     # vector of length T
  I0_draw,       # initial infected by age (vector of length A)
  R0,            # initial underlying immunity fraction by age (vector of length A)
  rho,           # detection probability
  gamma,         # recovery rate
  
  # Vaccination parameters
  delay,                 # the week index to start physical vaccination
  coverage_threshold,    # not used in this version
  target_age,            # 0/1 vector of length A (1 = target group)
  total_coverage,        # desired fraction of each target group eventually vaccinated (e.g., 0.6)
  weekly_delivery_speed, # fraction of the remaining global stock delivered per week
  VE_block,              # vaccine efficacy
  
  # Toggle for switching (here we implement reallocation for equal per-capita coverage)
  use_coverage_switch = TRUE
){
  # ------------------------------------------------------------------------
  # 1) Initialize compartments and data structures (all in absolute numbers)
  # ------------------------------------------------------------------------
  S  <- matrix(0, nrow = A, ncol = T)   # Susceptible individuals
  I  <- matrix(0, nrow = A, ncol = T)   # Infected
  R_ <- matrix(0, nrow = A, ncol = T)   # Recovered
  V  <- matrix(0, nrow = A, ncol = T)   # Vaccinated
  
  # Vaccination delay: doses allocated at time t become effective 2 weeks later.
  # Here we record the allocated (pending) doses.
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # Reporting variables
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)
  coverage_frac        <- matrix(0, nrow = A, ncol = T)
  coverage_target      <- numeric(T)
  phi_vec              <- numeric(T)
  
  # Identify target groups
  target_indices <- which(target_age == 1)
  target_pop     <- sum(N[target_indices])
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # ------------------------------------------------------------------------
  # 2) Define Desired Doses per Age Group (for equal per-capita coverage)
  # ------------------------------------------------------------------------
  # We want each target age group to eventually have:
  # desired_doses[a] = total_coverage * N[a]
  desired_doses <- total_coverage * N
  
  # ------------------------------------------------------------------------
  # 3) Initialize Global Stockpile
  # ------------------------------------------------------------------------
  total_supply <- sum(N) * total_coverage  # total doses available (for the entire population)
  global_stock <- total_supply             # initial stock (finite, non-replenished)
  
  global_stock_history <- numeric(T)
  global_stock_history[1] <- global_stock
  
  # ------------------------------------------------------------------------
  # 4) Initialize Initial Allocation for Target Group
  # ------------------------------------------------------------------------
  # For t < delay, no doses are delivered.
  weekly_supply_age <- numeric(A)
  if(!no_target) {
    sum_target_pop <- sum(N[target_indices])
    weekly_dose_total <- global_stock * weekly_delivery_speed
    for(a in target_indices) {
      # Initially, allocate proportional to the target population
      weekly_supply_age[a] <- weekly_dose_total * (N[a] / sum_target_pop)
    }
    for(a in setdiff(seq_len(A), target_indices)) {
      weekly_supply_age[a] <- 0
    }
  } else {
    weekly_supply_age <- rep(0, A)
  }
  
  # Also track per-age maximum available (initial allocation)
  total_available_supply_age <- numeric(A)
  if(!no_target) {
    for(a in target_indices) {
      total_available_supply_age[a] <- global_stock * (N[a] / sum_target_pop)
    }
  } else {
    total_available_supply_age <- rep(0, A)
  }
  
  total_supply_age  <- numeric(A)   # cumulative doses allocated per age group
  unused_supply_age <- total_available_supply_age  # initially, unused equals the allocated amount
  
  # We ignore switching beyond Phase 1.
  coverage_switch <- FALSE
  switch_day <- NA
  
  # ------------------------------------------------------------------------
  # 5) Vaccination Duration Tracking
  # ------------------------------------------------------------------------
  vacc_start_week <- rep(NA, A)
  vacc_end_week   <- rep(NA, A)
  
  # ------------------------------------------------------------------------
  # 6) Set initial conditions at t = 1
  # ------------------------------------------------------------------------
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a,1] - R_[a,1]
    V[a, 1]  <- 0
    age_stratified_cases[a, 1] <- rho * I[a,1]
    coverage_frac[a, 1] <- V[a,1] / N[a]
  }
  coverage_target[1] <- 0
  
  # ------------------------------------------------------------------------
  # 7) Main Loop, t = 2 .. T (Finite Stockpile, with Reallocation for Equal Coverage)
  # ------------------------------------------------------------------------
  weekly_dose_total_history <- numeric(T)
  
  for(t_ in 2:T) {
    # ---- A) Apply Vaccination Effect: doses allocated at (t_-2) become effective now.
    if(t_ - 2 >= 1) {
      for(a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t_ - 2]
        # Now, at time t_ - 1, subtract these doses from S and add them to V.
        actual_move <- min(S[a, t_ - 1], newly_effective)
        S[a, t_ - 1] <- S[a, t_ - 1] - actual_move
        V[a, t_ - 1] <- V[a, t_ - 1] + actual_move
      }
    }
    
    # ---- B) Compute Force of Infection for week t_
    infected_prev <- sum(I[, t_ - 1])
    region_pop <- sum(N)
    phi <- base_beta[t_ - 1] * (infected_prev / region_pop)
    phi_vec[t_] <- phi
    
    # ---- C) Update Weekly Dose Total Based on Remaining Global Stock
    if(t_ < delay) {
      weekly_dose_total <- 0
    } else {
      weekly_dose_total <- global_stock * weekly_delivery_speed
    }
    weekly_dose_total_history[t_] <- weekly_dose_total
    
    # ---- D) Reallocate Weekly Doses Among Target Groups Based on Remaining Need
    if(!no_target) {
      # Compute remaining need (absolute number) for each target group:
      remaining_need <- numeric(A)
      for(a in target_indices) {
        remaining_need[a] <- max(0, desired_doses[a] - total_supply_age[a])
      }
      sum_remaining_need <- sum(remaining_need[target_indices])
      if(sum_remaining_need > 0) {
        # Limit the weekly dose total to the total remaining need if smaller.
        weekly_dose_target <- min(weekly_dose_total, sum_remaining_need)
        for(a in target_indices) {
          weekly_supply_age[a] <- weekly_dose_target * (remaining_need[a] / sum_remaining_need)
        }
      } else {
        for(a in target_indices) {
          weekly_supply_age[a] <- 0
        }
      }
      for(a in setdiff(seq_len(A), target_indices)) {
        weekly_supply_age[a] <- 0
      }
    } else {
      weekly_supply_age <- rep(0, A)
    }
    
    # ---- E) Allocate Vaccines This Week & (IMPORTANTLY) Do NOT Remove from S Immediately
    # Instead, we only record the allocation and update the per-age cumulative totals.
    # The removal from S (and addition to V) will occur in step A after the delay.
    weekly_used <- 0
    for(a in seq_len(A)) {
      if(t_ < delay) {
        doses_planned <- 0
      } else {
        doses_planned <- weekly_supply_age[a]
      }
      susceptibles_now <- S[a, t_ - 1]
      remaining_supply <- unused_supply_age[a]
      
      feasible <- min(doses_planned, remaining_supply, susceptibles_now)
      feasible <- max(feasible, 0)
      
      raw_allocation_age[a, t_] <- feasible
      
      if(feasible > 0) {
        if(is.na(vacc_start_week[a])) {
          vacc_start_week[a] <- t_
        }
        vacc_end_week[a] <- t_
        
        total_supply_age[a] <- total_supply_age[a] + feasible
        unused_supply_age[a] <- total_available_supply_age[a] - total_supply_age[a]
        
        # **Key change:** DO NOT subtract feasible from S here.
        # Instead, record the allocated doses in vacc_delayed for later application.
        vacc_delayed[a, t_] <- feasible
        
        weekly_used <- weekly_used + feasible
      }
    }
    
    # ---- F) Update Global Stockpile by Subtracting the Doses Allocated This Week
    global_stock <- global_stock - weekly_used
    if(global_stock < 0) { global_stock <- 0 }
    global_stock_history[t_] <- global_stock
    
    # ---- G) S-I-R Transitions for Week t_
    for(a in seq_len(A)) {
      if(a == 1) {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1]
      } else {
        S[a, t_] <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1] +
          r[a - 1] * S[a - 1, t_ - 1]
        I[a, t_] <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1] +
          r[a - 1] * I[a - 1, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1] +
          r[a - 1] * R_[a - 1, t_ - 1]
        V[a, t_] <- V[a, t_ - 1] - r[a] * V[a, t_ - 1] +
          r[a - 1] * V[a - 1, t_ - 1]
      }
      S[a, t_] <- pmax(S[a, t_], 0)
      I[a, t_] <- pmax(I[a, t_], 0)
      R_[a, t_] <- pmax(R_[a, t_], 0)
      V[a, t_] <- pmax(V[a, t_], 0)
    }
    
    # ---- H) Update Reporting for Week t_
    for(a in seq_len(A)) {
      coverage_frac[a, t_] <- V[a, t_] / N[a]
      age_stratified_cases[a, t_] <- rho * I[a, t_]
    }
    if(!no_target) {
      target_vacc_cumulative <- sum(raw_allocation_age[target_indices, 1:t_])
      coverage_target[t_] <- target_vacc_cumulative / target_pop
    } else {
      coverage_target[t_] <- 0
    }
    
  }  # End main loop
  
  # ------------------------------------------------------------------------
  # 8) Compute Vaccination Duration by Age Group
  # ------------------------------------------------------------------------
  real_vaccine_duration_age <- numeric(A)
  for(a in seq_len(A)) {
    if(!is.na(vacc_start_week[a]) && !is.na(vacc_end_week[a])) {
      real_vaccine_duration_age[a] <- vacc_end_week[a] - vacc_start_week[a] + 1
    } else {
      real_vaccine_duration_age[a] <- 0
    }
  }
  
  # ------------------------------------------------------------------------
  # 9) Return Results
  # ------------------------------------------------------------------------
  return(list(
    S = S,
    I = I,
    R = R_,
    V = V,
    vacc_delayed = vacc_delayed,
    total_supply = total_supply,           # initial total supply (for reference)
    total_stock_final = global_stock,        # remaining global stock at end
    global_stock_history = global_stock_history,
    weekly_dose_total_history = weekly_dose_total_history,
    total_supply_age = total_supply_age,
    unused_supply_age = unused_supply_age,
    coverage_frac = coverage_frac,
    coverage_target = coverage_target,
    age_stratified_cases = age_stratified_cases,
    phi = phi_vec,
    vacc_start_week = vacc_start_week,
    vacc_end_week = vacc_end_week,
    real_vaccine_duration_age = real_vaccine_duration_age,
    raw_allocation_age = raw_allocation_age
  ))
}

# allow for sequential distribution
sirv_sim_coverageSwitch_seq <- function(
    # Model dimensions
  T, A,           # T = number of weeks, A = number of age groups
  N,              # population by age (vector of length A)
  r,              # aging rates (vector length A)
  
  # Epidemic parameters
  base_beta,      # vector of length T (or function of time)
  I0_draw,        # initial infected by age (vector of length A)
  R0,             # initial underlying immunity fraction by age (vector of length A)
  rho,            # detection probability
  gamma,          # recovery rate
  
  # Vaccination parameters
  delay,                 # week index to START vaccination
  coverage_threshold,    # e.g., 0.7 means 70% coverage needed in a group before switching
  total_coverage,        # fraction of total population (or of a group) eventually vaccinated
  weekly_delivery_speed, # fraction of available supply used per week
  VE_block,              # vaccine efficacy
  
  # Option to toggle phase-2 reallocation (not used in this example)
  use_coverage_switch = FALSE  
) {
  # ------------------------------------------------------------------------
  # 1) Initialize compartments and data structures
  # ------------------------------------------------------------------------
  S  <- matrix(0, nrow = A, ncol = T)
  I  <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V  <- matrix(0, nrow = A, ncol = T)
  
  # For delayed vaccine effect (becomes effective after 2 weeks)
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # For reporting
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)  # weekly detected cases by age
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)  # doses allocated each week by age
  coverage_frac        <- matrix(0, nrow = A, ncol = T)  # fraction vaccinated in each age group
  coverage_target      <- numeric(T)                     # coverage in the currently targeted group
  phi_vec              <- numeric(T)                     # force of infection each week
  
  # ------------------------------------------------------------------------
  # 2) Set initial conditions at t = 1
  # ------------------------------------------------------------------------
  for(a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a, 1] - R_[a, 1]
    V[a, 1]  <- 0
    age_stratified_cases[a, 1] <- rho * I[a, 1]
    coverage_frac[a, 1]       <- V[a, 1] / N[a]
  }
  coverage_target[1] <- 0
  
  # ------------------------------------------------------------------------
  # 3) Set up sequential vaccination parameters
  # ------------------------------------------------------------------------
  # Priority: vaccinate one age group at a time from oldest (group 1) downward.
  current_target <- 1  # start with age group 1 (oldest)
  # For each age group, the total available doses are a fraction of that group's population.
  # (Here we assume we aim to vaccinate 'total_coverage' fraction of that group.)
  total_supply_age <- total_coverage * N  # vector of available doses per age group
  
  # We'll track doses already used in each group.
  used_supply_age <- rep(0, A)
  
  # Compute the fixed weekly dose amount for the group in focus:
  # (For simplicity, we allocate a fixed fraction of the available supply each week.)
  # You may wish to modify this logic.
  weekly_dose_fraction <- weekly_delivery_speed
  
  # ------------------------------------------------------------------------
  # 4) Main simulation loop, t = 2..T
  # ------------------------------------------------------------------------
  for(t in 2:T) {
    # -----------------------------
    # 4A) Apply delayed vaccine effect (2 weeks lag)
    # -----------------------------
    if(t - 2 >= 1) {
      for(a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t - 2]
        actual_move <- min(S[a, t - 1], newly_effective)
        S[a, t - 1] <- S[a, t - 1] - actual_move
        V[a, t - 1] <- V[a, t - 1] + actual_move
      }
    }
    
    # -----------------------------
    # 4B) Check if current target group has reached its threshold
    # -----------------------------
    if(current_target <= A) {
      current_cov <- V[current_target, t - 1] / N[current_target]
      if(current_cov >= coverage_threshold) {
        # Move to the next age group
        current_target <- current_target + 1
      }
    }
    
    # -----------------------------
    # 4C) Set vaccination supply for this week
    # -----------------------------
    # Allocate weekly supply only to the current target group (if any remain)
    weekly_supply_age <- rep(0, A)
    if(current_target <= A && t >= delay) {
      remaining <- total_supply_age[current_target] - used_supply_age[current_target]
      # Compute weekly doses as a fraction of the remaining supply (or a fixed amount)
      weekly_doses <- weekly_dose_fraction * total_supply_age[current_target]
      # Do not allocate more than the remaining susceptibles in that age group:
      weekly_supply_age[current_target] <- min(weekly_doses, remaining, S[current_target, t - 1])
    }
    
    # -----------------------------
    # 4D) Allocate vaccines (record in raw_allocation_age and vacc_delayed)
    # -----------------------------
    for(a in seq_len(A)) {
      doses_planned <- if(t >= delay) weekly_supply_age[a] else 0
      feasible <- min(doses_planned, S[a, t - 1])
      feasible <- max(feasible, 0)
      raw_allocation_age[a, t] <- feasible
      if(feasible > 0) {
        vacc_delayed[a, t] <- feasible  # these will become effective 2 weeks later
        if(a == current_target) {
          used_supply_age[a] <- used_supply_age[a] + feasible
        }
      }
    }
    
    # -----------------------------
    # 4E) Compute Force of Infection and update compartments
    # -----------------------------
    infected_prev <- sum(I[, t - 1])
    region_pop <- sum(N)
    phi <- base_beta[t - 1] * (infected_prev / region_pop)
    phi_vec[t] <- phi
    
    for(a in seq_len(A)) {
      if(a == 1) {
        S[a, t] <- S[a, t - 1] - phi * S[a, t - 1] - r[a] * S[a, t - 1]
        I[a, t] <- I[a, t - 1] + phi * S[a, t - 1] - gamma * I[a, t - 1] - r[a] * I[a, t - 1]
        R_[a, t] <- R_[a, t - 1] + gamma * I[a, t - 1] - r[a] * R_[a, t - 1]
        V[a, t] <- V[a, t - 1] - r[a] * V[a, t - 1]
      } else {
        S[a, t] <- S[a, t - 1] - phi * S[a, t - 1] - r[a] * S[a, t - 1] + r[a - 1] * S[a - 1, t - 1]
        I[a, t] <- I[a, t - 1] + phi * S[a, t - 1] - gamma * I[a, t - 1] - r[a] * I[a, t - 1] + r[a - 1] * I[a - 1, t - 1]
        R_[a, t] <- R_[a, t - 1] + gamma * I[a, t - 1] - r[a] * R_[a, t - 1] + r[a - 1] * R_[a - 1, t - 1]
        V[a, t] <- V[a, t - 1] - r[a] * V[a, t - 1] + r[a - 1] * V[a - 1, t - 1]
      }
      S[a, t] <- pmax(0, S[a, t])
      I[a, t] <- pmax(0, I[a, t])
      R_[a, t] <- pmax(0, R_[a, t])
      V[a, t] <- pmax(0, V[a, t])
    }
    
    # -----------------------------
    # 4F) Update coverage and record cases
    # -----------------------------
    for(a in seq_len(A)) {
      coverage_frac[a, t] <- V[a, t] / N[a]
      age_stratified_cases[a, t] <- rho * I[a, t]
    }
    if(current_target <= A) {
      coverage_target[t] <- V[current_target, t] / N[current_target]
    } else {
      coverage_target[t] <- 1
    }
  }  # End of main loop
  
  # ------------------------------------------------------------------------
  # 5) Return results
  # ------------------------------------------------------------------------
  return(list(
    S = S,
    I = I,
    R = R_,
    V = V,
    vacc_delayed = vacc_delayed,
    raw_allocation_age = raw_allocation_age,
    coverage_frac = coverage_frac,
    coverage_target = coverage_target,
    phi = phi_vec,
    age_stratified_cases = age_stratified_cases,
    used_supply_age = used_supply_age  # for debugging/inspection
  ))
}



