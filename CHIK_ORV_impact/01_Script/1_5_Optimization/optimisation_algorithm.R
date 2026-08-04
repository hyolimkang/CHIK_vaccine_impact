# finding grid

grid <- expand.grid(
  age1 = seq(0, 1, by = 0.05),
  age2 = seq(0, 1, by = 0.05),
  age3 = seq(0, 1, by = 0.05),
  age4 = seq(0, 1, by = 0.05)
)
grid <- grid[rowSums(grid) == 1, ]  # Ensure proportions sum to 1

sirv_sim_coverageSwitch_Fatal <- function(
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
  
  allocation_proportions,  # NEW: Proportions for vaccine allocation
  VE_block,                # vaccine efficacy (infection blocking leaky)
  VE_p,                    # vaccine efficacy (disease blocking)
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
  target_indices <- which(allocation_proportions > 0)  
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
      # If no_target and allocation_proportions = 0 => no vaccine
      if(t_ < delay || t_ > vaccine_end) {
        # Outside vaccination window
        raw_allocation[a, t_] <- 0
        vacc_eff <- 0
        
      } else {
        # We are within [delay, vaccine_end]
        if(!coverage_switch) {
          # PHASE 1: target only
          if(!no_target && allocation_proportions[a] > 0) {
            # Allocate proportionally
            raw_allocation_a <- allocation_proportions[a] * sum(N) * 0.1  # Proportional allocation (adjust factor as needed)
            raw_allocation[a, t_] <- raw_allocation_a
            
            # Convert to per-capita rate
            vacc_eff <- raw_allocation_a / N[a]
          } else {
            raw_allocation[a, t_] <- 0
            vacc_eff <- 0
          }
        } else {
          # PHASE 2: all ages
          raw_allocation_a <- allocation_proportions[a] * sum(N) * 0.1  # Proportional allocation (adjust factor as needed)
          raw_allocation[a, t_] <- raw_allocation_a
          
          vacc_eff <- raw_allocation_a / N[a]
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

fr_avg <- c(mean(fatal[1:4]),
            mean(fatal[5:8]),
            mean(fatal[9:12]),
            mean(fatal[13:18])
            )

results_opt <- lapply(1:nrow(grid), function(i) {
  proportions <- as.numeric(grid[i, ])
  sirv_sim_coverageSwitch_Fatal(
    T = 52, A = 4, N = N_agegr,
    r = c(0, 0, 0, 0), base_beta = base_beta,
    I0_draw = I0, rho = rho, gamma = gamma,
    delay = 1, vaccine_end = 25, coverage_threshold = 0.7,
    allocation_proportions = proportions,
    VE_block = 0.75, VE_p = 0.8, FR_infection = fr_avg
  )
})

metrics <- sapply(results_opt, function(res) {
  # Sum infections and mortality by week
  weekly_infections <- colSums(res$age_stratified_cases)  # Sum across rows (age groups)
  total_infections <- sum(weekly_infections)  # Total across all weeks
  
  weekly_mortality <- colSums(res$IFatal_detect + res$IVFatal_detect)  # Sum across age groups
  total_mortality <- sum(weekly_mortality)  # Total across all weeks
  
  c(infections = total_infections, mortality = total_mortality)
})


metrics <- as.data.frame(t(metrics))

# Add metrics to the grid
grid$infections <- metrics$infections
grid$mortality <- metrics$mortality
