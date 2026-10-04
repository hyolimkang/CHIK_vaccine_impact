## Final simulation functions for Ixchiq and Vimkunya incorporating waning dynamics

# for prevacc
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma, sigma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,      # disease-blocking efficacy
  VE_inf  = 0,   # infection-blocking efficacy
  coverage_threshold  = 1,
  use_coverage_switch = FALSE,
  p_sym = 0.524  # symptomatic probability; rho = p_sym * rho_sym (age_struc_bra_seir_m3.stan)
) {
  
  pmax0      <- function(x) pmax(0, x)
  total_pop  <- sum(N)
  R0_scalar  <- base_beta[1] / gamma
  R0_vec     <- base_beta / gamma
  R_eff_vec  <- numeric(T)
  
  # ───────────────────────────────────────
  # state matrices
  # ───────────────────────────────────────
  S  <- I  <- R_ <- V <- E <- matrix(0, A, T)
  V_covered <- matrix(0, A, T)
  vacc_delayed        <- raw_allocation_age <- wasted_dose <- matrix(0, A, T)
  unvaccinated        <- N
  
  # IMPORTANT:
  # - age_stratified_cases_raw: true infection incidence per time step (new infections)
  # - age_stratified_cases:     reported infection incidence per time step
  coverage_frac       <- age_stratified_cases <- age_stratified_cases_raw <- matrix(0, A, T)
  
  # IMPORTANT:
  # - true_symptomatic: reported symptomatic incidence per time step
  true_symptomatic    <- matrix(0, A, T)
  
  phi_vec             <- FOI_times_S <- coverage_target <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  # ───────────────────────────────────────
  # vaccine supply
  # ───────────────────────────────────────
  total_supply       <- target_pop * total_coverage
  weekly_dose_total  <- total_supply * weekly_delivery_speed
  
  total_avail_age <- unused_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  # ───────────────────────────────────────
  # initial state (t = 1)
  # ───────────────────────────────────────
  S[, 1] <- N - I0_draw - R0 * N
  E[, 1] <- 0
  I[, 1] <- I0_draw
  R_[, 1] <- R0 * N
  V[, 1] <- 0
  
  coverage_frac[, 1] <- 0
  
  # incidence at t=1 set to 0 (no "new infections" counted at initialisation)
  age_stratified_cases_raw[, 1] <- 0
  age_stratified_cases[,  1]    <- 0
  true_symptomatic[, 1]         <- 0
  
  coverage_target[1]   <- 0
  eff_S0               <- sum(S[, 1])
  R_eff_vec[1]         <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  # ───────────────────────────────────────
  # main loop (t = 2 … T)
  # ───────────────────────────────────────
  for (t_ in 2:T) {
    
    if (t_ - 2 >= 1) {
      effective_dose      <- vacc_delayed[, t_ - 2]             # vaccinated 2 weeks before
      immunized           <- round(VE_inf * effective_dose)     # successfully immunised
      V[, t_]             <- V[, t_ - 1] + immunized
      V_covered[, t_]     <- V_covered[, t_ - 1] + effective_dose
    } else {
      immunized       <- rep(0, A)
      V[, t_] <- V[, t_ - 1]
      V_covered[, t_] <- V_covered[, t_ - 1]
    }
    
    coverage_frac[, t_] <- V_covered[, t_] / N
    
    # FOI (hazard)
    eff_S_prev_vec <- pmax0(S[, t_-1] - immunized)
    phi            <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / total_pop)
    phi_vec[t_]    <- phi
    FOI_times_S[t_] <- phi * sum(eff_S_prev_vec)
    
    # after delay: allocate doses
    if (t_ >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(
          ceiling(weekly_dose_total * (N[a] / target_pop)),
          rem,
          unvaccinated[a],
          total_avail_age[a] - total_used_age[a]
        )
        if (alloc > 0) {
          prop_S      <- S[a, t_ - 1] / N[a]
          vacc_to_S   <- round(alloc * prop_S)
          wasted      <- alloc - vacc_to_S
          
          raw_allocation_age[a, t_] <- alloc
          vacc_delayed[a,       t_] <- vacc_to_S
          wasted_dose[a,        t_] <- wasted
          
          total_used_age[a] <- total_used_age[a] + alloc
          unvaccinated[a]   <- unvaccinated[a] - alloc
          rem               <- rem - alloc
          
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    for (a in seq_len(A)) {
      
      S_prev <- S[a, t_ - 1];  E_prev <- E[a, t_ - 1]
      I_prev <- I[a, t_ - 1];  R_prev <- R_[a, t_ - 1]
      
      immunized_a <- immunized[a]
      S_prev_eff  <- pmax0(S_prev - immunized_a)
      
      # new infections (incidence of infection event)
      # Match Stan's p_infection = 1-exp(-lambda), not a linear phi
      # approximation -- diverges from Stan for non-small phi (see chat
      # record).
      p_inf <- 1 - exp(-phi)
      new_e <- p_inf * S_prev_eff
      
      # latent -> infectious; recovery
      new_i <- sigma * E_prev
      recov <- gamma * I_prev
      
      # update compartments
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # REDEFINED OUTPUTS (names unchanged):
      # - infections are incidence (new_e), not prevalence (I)
      age_stratified_cases_raw[a, t_] <- new_e
      age_stratified_cases[a,     t_] <- rho * new_e
      
      # true (unreported-rate-adjusted) symptomatic incidence from infection
      # incidence -- downstream (summarise_presim_ui) multiplies this
      # directly by a clinical hospitalisation rate, so it must be the TRUE
      # symptomatic burden, not rho-discounted to reported-case scale.
      # Correction (see chat record): dropped the erroneous extra `rho`
      # factor, which would have applied a reporting-rate discount here.
      # REPORTED symptomatic incidence: hosp_rate downstream (summarise_presim_ui)
      # is calibrated against reported cases, not true symptomatic burden, so
      # this must stay on the reported-case scale. rho = p_sym * rho_sym
      # already includes p_sym (age_struc_bra_seir_m3.stan), so rho * new_e
      # alone is "reported cases" -- matches Stan's own
      # expected_reported_cases = rho * new_infectious. (Reverted back from
      # the p_sym-only version -- see chat record, 2026-08-20.)
      # TRUE symptomatic incidence (p_sym-scaled, NOT rho-discounted):
      # benefit/risk quantities (symptomatic/hosp/death/DALY averted) must
      # reflect true population impact, not the surveillance-detected
      # fraction -- rho would make this "reported", not "true". Final
      true_symptomatic[a, t_] <-
        new_e *
        p_sym *
        (1 - VE_block * coverage_frac[a, t_])

      # aging
      if (a == 1) {
        S[a, t_] <- pmax0(S[a, t_] - r[a] * S[a, t_])
        E[a, t_] <- pmax0(E[a, t_] - r[a] * E[a, t_])
        I[a, t_] <- pmax0(I[a, t_] - r[a] * I[a, t_])
        R_[a,t_] <- pmax0(R_[a,t_] - r[a] * R_[a,t_])
        V[a, t_] <- pmax0(V[a, t_] - r[a] * V[a, t_])
      } else {
        S[a, t_] <- pmax0(S[a, t_] - r[a] * S[a, t_] + r[a-1] * S[a-1, t_-1])
        E[a, t_] <- pmax0(E[a, t_] - r[a] * E[a, t_] + r[a-1] * E[a-1, t_-1])
        I[a, t_] <- pmax0(I[a, t_] - r[a] * I[a, t_] + r[a-1] * I[a-1, t_-1])
        R_[a,t_] <- pmax0(R_[a,t_] - r[a] * R_[a,t_] + r[a-1] * R_[a-1, t_-1])
        V[a, t_] <- pmax0(V[a, t_] - r[a] * V[a, t_] + r[a-1] * V[a-1, t_-1])
      }
    }
    
    eff_S         <- sum(S[, t_])
    R_eff_vec[t_] <- (base_beta[t_] / gamma) * (eff_S / total_pop)
  }
  
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end),
                         0,
                         vacc_end - vacc_start + 1)
  
  list(
    S = S, E = E, I = I, R = R_, V = V,
    vacc_delayed = vacc_delayed,
    raw_allocation_age = raw_allocation_age,
    wasted_dose  = wasted_dose,
    total_supply = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,
    unused_supply_age          = unused_age,
    coverage_target = coverage_target,
    coverage_frac   = coverage_frac,
    phi             = phi_vec,
    FOI_times_S     = FOI_times_S,
    R0_basic        = R0_scalar,
    R0_vec          = R0_vec,
    R_eff_vec       = R_eff_vec,
    age_stratified_cases      = age_stratified_cases,
    age_stratified_cases_raw  = age_stratified_cases_raw,
    true_symptomatic          = true_symptomatic,
    vacc_start_week           = vacc_start,
    vacc_end_week             = vacc_end,
    real_vaccine_duration_age = real_vac_dur
  )
}


# updated final seir 
sirv_sim_coverageSwitch_ixchiq <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma, sigma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,      # disease-blocking efficacy
  VE_inf  = 0,   # infection-blocking efficacy
  coverage_threshold  = 1,
  use_coverage_switch = FALSE,
  p_sym = 0.524  # symptomatic probability; rho = p_sym * rho_sym (age_struc_bra_seir_m3.stan)
) {
  
  pmax0      <- function(x) pmax(0, x)
  total_pop  <- sum(N)
  R0_scalar  <- base_beta[1] / gamma
  R0_vec     <- base_beta / gamma
  R_eff_vec  <- numeric(T)
  
  # ───────────────────────────────────────
  # initial matrix
  # ───────────────────────────────────────
  S  <- I  <- R_ <- V <- E <- matrix(0, A, T)
  V_covered <- matrix(0, A, T)
  vacc_delayed        <- raw_allocation_age <- wasted_dose <- matrix(0, A, T)
  unvaccinated        <- N
  
  coverage_frac       <- age_stratified_cases <- age_stratified_cases_raw <- matrix(0, A, T)
  true_symptomatic    <- matrix(0, A, T)
  
  phi_vec             <- FOI_times_S <- coverage_target <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  # ───────────────────────────────────────
  # vaccine supply 
  # ───────────────────────────────────────
  total_supply       <- target_pop * total_coverage
  weekly_dose_total  <- total_supply * weekly_delivery_speed
  
  total_avail_age <- unused_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  # ───────────────────────────────────────
  # initial state (t = 1)
  # ───────────────────────────────────────
  S[, 1] <- N - I0_draw - R0 * N
  E[, 1] <- 0
  I[, 1] <- I0_draw
  R_[, 1] <- R0 * N
  V[, 1] <- 0
  
  coverage_frac[, 1] <- 0
  # incidence at t=1 set to 0 (no "new infections" counted at initialisation)
  age_stratified_cases_raw[, 1] <- 0
  age_stratified_cases[,  1]    <- 0
  true_symptomatic[, 1] <- 0
  coverage_target[1]   <- 0
  eff_S0               <- sum(S[, 1])    
  R_eff_vec[1]         <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  # ───────────────────────────────────────
  # main loop (t = 2 … T)
  # ───────────────────────────────────────
  for (t_ in 2:T) {
    
    if (t_ - 2 >= 1) {
      effective_dose      <- vacc_delayed[, t_ - 2]              # vaccinated 2 weeks before
      immunized           <- round(VE_inf * effective_dose)      # dose * VE = successfully vaccinated
      V[, t_]             <- V[, t_ - 1] + immunized             # V cumulative
      V_covered[, t_]     <- V_covered[, t_ - 1] + effective_dose # VE block track
      #S[, t_ - 1]         <- pmax0(S[, t_ - 1] - immunized)      
    } else {
      immunized       <- rep(0, A)
      V[, t_] <- V[, t_ - 1]
      V_covered[, t_] <- V_covered[, t_ - 1]
    }
    
    coverage_frac[, t_] <- V_covered[, t_] / N                           
    
    # FOI 
    eff_S_prev_vec <- pmax0(S[, t_-1] - immunized) 
    phi             <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / total_pop)
    phi_vec[t_]     <- phi
    FOI_times_S[t_] <- phi * sum(eff_S_prev_vec)
    
    # after delay
    if (t_ >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(
          ceiling(weekly_dose_total * (N[a] / target_pop)),
          rem, # remaining doses at this week 
          unvaccinated[a],
          total_avail_age[a] - total_used_age[a]
        )
        if (alloc > 0) {
          prop_S      <- S[a, t_ - 1] / N[a]
          vacc_to_S   <- round(alloc * prop_S)        
          wasted      <- alloc - vacc_to_S            
          raw_allocation_age[a, t_] <- alloc
          vacc_delayed[a,       t_] <- vacc_to_S
          wasted_dose[a,        t_] <- wasted
          
          total_used_age[a] <- total_used_age[a] + alloc
          unvaccinated[a]   <- unvaccinated[a] - alloc
          rem               <- rem - alloc
          
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    for (a in seq_len(A)) {
      S_prev <- S[a, t_ - 1];  E_prev <- E[a, t_ - 1]
      I_prev <- I[a, t_ - 1];  R_prev <- R_[a, t_ - 1]
      
      immunized_a <- immunized[a]                        
      S_prev_eff <- pmax0(S_prev - immunized_a) 
      
      # Match Stan's p_infection = 1-exp(-lambda), not a linear phi
      # approximation (see chat record).
      p_inf  <- 1 - exp(-phi)
      new_e  <- p_inf * S_prev_eff
      new_i  <- sigma * E_prev
      recov  <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # true (unreported-rate-adjusted) symptomatic incidence -- downstream
      # multiplies this directly by a clinical hospitalisation rate, so it
      # must be the TRUE symptomatic burden, not rho-discounted. Correction
      # (see chat record): dropped the erroneous `rho` factor.
      # REPORTED symptomatic incidence -- rho = p_sym * rho_sym already
      # includes p_sym (age_struc_bra_seir_m3.stan), so rho * new_e alone is
      # "reported cases" (matches Stan's expected_reported_cases = rho *
      # new_infectious). hosp_rate downstream is reported-case-calibrated.
      # (Reverted back from the p_sym-only version -- see chat record,
      # 2026-08-20.)
      # TRUE symptomatic incidence (p_sym-scaled, NOT rho-discounted) --
      # benefit/risk quantities must reflect true population impact, not the
      # surveillance-detected fraction. Final (chat record, 2026-08-20).
      true_symptomatic[a, t_] <-
        new_e *
        p_sym *
        (1 - VE_block * coverage_frac[a, t_])

      # ── 4-1) Aging ──
      if (a == 1) {
        S[a, t_] <- pmax0(S[a, t_] - r[a] * S[a, t_])
        E[a, t_] <- pmax0(E[a, t_] - r[a] * E[a, t_])
        I[a, t_] <- pmax0(I[a, t_] - r[a] * I[a, t_])
        R_[a,t_] <- pmax0(R_[a,t_] - r[a] * R_[a,t_])
        V[a, t_] <- pmax0(V[a, t_] - r[a] * V[a, t_])
      } else {
        S[a, t_] <- pmax0(S[a, t_] - r[a] * S[a, t_] + r[a-1] * S[a-1, t_-1])
        E[a, t_] <- pmax0(E[a, t_] - r[a] * E[a, t_] + r[a-1] * E[a-1, t_-1])
        I[a, t_] <- pmax0(I[a, t_] - r[a] * I[a, t_] + r[a-1] * I[a-1, t_-1])
        R_[a,t_] <- pmax0(R_[a,t_] - r[a] * R_[a,t_] + r[a-1] * R_[a-1, t_-1])
        V[a, t_] <- pmax0(V[a, t_] - r[a] * V[a, t_] + r[a-1] * V[a-1, t_-1])
      }
      
      age_stratified_cases_raw[a, t_] <- new_e
      age_stratified_cases[a,     t_] <- rho * new_e
    }
    
    # effective S
    eff_S              <- sum(S[, t_])          
    R_eff_vec[t_]      <- (base_beta[t_] / gamma) * (eff_S / total_pop)
  }
  
  # ───────────────────────────────────────
  # output
  # ───────────────────────────────────────
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end),
                         0,
                         vacc_end - vacc_start + 1)
  
  list(
    S = S, E = E, I = I, R = R_, V = V,
    vacc_delayed = vacc_delayed,
    raw_allocation_age = raw_allocation_age,
    wasted_dose  = wasted_dose,
    total_supply = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,
    unused_supply_age          = unused_age,
    coverage_target = coverage_target,
    coverage_frac   = coverage_frac,
    phi             = phi_vec,
    FOI_times_S     = FOI_times_S,
    R0_basic        = R0_scalar,
    R0_vec          = R0_vec,
    R_eff_vec       = R_eff_vec,
    age_stratified_cases      = age_stratified_cases,
    age_stratified_cases_raw  = age_stratified_cases_raw,
    true_symptomatic          = true_symptomatic,
    vacc_start_week           = vacc_start,
    vacc_end_week             = vacc_end,
    real_vaccine_duration_age = real_vac_dur
  )
}


## time-varying waning dyanmics (linear decrease)
sirv_sim_coverageSwitch_vimkun <- function(
    ## 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma, sigma,
  ## 2) Vaccine parameters
  delay                = 2,   
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,                    # disease‑blocking efficacy
  VE_inf          = 0.978,     # infection‑blocking efficacy 
  VE_inf_postwane = 0.855,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE,
  p_sym = 0.524  # symptomatic probability; rho = p_sym * rho_sym (age_struc_bra_seir_m3.stan)
) {
  
  if (length(VE_block) == 1) VE_block <- rep(VE_block, T)
  if (length(VE_inf)   == 1) VE_inf   <- rep(VE_inf,   T)
  
  immun_delay  <- 3    
  waning_weeks <- 26   
  
  pmax0      <- function(x) pmax(0, x)
  total_pop  <- sum(N)
  R0_scalar  <- base_beta[1] / gamma
  R0_vec     <- base_beta / gamma
  
  S <- I <- R_ <- V <- E <- matrix(0, A, T)
  V_covered <- matrix(0, A, T)
  vacc_delayed <- raw_allocation_age <- wasted_dose <- matrix(0, A, T)
  vacc_success <- matrix(0, A, T)  
  unvaccinated <- N
  
  coverage_frac  <- age_stratified_cases <- age_stratified_cases_raw <- matrix(0, A, T)
  true_symptomatic <- matrix(0, A, T)
  phi_vec <- FOI_times_S <- coverage_target <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  #total_supply      <- target_pop * total_coverage
  total_supply <- target_pop * total_coverage  # cepi doses (absolute dose model)
  weekly_dose_total <- total_supply * weekly_delivery_speed
  total_avail_age   <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  S[, 1] <- N - I0_draw - R0 * N
  I[, 1] <- I0_draw
  R_[, 1] <- R0 * N
  coverage_frac[, 1] <- 0
  age_stratified_cases_raw[, 1] <- I[, 1]
  age_stratified_cases[, 1]     <- rho * I[, 1]
  true_symptomatic[, 1]         <- rho * I[, 1]
  coverage_target[1] <- 0
  eff_S0 <- sum(S[, 1])
  R_eff_vec <- numeric(T)
  R_eff_vec[1] <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  for (t_ in 2:T) {
    immunized_success <- rep(0, A)  
    if (t_ - immun_delay >= 1) {
      effective_dose    <- vacc_delayed[, t_ - immun_delay]  
      immunized_success <- round(VE_inf[t_] * effective_dose)
      
      vacc_success[, t_] <- immunized_success           
      V[, t_]         <- V[, t_-1] + immunized_success
      V_covered[, t_] <- V_covered[, t_-1] + effective_dose
    } else {
      V[, t_]         <- V[, t_-1]
      V_covered[, t_] <- V_covered[, t_-1]
    }
    
    reinstate <- rep(0, A)
    if (t_ - (waning_weeks - immun_delay) >= 1) {   
      idx_wane   <- t_ - (waning_weeks - immun_delay)  
      waning_group <- vacc_success[, idx_wane]
      
      VE_inf_prev <- VE_inf[idx_wane]
      if (VE_inf_prev > 0) {
        #VE_inf_postwane <- 0.855
        current_reinstate <- round((1 - VE_inf_postwane / VE_inf_prev) * waning_group)
        current_reinstate[is.na(current_reinstate) | current_reinstate < 0] <- 0
        reinstate <- current_reinstate
      }
      
      if (sum(reinstate) > 0) {
        V[, t_]         <- pmax0(V[, t_] - reinstate)
        V_covered[, t_] <- pmax0(V_covered[, t_] - reinstate)
      }
    }
    
    phi             <- base_beta[t_-1] * (sum(I[, t_-1]) / total_pop)
    phi_vec[t_]     <- phi
    FOI_times_S[t_] <- phi * sum(S[, t_-1])
    
    if (t_ >= delay && target_pop > 0) {
      rem <- weekly_dose_total
      for (a in target_idx) {
        alloc <- min(
          ceiling(weekly_dose_total * (N[a] / target_pop)),
          rem,
          unvaccinated[a],
          total_avail_age[a] - total_used_age[a]
        )
        if (alloc > 0) {
          prop_S    <- S[a, t_-1] / N[a]
          vacc_to_S <- round(alloc * prop_S)
          wasted    <- alloc - vacc_to_S
          
          raw_allocation_age[a, t_] <- alloc
          vacc_delayed[a,       t_] <- vacc_to_S
          wasted_dose[a,        t_] <- wasted
          
          total_used_age[a] <- total_used_age[a] + alloc
          unvaccinated[a]   <- unvaccinated[a] - alloc
          rem               <- rem - alloc
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_-1]
    }
    
    coverage_frac[, t_] <- V_covered[, t_] / N
    
    ## 4) SEIRV 
    for (a in seq_len(A)) {
      S_prev <- S[a, t_-1]; E_prev <- E[a, t_-1]
      I_prev <- I[a, t_-1]; R_prev <- R_[a, t_-1]
      
      immunized_a <- immunized_success[a]
      S_prev_eff  <- pmax0(S_prev - immunized_a + reinstate[a])

      # Match Stan's p_infection = 1-exp(-lambda), not a linear phi
      # approximation (see chat record).
      p_inf <- 1 - exp(-phi)
      new_e <- p_inf * S_prev_eff
      new_i <- sigma * E_prev
      recov <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # true (unreported-rate-adjusted) symptomatic prevalence -- downstream
      # multiplies this directly by a clinical hospitalisation rate, so it
      # must be the TRUE symptomatic burden, not rho-discounted. Correction
      # (see chat record): dropped the erroneous `rho` factor.
      # REPORTED symptomatic prevalence -- rho = p_sym * rho_sym already
      # includes p_sym, so rho * I alone is "reported cases". (Reverted back
      # from the p_sym-only version -- see chat record, 2026-08-20.)
      # TRUE symptomatic prevalence (p_sym-scaled, NOT rho-discounted) --
      true_symptomatic[a, t_] <-
        I[a, t_] *
        p_sym *
        (1 - VE_block[t_] * coverage_frac[a, t_])
    }
    
    age_stratified_cases_raw[, t_] <- I[, t_]
    
    coverage_frac[, t_] <- V_covered[, t_] / N
    eff_S          <- sum(S[, t_])
    R_eff_vec[t_]  <- (base_beta[t_] / gamma) * (eff_S / total_pop)
  }
  
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end),
                         0, vacc_end - vacc_start + 1)
  
  list(
    S = S, E = E, I = I, R = R_, V = V,
    vacc_delayed = vacc_delayed,
    raw_allocation_age = raw_allocation_age,
    wasted_dose  = wasted_dose,
    total_supply = total_supply,
    coverage_target = coverage_target,
    coverage_frac   = coverage_frac,
    phi             = phi_vec,
    FOI_times_S     = FOI_times_S,
    R0_basic        = R0_scalar,
    R0_vec          = R0_vec,
    R_eff_vec       = R_eff_vec,
    age_stratified_cases      = age_stratified_cases,
    age_stratified_cases_raw  = age_stratified_cases_raw,
    true_symptomatic          = true_symptomatic,
    vacc_start_week           = vacc_start,
    vacc_end_week             = vacc_end,
    real_vaccine_duration_age = real_vac_dur
  )
}
