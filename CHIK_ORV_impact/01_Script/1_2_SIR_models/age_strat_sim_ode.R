sir_age_model_daily <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    n <- length(r)  # Number of age groups
    
    # Extract compartments
    S  <- state[1:n]                   # Susceptible
    I  <- state[(n + 1):(2 * n)]       # Infectious
    R  <- state[(2 * n + 1):(3 * n)]   # Recovered
    V  <- state[(3 * n + 1):(4 * n)]   # Vaccinated
    IV <- state[(4 * n + 1):(5 * n)]   # Vaccinated and infected
    RV <- state[(5 * n + 1):(6 * n)]   # Vaccinated and recovered from infection
    
    # Per-age total
    N_vec <- S + I + R + V + IV + RV
    # Region-wide population
    total_pop <- sum(N_vec)
    
    N_target <- sum(N_vec[target_age == 1])
    coverage_target_value <- 0
    if (N_target > 0) {
      coverage_target_value <- sum(V[target_age == 1] + RV[target_age == 1]) / N_target
    }
    
    coverage_switch_bool <- (coverage_target_value >= coverage_threshold)
    
    # Time-varying params
    t_index <- floor(t)  # or ceiling, or a spline interpolation, etc.
    if (t_index < 1) t_index <- 1
    if (t_index > length(base_beta)) t_index <- length(base_beta)
    beta_t <- base_beta[t_index]
    # FOI
    infected_all <- sum(I + IV)
    phi <- beta_t * (infected_all / total_pop)
    
    daily_supply_phase1 <- weekly_supply_phase1 
    daily_supply_phase2 <- weekly_supply_phase2 
    
    vacc_eff <- rep(0, n)
    
    if (t >= delay && t <= vaccine_end) {
      if (!coverage_switch_bool) {
        # PHASE 1: target only
        #  1) total target pop = sum(N_vec[target_age==1])
        #  2) fraction of supply allocated to age a = N[a]/(target_pop)
        #  3) raw_allocation[a] = daily_supply_phase1 * pop_frac_a
        #  4) vacc_eff[a] = raw_allocation[a] / N[a]
        if (N_target > 0) {
          for (a in seq_len(n)) {
            if (target_age[a] == 1) {
              pop_frac_a  <- N_vec[a] / N_target
              allocation_a <- daily_supply_phase1 * pop_frac_a
              if (N_vec[a] > 0) {
                vacc_eff[a] <- allocation_a / N_vec[a]
              } else {
                vacc_eff[a] <- 0
              }
            }
          }
        }
      } else {
        # PHASE 2: entire population
        for (a in seq_len(n)) {
          pop_frac_a  <- N_vec[a] / total_pop
          allocation_a <- daily_supply_phase2 * pop_frac_a
          if (N_vec[a] > 0) {
            vacc_eff[a] <- allocation_a / N_vec[a]
          } else {
            vacc_eff[a] <- 0
          }
        }
      }
    } # otherwise vacc_eff = 0 if outside [delay, vaccine_end]
    
    # Initialize derivatives
    dS  <- numeric(n)
    dI  <- numeric(n)
    dR  <- numeric(n)
    dV  <- numeric(n)
    dIV <- numeric(n)
    dRV <- numeric(n)
    
    # Equations for the first age group
    
    # For each age group a:
    for (a in seq_len(n)) {
      
      # Age group "a" (1-based)
      # If you have age-specific aging out rates r[a], 
      #   you might also have inflow from group a-1, outflow to a+1, etc.
      #
      # We'll assume r[a] is your "aging out" from group a. 
      # Then typically there's an inflow r[a-1]*S[a-1], etc. 
      # (like in your discrete code).
      #
      # Let's define r_in[a] = r[a-1] * S[a-1] except a=1 has none,
      #    and r_out[a] = r[a] * S[a].
      #
      # For simplicity, let's just replicate your discrete logic:
      
      #---------------------
      #   S derivative
      #---------------------
      # Infectious force & vaccination
      loss_S <- phi * S[a] + vacc_eff[a] * S[a] + r[a]*S[a]
      inflow_S <- 0
      if (a > 1) {
        inflow_S <- r[a-1] * S[a-1]
      }
      dS[a] <- -loss_S + inflow_S
      
      #---------------------
      #   I derivative
      #---------------------
      gain_I <- phi * S[a]
      loss_I <- gamma*I[a] + r[a]*I[a]
      inflow_I <- 0
      if (a > 1) {
        inflow_I <- r[a-1]*I[a-1]
      }
      dI[a] <- +gain_I - loss_I + inflow_I
      
      #---------------------
      #   R derivative
      #---------------------
      gain_R <- gamma*I[a]
      loss_R <- r[a]*R[a]
      inflow_R <- 0
      if (a > 1) {
        inflow_R <- r[a-1]*R[a-1]
      }
      dR[a] <- +gain_R - loss_R + inflow_R
      
      #---------------------
      #   V derivative
      #---------------------
      gain_V <- vacc_eff[a]*S[a]
      loss_V <- (1 - VE_block)*phi*V[a] + r[a]*V[a]
      inflow_V <- 0
      if (a > 1) {
        inflow_V <- r[a-1]*V[a-1]
      }
      dV[a] <- +gain_V - loss_V + inflow_V
      
      #---------------------
      #   IV derivative
      #---------------------
      gain_IV <- (1 - VE_block)*phi*V[a]
      loss_IV <- gamma*IV[a] + r[a]*IV[a]
      inflow_IV <- 0
      if (a > 1) {
        inflow_IV <- r[a-1]*IV[a-1]
      }
      dIV[a] <- +gain_IV - loss_IV + inflow_IV
      
      #---------------------
      #   RV derivative
      #---------------------
      gain_RV <- gamma*IV[a]
      loss_RV <- r[a]*RV[a]
      inflow_RV <- 0
      if (a > 1) {
        inflow_RV <- r[a-1]*RV[a-1]
      }
      dRV[a] <- +gain_RV - loss_RV + inflow_RV
    }
    
    age_strat_cases <- rho * (I + IV)
    
    # Combine results into a single vector
    list(
      c(dS, dI, dR, dV, dIV, dRV),  # Flatten all derivatives into a single vector
      age_strat_cases = age_strat_cases
    )
    
  })
}

sirv_ode_coverageSwitch <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #---------------------------------------------------------------
    # We assume time t is an integer week (t = 1, 2, …, T).
    # Number of age groups:
    n <- A
    
    #---------------------------------------------------------------
    # (1) Extract state compartments.
    # The state vector order is: 
    #    S[1:n], I[1:n], R[1:n], V[1:n]
    S_current <- state[1:n]
    I_current <- state[(n + 1):(2 * n)]
    R_current <- state[(2 * n + 1):(3 * n)]
    V_current <- state[(3 * n + 1):(4 * n)]
    
    #---------------------------------------------------------------
    # (2) Region-wide and target population.
    region_pop <- sum(N)
    target_indices <- which(target_age == 1)
    if (length(target_indices) > 0) {
      target_pop <- sum(N[target_indices])
    } else {
      target_pop <- 0
    }
    
    # Current target coverage (based on vaccinated V in target groups).
    current_target_coverage <- if (target_pop > 0) sum(V_current[target_indices]) / target_pop else 0
    coverage_switch_bool <- (current_target_coverage >= coverage_threshold)
    
    #---------------------------------------------------------------
    # (3) Force of Infection and Beta Inference.
    # Use FOI_vec[t] (t in weeks) as the target FOI.
    FOI_target <- FOI_vec[t]
    I_total <- sum(I_current)
    if (I_total > 0) {
      Beta_t <- FOI_target * region_pop / I_total
    } else {
      Beta_t <- 0
    }
    # The shared hazard is then:
    phi_shared <- Beta_t * (I_total / region_pop)  # equals FOI_target if I_total > 0
    
    #---------------------------------------------------------------
    # (4) Vaccination Allocation.
    # Compute total supply (for target) and weekly dose.
    total_supply <- target_pop * total_coverage
    weekly_dose_total <- total_supply * weekly_delivery_speed
    
    # For simplicity, we assume that if t < delay, no vaccines are allocated.
    vacc_allocation <- rep(0, n)
    if (t >= delay) {
      if (!coverage_switch_bool) {
        # Phase 1: only target groups receive vaccine.
        if (target_pop > 0) {
          for (a in 1:n) {
            if (target_age[a] == 1) {
              vacc_allocation[a] <- weekly_dose_total * (N[a] / target_pop)
            }
          }
        }
      } else {
        # Phase 2: entire population receives vaccine.
        for (a in 1:n) {
          vacc_allocation[a] <- weekly_dose_total * (N[a] / region_pop)
        }
      }
    }
    # In the discrete simulation, vaccine doses become effective after a 2-week delay.
    # Here we mimic that effect only if t >= 3; otherwise, effective vaccination is zero.
    effective_vaccination <- rep(0, n)
    if (t >= 3) {
      # (Ideally one would use the allocation from t-2; here we approximate by using the current allocation.)
      effective_vaccination <- VE_block * vacc_allocation
    }
    
    #---------------------------------------------------------------
    # (5) Update the compartments following the discrete S–I–R–V logic with aging.
    S_next <- numeric(n)
    I_next <- numeric(n)
    R_next <- numeric(n)
    V_next <- numeric(n)
    
    for (a in 1:n) {
      if (a == 1) {
        S_next[a] <- S_current[a] - phi_shared * S_current[a] - r[a] * S_current[a]
        I_next[a] <- I_current[a] + phi_shared * S_current[a] - gamma * I_current[a] - r[a] * I_current[a]
        R_next[a] <- R_current[a] + gamma * I_current[a] - r[a] * R_current[a]
        V_next[a] <- V_current[a] + effective_vaccination[a] - r[a] * V_current[a]
      } else {
        S_next[a] <- S_current[a] - phi_shared * S_current[a] - r[a] * S_current[a] + r[a - 1] * S_current[a - 1]
        I_next[a] <- I_current[a] + phi_shared * S_current[a] - gamma * I_current[a] - r[a] * I_current[a] + r[a - 1] * I_current[a - 1]
        R_next[a] <- R_current[a] + gamma * I_current[a] - r[a] * R_current[a] + r[a - 1] * R_current[a - 1]
        V_next[a] <- V_current[a] + effective_vaccination[a] - r[a] * V_current[a] + r[a - 1] * V_current[a - 1]
      }
      # Prevent negatives:
      S_next[a] <- max(0, S_next[a])
      I_next[a] <- max(0, I_next[a])
      R_next[a] <- max(0, R_next[a])
      V_next[a] <- max(0, V_next[a])
    }
    
    #---------------------------------------------------------------
    # (6) Compute the “difference” (update) as the derivative approximation over one week.
    dS <- S_next - S_current
    dI <- I_next - I_current
    dR <- R_next - R_current
    dV <- V_next - V_current
    
    #---------------------------------------------------------------
    # (7) Compute additional reporting / debugging outputs.
    age_strat_cases <- rho * I_next
    coverage_frac   <- V_next / N  # per-age vaccination coverage
    cumulative_target_vacc <- if (target_pop > 0) sum(V_next[target_indices]) else 0
    coverage_target_value <- if (target_pop > 0) cumulative_target_vacc / target_pop else 0
    
    # (For the variables that require tracking across weeks, we return NA here.)
    vacc_delayed                <- rep(NA, n)
    total_available_supply_age  <- rep(NA, n)
    total_supply_age            <- rep(NA, n)
    unused_supply_age           <- rep(NA, n)
    vacc_start_week             <- rep(NA, n)
    vacc_end_week               <- rep(NA, n)
    real_vaccine_duration_age   <- rep(NA, n)
    raw_allocation_age          <- matrix(vacc_allocation, nrow = n, ncol = 1)
    
    # Record Beta and FOI used.
    Beta_store <- Beta_t
    FOI_used   <- FOI_target
    
    #---------------------------------------------------------------
    # (8) Return the “derivative” vector (the difference over one week) 
    #     as well as all the extra outputs.
    list(
      c(dS, dI, dR, dV),
      S = S_next,
      I = I_next,
      R = R_next,
      V = V_next,
      vacc_delayed = vacc_delayed,
      total_supply = total_supply,
      total_available_supply_age = total_available_supply_age,
      total_supply_age = total_supply_age,
      unused_supply_age = unused_supply_age,
      coverage_switch = coverage_switch_bool,
      switch_day = NA,
      age_strat_cases = age_strat_cases,
      coverage_frac = coverage_frac,
      coverage_target = coverage_target_value,
      Beta_store = Beta_store,
      FOI_used = FOI_used,
      FOI_input = FOI_vec,
      vacc_start_week = vacc_start_week,
      vacc_end_week = vacc_end_week,
      real_vaccine_duration_age = real_vaccine_duration_age,
      raw_allocation_age = raw_allocation_age
    )
  })
}

params <- list(
  T = 52*n_year, 
  A = 18,
  N = N_2023,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  FOI_vec = FOI_weekly,
  I0_draw = rep(100, 18),
  R0  = 1 - exp(-0.01 * age_groups),
  rho = 0.3,
  gamma = 1/7, # 7 days 
  delay = 6,                   # Vaccination starts from Week x
  VE_block = 0.75,                 # Vaccine efficacy
  target_age = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0.6,
  weekly_delivery_speed = 0.1
)

params_novacc <- list(
    T = 52*n_year, 
    A = 18,
    N = N_2023,
    r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
    
    FOI_vec = FOI_weekly,
    I0_draw = rep(100, 18),
    R0  = 1 - exp(-0.01 * age_groups),
    rho = 0.3,
    gamma = 1/7, # 7 days 
    delay = 53,                   # Vaccination starts from Week x
    VE_block = 0,                 # Vaccine efficacy
    target_age = c(rep(0,18)),    # Targeting specific age groups
    coverage_threshold = 0,
    total_coverage = 0,
    weekly_delivery_speed = 0
)


I0  <- I0
S0  <- N_2023 - I0
R0  <- c(rep(0,18))
V0  <- c(rep(0,18))

state0 <- c(S0, I0, R0, V0)
names(state0) <- c(
  paste0("S",  seq_len(18)),  # e.g. "S1"..."S18"
  paste0("I",  seq_len(18)),
  paste0("R",  seq_len(18)),
  paste0("V",  seq_len(18))
)
# Times
times <- seq(1, 52*3, by=1)  

# Solve with deSolve
library(deSolve)
out <- ode(y = state0, times = times,
           func = sirv_ode_coverageSwitch,
           parms = params,
           method = "euler")
out <- as.data.frame(out)


out_novacc <- ode(y = state0, times = times,
           func = sirv_ode_coverageSwitch,
           parms = params_novacc,
           method = "euler")
out_novacc <- as.data.frame(out_novacc)

# indexing age_strat_cases
age_strat_cases <- out %>% dplyr::select(starts_with("age_strat_"))
age_strat_cases_novacc <- out_novacc %>% dplyr::select(starts_with("age_strat_"))

# reshape
age_strat_long <- age_strat_cases %>% pivot_longer(
  cols = 1:18,
  names_to = "Age",
  values_to = "cases"
)

age_strat_long$age <- rep(c(1:18), 52 * n_year)
age_strat_long$week <- rep(1: (52 * n_year), each = 18)
age_strat_long$scenario <- "Post-vaccination"

# reshape
age_strat_long_novacc <- age_strat_cases_novacc %>% pivot_longer(
  cols = 1:18,
  names_to = "Age",
  values_to = "cases"
)

age_strat_long_novacc$age <- rep(c(1:18), 52 * n_year)
age_strat_long_novacc$week <- rep(1: (52* n_year), each = 18)
age_strat_long_novacc$scenario <- "Pre-vaccination"

ggplot()+
  geom_line(data = age_strat_long, aes(x = week, y = cases, color = scenario))+
  geom_line(data = age_strat_long_novacc, aes(x = week, y = cases, color = scenario))+
  facet_wrap(~age)

