## old sim_function_final

## for owsa
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
  sigma,
  # Vaccination parameters
  delay,                 # the week index to START physical vaccination
  coverage_threshold,    # e.g., 0.7 => 70%
  target_age,            # 0/1 vector of length A
  total_coverage,        # fraction of total population eventually vaccinated
  weekly_delivery_speed, # fraction of total supply used per week
  VE_block,              # vaccine efficacy
  VE_inf,
  time_until_immunity, 
  # NEW: toggle coverage switch
  use_coverage_switch = FALSE  # set FALSE to disable Phase 2 reallocation
){
  # ───────────────────────────────────────
  # 도움 함수 및 상수
  # ───────────────────────────────────────
  pmax0      <- function(x) pmax(0, x)
  total_pop  <- sum(N)
  R0_scalar  <- base_beta[1] / gamma
  R0_vec     <- base_beta / gamma
  R_eff_vec  <- numeric(T)
  
  # ───────────────────────────────────────
  # 상태 행렬 초기화
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
  # 백신 공급량 계산
  # ───────────────────────────────────────
  total_supply       <- target_pop * total_coverage
  weekly_dose_total  <- total_supply * weekly_delivery_speed
  
  total_avail_age <- unused_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  # ───────────────────────────────────────
  # 초기 상태 (t = 1)
  # ───────────────────────────────────────
  S[, 1] <- N - I0_draw - R0 * N
  E[, 1] <- 0
  I[, 1] <- I0_draw
  R_[, 1] <- R0 * N
  V[, 1] <- 0
  
  coverage_frac[, 1] <- 0
  age_stratified_cases_raw[, 1] <- I[, 1]
  age_stratified_cases[,  1]    <- rho * I[, 1]
  true_symptomatic[, 1] <- rho * I[, 1]   # t=1에는 아직 백신효과 X
  coverage_target[1]   <- 0
  eff_S0               <- sum(S[, 1])     # ★ 완전 차단형이므로 그대로 S 합계
  R_eff_vec[1]         <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  # ───────────────────────────────────────
  # 메인 루프 (t = 2 … T)
  # ───────────────────────────────────────
  for (t_ in 2:T) {
    
    # ───── 1) 2주 전 접종분 → 면역 발현 ─────
    if (t_ - 2 >= 1) {
      effective_dose      <- vacc_delayed[, t_ - 2]              # 2주 전 S에 투입된 백신
      immunized           <- round(VE_inf * effective_dose)      # ★ 완전 면역 성공자
      V[, t_]             <- V[, t_ - 1] + immunized             # ★ V 누적
      V_covered[, t_]     <- V_covered[, t_ - 1] + effective_dose # ★ VE block track
      #S[, t_ - 1]         <- pmax0(S[, t_ - 1] - immunized)      # ★ S에서 제거
    } else {
      immunized       <- rep(0, A)
      V[, t_] <- V[, t_ - 1]
      V_covered[, t_] <- V_covered[, t_ - 1]
    }
    
    coverage_frac[, t_] <- V_covered[, t_] / N                           # ★ 완전 면역 비율
    
    # ───── 2) FOI 계산 ─────
    eff_S_prev_vec <- pmax0(S[, t_-1] - immunized) 
    phi             <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / total_pop)
    phi_vec[t_]     <- phi
    FOI_times_S[t_] <- phi * sum(eff_S_prev_vec)
    
    # ───── 3) 주간 할당 (delay 이후) ─────
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
          vacc_to_S   <- round(alloc * prop_S)        # S에게만 투입
          wasted      <- alloc - vacc_to_S            # I·R분은 폐기
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
    
    # ───── 4) SEIRV 업데이트 ─────
    for (a in seq_len(A)) {
      S_prev <- S[a, t_ - 1];  E_prev <- E[a, t_ - 1]
      I_prev <- I[a, t_ - 1];  R_prev <- R_[a, t_ - 1]
      
      immunized_a <- immunized[a]                        # t_ 시점에 면역 발현
      S_prev_eff <- pmax0(S_prev - immunized_a) 
      
      new_e  <- phi * S_prev_eff                     # ★ VE_inf 곱 제거 (완전 차단)
      new_i  <- sigma * E_prev
      recov  <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # 증상자 
      true_symptomatic[a, t_] <- 0.5242478 *
        I[a, t_] *
        (1 - VE_block * coverage_frac[a, t_]) *
        rho
      
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
      
      age_stratified_cases_raw[a, t_] <- I[a, t_]
      age_stratified_cases[a,     t_] <- rho * I[a, t_]
    }
    
    # ───── 5) 유효 감수자 / R_eff ─────
    eff_S              <- sum(S[, t_])          # ★ 완전 차단형 → 그냥 S 합계
    R_eff_vec[t_]      <- (base_beta[t_] / gamma) * (eff_S / total_pop)
  }
  
  # ───────────────────────────────────────
  # 출력
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

## same dose model
sirv_sim_coverageSwitch <- function(
    #─────────────────────────────────────────────────────────────────────────#
  # 1) Model dimensions and transmission parameters
  #─────────────────────────────────────────────────────────────────────────#
  T, A,            # T = number of weeks, A = number of age groups (e.g., 18)
  N, r,            # age‐specific population vector (length A), aging rates vector (length A)
  
  base_beta,       # transmission rate time series (length T)
  I0_draw, R0,     # initial infected (vector length A), initial immune proportion (vector length A)
  rho, gamma,      # detection probability, recovery rate
  
  #─────────────────────────────────────────────────────────────────────────#
  # 2) Vaccine parameters
  #─────────────────────────────────────────────────────────────────────────#
  delay,                 # vaccine rollout start week (integer)
  coverage_threshold,    # coverage threshold for Phase 2 switch (e.g., 0.7)
  target_age,            # 0/1 vector (length A): which age groups to target
  total_coverage,        # final vaccination proportion of total population (e.g., 0.3 = 30% of population doses)
  weekly_delivery_speed, # fraction of total supply delivered per week (e.g., 0.05 = 5% per week)
  VE_block,              # vaccine efficacy (conversion to immune, e.g., 0.988)
  
  #─────────────────────────────────────────────────────────────────────────#
  # 3) Phase 2 toggle
  #─────────────────────────────────────────────────────────────────────────#
  use_coverage_switch = FALSE  # if FALSE, skip Phase 2 (do not reallocate to non‐targets)
){
  # ───────────────────────────────────────────────────────────────────────#
  # 1) Initialize
  # ───────────────────────────────────────────────────────────────────────#
  S  <- matrix(0, nrow = A, ncol = T)  # Susceptible[a, t]
  I  <- matrix(0, nrow = A, ncol = T)  # Infected[a, t]
  R_ <- matrix(0, nrow = A, ncol = T)  # Recovered[a, t]
  V  <- matrix(0, nrow = A, ncol = T)  # Vaccinated (converted to immune)[a, t]
  
  # “Scheduled” vaccine: dose given at week t → moves from S to V at t+2
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  
  # For recording outputs
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)  # age‐specific detected cases
  age_stratified_cases_raw <- matrix(0, nrow = A, ncol = T) # raw number of infection
  true_symptomatic        <- matrix(0, nrow = A, ncol = T) # true number of symp
  raw_allocation_age   <- matrix(0, nrow = A, ncol = T)  # age‐specific weekly allocations
  coverage_frac        <- matrix(0, nrow = A, ncol = T)  # age‐specific coverage fraction
  coverage_target      <- numeric(T)                     # target group coverage
  phi_vec              <- numeric(T)                     # force of infection each week
  
  # Identify target age indices
  target_indices <- which(target_age == 1)
  target_pop     <- sum(N[target_indices])
  no_target      <- (length(target_indices) == 0 || target_pop == 0)
  
  # ───────────────────────────────────────────────────────────────────────#
  # 1A) Compute total vaccine supply (same across scenarios)
  # ───────────────────────────────────────────────────────────────────────#
  total_supply <- sum(N) * total_coverage

  # ───────────────────────────────────────────────────────────────────────#
  # 1B) Compute weekly Phase 1 allocation
  # ───────────────────────────────────────────────────────────────────────#
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  # “Equally reserved” supply for each target age group
  total_available_supply_age <- numeric(A)
  if (!no_target) {
    for (a in target_indices) {
      total_available_supply_age[a] <- total_supply / length(target_indices)
    }
  }
  total_supply_age  <- numeric(A)  # cumulative doses given per age
  unused_supply_age <- total_available_supply_age  # remaining supply per age
  
  # Phase 2 toggle & switch week record
  coverage_switch <- FALSE
  switch_day      <- NA
  
  # ───────────────────────────────────────────────────────────────────────#
  # 1C) Initial state at t = 1
  # ───────────────────────────────────────────────────────────────────────#
  vacc_start_week <- rep(NA, A)
  vacc_end_week   <- rep(NA, A)
  
  for (a in seq_len(A)) {
    R_[a, 1] <- R0[a] * N[a]
    I[a, 1]  <- I0_draw[a]
    S[a, 1]  <- N[a] - I[a, 1] - R_[a, 1]
    V[a, 1]  <- 0
    
    age_stratified_cases[a, 1] <- rho * I[a, 1]
    coverage_frac[a, 1]       <- 0
  }
  coverage_target[1] <- 0
  
  # ───────────────────────────────────────────────────────────────────────#
  # 2) Main loop: t = 2..T
  # ───────────────────────────────────────────────────────────────────────#
  for (t_ in 2:T) {
    # ============================================================
    # 2A) Apply vaccine effect delayed by 2 weeks
    # ============================================================
    if (t_ - 2 >= 1) {
      for (a in seq_len(A)) {
        newly_effective <- VE_block * vacc_delayed[a, t_ - 2]
        actual_move     <- min(S[a, t_ - 1], newly_effective)
        S[a, t_ - 1]    <- S[a, t_ - 1] - actual_move
        V[a, t_ - 1]    <- V[a, t_ - 1] + actual_move
      }
    }
    
    # ============================================================
    # 2B) Update cumulative coverage in target group
    # ============================================================
    if (!no_target) {
      cumulative_vacc       <- sum(raw_allocation_age[target_indices, 1:(t_ - 1)])
      coverage_target[t_]   <- cumulative_vacc / target_pop
    } else {
      coverage_target[t_] <- 0
    }
    
    # ============================================================
    # 2C) Compute force of infection
    # ============================================================
    infected_prev <- sum(I[, t_ - 1])
    region_pop    <- sum(N)
    phi           <- base_beta[t_ - 1] * (infected_prev / region_pop)
    phi_vec[t_]   <- phi
    
    # ============================================================
    # 2D) Vaccine allocation per week: Phase 1 → Phase 2 logic
    #     (Modify so that t_ == delay 시점부터 투여 시작)
    # ============================================================
    if (t_ >= delay) {
      weekly_remaining <- weekly_dose_total
      
      # “Compute how many scheduled doses are not yet effective”
      pending_matrix <- vacc_delayed[, 1:(t_ - 1), drop = FALSE]
      if (t_ - 2 >= 1) {
        u_max <- t_ - 2
        effective_cumulative <- rowSums(vacc_delayed[, 1:u_max, drop = FALSE]) * VE_block
      } else {
        effective_cumulative <- rep(0, A)
      }
      pending_not_effective <- rowSums(pending_matrix) - effective_cumulative
      pending_not_effective <- pmax(0, pending_not_effective)
      
      # ─── (1) Phase 1: if coverage_switch == FALSE, allocate only to target ages
      if (!coverage_switch) {
        rem_targets <- target_indices[ unused_supply_age[target_indices] > 0 ]
        
        if (length(rem_targets) > 0) {
          # (a) Equal share among remaining target ages
          equal_share <- floor(weekly_remaining / length(rem_targets))
          allocated_this_round <- integer(length(rem_targets))
          
          for (i in seq_along(rem_targets)) {
            a <- rem_targets[i]
            sus_true       <- pmax(0, S[a, t_ - 1] - pending_not_effective[a])
            rem_supply_age <- unused_supply_age[a]
            feasible <- min(equal_share, sus_true, rem_supply_age)
            if (feasible > 0) {
              raw_allocation_age[a, t_] <- feasible
              vacc_delayed[a, t_]       <- feasible
              total_supply_age[a]       <- total_supply_age[a] + feasible
              unused_supply_age[a]      <- unused_supply_age[a] - feasible
              weekly_remaining          <- weekly_remaining - feasible
              allocated_this_round[i]   <- feasible
              if (is.na(vacc_start_week[a])) {
                vacc_start_week[a] <- t_
              }
              vacc_end_week[a] <- t_
            }
          }
          
          # (b) If leftover remains, redistribute 1 dose at a time among remaining ages
          if (weekly_remaining > 0) {
            rem_targets2 <- rem_targets[
              (unused_supply_age[rem_targets] > 0) &
                ((S[rem_targets, t_ - 1] - pending_not_effective[rem_targets]) > 0)
            ]
            j <- 1
            while (weekly_remaining > 0 && length(rem_targets2) > 0) {
              a <- rem_targets2[j]
              sus_true       <- pmax(0, S[a, t_ - 1] - pending_not_effective[a])
              rem_supply_age <- unused_supply_age[a]
              feasible <- min(1, sus_true, rem_supply_age, weekly_remaining)
              if (feasible > 0) {
                raw_allocation_age[a, t_] <- raw_allocation_age[a, t_] + feasible
                vacc_delayed[a, t_]       <- vacc_delayed[a, t_] + feasible
                total_supply_age[a]       <- total_supply_age[a] + feasible
                unused_supply_age[a]      <- unused_supply_age[a] - feasible
                weekly_remaining          <- weekly_remaining - feasible
                if (is.na(vacc_start_week[a])) {
                  vacc_start_week[a] <- t_
                }
                vacc_end_week[a] <- t_
              }
              j <- j + 1
              if (j > length(rem_targets2)) {
                rem_targets2 <- rem_targets2[
                  (unused_supply_age[rem_targets2] > 0) &
                    ((S[rem_targets2, t_ - 1] - pending_not_effective[rem_targets2]) > 0)
                ]
                if (length(rem_targets2) == 0) break
                j <- 1
              }
            }
          }
        }
        
      } else {
        # ─── (2) Phase 2: if coverage_switch == TRUE, allocate to all ages with unused supply
        rem_all <- which(unused_supply_age > 0)
        if (length(rem_all) > 0) {
          equal_share2 <- floor(weekly_remaining / length(rem_all))
          for (i in seq_along(rem_all)) {
            a <- rem_all[i]
            sus_true       <- pmax(0, S[a, t_ - 1] - pending_not_effective[a])
            rem_supply_age <- unused_supply_age[a]
            feasible <- min(equal_share2, sus_true, rem_supply_age)
            if (feasible > 0) {
              raw_allocation_age[a, t_] <- feasible
              vacc_delayed[a, t_]       <- feasible
              total_supply_age[a]       <- total_supply_age[a] + feasible
              unused_supply_age[a]      <- unused_supply_age[a] - feasible
              weekly_remaining          <- weekly_remaining - feasible
              if (is.na(vacc_start_week[a])) {
                vacc_start_week[a] <- t_
              }
              vacc_end_week[a] <- t_
            }
          }
          # Redistribute leftover
          if (weekly_remaining > 0) {
            rem_all2 <- rem_all[
              (unused_supply_age[rem_all] > 0) &
                ((S[rem_all, t_ - 1] - pending_not_effective[rem_all]) > 0)
            ]
            k <- 1
            while (weekly_remaining > 0 && length(rem_all2) > 0) {
              a <- rem_all2[k]
              sus_true       <- pmax(0, S[a, t_ - 1] - pending_not_effective[a])
              rem_supply_age <- unused_supply_age[a]
              feasible <- min(1, sus_true, rem_supply_age, weekly_remaining)
              if (feasible > 0) {
                raw_allocation_age[a, t_] <- raw_allocation_age[a, t_] + feasible
                vacc_delayed[a, t_]       <- vacc_delayed[a, t_] + feasible
                total_supply_age[a]       <- total_supply_age[a] + feasible
                unused_supply_age[a]      <- unused_supply_age[a] - feasible
                weekly_remaining          <- weekly_remaining - feasible
                if (is.na(vacc_start_week[a])) {
                  vacc_start_week[a] <- t_
                }
                vacc_end_week[a] <- t_
              }
              k <- k + 1
              if (k > length(rem_all2)) {
                rem_all2 <- rem_all2[
                  (unused_supply_age[rem_all2] > 0) &
                    ((S[rem_all2, t_ - 1] - pending_not_effective[rem_all2]) > 0)
                ]
                if (length(rem_all2) == 0) break
                k <- 1
              }
            }
          }
        }
      }
    }
    # (t_ < delay)인 경우 raw_allocation_age와 vacc_delayed가 모두 0인 상태로 유지됩니다.
    
    # ============================================================
    # 2E) Update S‐I‐R‐V states at time t_
    # ============================================================
    for (a in seq_len(A)) {
      if (a == 1) {
        S[a, t_]  <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1]
        I[a, t_]  <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1]
        V[a, t_]  <- V[a, t_ - 1] - r[a] * V[a, t_ - 1]
      } else {
        S[a, t_]  <- S[a, t_ - 1] - phi * S[a, t_ - 1] - r[a] * S[a, t_ - 1] + r[a - 1] * S[a - 1, t_ - 1]
        I[a, t_]  <- I[a, t_ - 1] + phi * S[a, t_ - 1] - gamma * I[a, t_ - 1] - r[a] * I[a, t_ - 1] + r[a - 1] * I[a - 1, t_ - 1]
        R_[a, t_] <- R_[a, t_ - 1] + gamma * I[a, t_ - 1] - r[a] * R_[a, t_ - 1] + r[a - 1] * R_[a - 1, t_ - 1]
        V[a, t_]  <- V[a, t_ - 1] - r[a] * V[a, t_ - 1] + r[a - 1] * V[a - 1, t_ - 1]
      }
      S[a, t_]  <- pmax(0, S[a, t_])
      I[a, t_]  <- pmax(0, I[a, t_])
      R_[a, t_] <- pmax(0, R_[a, t_])
      V[a, t_]  <- pmax(0, V[a, t_])
    }
    
    # ============================================================
    # 2F) Record coverage & cases
    # ============================================================
    for (a in seq_len(A)) {
      coverage_frac[a, t_]       <- V[a, t_] / N[a]
      age_stratified_cases[a, t_] <- rho * I[a, t_]
      age_stratified_cases_raw[a, t_] <-  I[a, t_]
      true_symptomatic[a, t_] <-  age_stratified_cases_raw[a, t_] * 0.5242478
    }
    
    # ============================================================
    # 2G) Determine Phase 2 entry
    # ============================================================
    if (use_coverage_switch && !coverage_switch && !no_target) {
      if (coverage_target[t_] >= coverage_threshold) {
        coverage_switch <- TRUE
        switch_day      <- t_
      }
    }
  }  # end of main loop
  
  # ───────────────────────────────────────────────────────────────────────#
  # 3) Compute actual vaccination duration per age (start → end week)
  # ───────────────────────────────────────────────────────────────────────#
  real_vaccine_duration_age <- numeric(A)
  for (a in seq_len(A)) {
    if (!is.na(vacc_start_week[a]) && !is.na(vacc_end_week[a])) {
      real_vaccine_duration_age[a] <- vacc_end_week[a] - vacc_start_week[a] + 1
    } else {
      real_vaccine_duration_age[a] <- 0
    }
  }
  
  # ───────────────────────────────────────────────────────────────────────#
  # 4) Return results
  # ───────────────────────────────────────────────────────────────────────#
  return(list(
    S                          = S,
    I                          = I,
    R                          = R_,
    V                          = V,
    vacc_delayed               = vacc_delayed,
    total_supply               = total_supply,
    total_available_supply_age = total_available_supply_age,
    total_supply_age           = total_supply_age,
    unused_supply_age          = unused_supply_age,
    coverage_switch            = coverage_switch,
    switch_day                 = switch_day,
    coverage_target            = coverage_target,
    age_stratified_cases       = age_stratified_cases,
    age_stratified_cases_raw   = age_stratified_cases_raw,
    true_symptomatic           = true_symptomatic,
    coverage_frac              = coverage_frac,
    phi                        = phi_vec,
    vacc_start_week            = vacc_start_week,
    vacc_end_week              = vacc_end_week,
    real_vaccine_duration_age  = real_vaccine_duration_age,
    raw_allocation_age         = raw_allocation_age
  ))
}

# with immune delay 
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
){
  pmax0 <- function(x) pmax(0, x)
  
  ## 1) 초기화
  S <- I <- R_ <- V <- matrix(0, A, T)
  vacc_delayed   <- matrix(0, A, T)
  raw_allocation_age <- matrix(0, A, T)
  wasted_dose        <- matrix(0, A, T)
  unvaccinated       <- N
  coverage_frac <- matrix(0, A, T)
  age_stratified_cases     <- matrix(0, A, T)
  age_stratified_cases_raw <- matrix(0, A, T)
  true_symptomatic         <- matrix(0, A, T)
  phi_vec         <- numeric(T)
  FOI_times_S <- numeric(T)
  coverage_target <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  ## ↓ 수정된 부분: “시나리오 간 동일 도즈”를 위해 전체 인구 기준으로 총도즈 계산
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  ## 여전히 “그 총도즈를 타깃 연령대별 인구 비율로 배분”하되
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  ## 초기 상태 (t = 1)
  S[,1] <- N - I0_draw - R0*N
  I[,1] <- I0_draw
  R_[,1] <- R0*N
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1]         <- I[,1] * 0.5242478
  coverage_frac[,1] <- 0
  coverage_target[1] <- 0
  
  ## 2) 메인 루프
  for(t_ in 2:T){
    ## 2A. 백신 효과(2주 지연)
    if(t_ - 2 >= 1){
      for(a in seq_len(A)){
        eff  <- VE_block * vacc_delayed[a, t_-2]
        move <- min(S[a,t_-1], eff)
        S[a,t_-1] <- S[a,t_-1] - move
        V[a,t_-1] <- V[a,t_-1] + move
      }
    }
    
    
    ## 2B. FOI
    phi <- base_beta[t_-1] * (sum(I[,t_-1]) / sum(N))
    phi_vec[t_] <- phi
    total_S      <- sum(S[, t_ - 1])
    FOI_times_S[t_] <- phi * total_S
    
    ## 2C. 백신 분배 (중복 금지 + wasted 계산)
    if(t_ >= delay && target_pop > 0){
      remaining <- weekly_dose_total
      for(a in target_idx){
        base_alloc <- ceiling(weekly_dose_total * N[a] / target_pop)
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        
        if(alloc > 0){
          avail_S   <- S[a, t_-1]
          effective <- min(alloc, avail_S)
          wasted    <- alloc - effective
          
          raw_allocation_age[a, t_] <- alloc
          wasted_dose[a, t_]        <- wasted
          if(effective > 0) vacc_delayed[a, t_] <- effective
          
          total_used_age[a] <- total_used_age[a] + alloc
          unused_age[a]     <- total_avail_age[a] - total_used_age[a]
          unvaccinated[a]   <- unvaccinated[a] - alloc
          remaining         <- remaining - alloc
          
          if(is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_-1]
    }
    
    ## 2D. S-I-R-V 업데이트
    for(a in seq_len(A)){
      S_prev <- S[a,t_-1]; I_prev <- I[a,t_-1]; R_prev <- R_[a,t_-1]; V_prev <- V[a,t_-1]
      new_inf <- phi * S_prev; recov <- gamma * I_prev
      S_next <- S_prev - new_inf; I_next <- I_prev + new_inf - recov
      R_next <- R_prev + recov;  V_next <- V_prev
      
      if(a==1){
        S[a,t_] <- S_next - r[a]*S_next
        I[a,t_] <- I_next - r[a]*I_next
        R_[a,t_] <- R_next - r[a]*R_next
        V[a,t_] <- V_next - r[a]*V_next
      } else {
        S[a,t_] <- S_next - r[a]*S_next + r[a-1]*S[a-1,t_-1]
        I[a,t_] <- I_next - r[a]*I_next + r[a-1]*I[a-1,t_-1]
        R_[a,t_] <- R_next - r[a]*R_next + r[a-1]*R_[a-1,t_-1]
        V[a,t_] <- V_next - r[a]*V_next + r[a-1]*V[a-1,t_-1]
      }
      S[a,t_] <- pmax0(S[a,t_]); I[a,t_] <- pmax0(I[a,t_])
      R_[a,t_] <- pmax0(R_[a,t_]); V[a,t_] <- pmax0(V[a,t_])
    }
    
    ## 2E. 기록
    coverage_frac[,t_]            <- V[,t_] / N
    age_stratified_cases_raw[,t_] <- I[,t_]
    age_stratified_cases[,t_]     <- rho * I[,t_]
    true_symptomatic[,t_]         <- I[,t_] * 0.5242478
  }
  
  ## 3) 접종 기간
  real_vac_dur <- ifelse(is.na(vacc_start)|is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) 반환
  list(
    S = S, I = I, R = R_, V = V,
    vacc_delayed            = vacc_delayed,
    raw_allocation_age      = raw_allocation_age,
    wasted_dose             = wasted_dose,
    total_supply            = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age        = total_used_age,   # 효과 유무 상관없이 배정된 모든 도즈
    unused_supply_age       = unused_age,
    coverage_target         = coverage_target,
    coverage_frac           = coverage_frac,
    phi                     = phi_vec,
    FOI_times_S             = FOI_times_S,
    age_stratified_cases    = age_stratified_cases,
    age_stratified_cases_raw= age_stratified_cases_raw,
    true_symptomatic        = true_symptomatic,
    vacc_start_week         = vacc_start,
    vacc_end_week           = vacc_end,
    real_vaccine_duration_age = real_vac_dur
  )
}


# no delay in immunity
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  pmax0 <- function(x) pmax(0, x)
  
  ## 1) 초기화
  S <- I <- R_ <- V <- matrix(0, A, T)
  vacc_delayed       <- matrix(0, A, T)
  raw_allocation_age <- matrix(0, A, T)
  wasted_dose        <- matrix(0, A, T)
  effective_dose     <- matrix(0, A, T)
  unvaccinated       <- N
  coverage_frac      <- matrix(0, A, T)
  age_stratified_cases     <- matrix(0, A, T)
  age_stratified_cases_raw <- matrix(0, A, T)
  true_symptomatic         <- matrix(0, A, T)
  phi_vec         <- numeric(T)
  coverage_target  <- numeric(T)
  
  target_idx  <- which(target_age == 1)
  target_pop  <- sum(N[target_idx])
  
  ## ↓ “시나리오 간 동일 도즈”를 위해 전체 인구 기준으로 총도즈 계산
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  ## “그 총도즈를 타깃 연령대별 인구 비율로 배분”하되
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start     <- vacc_end <- rep(NA_integer_, A)
  
  ## 초기 상태 (t = 1)
  S[,1] <- N - I0_draw - R0 * N
  I[,1] <- I0_draw
  R_[,1] <- R0 * N
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1]         <- I[,1] * 0.5242478
  coverage_frac[,1] <- 0
  coverage_target[1] <- 0
  
  ## 2) 메인 루프
  for (t_ in 2:T) {
    ## ───────────────────────────────────────────────
    ## 2A. 백신 분배 및 즉시 효과 적용
    ## ───────────────────────────────────────────────
    if (t_ >= delay && target_pop > 0) {
      remaining <- weekly_dose_total
      
      for (a in target_idx) {
        base_alloc <- ceiling(weekly_dose_total * N[a] / target_pop)
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        
        if (alloc > 0) {
          S_prev  <- S[a, t_-1]
          effective <- min(alloc, S_prev)
          wasted    <- alloc - effective
          
          # 배분 기록
          raw_allocation_age[a, t_] <- alloc
          wasted_dose[a, t_]        <- wasted
          effective_dose[a, t_]     <- effective
          # 즉시 효과: S → V
          S[a, t_-1] <- S_prev - effective
          V[a, t_-1] <- V[a, t_-1] + effective
          
          # 누적 사용량, 미사용량, 미접종 인구 업데이트
          total_used_age[a] <- total_used_age[a] + alloc
          unused_age[a]     <- total_avail_age[a] - total_used_age[a]
          unvaccinated[a]   <- unvaccinated[a] - alloc
          remaining         <- remaining - alloc
          
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_-1]
    }
    
    ## ───────────────────────────────────────────────
    ## 2B. FOI (Force of Infection) 계산
    ## ───────────────────────────────────────────────
    phi <- base_beta[t_-1] * (sum(I[, t_-1]) / sum(N))
    phi_vec[t_] <- phi
    
    ## ───────────────────────────────────────────────
    ## 2C. S-I-R-V 업데이트 (전염·회복·노령화)
    ## ───────────────────────────────────────────────
    for (a in seq_len(A)) {
      S_prev <- S[a, t_-1]
      I_prev <- I[a, t_-1]
      R_prev <- R_[a, t_-1]
      V_prev <- V[a, t_-1]
      
      new_inf <- phi * S_prev
      recov   <- gamma * I_prev
      
      S_next <- S_prev - new_inf
      I_next <- I_prev + new_inf - recov
      R_next <- R_prev + recov
      V_next <- V_prev  # 백신 효과는 이미 위에서 반영됨
      
      if (a == 1) {
        S[a, t_]  <- S_next - r[a] * S_next
        I[a, t_]  <- I_next - r[a] * I_next
        R_[a, t_] <- R_next - r[a] * R_next
        V[a, t_]  <- V_next - r[a] * V_next
      } else {
        S[a, t_]  <- S_next - r[a] * S_next + r[a-1] * S[a-1, t_-1]
        I[a, t_]  <- I_next - r[a] * I_next + r[a-1] * I[a-1, t_-1]
        R_[a, t_] <- R_next - r[a] * R_next + r[a-1] * R_[a-1, t_-1]
        V[a, t_]  <- V_next - r[a] * V_next + r[a-1] * V[a-1, t_-1]
      }
      
      S[a, t_]  <- pmax0(S[a, t_])
      I[a, t_]  <- pmax0(I[a, t_])
      R_[a, t_] <- pmax0(R_[a, t_])
      V[a, t_]  <- pmax0(V[a, t_])
    }
    
    ## ───────────────────────────────────────────────
    ## 2D. 기록
    ## ───────────────────────────────────────────────
    coverage_frac[, t_]            <- V[, t_] / N
    age_stratified_cases_raw[, t_] <- I[, t_]
    age_stratified_cases[, t_]     <- rho * I[, t_]
    true_symptomatic[, t_]         <- I[, t_] * 0.5242478
  }
  
  ## 3) 접종 기간 계산
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) 반환
  list(
    S                         = S,
    I                         = I,
    R                         = R_,
    V                         = V,
    vacc_delayed              = vacc_delayed,
    raw_allocation_age        = raw_allocation_age,
    wasted_dose               = wasted_dose,
    effective_dose            = effective_dose,
    total_supply              = total_supply,
    total_available_supply_age= total_avail_age,
    total_supply_age          = total_used_age,
    unused_supply_age         = unused_age,
    coverage_target           = coverage_target,
    coverage_frac             = coverage_frac,
    phi                       = phi_vec,
    age_stratified_cases      = age_stratified_cases,
    age_stratified_cases_raw  = age_stratified_cases_raw,
    true_symptomatic          = true_symptomatic,
    vacc_start_week           = vacc_start,
    vacc_end_week             = vacc_end,
    real_vaccine_duration_age = real_vac_dur
  )
}


# final version.
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
){
  pmax0 <- function(x) pmax(0, x)
  
  ## 1) initialisation
  S  <- I  <- R_ <- V   <- matrix(0, A, T)
  vacc_delayed     <- matrix(0, A, T)
  raw_allocation_age <- matrix(0, A, T)
  wasted_dose        <- matrix(0, A, T) ## ineffective doses 
  unvaccinated       <- N
  coverage_frac      <- matrix(0, A, T)
  age_stratified_cases     <- matrix(0, A, T)
  age_stratified_cases_raw <- matrix(0, A, T)
  true_symptomatic         <- matrix(0, A, T)
  phi_vec         <- numeric(T)
  FOI_times_S     <- numeric(T)
  coverage_target <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  ## ↓ 수정된 부분: “시나리오 간 동일 도즈”를 위해 전체 인구 기준으로 총도즈 계산
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  ## “그 총도즈를 타깃 연령대별 인구 비율로 배분”
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  ## 초기 상태 (t = 1)
  S[,1] <- N - I0_draw - R0 * N
  I[,1] <- I0_draw
  R_[,1] <- R0 * N
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1]         <- I[,1] * 0.5242478
  coverage_frac[,1] <- 0
  coverage_target[1] <- 0
  
  ## 2) 메인 루프
  for(t_ in 2:T){
    ## 2A. 백신 효과(2주 지연)
    if(t_ - 4 >= 1){
      for(a in seq_len(A)){
        eff  <- VE_block * vacc_delayed[a, t_ - 4]
        move <- min(S[a, t_ - 1], eff)
        S[a, t_ - 1] <- S[a, t_ - 1] - move
        V[a, t_ - 1] <- V[a, t_ - 1] + move
      }
    }
    
    ## 2B. FOI 계산
    phi <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / sum(N))
    phi_vec[t_] <- phi
    total_S      <- sum(S[, t_ - 1])
    FOI_times_S[t_] <- phi * total_S
    
    ## 2C. 백신 분배 (무작위 배분 → S, I, R 비율에 맞춰 할당 + wasted 계산)
    if(t_ >= delay && target_pop > 0){
      remaining <- weekly_dose_total
      
      for(a in target_idx){
        # 1) 연령 a에 배정되어야 할 기본 할당량 (round 또는 ceiling/ floor 취향대로 조절)
        base_alloc <- ceiling(weekly_dose_total * (N[a] / target_pop))
        
        # 2) 실제 할당가능량을 남은 unvaccinated, unused_age, remaining을 고려해 제한
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        
        if(alloc > 0){
          ## ─┬─ 인구 전체 중에서 비율대로 분배 ───────────────
          S_prev <- S[a, t_ - 1]
          I_prev <- I[a, t_ - 1]
          R_prev <- R_[a, t_ - 1]
          N_a    <- N[a]
          
          # (1) S, I, R 비율 계산
          prop_S <- ifelse(N_a > 0, S_prev / N_a, 0)
          prop_I <- ifelse(N_a > 0, I_prev / N_a, 0)
          prop_R <- ifelse(N_a > 0, R_prev / N_a, 0)
          
          # (2) 할당량을 비율대로 배정 (round → 적절히 조절 가능)
          vacc_to_S <- round(alloc * prop_S)
          vacc_to_I <- round(alloc * prop_I)
          #  나머지는 R에게 할당 (납부 검증: 숫자 오차로 합이 alloc와 다를 수 있어 조정)
          vacc_to_R <- alloc - vacc_to_S - vacc_to_I
          
          # (3) 감수성자 내에서 가능한 만큼만 면역 후보로 둔다
          #     (S_prev보다 더 많이 배정될 수도 있으므로 제한)
          if(vacc_to_S > S_prev){
            vacc_to_S   <- S_prev
            # 남은 alloc에서 I와 R로 자동 재분배 (단순히 낭비 처리)
            vacc_to_I <- alloc - vacc_to_S - vacc_to_R
            # (round 로 인해 I, R의 합산이 바뀌었다면 추가 조정해 주세요)
          }
          
          # (4) “감염자·회복자에게 할당된 부분” = 낭비량
          wasted_IR <- vacc_to_I + vacc_to_R
          
          # └───────────────────────────────────────────────
          
          # 3) 기록
          raw_allocation_age[a, t_] <- vacc_to_S + vacc_to_I + vacc_to_R  # = alloc
          wasted_dose[a, t_]        <- wasted_IR
          
          # 4) 실제 면역 후보(감수성자 중)만 vacc_delayed로 이관
          if(vacc_to_S > 0){
            vacc_delayed[a, t_] <- vacc_to_S
          }
          
          # 5) 연령군 단위 재고 관리
          total_used_age[a] <- total_used_age[a] + alloc
          unused_age[a]     <- total_avail_age[a] - total_used_age[a]
          unvaccinated[a]   <- unvaccinated[a] - alloc
          remaining         <- remaining - alloc
          
          # 6) 접종 시작/종료 시점 기록
          if(is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      } # for(a)
      
      # 목표연령 커버리지 계산
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    ## 2D. S-I-R-V 상태 업데이트
    for(a in seq_len(A)){
      S_prev <- S[a, t_ - 1]
      I_prev <- I[a, t_ - 1]
      R_prev <- R_[a, t_ - 1]
      V_prev <- V[a, t_ - 1]
      
      new_inf <- phi * S_prev
      recov   <- gamma * I_prev
      
      # 기본 SIR 업데이트 (백신으로 인한 이동은 2주 지연된 vacc_delayed에서 이미 처리됨)
      S_next <- S_prev - new_inf
      I_next <- I_prev + new_inf - recov
      R_next <- R_prev + recov
      V_next <- V_prev  # 백신 이동(면역획득)은 2주 지연 블록에서 이미 반영됨
      
      # 연령 이동 (aging/death) 적용
      if(a == 1){
        S[a, t_]  <- S_next - r[a] * S_next
        I[a, t_]  <- I_next - r[a] * I_next
        R_[a, t_] <- R_next - r[a] * R_next
        V[a, t_]  <- V_next - r[a] * V_next
      } else {
        S[a, t_]  <- S_next - r[a] * S_next + r[a - 1] * S[a - 1, t_ - 1]
        I[a, t_]  <- I_next - r[a] * I_next + r[a - 1] * I[a - 1, t_ - 1]
        R_[a, t_] <- R_next - r[a] * R_next + r[a - 1] * R_[a - 1, t_ - 1]
        V[a, t_]  <- V_next - r[a] * V_next + r[a - 1] * V[a - 1, t_ - 1]
      }
      
      # 음수 방지
      S[a, t_]  <- pmax0(S[a, t_])
      I[a, t_]  <- pmax0(I[a, t_])
      R_[a, t_] <- pmax0(R_[a, t_])
      V[a, t_]  <- pmax0(V[a, t_])
    }
    
    ## 2E. 기록
    coverage_frac[, t_]            <- V[, t_] / N
    age_stratified_cases_raw[, t_] <- I[, t_]
    age_stratified_cases[, t_]     <- rho * I[, t_]
    true_symptomatic[, t_]         <- I[, t_] * 0.5242478
  } # for t_
  
  ## 3) 접종 기간
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) 결과 반환
  list(
    S  = S,  I  = I,  R  = R_,  V  = V,
    vacc_delayed               = vacc_delayed,
    raw_allocation_age         = raw_allocation_age,
    wasted_dose                = wasted_dose,
    total_supply               = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,   # 연령별로 할당된 총 도즈
    unused_supply_age          = unused_age,
    coverage_target            = coverage_target,
    coverage_frac              = coverage_frac,
    phi                        = phi_vec,
    FOI_times_S                = FOI_times_S,
    age_stratified_cases       = age_stratified_cases,
    age_stratified_cases_raw   = age_stratified_cases_raw,
    true_symptomatic           = true_symptomatic,
    vacc_start_week            = vacc_start,
    vacc_end_week              = vacc_end,
    real_vaccine_duration_age  = real_vac_dur
  )
}

# disease blocking 
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,               # disease-blocking efficacy
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  pmax0 <- function(x) pmax(0, x)
  
  ## 1) Initialization
  S  <- I  <- R_ <- V   <- matrix(0, A, T)
  vacc_delayed        <- matrix(0, A, T)
  raw_allocation_age  <- matrix(0, A, T)
  wasted_dose         <- matrix(0, A, T)
  unvaccinated        <- N
  coverage_frac       <- matrix(0, A, T)
  age_stratified_cases        <- matrix(0, A, T)
  age_stratified_cases_raw    <- matrix(0, A, T)
  true_symptomatic    <- matrix(0, A, T)
  phi_vec             <- numeric(T)
  FOI_times_S         <- numeric(T)
  coverage_target     <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  ## Compute total supply and weekly doses
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start     <- vacc_end <- rep(NA_integer_, A)
  
  ## Initial state at t = 1
  S[,1] <- N - I0_draw - R0 * N
  I[,1] <- I0_draw
  R_[,1] <- R0 * N
  V[,1] <- 0
  coverage_frac[,1] <- 0
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1] <- rho * I[,1] * (1 - VE_block)
  coverage_target[1] <- 0
  
  ## 2) Main loop
  for (t_ in 2:T) {
    ## 2A. Accumulate vaccine doses (for symptomatic blocking)
    if (t_ - 4 >= 1) {
      V[, t_] <- V[, t_ - 1] + vacc_delayed[, t_ - 4]
    } else {
      V[, t_] <- V[, t_ - 1]
    }
    coverage_frac[, t_] <- ifelse(N > 0, V[, t_] / N, 0)
    
    ## 2B. Force-of-Infection
    phi <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / sum(N))
    phi_vec[t_] <- phi
    FOI_times_S[t_] <- phi * sum(S[, t_ - 1])
    
    ## 2C. Vaccine allocation (same as original)
    if (t_ >= delay && target_pop > 0) {
      remaining <- weekly_dose_total
      for (a in target_idx) {
        base_alloc <- ceiling(weekly_dose_total * (N[a] / target_pop))
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        if (alloc > 0) {
          raw_allocation_age[a, t_] <- alloc
          S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
          prop_S <- ifelse(N[a] > 0, S_prev / N[a], 0)
          prop_I <- ifelse(N[a] > 0, I_prev / N[a], 0)
          vacc_to_S <- round(alloc * prop_S)
          vacc_to_I <- round(alloc * prop_I)
          vacc_to_R <- alloc - vacc_to_S - vacc_to_I
          wasted_IR <- vacc_to_I + vacc_to_R
          wasted_dose[a, t_]  <- wasted_IR
          vacc_delayed[a, t_] <- vacc_to_S
          total_used_age[a]   <- total_used_age[a] + alloc
          unused_age[a]       <- unused_age[a] - alloc
          unvaccinated[a]     <- unvaccinated[a] - alloc
          remaining           <- remaining - alloc
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    ## 2D. S-I-R-V update with disease-blocking
    for (a in seq_len(A)) {
      S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
      # infection unchanged
      new_inf <- phi * S_prev
      recov <- gamma * I_prev
      S[a, t_]  <- pmax0(S_prev - new_inf)
      I[a, t_]  <- pmax0(I_prev + new_inf - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # symptomatic blocking (true symp means initial Infection that goes to biologically symptomatic before reporting)
      true_symptomatic[a, t_] <- 0.5242478 * I[a, t_] * (1 - VE_block * coverage_frac[a, t_])
      
      # aging/death
      if (a == 1) {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_])
      } else {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_]  + r[a-1] * S[a-1, t_-1])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_]  + r[a-1] * I[a-1, t_-1])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_] + r[a-1] * R_[a-1, t_-1])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_]  + r[a-1] * V[a-1, t_-1])
      }
    }
    
    # track cases
    age_stratified_cases_raw[, t_] <- I[, t_]
    age_stratified_cases[, t_]     <- rho * I[, t_]
  }
  
  ## 3) vaccination duration
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) return
  list(
    S  = S, I  = I, R  = R_, V  = V,
    vacc_delayed               = vacc_delayed,
    raw_allocation_age         = raw_allocation_age,
    wasted_dose                = wasted_dose,
    total_supply               = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,
    unused_supply_age          = unused_age,
    coverage_target            = coverage_target,
    coverage_frac              = coverage_frac,
    phi                        = phi_vec,
    FOI_times_S                = FOI_times_S,
    age_stratified_cases       = age_stratified_cases,
    age_stratified_cases_raw   = age_stratified_cases_raw,
    true_symptomatic           = true_symptomatic,
    vacc_start_week            = vacc_start,
    vacc_end_week              = vacc_end,
    real_vaccine_duration_age  = real_vac_dur
  )
}


# both infection and disease block
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,      # disease-blocking efficacy
  VE_inf = 0,    # infection-blocking efficacy (default 0)
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  pmax0 <- function(x) pmax(0, x)
  
  ## 1) Initialization
  S  <- I  <- R_ <- V   <- matrix(0, A, T)
  vacc_delayed        <- matrix(0, A, T)
  raw_allocation_age  <- matrix(0, A, T)
  wasted_dose         <- matrix(0, A, T)
  unvaccinated        <- N
  coverage_frac       <- matrix(0, A, T)
  age_stratified_cases        <- matrix(0, A, T)
  age_stratified_cases_raw    <- matrix(0, A, T)
  true_symptomatic    <- matrix(0, A, T)
  phi_vec             <- numeric(T)
  FOI_times_S         <- numeric(T)
  coverage_target     <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  ## Compute total supply and weekly doses
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start     <- vacc_end <- rep(NA_integer_, A)
  
  ## Initial state at t = 1
  S[,1] <- N - I0_draw - R0 * N
  I[,1] <- I0_draw
  R_[,1] <- R0 * N
  V[,1] <- 0
  coverage_frac[,1] <- 0
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1] <- rho * I[,1] * (1 - VE_block * coverage_frac[,1])
  coverage_target[1] <- 0
  
  ## 2) Main loop
  for (t_ in 2:T) {
    ## 2A. Accumulate vaccine doses (for symptomatic blocking)
    if (t_ - 2 >= 1) {
      V[, t_] <- V[, t_ - 1] + vacc_delayed[, t_ - 2]
    } else {
      V[, t_] <- V[, t_ - 1]
    }
    coverage_frac[, t_] <- ifelse(N > 0, V[, t_] / N, 0)
    
    ## 2B. Force-of-Infection
    phi <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / sum(N))
    phi_vec[t_] <- phi
    FOI_times_S[t_] <- phi * sum(S[, t_ - 1])
    
    ## 2C. Vaccine allocation (same as original)
    if (t_ >= delay && target_pop > 0) {
      remaining <- weekly_dose_total
      for (a in target_idx) {
        base_alloc <- ceiling(weekly_dose_total * (N[a] / target_pop))
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        if (alloc > 0) {
          raw_allocation_age[a, t_] <- alloc
          S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
          prop_S <- ifelse(N[a] > 0, S_prev / N[a], 0)
          prop_I <- ifelse(N[a] > 0, I_prev / N[a], 0)
          vacc_to_S <- round(alloc * prop_S)
          vacc_to_I <- round(alloc * prop_I)
          vacc_to_R <- alloc - vacc_to_S - vacc_to_I
          wasted_IR <- vacc_to_I + vacc_to_R
          wasted_dose[a, t_]  <- wasted_IR
          vacc_delayed[a, t_] <- vacc_to_S
          total_used_age[a]   <- total_used_age[a] + alloc
          unused_age[a]       <- unused_age[a] - alloc
          unvaccinated[a]     <- unvaccinated[a] - alloc
          remaining           <- remaining - alloc
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    ## 2D. S-I-R-V update with combined blocking
    for (a in seq_len(A)) {
      S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
      # infection with infection-blocking
      new_inf <- phi * S_prev * (1 - VE_inf * coverage_frac[a, t_])
      recov   <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev - new_inf)
      I[a, t_]  <- pmax0(I_prev + new_inf - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # symptomatic blocking
      true_symptomatic[a, t_] <-
        0.5242478 * I[a, t_] * (1 - VE_block * coverage_frac[a, t_])
      
      # aging/death
      if (a == 1) {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_])
      } else {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_]  + r[a-1] * S[a-1, t_-1])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_]  + r[a-1] * I[a-1, t_-1])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_] + r[a-1] * R_[a-1, t_-1])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_]  + r[a-1] * V[a-1, t_-1])
      }
      
      # track cases
      age_stratified_cases_raw[a, t_] <- I[a, t_]
      age_stratified_cases[a, t_]     <- rho * I[a, t_]
    }
  }
  
  ## 3) vaccination duration
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) return
  list(
    S  = S, I  = I, R  = R_, V  = V,
    vacc_delayed               = vacc_delayed,
    raw_allocation_age         = raw_allocation_age,
    wasted_dose                = wasted_dose,
    total_supply               = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,
    unused_supply_age          = unused_age,
    coverage_target            = coverage_target,
    coverage_frac              = coverage_frac,
    phi                        = phi_vec,
    FOI_times_S                = FOI_times_S,
    age_stratified_cases       = age_stratified_cases,
    age_stratified_cases_raw   = age_stratified_cases_raw,
    true_symptomatic           = true_symptomatic,
    vacc_start_week            = vacc_start,
    vacc_end_week              = vacc_end,
    real_vaccine_duration_age  = real_vac_dur
  )
}

## scaling rho 
sirv_sim_coverageSwitch <- function(
    # 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma,
  # 2) Vaccine parameters
  delay,
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,      # disease-blocking efficacy
  VE_inf = 0,    # infection-blocking efficacy (default 0)
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  pmax0 <- function(x) pmax(0, x)
  total_pop   <- sum(N)
  R0_scalar   <- base_beta[1] / gamma      # classical basic R0 (assumes β[1] is baseline)
  R_eff_vec   <- numeric(T)  
  R0_vec      <- base_beta / gamma         # time-varying “potential R0(t)”
  R0_max      <- max(R0_vec, na.rm = TRUE)

  
  R_eff_vec <- numeric(T)     # container for effective Rt
  
  ## 1) Initialization
  S  <- I  <- R_ <- V   <- matrix(0, A, T)
  vacc_delayed        <- matrix(0, A, T)
  raw_allocation_age  <- matrix(0, A, T)
  wasted_dose         <- matrix(0, A, T)
  unvaccinated        <- N
  coverage_frac       <- matrix(0, A, T)
  age_stratified_cases        <- matrix(0, A, T)
  age_stratified_cases_raw    <- matrix(0, A, T)
  true_symptomatic    <- matrix(0, A, T)
  phi_vec             <- numeric(T)
  FOI_times_S         <- numeric(T)
  coverage_target     <- numeric(T)
  
  target_idx <- which(target_age == 1)
  target_pop <- sum(N[target_idx])
  
  ## Compute total supply and weekly doses
  total_supply      <- sum(N) * total_coverage
  weekly_dose_total <- total_supply * weekly_delivery_speed
  
  total_avail_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  unused_age     <- total_avail_age
  vacc_start     <- vacc_end <- rep(NA_integer_, A)

  ## Initial state at t = 1
  S[,1] <- N - I0_draw - R0 * N
  I[,1] <- I0_draw
  R_[,1] <- R0 * N
  V[,1] <- 0
  coverage_frac[,1] <- 0
  age_stratified_cases_raw[,1] <- I[,1]
  age_stratified_cases[,1]     <- rho * I[,1]
  true_symptomatic[,1] <- rho * I[,1] * (1 - VE_block * coverage_frac[,1])
  coverage_target[1] <- 0
  
  eff_S0         <- sum(S[, 1] * (1 - VE_inf * coverage_frac[, 1]))
  R_eff_vec[1]   <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  ## 2) Main loop
  for (t_ in 2:T) {
    ## 2A. Accumulate vaccine doses (for symptomatic blocking)
    if (t_ - 2 >= 1) {
      V[, t_] <- V[, t_ - 1] + vacc_delayed[, t_ - 2]
    } else {
      V[, t_] <- V[, t_ - 1]
    }
    coverage_frac[, t_] <- ifelse(N > 0, V[, t_] / N, 0)
    
    ## 2B. Force-of-Infection
    phi <- base_beta[t_ - 1] * (sum(I[, t_ - 1]) / sum(N))
    phi_vec[t_] <- phi
    FOI_times_S[t_] <- phi * sum(S[, t_ - 1])
    
    ## 2C. Vaccine allocation (same as original)
    if (t_ >= delay && target_pop > 0) {
      remaining <- weekly_dose_total
      for (a in target_idx) {
        base_alloc <- ceiling(weekly_dose_total * (N[a] / target_pop))
        max_people <- min(unvaccinated[a], remaining, unused_age[a])
        alloc      <- min(base_alloc, max_people)
        if (alloc > 0) {
          raw_allocation_age[a, t_] <- alloc
          S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
          prop_S <- ifelse(N[a] > 0, S_prev / N[a], 0)
          prop_I <- ifelse(N[a] > 0, I_prev / N[a], 0)
          vacc_to_S <- round(alloc * prop_S)
          vacc_to_I <- round(alloc * prop_I)
          vacc_to_R <- alloc - vacc_to_S - vacc_to_I
          wasted_IR <- vacc_to_I + vacc_to_R
          wasted_dose[a, t_]  <- wasted_IR
          vacc_delayed[a, t_] <- vacc_to_S
          total_used_age[a]   <- total_used_age[a] + alloc
          unused_age[a]       <- unused_age[a] - alloc
          unvaccinated[a]     <- unvaccinated[a] - alloc
          remaining           <- remaining - alloc
          if (is.na(vacc_start[a])) vacc_start[a] <- t_
          vacc_end[a] <- t_
        }
      }
      coverage_target[t_] <- sum(total_used_age[target_idx]) / target_pop
    } else {
      coverage_target[t_] <- coverage_target[t_ - 1]
    }
    
    ## 2D. S-I-R-V update with combined blocking
    for (a in seq_len(A)) {
      S_prev <- S[a, t_ - 1]; I_prev <- I[a, t_ - 1]; R_prev <- R_[a, t_ - 1]
      # infection with infection-blocking
      new_inf <- phi * S_prev * (1 - VE_inf * coverage_frac[a, t_])
      recov   <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev - new_inf)
      I[a, t_]  <- pmax0(I_prev + new_inf - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # symptomatic blocking
      true_symptomatic[a, t_] <-
        0.5242478 * I[a, t_] * (1 - VE_block * coverage_frac[a, t_]) * rho
      
      # aging/death
      if (a == 1) {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_])
      } else {
        S[a, t_]  <- pmax0(S[a, t_]  - r[a] * S[a, t_]  + r[a-1] * S[a-1, t_-1])
        I[a, t_]  <- pmax0(I[a, t_]  - r[a] * I[a, t_]  + r[a-1] * I[a-1, t_-1])
        R_[a, t_] <- pmax0(R_[a, t_] - r[a] * R_[a, t_] + r[a-1] * R_[a-1, t_-1])
        V[a, t_]  <- pmax0(V[a, t_]  - r[a] * V[a, t_]  + r[a-1] * V[a-1, t_-1])
      }
      
      # track cases
      age_stratified_cases_raw[a, t_] <- I[a, t_]
      age_stratified_cases[a, t_]     <- rho * I[a, t_]
    }
    eff_S_t        <- sum(S[, t_] * (1 - VE_inf * coverage_frac[, t_]))
    R_eff_vec[t_]  <- (base_beta[t_] / gamma) * (eff_S_t / total_pop)
  }
  
  ## 3) vaccination duration
  real_vac_dur <- ifelse(is.na(vacc_start) | is.na(vacc_end), 0,
                         vacc_end - vacc_start + 1)
  
  ## 4) return
  list(
    S  = S, I  = I, R  = R_, V  = V,
    vacc_delayed               = vacc_delayed,
    raw_allocation_age         = raw_allocation_age,
    wasted_dose                = wasted_dose,
    total_supply               = total_supply,
    total_available_supply_age = total_avail_age,
    total_supply_age           = total_used_age,
    unused_supply_age          = unused_age,
    coverage_target            = coverage_target,
    coverage_frac              = coverage_frac,
    phi                        = phi_vec,
    FOI_times_S                = FOI_times_S,
    R0_basic                   = R0_scalar,   # ← scalar classical R0
    R0_vec                     = R0_vec,
    R_eff_vec                  = R_eff_vec,
    age_stratified_cases       = age_stratified_cases,
    age_stratified_cases_raw   = age_stratified_cases_raw,
    true_symptomatic           = true_symptomatic,
    vacc_start_week            = vacc_start,
    vacc_end_week              = vacc_end,
    real_vaccine_duration_age  = real_vac_dur
  )
}


## seir model (coverage model)
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
  VE_inf  = 0,   # infection-blocking efficacy (all-or-nothing 성공 확률)
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  
  # ───────────────────────────────────────
  # 도움 함수 및 상수
  # ───────────────────────────────────────
  pmax0      <- function(x) pmax(0, x)
  total_pop  <- sum(N)
  R0_scalar  <- base_beta[1] / gamma
  R0_vec     <- base_beta / gamma
  R_eff_vec  <- numeric(T)
  
  # ───────────────────────────────────────
  # 상태 행렬 초기화
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
  # 백신 공급량 계산
  # ───────────────────────────────────────
  total_supply       <- target_pop * total_coverage
  weekly_dose_total  <- total_supply * weekly_delivery_speed
  
  total_avail_age <- unused_age <- rep(0, A)
  total_avail_age[target_idx] <- total_supply * (N[target_idx] / target_pop)
  total_used_age <- rep(0, A)
  
  vacc_start <- vacc_end <- rep(NA_integer_, A)
  
  # ───────────────────────────────────────
  # 초기 상태 (t = 1)
  # ───────────────────────────────────────
  S[, 1] <- N - I0_draw - R0 * N
  E[, 1] <- 0
  I[, 1] <- I0_draw
  R_[, 1] <- R0 * N
  V[, 1] <- 0
  
  coverage_frac[, 1] <- 0
  age_stratified_cases_raw[, 1] <- I[, 1]
  age_stratified_cases[,  1]    <- rho * I[, 1]
  true_symptomatic[, 1] <- rho * I[, 1]   # t=1에는 아직 백신효과 X
  coverage_target[1]   <- 0
  eff_S0               <- sum(S[, 1])     # ★ 완전 차단형이므로 그대로 S 합계
  R_eff_vec[1]         <- (base_beta[1] / gamma) * (eff_S0 / total_pop)
  
  # ───────────────────────────────────────
  # 메인 루프 (t = 2 … T)
  # ───────────────────────────────────────
  for (t_ in 2:T) {
    
    # ───── 1) 2주 전 접종분 → 면역 발현 ─────
    if (t_ - 2 >= 1) {
      effective_dose      <- vacc_delayed[, t_ - 2]              # 2주 전 S에 투입된 백신
      immunized           <- round(VE_inf * effective_dose)      # ★ 완전 면역 성공자
      V[, t_]             <- V[, t_ - 1] + immunized             # ★ V 누적
      V_covered[, t_]     <- V_covered[, t_ - 1] + effective_dose # ★ VE block track
      #S[, t_ - 1]         <- pmax0(S[, t_ - 1] - immunized)      # ★ S에서 제거
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
      S_prev_eff <- pmax0(S_prev - immunized_a) 
      
      new_e  <- phi * S_prev_eff                    
      new_i  <- sigma * E_prev
      recov  <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      # 증상자 
      true_symptomatic[a, t_] <- 0.5242478 *
        I[a, t_] *
        (1 - VE_block * coverage_frac[a, t_]) *
        rho
      
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
      
      age_stratified_cases_raw[a, t_] <- I[a, t_]
      age_stratified_cases[a,     t_] <- rho * I[a, t_]
    }
    
    # ───── 5) 유효 감수자 / R_eff ─────
    eff_S              <- sum(S[, t_])          
    R_eff_vec[t_]      <- (base_beta[t_] / gamma) * (eff_S / total_pop)
  }
  
  # ───────────────────────────────────────
  # 출력
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

## vimkunya
sirv_sim_coverageSwitch <- function(
    ## 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma, sigma,
  ## 2) Vaccine parameters
  delay                = 2,   # 공급(roll‑out) 지연 주수
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,                    # disease‑blocking efficacy
  VE_inf          = 0.978,     # infection‑blocking efficacy (immunogenic 성공 확률)
  VE_inf_postwane = 0.855,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
) {
  
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
      immunized_success <- round(VE_inf * effective_dose)
      
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
      
      if (VE_inf > 0) {
        #VE_inf_postwane <- 0.855
        current_reinstate <- round((1 - VE_inf_postwane / VE_inf) * waning_group)
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
    
    ## 4) SEIRV 업데이트
    for (a in seq_len(A)) {
      S_prev <- S[a, t_-1]; E_prev <- E[a, t_-1]
      I_prev <- I[a, t_-1]; R_prev <- R_[a, t_-1]
      
      immunized_a <- immunized_success[a]
      S_prev_eff  <- pmax0(S_prev - immunized_a + reinstate[a])
      
      new_e <- phi * S_prev_eff
      new_i <- sigma * E_prev
      recov <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      true_symptomatic[a, t_] <- 0.5242478 *
        I[a, t_] *
        (1 - VE_block * coverage_frac[a, t_]) *
        rho
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

## time-varying waning dyanmics (linear decrease)
sirv_sim_coverageSwitch_vimkun <- function(
    ## 1) Model dimensions & epidemic parameters
  T, A, N, r,
  base_beta, I0_draw, R0,
  rho, gamma, sigma,
  ## 2) Vaccine parameters
  delay                = 2,   # 공급(roll‑out) 지연 주수
  target_age,
  total_coverage,
  weekly_delivery_speed,
  VE_block,                    # disease‑blocking efficacy
  VE_inf          = 0.978,     # infection‑blocking efficacy (immunogenic 성공 확률)
  VE_inf_postwane = 0.855,
  coverage_threshold  = 1,
  use_coverage_switch = FALSE
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
    
    ## 4) SEIRV 업데이트
    for (a in seq_len(A)) {
      S_prev <- S[a, t_-1]; E_prev <- E[a, t_-1]
      I_prev <- I[a, t_-1]; R_prev <- R_[a, t_-1]
      
      immunized_a <- immunized_success[a]
      S_prev_eff  <- pmax0(S_prev - immunized_a + reinstate[a])
      
      new_e <- phi * S_prev_eff
      new_i <- sigma * E_prev
      recov <- gamma * I_prev
      
      S[a, t_]  <- pmax0(S_prev_eff - new_e)
      E[a, t_]  <- pmax0(E_prev + new_e - new_i)
      I[a, t_]  <- pmax0(I_prev + new_i - recov)
      R_[a, t_] <- pmax0(R_prev + recov)
      
      true_symptomatic[a, t_] <- 0.5242478 *
        I[a, t_] *
        (1 - VE_block[t_] * coverage_frac[a, t_]) *
        rho
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

