### one-way sensitivity analysis ----------------------------------------------

# ----- 1. Define baseline values and sensitivity ranges -----
make_params <- function(region, bra_foi_state_summ, posterior, lhs_sample_young, hosp, fatal, nh_fatal) {
  # Create baseline parameters based on the specified region
  baseline_params <- list(
    foi          = bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region],
    hosp         = hosp,
    fatal        = fatal,
    nh_fatal     = nh_fatal,
    VE_block     = 0.989,
    VE_inf       = 0.989,
    weekly_delivery_speed = 0.1,
    delay = 2,
    total_coverage = 0.5,
    time_until_immunity = 2
  )
  
  # Define sensitivity ranges for each parameter. (Ensure these ranges are wide enough to affect outcomes.)
  sensitivity_ranges <- list(
    foi          = c(bra_foi_state_summ$foi_lo[bra_foi_state_summ$NAME_1 == region],
                     bra_foi_state_summ$foi_hi[bra_foi_state_summ$NAME_1 == region]),
    VE_block    = c(0.967, 0.998),
    VE_inf      = c(0.967, 0.998),
    weekly_delivery_speed = c(0.09, 0.11),
    delay = c(1,3),
    total_coverage = c(0.4, 0.6),
    time_until_immunity = c(1, 3)
  )
  
  # Return both sets of parameters in a list
  return(list(baseline_params = baseline_params, sensitivity_ranges = sensitivity_ranges))
}


# --- 3. modify sim functions ---
simulate_pre_ui_age_owsa <- function(
    posterior,
    bra_foi_state_summ,
    age_groups,
    N,
    region,
    observed,
    A = 20, 
    delay = 53,
    VE_block = 0,
    VE_inf   = 0,
    target_age = rep(0, A),
    coverage_threshold = 0,
    total_coverage = 0,
    total_supply = 0,
    weekly_delivery_speed = 0,
    time_until_immunity = 0,
    
    # New argument: FOI override (default = region's avg_foi)
    foi_value = bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region]
) {
  # Number of weeks = number of rows in observed
  T <- nrow(observed)
  
  # Number of posterior draws
  n_draws <- length(posterior$gamma)
  
  # Prepare list to store results from each draw
  sim_results_list <- vector("list", n_draws)
  
  # Loop over posterior draws
  for (i in seq_len(n_draws)) {
    base_beta_draw <- posterior$base_beta[i, ]
    I0_draw        <- posterior$I0[i, ]
    gamma_draw     <- posterior$gamma[i]
    rho_draw       <- posterior$rho[i]
    sigma_draw     <- posterior$sigma[i]
    
    # Use the user-specified FOI value
    R0_vec <- 1 - exp(-foi_value * age_groups)
    
    # Run your custom SIRV simulation function
   
    sim_out <- sirv_sim_coverageSwitch(
      T = T,
      A = length(age_gr_levels),
      N = N,
      r = rep(0, length(age_gr_levels)),
      base_beta = base_beta_draw,
      I0_draw = I0_draw,
      R0 = R0_vec,
      rho = rho_draw,
      gamma = gamma_draw,
      sigma = sigma_draw,
      delay = delay,
      VE_block = VE_block,
      VE_inf = VE_inf,
      #coverage_threshold = 1,
      target_age = target_age,
      #total_coverage = 2e6 * 0.2 / sum(N),
      #total_coverage = total_coverage_frac,
      total_coverage = total_coverage,
      #total_coverage = 0.01,
      weekly_delivery_speed = weekly_delivery_speed,
      time_until_immunity = time_until_immunity
    )
    # Each draw returns an age-stratified cases matrix [A x T]
    sim_results_list[[i]] <- sim_out$true_symptomatic
  }
  
  # Combine into 3D array: [A, T, n_draws]
  age_array <- array(unlist(sim_results_list), dim = c(A, T, n_draws))
  
  # Prepare empty matrices to store median and quantiles by age x week
  median_by_age <- matrix(0, nrow = A, ncol = T)
  low95_by_age  <- matrix(0, nrow = A, ncol = T)
  hi95_by_age   <- matrix(0, nrow = A, ncol = T)
  
  # For each age group and week, compute median, low 95%, and high 95% across draws
  for (a in seq_len(A)) {
    for (t in seq_len(T)) {
      draws <- age_array[a, t, ]
      median_by_age[a, t] <- median(draws)
      low95_by_age[a, t]  <- quantile(draws, 0.025)
      hi95_by_age[a, t]   <- quantile(draws, 0.975)
    }
  }
  
  # Weekly totals across age groups for each draw => matrix [T x n_draws]
  weekly_totals <- apply(age_array, c(2, 3), sum)
  
  weekly_cases_median <- apply(weekly_totals, 1, median)
  weekly_cases_low95  <- apply(weekly_totals, 1, quantile, probs = 0.025)
  weekly_cases_hi95   <- apply(weekly_totals, 1, quantile, probs = 0.975)
  
  # Create a summary data frame
  df <- data.frame(
    week   = seq_len(T),
    median = weekly_cases_median,
    low95  = weekly_cases_low95,
    hi95   = weekly_cases_hi95
  )
  
  return(list(
    sim_results_list    = sim_results_list,
    age_array           = age_array,
    median_by_age       = median_by_age,
    low95_by_age        = low95_by_age,
    hi95_by_age         = hi95_by_age,
    weekly_totals       = weekly_totals,
    weekly_cases_median = weekly_cases_median,
    weekly_cases_low95  = weekly_cases_low95,
    weekly_cases_hi95   = weekly_cases_hi95,
    df                  = df
  ))
}


summarise_presim_ui_owsa <- function(
    sim_result,        # output from simulate_pre_ui_age (or simulate_pre_ui_age_owsa)
    observed,
    age_gr_levels = c(
      "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
      "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
      "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
    ),
    lhs_sample_young,  # data for younger ages (if needed for reference)
    lhs_old,           # data for older ages (if needed for reference)
    le_sample,         # life-expectancy sample (if needed)
    hosp,              # hospitalization rate (scalar or vector)
    fatal,             # fatality rate for hospitalized
    nh_fatal,          # fatality rate for non-hospitalized
    region
) {
  # Number of weeks
  T <- nrow(observed)
  
  # 1) Prepare data frames from sim_result
  default_age_vector <-  c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  # Age group repeated by T weeks (for merging later)
  age_gr <- rep(default_age_vector, each = T)
  
  median_df <- as.data.frame.table(sim_result$median_by_age, responseName = "Median")
  colnames(median_df) <- c("AgeGroup", "Week", "Median")
  median_df <- median_df %>%
    dplyr::mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week = as.numeric(Week)
    )
  
  low95_df <- as.data.frame.table(sim_result$low95_by_age, responseName = "low95")
  hi95_df  <- as.data.frame.table(sim_result$hi95_by_age, responseName = "hi95")
  
  summary_cases_pre <- median_df %>%
    dplyr::mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]
    )
  summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)
  
  # Weekly summary data frame (overall time series)
  weekly_df <- data.frame(
    Week = seq_along(sim_result$weekly_cases_median),
    weekly_median = sim_result$weekly_cases_median,
    weekly_low95  = sim_result$weekly_cases_low95,
    weekly_hi95   = sim_result$weekly_cases_hi95
  )
  
  # 2) Summarize across ages for each week -> summary_cases_pre_all
  summary_cases_pre_all <- summary_cases_pre %>%
    dplyr::group_by(Week) %>%
    dplyr::summarise(Median = sum(Median, na.rm = TRUE)) %>%
    dplyr::mutate(Scenario = "Pre-vaccination") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      lo95  = weekly_low95,
      hi95  = weekly_hi95
    )
  
  # 3) Add columns for hospitalization/fatalities
  summary_cases_pre <- summary_cases_pre %>%
    dplyr::mutate(
      hosp_rate        = rep(hosp, T),
      hospitalised     = Median * hosp_rate,
      hospitalised_lo  = low95 * hosp_rate,
      hospitalised_hi  = hi95  * hosp_rate,
      
      non_hospitalised    = Median - hospitalised,
      non_hospitalised_lo = low95 - hospitalised_lo,
      non_hospitalised_hi = hi95  - hospitalised_hi,
      
      fatality    = rep(fatal, T),
      nh_fatality = rep(nh_fatal, T),
      
      fatal    = hospitalised     * fatality   + non_hospitalised     * nh_fatality,
      fatal_lo = hospitalised_lo  * fatality   + non_hospitalised_lo  * nh_fatality,
      fatal_hi = hospitalised_hi  * fatality   + non_hospitalised_hi  * nh_fatality
    ) %>%
    dplyr::arrange(Week, AgeGroup) %>%
    dplyr::group_by(Week, AgeGroup) %>%
    dplyr::mutate(
      cum_fatal = cumsum(fatal),
      cum_hosp  = cumsum(hospitalised)
    ) %>%
    dplyr::ungroup()
  
  # 4) Age-group summary (aggregates all weeks for each age)
  summary_cases_pre_age <- summary_cases_pre %>% 
    dplyr::group_by(AgeGroup) %>%
    dplyr::summarise(
      Median         = sum(Median, na.rm = TRUE),
      hospitalised   = sum(hospitalised, na.rm = TRUE),
      hospitalised_lo= sum(hospitalised_lo, na.rm = TRUE),
      hospitalised_hi= sum(hospitalised_hi, na.rm = TRUE),
      fatal          = sum(fatal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Scenario = "Pre-vaccination")
  
  summary_cases_pre_age$age_gr <- age_gr[1:20]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # 5) Summarize by Week for the entire population (detailed)
  summary_cases_pre_all <- summary_cases_pre %>%
    dplyr::group_by(Week) %>%
    dplyr::summarise(
      Median       = sum(Median, na.rm = TRUE),
      hospitalised = sum(hospitalised, na.rm = TRUE),
      fatal        = sum(fatal, na.rm = TRUE),
      fatal_lo     = sum(fatal_lo, na.rm = TRUE),
      fatal_hi     = sum(fatal_hi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal  = cumsum(fatal),
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      cum_hosp   = cumsum(hospitalised)
    ) %>%
    dplyr::left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      pre_weekly_low95  = weekly_low95,
      pre_weekly_hi95   = weekly_hi95
    )
  
  # [Removed the DALY calculations and arguments that overridden dw_hosp, dur_acute, etc.]
  
  # Add region info
  summary_cases_pre$region      <- region
  summary_cases_pre_all$region  <- region
  summary_cases_pre_age$region  <- region
  
  # Return a list of data frames
  return(list(
    summary_cases_pre     = summary_cases_pre,
    summary_cases_pre_all = summary_cases_pre_all,
    summary_cases_pre_age = summary_cases_pre_age
  ))
}


run_simulation_scenarios_ui_owsa <- function(
    target_age_list,
    observed,
    N,
    bra_foi_state_summ,
    age_groups,
    region_name,
    
    # Hospitalization & fatality rates
    hosp       = NULL,
    fatal      = NULL,
    nh_fatal   = NULL,
    
    # Young/old samples (for default fallback)
    lhs_sample_young,
    lhs_old,
    le_sample,
    age_gr_levels,
    prevacc_ui = NULL,
    
    # Posterior draws
    posterior,
    total_coverage = NULL,
    delay = 2,
    # Remaining possible overrides
    foi                    = NULL,
    VE_block              = NULL,
    VE_inf                = NULL,
    weekly_delivery_speed = NULL,
    time_until_immunity = NULL
    
) {
  # Number of scenarios = length of target_age_list
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  
  # Number of posterior draws
  n_draws <- length(posterior$gamma)
  
  # Number of weeks = rows in observed data
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <-c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  # If user doesn't specify foi, use the region-specific default from bra_foi_state_summ
  if (is.null(foi)) {
    foi <- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name]
  }
  if (is.null(hosp))     hosp     <- hosp
  if (is.null(fatal))    fatal    <- fatal
  if (is.null(nh_fatal)) nh_fatal <- nh_fatal
  
  # Vaccine overrides
  if (is.null(VE_block))              VE_block <- 0.989
  if (is.null(weekly_delivery_speed)) weekly_delivery_speed <- 0.1
  if (is.null(VE_inf)) VE_inf <- 0.989
  # Loop over each scenario
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    draw_results <- vector("list", n_draws)
    
    for (d in seq_len(n_draws)) {
      # Extract posterior draws
      base_beta_draw <- posterior$base_beta[d, ]
      I0_draw        <- posterior$I0[d, ]
      gamma_draw     <- posterior$gamma[d]
      rho_draw       <- posterior$rho[d]
      sigma_draw     <- posterior$sigma[d]
      
      # Compute R0 using the (possibly overridden) foi
      R0_vec <- 1 - exp(-foi * age_groups)
      
      sim_out <- sirv_sim_coverageSwitch(
        T = T,
        A = length(age_gr_levels),
        N = N,
        r = rep(0, length(age_gr_levels)),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        R0 = R0_vec,
        rho = rho_draw,
        gamma = gamma_draw,
        sigma = sigma_draw,
        delay = delay,
        VE_block = VE_block,
        VE_inf = VE_inf,
        #VE_inf_postwane   = this_postwane, 
        coverage_threshold = 1,
        target_age = target,
        #total_coverage = 2e6 * 0.2 / sum(N),
        #total_coverage = total_coverage_frac,
        total_coverage = total_coverage,
        #total_coverage = 0.01,
        weekly_delivery_speed = weekly_delivery_speed,
        time_until_immunity = time_until_immunity
      )
      draw_results[[d]] <- sim_out$true_symptomatic  # A x T matrix
    }
    
    # Combine draws into a 3D array [18, T, n_draws]
    age_array <- array(unlist(draw_results), dim = c(20, T, n_draws))
    
    # Compute quantiles across draws
    median_by_age <- apply(age_array, c(1, 2), median)
    low95_by_age  <- apply(age_array, c(1, 2), quantile, probs = 0.025)
    hi95_by_age   <- apply(age_array, c(1, 2), quantile, probs = 0.975)
    
    # Weekly totals
    weekly_totals <- apply(age_array, c(2, 3), sum)
    weekly_cases_median <- apply(weekly_totals, 1, median)
    weekly_cases_low95  <- apply(weekly_totals, 1, quantile, probs = 0.025)
    weekly_cases_hi95   <- apply(weekly_totals, 1, quantile, probs = 0.975)
    
    # Create simulation data frame from median_by_age
    sim_df <- as.data.frame.table(median_by_age, responseName = "Cases")
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      dplyr::mutate(
        Scenario = s,
        AgeGroup = as.numeric(AgeGroup),
        Week     = as.numeric(Week)
      )
    
    # Attach UI for cases
    low95_df <- as.data.frame.table(low95_by_age, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age, responseName = "post_hi95")
    
    sim_df <- sim_df %>%
      dplyr::mutate(
        post_low95 = low95_df$post_low95,
        post_hi95  = hi95_df$post_hi95
      )
    
    # Weekly summary data frame
    weekly_df <- data.frame(
      Week = seq_len(T),
      weekly_median = weekly_cases_median,
      weekly_low95  = weekly_cases_low95,
      weekly_hi95   = weekly_cases_hi95
    )
    
    # Attach age group labels
    sim_df <- sim_df %>%
      dplyr::mutate(
        age_gr = rep(default_age_vector, T)
      )
    sim_df$age_gr <- factor(sim_df$age_gr, levels = age_gr_levels)
    
    # Calculate hospitalisations, fatalities, etc.
    sim_df <- sim_df %>%
      dplyr::mutate(
        hosp_rate       = rep(hosp, T),
        hospitalised    = Cases * hosp_rate,
        hospitalised_lo = post_low95 * hosp_rate,
        hospitalised_hi = post_hi95 * hosp_rate,
        
        non_hospitalised    = Cases - hospitalised,
        non_hospitalised_lo = post_low95 - hospitalised_hi,
        non_hospitalised_hi = post_hi95 - hospitalised_lo,
        
        fatality    = rep(fatal, T),
        nh_fatality = rep(nh_fatal, T),
        
        fatal    = hospitalised * fatality + non_hospitalised * nh_fatality,
        fatal_lo = hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality,
        fatal_hi = hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality
      ) %>%
      dplyr::arrange(Week, AgeGroup) %>%
      dplyr::group_by(Week, AgeGroup) %>%
      dplyr::mutate(
        cum_fatal = cumsum(fatal),
        cum_hosp  = cumsum(hospitalised)
      ) %>%
      dplyr::ungroup()
    
    sim_df <- sim_df %>%
      mutate(
        age_numeric = case_when(
          age_gr == "<1" ~ 0,
          age_gr == "1-4" ~ 2.5,
          age_gr == "5–9" ~ 7,
          age_gr == "10-11" ~ 10.5,
          age_gr == "12-17" ~ 14.5,
          age_gr == "18–19" ~ 18.5,
          age_gr == "20–24" ~ 22,
          age_gr == "25–29" ~ 27,
          age_gr == "30–34" ~ 32,
          age_gr == "35–39" ~ 37,
          age_gr == "40–44" ~ 42,
          age_gr == "45–49" ~ 47,
          age_gr == "50–54" ~ 52,
          age_gr == "55–59" ~ 57,
          age_gr == "60–64" ~ 62,
          age_gr == "65–69" ~ 67,
          age_gr == "70–74" ~ 72,
          age_gr == "75–79" ~ 77,
          age_gr == "80–84" ~ 82,
          age_gr == "85+" ~ 87.5,
          TRUE ~ NA_real_
        )
      )
    
    # Add YLD parameters (age independent)
    sim_df <- sim_df %>%
      mutate(
        dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Age-dependent parameters:
      mutate(
        subac_prop = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)
        ),
        chr_6m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr6m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr6m, 0.5)
        ),
        chr_12m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr12m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr12m, 0.5)
        ),
        chr_30m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr30m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr30m, 0.5)
        ),
        chr_prop = chr_6m + chr_12m + chr_30m
      ) %>%
      # Calculate YLD with uncertainty bounds
      mutate(
        yld_acute = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_acute_lo = (hospitalised_lo * dw_hosp * dur_acute) +
          ((post_low95 - hospitalised_hi) * dw_nonhosp * dur_acute),
        yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
          ((post_hi95 - hospitalised_lo) * dw_nonhosp * dur_acute),
        
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
          ((post_low95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
          ((post_hi95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
        
        yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
          ((post_low95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
          ((post_hi95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
        
        yld_total = yld_acute + yld_subacute + yld_chronic,
        yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
        yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
      ) %>%
      # Determine remaining life expectancy by age and calculate YLL and DALY
      mutate(
        le_left = case_when(
          age_numeric <= 1 ~ quantile(le_sample$le_1, 0.5),
          age_numeric > 1 & age_numeric < 10 ~ quantile(le_sample$le_1, 0.5),
          age_numeric >= 10 & age_numeric < 20 ~ quantile(le_sample$le_2, 0.5),
          age_numeric >= 20 & age_numeric < 30 ~ quantile(le_sample$le_2, 0.5),
          age_numeric >= 30 & age_numeric < 40 ~ quantile(le_sample$le_3, 0.5),
          age_numeric >= 40 & age_numeric < 50 ~ quantile(le_sample$le_4, 0.5),
          age_numeric >= 50 & age_numeric < 60 ~ quantile(le_sample$le_5, 0.5),
          age_numeric >= 60 & age_numeric < 70 ~ quantile(le_sample$le_6, 0.5),
          age_numeric >= 70 & age_numeric < 80 ~ quantile(le_sample$le_7, 0.5),
          age_numeric >= 80 ~ quantile(le_sample$le_8, 0.5),
          TRUE ~ NA_real_
        )
      ) %>%
      mutate(
        yll = fatal * le_left,
        yll_lo = fatal_lo * le_left,
        yll_hi = fatal_hi * le_left,
        daly_tot = yld_total + yll,
        daly_tot_lo = yld_total_lo + yll_lo,
        daly_tot_hi = yld_total_hi + yll_hi,
        cum_daly = cumsum(daly_tot)
      )
    
    # If pre-vacc UI is provided, attach it
    if (!is.null(prevacc_ui)) {
      sim_df <- dplyr::left_join(
        sim_df,
        dplyr::select(prevacc_ui, AgeGroup, Week, 
                      pre_Median = Median, pre_low95 = low95, pre_hi95 = hi95),
        by = c("AgeGroup", "Week")
      )
    }
    
    # Collect results for this scenario
    scenario_result[[s]] <- list(
      sim_result = list(
        age_array            = age_array,
        weekly_cases_median  = weekly_cases_median,
        weekly_cases_low95   = weekly_cases_low95,
        weekly_cases_hi95    = weekly_cases_hi95
      ),
      sim_out   = NULL,  # or store the entire sim_out if desired
      sim_df    = sim_df,
      weekly_df = weekly_df
    )
  }
  
  return(scenario_result)
}

postsim_all_ui_owsa <- function(
    scenario_result,   # list of scenario outputs: each with $sim_df and $weekly_df
    observed,
    age_gr_levels,
    pre_summary_cases_age,   # aggregated pre-vacc by age
    pre_summary_cases,       # pre-vacc data at age-week level (if needed)
    pre_summary_cases_all,   # aggregated pre-vacc weekly UI data (with Week, Median, low95, hi95)
    region,
    
    # Keep these overrides if you like
    foi                    = NULL,
    hosp                  = NULL,
    fatal                 = NULL,
    nh_fatal              = NULL,
    VE_block              = NULL,
    VE_inf                = NULL, 
    weekly_delivery_speed = NULL,
    delay = NULL,
    total_coverage = NULL
    #base_beta,
    #I0,
    #gamma,
    #rho
) {
  # Number of weeks
  T <- nrow(observed)
  
  # Default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  # Create an age group vector for age-week summary
  age_gr <- rep(default_age_vector, T)
  
  # Number of scenarios
  n_scenarios <- length(scenario_result)
  
  # Extract the age-week simulation data for each scenario
  scenario_df <- lapply(scenario_result, function(x) x$sim_df)
  
  # Build a summary list comparing pre-vacc vs post-vacc for each scenario
  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    # Attach aggregated pre-vaccination values (assumed same resolution)
    df$pre_vacc  <- pre_summary_cases$Median  
    df$pre_fatal <- pre_summary_cases$fatal
    
    # Just do differences in Cases & Fatal
    df <- df %>%
      dplyr::mutate(
        diff   = pre_vacc - Cases,
        impact = diff / pre_vacc * 100
      )
    return(df)
  })
  
  # Combine them into one DataFrame
  summary_list_df <- do.call(rbind, summary_list)
  summary_list_df$age_gr <- rep(age_gr, length.out = nrow(summary_list_df))
  summary_list_df$age_gr <- factor(summary_list_df$age_gr, levels = age_gr_levels)
  
  # Optionally merge aggregated pre-vacc UI (weekly) from pre_summary_cases_all
  if (all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      dplyr::distinct(Week, .keep_all = TRUE)
    summary_list_df <- dplyr::left_join(summary_list_df, pre_weekly_df, by = "Week")
  }
  
  # Summarize by AgeGroup -> sum over weeks
  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%  
      summarise(
        Median           = sum(Cases, na.rm = TRUE),
        low95            = sum(post_low95, na.rm = TRUE),
        hi95             = sum(post_hi95, na.rm = TRUE),
        fatal            = sum(fatal, na.rm = TRUE),
        fatal_lo         = sum(fatal_lo, na.rm = TRUE),
        fatal_hi         = sum(fatal_hi, na.rm = TRUE)
      ) %>% 
      mutate(
        pre_vacc              = pre_summary_cases_age$Median,
        pre_vacc_low95        = pre_summary_cases_age$low95,
        pre_vacc_hi95         = pre_summary_cases_age$hi95,

        pre_fatal             = pre_summary_cases_age$fatal,
        pre_fatal_low95       = pre_summary_cases_age$fatal_lo,
        pre_fatal_hi          = pre_summary_cases_age$fatal_hi,
        
        # Differences for Cases:
        diff       = pre_vacc - Median,
        diff_low   = pre_vacc_low95 - low95,
        diff_hi    = pre_vacc_hi95 - hi95,
        impact     = diff / pre_vacc * 100,
        impact_low = diff_low / pre_vacc_low95 * 100,
        impact_hi  = diff_hi / pre_vacc_hi95 * 100,
        
        # Differences for Fatal outcomes:
        diff_fatal       = pre_fatal - fatal,
        diff_fatal_low   = pre_fatal_low95 - fatal_lo,
        diff_fatal_hi    = pre_fatal_hi - fatal_hi,
        impact_fatal     = diff_fatal / pre_fatal * 100,
        impact_fatal_low = diff_fatal_low / pre_fatal_low95 * 100,
        impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,
      )
  })
  
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_cases = sum(Cases, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        Scenario   = first(Scenario),
        .groups = "drop"
      )
  })
  
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_cases = sum(Cases, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        Scenario   = first(Scenario),
        .groups = "drop"
      )
  })
  
  summary_week_df <- do.call(rbind, summary_week)
  
  # Incorporate the weekly UI from each scenario.
  for (i in seq_len(n_scenarios)) {
    scenario_ui <- scenario_result[[i]]$weekly_df %>%
      dplyr::rename(
        post_weekly_median = weekly_median,
        post_weekly_low95  = weekly_low95,
        post_weekly_hi95   = weekly_hi95
      ) %>%
      distinct(Week, .keep_all = TRUE)
    summary_week[[i]] <- left_join(summary_week[[i]], scenario_ui, by = "Week", relationship = "many-to-many")
  }
  
  summary_week_df <- do.call(rbind, summary_week)
  
  # Merge aggregated pre-vacc UI into the weekly summary.
  if(all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_week_df <- left_join(summary_week_df, pre_weekly_df, by = "Week", relationship = "many-to-many")
  }
  
  summary_week_df$region <- region
  
  summary_week_df <- summary_week_df %>%
    left_join(rho_df, by = "region") %>%
    mutate(
      across(
        c(post_cases, pre_cases, diff, post_fatal, pre_fatal, diff_fatal,
          post_weekly_median, post_weekly_low95,
          post_weekly_hi95, lo95, hi95),
        ~ .x / rho_p50
      ) 
    )
  
  # Compute global weekly impact
  global_impact_week <- summary_week_df %>%
    dplyr::group_by(Scenario) %>%
    dplyr::summarise(
      total_post_cases = sum(post_cases, na.rm = TRUE),
      total_pre_case   = sum(pre_cases, na.rm = TRUE),
      total_post_fatal = sum(post_fatal, na.rm = TRUE),
      total_pre_fatal  = sum(pre_fatal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      diff         = total_pre_case - total_post_cases,
      impact       = diff / total_pre_case * 100,
      diff_fatal   = total_pre_fatal - total_post_fatal,
      impact_fatal = diff_fatal / total_pre_fatal * 100
    )
  
  # Simple textual annotation
  annotation_text <- paste0(
    global_impact_week$Scenario, ": ", 
    round(global_impact_week$impact, 1), "%",
    collapse = "\n"
  )
  
  global_impact_week$region <- region
  
  # Example references to vaccination weeks (if the data exist in scenario_result)
  # Adjust or remove as needed
  vacc_start_week_s1 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s1   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  
  vacc_start_week_s2 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s2   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  
  vacc_start_week_s3 <- scenario_result[[3]]$sim_out$vacc_start_week[13]
  vacc_end_week_s3   <- scenario_result[[3]]$sim_out$vacc_end_week[13]
  
  return(list(
    scenario_result     = scenario_result,
    scenario_data       = lapply(scenario_result, function(x) x$sim_result),
    scenario_df         = scenario_df,
    summary_list        = summary_list,
    summary_list_df     = summary_list_df,
    final_summ          = final_summ,
    summary_week        = summary_week,
    summary_week_df     = summary_week_df,
    global_impact_week  = global_impact_week,
    annotation_text     = annotation_text,
    vacc_weeks = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
    ),
    
    # Keep only these param overrides in metadata (optional)
    param_metadata = list(
      foi                    = foi,
      hosp                  = hosp,
      fatal                 = fatal,
      nh_fatal              = nh_fatal,
      VE_block              = VE_block,
      VE_inf                = VE_inf,
      weekly_delivery_speed = weekly_delivery_speed,
      base_beta             = base_beta,
      gamma                 = gamma,
      I0                    = I0,
      rho                   = rho,
      delay = delay
    )
  ))
}


# --- 2. Create a wrapper function to run the simulation with modified parameters ---


run_simulation_with_params <- function(
    param_overrides,
    baseline_params,
    posterior,
    bra_foi_state_summ,
    age_groups,
    N,
    region,
    observed,
    age_gr_levels,
    lhs_sample_young,
    lhs_old,
    le_sample,
    
    # New arguments for the scenario steps:
    target_age_list = target_age_list,        # list of post-vacc target age vectors
    total_coverage = 0.5,
    weekly_delivery_speed = 0.1,
    VE_block = 0.989,
    VE_inf = 0.989,
    delay = 2
    
    #base_beta,
    #I0,
    #gamma,
    #rho
) {
  # ---- 0. Combine baseline + overrides ----
  params <- modifyList(baseline_params, param_overrides)
  
  total_coverage <- if (!is.null(params$total_coverage)) params$total_coverage else total_coverage
  delay          <- if (!is.null(params$delay)) params$delay else delay
  VE_block           <- if (!is.null(params$VE_block))           params$VE_block else VE_block
  VE_inf             <- if (!is.null(params$VE_inf))             params$VE_inf else VE_inf
  weekly_delivery_speed <- if (!is.null(params$weekly_delivery_speed)) params$weekly_delivery_speed else weekly_delivery_speed
  time_until_immunity   <- if (!is.null(params$time_until_immunity))   params$time_until_immunity else time_until_immunity
  
  # ---- 1. PRE-VACC SIMULATION: mid 
  # This uses simulate_pre_ui_age_owsa() with FOI override, etc.
  sim_res_pre <- simulate_pre_ui_age_owsa(
    posterior             = posterior,
    bra_foi_state_summ    = bra_foi_state_summ,
    age_groups            = age_groups,
    N                     = N,
    region                = region,
    observed              = observed,
    delay                 = 53,
    VE_block              = 0,
    VE_inf                = 0,
    target_age            = rep(0, 20),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    time_until_immunity   = 0, 
    
    # FOI override and fixed parameters from params:
    foi_value = baseline_params$foi
    #base_beta = baseline_params$base_beta,
    #I0        = baseline_params$I0,
    #gamma     = baseline_params$gamma,
    #rho       = baseline_params$rho
  )
  
  sim_res_pre_owsa <- simulate_pre_ui_age_owsa(
    posterior             = posterior,
    bra_foi_state_summ    = bra_foi_state_summ,
    age_groups            = age_groups,
    N                     = N,
    region                = region,
    observed              = observed,
    delay                 = 53,
    VE_block              = 0,
    VE_inf                = 0,
    target_age            = rep(0, 20),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    time_until_immunity   = 0,
    # FOI override and fixed parameters from params:
    foi_value = params$foi
    #base_beta = baseline_params$base_beta,
    #I0        = baseline_params$I0,
    #gamma     = baseline_params$gamma,
    #rho       = baseline_params$rho
  )
  
  # ---- 2. SUMMARISE PRE-VACC ----
  summ_res_pre <- summarise_presim_ui_owsa(
    sim_result       = sim_res_pre,
    observed         = observed,
    age_gr_levels    = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = params$hosp,
    fatal            = params$fatal,
    nh_fatal         = params$nh_fatal,
    region           = region
  )
  
  summ_res_pre_owsa <- summarise_presim_ui_owsa(
    sim_result       = sim_res_pre_owsa,
    observed         = observed,
    age_gr_levels    = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = params$hosp,
    fatal            = params$fatal,
    nh_fatal         = params$nh_fatal,
    region           = region
  )
  
  # mid results 
  pre_summary_cases      <- summ_res_pre$summary_cases_pre      # age x week detail
  pre_summary_cases_all  <- summ_res_pre$summary_cases_pre_all  # weekly aggregated
  pre_summary_cases_age  <- summ_res_pre$summary_cases_pre_age  # age-aggregated
  
  # If we want a single "outcome metric" from pre-vacc, e.g. total cases:
  pre_tot_cases <- sum(pre_summary_cases_all$Median, na.rm = TRUE)
  
  # scenario results 
  # mid results 
  pre_summary_cases_owsa      <- summ_res_pre_owsa$summary_cases_pre      # age x week detail
  pre_summary_cases_all_owsa  <- summ_res_pre_owsa$summary_cases_pre_all  # weekly aggregated
  pre_summary_cases_age_owsa  <- summ_res_pre_owsa$summary_cases_pre_age  # age-aggregated
  
  pre_agecat_owsa <- pre_summary_cases_age_owsa %>%
    mutate(
      age_cat = case_when(
        AgeGroup %in% 2:4   ~ "2-4",
        AgeGroup %in% 5     ~ "5",
        AgeGroup %in% 6:15 ~ "6-15",
        AgeGroup %in% 16:20 ~ "16-20",
        TRUE                ~ NA_character_
      )
    )
  
  pre_tot_by_agecat_owsa <- pre_agecat_owsa %>%
    group_by(age_cat) %>%
    summarise(
      tot_pre_cases = sum(Median,      na.rm = TRUE),
      tot_pre_fatal = sum(fatal,      na.rm = TRUE),
      .groups = "drop"
    )
  
  pre_tot_cases_agegr1_owsa <- pre_tot_by_agecat_owsa$tot_pre_cases[ pre_tot_by_agecat_owsa$age_cat=="2-4" ]
  pre_tot_cases_agegr2_owsa <- pre_tot_by_agecat_owsa$tot_pre_cases[ pre_tot_by_agecat_owsa$age_cat=="5" ]
  pre_tot_cases_agegr3_owsa <- pre_tot_by_agecat_owsa$tot_pre_cases[ pre_tot_by_agecat_owsa$age_cat=="6-15" ]
  pre_tot_cases_agegr4_owsa <- pre_tot_by_agecat_owsa$tot_pre_cases[ pre_tot_by_agecat_owsa$age_cat=="16-20" ]
  
  pre_summary_cases_age_mid <- summ_res_pre$summary_cases_pre_age  
  
  pre_agecat_mid <- pre_summary_cases_age_mid %>%
    mutate(
      age_cat = case_when(
        AgeGroup %in% 2:4   ~ "2-4",
        AgeGroup %in% 5     ~ "5",
        AgeGroup %in% 6:15  ~ "6-15",
        AgeGroup %in% 16:20 ~ "16-20",
        TRUE                ~ NA_character_
      )
    )
  
  pre_tot_by_agecat_mid <- pre_agecat_mid %>%
    group_by(age_cat) %>%
    summarise(
      tot_pre_cases_mid = sum(Median, na.rm = TRUE),
      .groups = "drop"
    )
  
  pre_tot_cases_agegr1_mid <- pre_tot_by_agecat_mid$tot_pre_cases_mid[ pre_tot_by_agecat_mid$age_cat=="2-4" ]
  pre_tot_cases_agegr2_mid <- pre_tot_by_agecat_mid$tot_pre_cases_mid[ pre_tot_by_agecat_mid$age_cat=="5" ]
  pre_tot_cases_agegr3_mid <- pre_tot_by_agecat_mid$tot_pre_cases_mid[ pre_tot_by_agecat_mid$age_cat=="6-15" ]
  pre_tot_cases_agegr4_mid <- pre_tot_by_agecat_mid$tot_pre_cases_mid[ pre_tot_by_agecat_mid$age_cat=="16-20" ]
  
  # If we want a single "outcome metric" from pre-vacc, e.g. total cases:
  pre_tot_cases_owsa <- sum(pre_summary_cases_all_owsa$Median, na.rm = TRUE)
  
  
  # ---- 3. RUN POST-VACC SCENARIOS ----
  # We use run_simulation_scenarios_ui_owsa, passing in overrides from params
  scenario_result <- run_simulation_scenarios_ui_owsa(
    target_age_list    ,
    observed           = observed,
    N                  = N,
    bra_foi_state_summ = bra_foi_state_summ,
    age_groups         = age_groups,
    region_name        = region,
    
    # Pass rates
    hosp               = params$hosp,
    fatal              = params$fatal,
    nh_fatal           = params$nh_fatal,
    
    lhs_sample_young   = lhs_sample_young,
    lhs_old            = lhs_old,
    le_sample          = le_sample,
    age_gr_levels      = age_gr_levels,
    #prevacc_ui         = pre_summary_cases,   # optionally pass the pre-vacc detail
    posterior          = posterior,
    
    total_coverage     = params$total_coverage,
    delay = params$delay,
    # Overridden FOI, if you want post-vacc to also use a custom FOI:
    foi = params$foi,
    
    # Vaccine coverage parameters, possibly overridden
    VE_block           = params$VE_block,
    VE_inf             = params$VE_inf,
    weekly_delivery_speed = params$weekly_delivery_speed,
    time_until_immunity  = params$time_until_immunity,
    
    #base_beta        = baseline_params$base_beta,
    #I0               = baseline_params$$I0,
    #gamma            = baseline_params$$gamma,
    #rho              = baseline_params$$rho
    
  )
  
  scenario_week_sums <- lapply(scenario_result, function(one_scenario) {
    # For example, sum of weekly_median
    sum(one_scenario$weekly_df$weekly_median, na.rm = TRUE)
  })
  
  scenario_week_sums <- unlist(scenario_week_sums)
  
  #pre_cases_sum <- sum(sim_res_pre$weekly_cases_median, na.rm = TRUE)
  post_vacc_age_by_week <- lapply(scenario_result, function(one_scenario) {
    one_scenario$sim_df %>%
      summarise(
        # AgeGroup 1–4
        tot_post_cases_agegr1  = sum(Cases[AgeGroup %in% 2:4],    na.rm = TRUE),
        tot_post_fatal_agegr1  = sum(fatal[AgeGroup %in% 2:4],    na.rm = TRUE),
        # AgeGroup 5–12
        tot_post_cases_agegr2  = sum(Cases[AgeGroup %in% 5],   na.rm = TRUE),
        tot_post_fatal_agegr2  = sum(fatal[AgeGroup %in% 5],   na.rm = TRUE),
        # AgeGroup 13–18
        tot_post_cases_agegr3  = sum(Cases[AgeGroup %in% 6:15],  na.rm = TRUE),
        tot_post_fatal_agegr3  = sum(fatal[AgeGroup %in% 6:15],  na.rm = TRUE),
        # AgeGroup 13–18
        tot_post_cases_agegr4  = sum(Cases[AgeGroup %in% 16:20],  na.rm = TRUE),
        tot_post_fatal_agegr4  = sum(fatal[AgeGroup %in% 16:20],  na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  
  # ---- 4. RUN POST-VACC SCENARIO (MID values from baseline_params) ----
  # This simulation uses the unmodified baseline values as the mid scenario
  mid_scenario_result <- run_simulation_scenarios_ui_owsa(
    target_age_list    = target_age_list,
    observed           = observed,
    N                  = N,
    bra_foi_state_summ = bra_foi_state_summ,
    age_groups         = age_groups,
    region_name        = region,
    hosp               = baseline_params$hosp,
    fatal              = baseline_params$fatal,
    nh_fatal           = baseline_params$nh_fatal,
    lhs_sample_young   = lhs_sample_young,
    lhs_old            = lhs_old,
    le_sample          = le_sample,
    age_gr_levels      = age_gr_levels,
    posterior          = posterior,
    total_coverage     = baseline_params$total_coverage,
    foi                = baseline_params$foi,
    VE_block           = baseline_params$VE_block,
    VE_inf             = baseline_params$VE_inf,
    weekly_delivery_speed = baseline_params$weekly_delivery_speed,
    delay = baseline_params$delay,
    time_until_immunity = baseline_params$time_until_immunity
    #base_beta        = baseline_params$base_beta,
    #I0               = baseline_params$I0,
    #gamma            = baseline_params$gamma,
    #rho              = baseline_params$rho
  )
  
 #mid_week_sums <- lapply(mid_scenario_result, function(one_scenario) {
#    sum(one_scenario$weekly_df$weekly_median, na.rm = TRUE)
#  })
 # mid_week_sums <- unlist(mid_week_sums)
  
  
  mid_post_vacc_age_by_week <- lapply(mid_scenario_result, function(one_scenario) {
    one_scenario$sim_df %>% 
      ungroup() %>%       
      summarise(
        # AgeGroup 1–4
        tot_post_cases_agegr1 = sum(Cases[AgeGroup %in% 2:4],    na.rm = TRUE),
        # AgeGroup 5–12
        tot_post_cases_agegr2 = sum(Cases[AgeGroup %in% 5],   na.rm = TRUE),
        # AgeGroup 13–18
        tot_post_cases_agegr3 = sum(Cases[AgeGroup %in% 6:15],  na.rm = TRUE),
        # AgeGroup 13–18
        tot_post_cases_agegr4 = sum(Cases[AgeGroup %in% 16:20],  na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  return( list(
    # Pre-vacc overall
    pre_tot_cases                 = pre_tot_cases,
    pre_tot_cases_owsa            = pre_tot_cases_owsa,
    
    # Pre-vacc by age group
    pre_tot_by_agecat_owsa        = pre_tot_by_agecat_owsa,
    
    pre_tot_cases_agegr1_mid      = pre_tot_cases_agegr1_mid,
    pre_tot_cases_agegr2_mid      = pre_tot_cases_agegr2_mid,
    pre_tot_cases_agegr3_mid      = pre_tot_cases_agegr3_mid,
    pre_tot_cases_agegr4_mid      = pre_tot_cases_agegr4_mid,
    
    pre_tot_cases_agegr1_owsa     = pre_tot_cases_agegr1_owsa,
    pre_tot_cases_agegr2_owsa     = pre_tot_cases_agegr2_owsa,
    pre_tot_cases_agegr3_owsa     = pre_tot_cases_agegr3_owsa,
    pre_tot_cases_agegr4_owsa     = pre_tot_cases_agegr4_owsa,
    
    # Post-vacc scenarios
    scenario_result               = scenario_result,
    scenario_week_sums            = scenario_week_sums,
    post_vacc_age_by_week         = post_vacc_age_by_week,
    
    # Mid (baseline) scenario
    mid_scenario_result           = mid_scenario_result,
    mid_post_vacc_age_by_week     = mid_post_vacc_age_by_week
  ))
}


## updated run


run_simulation_with_params <- function(
    param_overrides,
    baseline_params,
    posterior,
    bra_foi_state_summ,
    age_groups,
    N,
    region,
    observed,
    age_gr_levels,
    lhs_sample_young,
    lhs_old,
    le_sample,
    target_age_list,
    total_coverage = 0.5,
    weekly_delivery_speed = 0.1,
    VE_block = 0.989,
    VE_inf = 0.989,
    delay = 2,
    time_until_immunity = 0
) {
  ## ----- tiny local helpers -----
  collapse_age4_local <- function(df, value_col = "Median") {
    df %>%
      dplyr::filter(AgeGroup >= 2, AgeGroup <= 20) %>% 
      dplyr::mutate(
        age_cat = dplyr::case_when(
          AgeGroup %in%  2:4  ~ "agegr1",
          AgeGroup %in%  5    ~ "agegr2",
          AgeGroup %in%  6:15 ~ "agegr3",
          AgeGroup %in% 16:20 ~ "agegr4",
          TRUE                ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(age_cat)) %>%
      dplyr::group_by(age_cat) %>%
      dplyr::summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(
        age_cat = c("agegr1","agegr2","agegr3","agegr4"),
        fill = list(val = 0)
      ) %>%
      dplyr::arrange(match(age_cat, c("agegr1","agegr2","agegr3","agegr4")))
  }
  get_age <- function(df, agegr) {
    v <- df$val[df$age_cat == agegr]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    v[1]
  }
  
  ## ----- merge overrides -----
  params <- modifyList(baseline_params, param_overrides)
  p_total_coverage      <- if (!is.null(params$total_coverage))      params$total_coverage      else total_coverage
  p_delay               <- if (!is.null(params$delay))               params$delay               else delay
  p_VE_block            <- if (!is.null(params$VE_block))            params$VE_block            else VE_block
  p_VE_inf              <- if (!is.null(params$VE_inf))              params$VE_inf              else VE_inf
  p_weekly_delivery     <- if (!is.null(params$weekly_delivery_speed))params$weekly_delivery_speed else weekly_delivery_speed
  p_time_until_immunity <- if (!is.null(params$time_until_immunity)) params$time_until_immunity else time_until_immunity
  p_foi                 <- if (!is.null(params$foi))                 params$foi                 else baseline_params$foi
  p_hosp                <- if (!is.null(params$hosp))                params$hosp                else baseline_params$hosp
  p_fatal               <- if (!is.null(params$fatal))               params$fatal               else baseline_params$fatal
  p_nh_fatal            <- if (!is.null(params$nh_fatal))            params$nh_fatal            else baseline_params$nh_fatal
  
  ## ----- PRE sims -----
  sim_res_pre <- simulate_pre_ui_age_owsa(
    posterior             = posterior,
    bra_foi_state_summ    = bra_foi_state_summ,
    age_groups            = age_groups,
    N                     = N,
    region                = region,
    observed              = observed,
    delay                 = 53,
    VE_block              = 0,
    VE_inf                = 0,
    target_age            = rep(0, 20),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    time_until_immunity   = 0,
    foi_value             = baseline_params$foi
  )
  sim_res_pre_owsa <- simulate_pre_ui_age_owsa(
    posterior             = posterior,
    bra_foi_state_summ    = bra_foi_state_summ,
    age_groups            = age_groups,
    N                     = N,
    region                = region,
    observed              = observed,
    delay                 = 53,
    VE_block              = 0,
    VE_inf                = 0,
    target_age            = rep(0, 20),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    time_until_immunity   = 0,
    foi_value             = p_foi
  )
  
  ## ----- PRE summaries -----
  summ_res_pre <- summarise_presim_ui_owsa(
    sim_result       = sim_res_pre,
    observed         = observed,
    age_gr_levels    = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = p_hosp,
    fatal            = p_fatal,
    nh_fatal         = p_nh_fatal,
    region           = region
  )
  summ_res_pre_owsa <- summarise_presim_ui_owsa(
    sim_result       = sim_res_pre_owsa,
    observed         = observed,
    age_gr_levels    = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = p_hosp,
    fatal            = p_fatal,
    nh_fatal         = p_nh_fatal,
    region           = region
  )
  
  pre_tot_cases_mid  <- sum(summ_res_pre$summary_cases_pre_all$Median,      na.rm = TRUE)
  pre_tot_cases_owsa <- sum(summ_res_pre_owsa$summary_cases_pre_all$Median, na.rm = TRUE)
  
  ## collapse → agegr1~4 + val
  pre_tot_by_agecat_mid  <- collapse_age4_local(summ_res_pre$summary_cases_pre_age,       "Median")
  pre_tot_by_agecat_owsa <- collapse_age4_local(summ_res_pre_owsa$summary_cases_pre_age, "Median")
  
  ## scalar pulls
  pre_tot_cases_agegr1_mid  <- get_age(pre_tot_by_agecat_mid,  "agegr1")
  pre_tot_cases_agegr2_mid  <- get_age(pre_tot_by_agecat_mid,  "agegr2")
  pre_tot_cases_agegr3_mid  <- get_age(pre_tot_by_agecat_mid,  "agegr3")
  pre_tot_cases_agegr4_mid  <- get_age(pre_tot_by_agecat_mid,  "agegr4")
  pre_tot_cases_agegr1_owsa <- get_age(pre_tot_by_agecat_owsa, "agegr1")
  pre_tot_cases_agegr2_owsa <- get_age(pre_tot_by_agecat_owsa, "agegr2")
  pre_tot_cases_agegr3_owsa <- get_age(pre_tot_by_agecat_owsa, "agegr3")
  pre_tot_cases_agegr4_owsa <- get_age(pre_tot_by_agecat_owsa, "agegr4")
  
  ## ----- POST sims OWSA -----
  scenario_result <- run_simulation_scenarios_ui_owsa(
    target_age_list        = target_age_list,
    observed               = observed,
    N                      = N,
    bra_foi_state_summ     = bra_foi_state_summ,
    age_groups             = age_groups,
    region_name            = region,
    hosp                   = p_hosp,
    fatal                  = p_fatal,
    nh_fatal               = p_nh_fatal,
    lhs_sample_young       = lhs_sample_young,
    lhs_old                = lhs_old,
    le_sample              = le_sample,
    age_gr_levels          = age_gr_levels,
    posterior              = posterior,
    total_coverage         = p_total_coverage,
    delay                  = p_delay,
    foi                    = p_foi,
    VE_block               = p_VE_block,
    VE_inf                 = p_VE_inf,
    weekly_delivery_speed  = p_weekly_delivery,
    time_until_immunity    = p_time_until_immunity
  )
  post_vacc_age_by_week <- purrr::map(scenario_result, function(one_scenario) {
    one_scenario$sim_df %>%
      dplyr::ungroup() %>%
      dplyr::summarise(
        tot_post_cases_agegr1 = sum(Cases[AgeGroup %in%  2:4],  na.rm = TRUE),
        tot_post_cases_agegr2 = sum(Cases[AgeGroup %in%  5],    na.rm = TRUE),
        tot_post_cases_agegr3 = sum(Cases[AgeGroup %in%  6:15], na.rm = TRUE),
        tot_post_cases_agegr4 = sum(Cases[AgeGroup %in% 16:20], na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  ## ----- POST sims MID -----
  mid_scenario_result <- run_simulation_scenarios_ui_owsa(
    target_age_list        = target_age_list,
    observed               = observed,
    N                      = N,
    bra_foi_state_summ     = bra_foi_state_summ,
    age_groups             = age_groups,
    region_name            = region,
    hosp                   = baseline_params$hosp,
    fatal                  = baseline_params$fatal,
    nh_fatal               = baseline_params$nh_fatal,
    lhs_sample_young       = lhs_sample_young,
    lhs_old                = lhs_old,
    le_sample              = le_sample,
    age_gr_levels          = age_gr_levels,
    posterior              = posterior,
    total_coverage         = baseline_params$total_coverage,
    delay                  = baseline_params$delay,
    foi                    = baseline_params$foi,
    VE_block               = baseline_params$VE_block,
    VE_inf                 = baseline_params$VE_inf,
    weekly_delivery_speed  = baseline_params$weekly_delivery_speed,
    time_until_immunity    = baseline_params$time_until_immunity
  )
  mid_post_vacc_age_by_week <- purrr::map(mid_scenario_result, function(one_scenario) {
    one_scenario$sim_df %>%
      dplyr::ungroup() %>%
      dplyr::summarise(
        tot_post_cases_agegr1 = sum(Cases[AgeGroup %in%  2:4],  na.rm = TRUE),
        tot_post_cases_agegr2 = sum(Cases[AgeGroup %in%  5],    na.rm = TRUE),
        tot_post_cases_agegr3 = sum(Cases[AgeGroup %in%  6:15], na.rm = TRUE),
        tot_post_cases_agegr4 = sum(Cases[AgeGroup %in% 16:20], na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  ## ----- return -----
  list(
    pre_tot_cases_mid             = pre_tot_cases_mid,
    pre_tot_cases_owsa            = pre_tot_cases_owsa,
    pre_tot_cases_agegr1_mid      = pre_tot_cases_agegr1_mid,
    pre_tot_cases_agegr2_mid      = pre_tot_cases_agegr2_mid,
    pre_tot_cases_agegr3_mid      = pre_tot_cases_agegr3_mid,
    pre_tot_cases_agegr4_mid      = pre_tot_cases_agegr4_mid,
    pre_tot_cases_agegr1_owsa     = pre_tot_cases_agegr1_owsa,
    pre_tot_cases_agegr2_owsa     = pre_tot_cases_agegr2_owsa,
    pre_tot_cases_agegr3_owsa     = pre_tot_cases_agegr3_owsa,
    pre_tot_cases_agegr4_owsa     = pre_tot_cases_agegr4_owsa,
    pre_tot_by_agecat_mid         = pre_tot_by_agecat_mid,
    pre_tot_by_agecat_owsa        = pre_tot_by_agecat_owsa,
    post_vacc_age_by_week         = post_vacc_age_by_week,
    mid_post_vacc_age_by_week     = mid_post_vacc_age_by_week,
    scenario_result               = scenario_result,
    mid_scenario_result           = mid_scenario_result
  )
}



# --- 4. Run one-way sensitivity analysis for each parameter ---

# Set dynamic parameters
region <- "Bahia"
hosp <- hosp      # default hospitalisation parameter
fatal <- fatal    # default fatality parameter
nh_fatal <- nh_fatal # default non-hospital fatality parameter

params_list <- make_params(region, bra_foi_state_summ, posterior_bh, lhs_sample_young, hosp, fatal, nh_fatal)
baseline_params    <- params_list$baseline_params
sensitivity_ranges <- params_list$sensitivity_ranges

# Sensitivity analysis loop: iterate over each parameter and run simulation with each value

library(dplyr)
library(purrr)
library(tibble)

# updated sensitivity results
sensitivity_results <- list()

for (param in names(sensitivity_ranges)) {
  param_values <- sensitivity_ranges[[param]]
  
  # 각 파라미터값마다 run_simulation…() 호출 후 tibble(row) 생성
  df_list <- lapply(param_values, function(val) {
    res <- run_simulation_with_params(
      param_overrides      = setNames(list(val), param),
      baseline_params      = baseline_params,
      posterior            = posterior_bh,
      bra_foi_state_summ   = bra_foi_state_summ,
      age_groups           = age_groups,
      N                    = N_bahia$Bahia,
      region               = region,
      observed             = observed_bh,
      age_gr_levels        = age_gr_levels,
      lhs_sample_young     = lhs_sample_young,
      lhs_old              = lhs_old,
      le_sample            = le_sample,
      target_age_list      = target_age_list
    )
    
    # 한 파라미터값에 대한 한 행(row) tibble
    tibble(
      Parameter                    = param,
      Value                        = val,
      pre_tot_cases_agegr1_owsa    = res$pre_tot_cases_agegr1_owsa,
      pre_tot_cases_agegr2_owsa    = res$pre_tot_cases_agegr2_owsa,
      pre_tot_cases_agegr3_owsa    = res$pre_tot_cases_agegr3_owsa,
      pre_tot_cases_agegr4_owsa    = res$pre_tot_cases_agegr4_owsa,
      
      pre_tot_cases_agegr1_mid    = res$pre_tot_cases_agegr1_mid,
      pre_tot_cases_agegr2_mid    = res$pre_tot_cases_agegr2_mid,
      pre_tot_cases_agegr3_mid    = res$pre_tot_cases_agegr3_mid,
      pre_tot_cases_agegr4_mid    = res$pre_tot_cases_agegr4_mid,
      
      # list‐column으로 전체 mid_post_vacc_age_by_week 붙이기
      post_vacc_age_by_week     = list(res$post_vacc_age_by_week),
      mid_post_vacc_age_by_week    = list(res$mid_post_vacc_age_by_week)
    )
  })
  
  # 리스트 of tibble → 하나의 data.frame 으로 바인딩
  sensitivity_results[[param]] <- bind_rows(df_list)
}



extract_owsa <- lapply(sensitivity_results, function(x) {
  list(
    foi = x$foi,
    lo = list(
      scenario_week_sums = x$Outcome.1$scenario_week_sums,
      pre_tot_cases_owsa = x$Outcome.1$pre_tot_cases_owsa
    ),
    pre_tot_cases = x$Outcome.1$pre_tot_cases,
    scenario_week_sums = x$Outcome.1$scenario_week_sums,
    mid_week_sums = x$Outcome.1$mid_week_sums,
    pre_cases_sum = x$Outcome.1$pre_cases_sum,
    hi = list(
      pre_tot_cases = x$Outcome.2$pre_tot_cases,
      scenario_week_sums = x$Outcome.2$scenario_week_sums,
      mid_week_sums = x$Outcome.2$mid_week_sums,
      pre_cases_sum = x$Outcome.2$pre_cases_sum,
      pre_tot_cases_owsa = x$Outcome.2$pre_tot_cases_owsa
    ),
    VE_block = x$VE_block,
    weekly_delivery_speed = x$weekly_delivery_speed
  )
})

# updated extract owsa
extract_owsa <- lapply(sensitivity_results, function(df) {
  # Low=1, High=2
  scen_names <- c("lo","hi")[ seq_len(nrow(df)) ]
  out <- vector("list", length = length(scen_names))
  names(out) <- scen_names
  
  for(i in seq_along(scen_names)) {
    # 1) 파라미터값
    val <- df$Value[i]
    
    # 2) pre‐vacc OWSA & MID, 연령그룹1–3
    pre_owsa <- df[ i, paste0("pre_tot_cases_agegr", 1:4, "_owsa") ]
    pre_mid  <- df[ i, paste0("pre_tot_cases_agegr", 1:4, "_mid") ]
    
    # 3) post‐vacc 시나리오 결과 (list of 3 tibbles)
    post_list <- df$post_vacc_age_by_week[[i]]
    post_vals <- sapply(seq_along(post_list), function(j) {
      tb <- post_list[[j]]
      tb[[ paste0("tot_post_cases_agegr", j) ]]
    })
    
    # 4) mid‐post‐vacc 시나리오 결과
    mid_list  <- df$mid_post_vacc_age_by_week[[i]]
    mid_vals  <- sapply(seq_along(mid_list), function(j) {
      tb <- mid_list[[j]]
      tb[[ paste0("tot_post_cases_agegr", j) ]]
    })
    
    out[[i]] <- list(
      parameter   = df$Parameter[i],
      value       = val,
      pre_owsa    = setNames(as.numeric(pre_owsa), paste0("agegrp", 1:4)),
      pre_mid     = setNames(as.numeric(pre_mid),  paste0("agegrp", 1:4)),
      post_vacc   = setNames(post_vals,            paste0("agegrp", 1:4)),
      mid_vacc    = setNames(mid_vals,             paste0("agegrp", 1:4))
    )
  }
  
  out
})


# Extract desired outcome components from each sensitivity result data frame
tornado_data <- bind_rows(lapply(names(extract_owsa), function(param_name) {
  lo_vec <- extract_owsa[[param_name]]$lo$scenario_week_sums
  hi_vec <- extract_owsa[[param_name]]$hi$scenario_week_sums
  
  # mid_week_sums and pre_cases_sum are outside lo/hi in your original list
  mid_vec <- extract_owsa[[param_name]]$mid_week_sums
  pre_cases_sum <- extract_owsa[[param_name]]$pre_cases_sum
  pre_tot_cases_owsa_lo <- extract_owsa[[param_name]]$lo$pre_tot_cases_owsa
  pre_tot_cases_owsa_hi <- extract_owsa[[param_name]]$hi$pre_tot_cases_owsa
  n_scenarios <- length(mid_vec)
  
  tibble(
    Parameter = param_name,
    Scenario = paste0("Scenario_", seq_len(n_scenarios)),
    Low = lo_vec,
    High = hi_vec,
    Mid = mid_vec,
    Pre_Cases_Sum = rep(pre_cases_sum, n_scenarios),  # repeat to match row count
    Pre_Tot_Cases_OWSA_Lo = rep(pre_tot_cases_owsa_lo, n_scenarios),
    Pre_Tot_Cases_OWSA_Hi = rep(pre_tot_cases_owsa_hi, n_scenarios)
  )
})) %>%
  mutate(
    Mid_Impact = Pre_Cases_Sum - Mid,
    Lo_Impact = Pre_Tot_Cases_OWSA_Lo - Low,
    Hi_Impact = Pre_Tot_Cases_OWSA_Hi - High,
    
    # Percent impacts
    Mid_Impact_pct = (Mid_Impact / Pre_Cases_Sum) * 100,
    Lo_Impact_pct  = (Lo_Impact / Pre_Tot_Cases_OWSA_Lo) * 100,
    Hi_Impact_pct  = (Hi_Impact / Pre_Tot_Cases_OWSA_Hi) * 100
  )


# updated tornado data
tornado_data <- map_dfr(names(extract_owsa), function(param) {
  ex <- extract_owsa[[param]]
  # low/hi 각각의 named 벡터
  lo <- ex$lo
  hi <- ex$hi
  # agegrp1–3
  agegrps <- names(lo$pre_owsa)
  tibble(
    Parameter        = param,
    Scenario         = agegrps,
    pre_owsa_lo      = as.numeric(lo$pre_owsa),
    pre_owsa_hi      = as.numeric(hi$pre_owsa),
    pre_mid          = as.numeric(lo$pre_mid),     # mid는 lo/hi 같으므로 lo만 사용
    post_owsa_lo     = as.numeric(lo$post_vacc),
    post_owsa_hi     = as.numeric(hi$post_vacc),
    mid_vacc         = as.numeric(lo$mid_vacc),    # mid_vacc도 lo/hi 같음
    mid_diff         = pre_mid - mid_vacc,
    lo_diff          = pre_owsa_lo  - post_owsa_lo,
    hi_diff          = pre_owsa_hi  - post_owsa_hi
  )
})


tornado_data <- tornado_data %>%
  mutate(
    totvacc = case_when(
      Scenario == "agegrp1" ~ 1883349,
      Scenario == "agegrp2" ~ 4025158,
      Scenario == "agegrp3" ~ 1079639,
      TRUE                  ~ NA_real_
    )
  ) %>% mutate(
    mid_impact = mid_diff / totvacc * 1e6,
    lo_impact  = lo_diff / totvacc * 1e6,
    hi_impact  = hi_diff / totvacc * 1e6
    
  )

tornado_long <- tornado_data %>%
  pivot_longer(
    cols      = c(lo_impact, hi_impact),
    names_to  = "Direction",
    values_to = "Impact"
  ) %>%
  mutate(
    Direction = dplyr::recode(Direction,
                       lo_impact = "Lower",
                       hi_impact = "Upper")
  )

tornado_long_b <- tornado_long %>% filter(Parameter != "foi")

ggplot(tornado_long_b, aes(y = fct_rev(Parameter))) +
  # Lower/Upper 를 색깔별로
  geom_segment(aes(x = mid_impact, xend = Impact, yend = fct_rev(Parameter),
                   color = Direction),
               size = 5) +
  # mid 점
  geom_point(data = tornado_long_b,
             aes(x = mid_impact),
             shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(x = "Impact per million vaccinated", y = NULL,
       color = "Direction") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position  = "top"
  )


# Long format for plotting with Lo_Impact and Hi_Impact
tornado_long <- tornado_data %>%
  select(Parameter, Scenario, Mid_Impact, Lo_Impact, Hi_Impact) %>%
  pivot_longer(cols = c(Lo_Impact, Hi_Impact), names_to = "Type", values_to = "Impact") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact", "Lower", "Upper"),
    Scenario = as.factor(Scenario)
  )

tornado_long_pct <- tornado_data %>%
  select(Parameter, Scenario, Mid_Impact_pct, Lo_Impact_pct, Hi_Impact_pct) %>%
  pivot_longer(cols = c(Lo_Impact_pct, Hi_Impact_pct), names_to = "Type", values_to = "Impact") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact_pct", "Lower", "Upper"),
    Scenario = as.factor(Scenario)
  )

param_order <- tornado_long %>%
  group_by(Scenario, Parameter) %>%
  summarise(MaxImpact = max(abs(Impact)), .groups = "drop") %>%
  arrange(Scenario, desc(MaxImpact)) %>%
  group_by(Scenario) %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

param_order_pct <- tornado_long_pct %>%
  group_by(Scenario, Parameter) %>%
  summarise(MaxImpact = max(abs(Impact)), .groups = "drop") %>%
  arrange(Scenario, desc(MaxImpact)) %>%
  group_by(Scenario) %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

# Join back
tornado_long <- left_join(tornado_long, param_order, by = c("Scenario", "Parameter"))

tornado_long_pct <- left_join(tornado_long_pct, param_order_pct, by = c("Scenario", "Parameter"))

tornado_bar <- tornado_data %>%
  select(Parameter, Scenario, Mid_Impact, Lo_Impact, Hi_Impact) %>%
  pivot_longer(cols = c(Lo_Impact, Hi_Impact), names_to = "Type", values_to = "Value") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact", "Lower", "Upper"),
    Bar_start = Mid_Impact,
    Bar_end = Value,
    Scenario = as.factor(Scenario)
  )

tornado_bar_pct <- tornado_data %>%
  select(Parameter, Scenario, Mid_Impact_pct, Lo_Impact_pct, Hi_Impact_pct) %>%
  pivot_longer(cols = c(Lo_Impact_pct, Hi_Impact_pct), names_to = "Type", values_to = "Value") %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact_pct", "Lower", "Upper"),
    Bar_start = Mid_Impact_pct,
    Bar_end = Value,
    Scenario = as.factor(Scenario)
  )

# Plot
ggplot(tornado_bar, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter), color = Direction),
               linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "Averted cases",
    y = "Parameter",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )

ggplot(tornado_bar_pct, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter), color = Direction),
               linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "% Averted cases",
    y = "Parameter",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )


ggplot(tornado_bar, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter), color = Direction),
               linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "% Averted cases",
    y = "Parameter",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )


ggplot(tornado_bar, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter), color = Direction),
               linewidth = 5) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white", color = "black", size = 2.5) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_color_manual(values = c("Lower" = "#F8766D", "Upper" = "#00BFC4")) +
  labs(
    x = "% Averted cases",
    y = "Parameter",
    color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )



