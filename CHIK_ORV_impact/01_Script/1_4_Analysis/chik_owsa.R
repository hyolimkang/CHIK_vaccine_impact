### one-way sensitivity analysis ----------------------------------------------

# ----- 1. Define baseline values and sensitivity ranges -----
make_params <- function(region, bra_foi_state_summ, posterior, lhs_sample_young, hosp, fatal, nh_fatal) {
  # Create baseline parameters based on the specified region
  baseline_params <- list(
    foi          = bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region],
    hosp         = hosp,
    fatal        = fatal,
    nh_fatal     = nh_fatal,
    VE_block     = 0.75,
    weekly_delivery_speed = 0.1
  )
  
  # Define sensitivity ranges for each parameter. (Ensure these ranges are wide enough to affect outcomes.)
  sensitivity_ranges <- list(
    foi          = c(bra_foi_state_summ$foi_lo[bra_foi_state_summ$NAME_1 == region],
                     bra_foi_state_summ$foi_hi[bra_foi_state_summ$NAME_1 == region]),
    VE_block    = c(0.5, 0.9),
    weekly_delivery_speed = c(0.09, 0.11)
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
    A = 18, 
    delay = 53,
    VE_block = 0,
    target_age = rep(0, A),
    coverage_threshold = 0,
    total_coverage = 0,
    total_supply = 0,
    weekly_delivery_speed = 0,
    
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
    
    # Use the user-specified FOI value
    R0_vec <- 1 - exp(-foi_value * age_groups)
    
    # Run your custom SIRV simulation function
    sim_out <- sirv_sim_coverageSwitch(
      T = T,
      A = A,
      N = N,
      r = rep(0, A),
      base_beta = base_beta_draw,
      I0_draw   = I0_draw,
      R0        = R0_vec,
      rho       = rho_draw,
      gamma     = gamma_draw,
      delay     = delay,
      VE_block  = VE_block,
      target_age = target_age,
      coverage_threshold = coverage_threshold,
      total_coverage = total_coverage,
      weekly_delivery_speed = weekly_delivery_speed
    )
    
    # Each draw returns an age-stratified cases matrix [A x T]
    sim_results_list[[i]] <- sim_out$age_stratified_cases
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
    age_gr_levels = c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                      "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                      "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                      "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                      "80-84 years", "85-89 years"),
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
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
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
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
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
    total_coverage = 0.5,
    
    # Remaining possible overrides
    foi                    = NULL,
    VE_block              = NULL,
    weekly_delivery_speed = NULL
    
) {
  # Number of scenarios = length of target_age_list
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  
  # Number of posterior draws
  n_draws <- length(posterior$gamma)
  
  # Number of weeks = rows in observed data
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<5 years", "5-9 years", "10-14 years", "15-19 years",
    "20-24 years", "25-29 years", "30-34 years", "35-39 years",
    "40-44 years", "45-49 years", "50-54 years", "55-59 years",
    "60-64 years", "65-69 years", "70-74 years", "75-79 years",
    "80-84 years", "85-89 years"
  )
  
  # If user doesn't specify foi, use the region-specific default from bra_foi_state_summ
  if (is.null(foi)) {
    foi <- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name]
  }
  if (is.null(hosp))     hosp     <- hosp
  if (is.null(fatal))    fatal    <- fatal
  if (is.null(nh_fatal)) nh_fatal <- nh_fatal
  
  # Vaccine overrides
  if (is.null(VE_block))              VE_block <- 0.75
  if (is.null(weekly_delivery_speed)) weekly_delivery_speed <- 0.1
  
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
      
      # Compute R0 using the (possibly overridden) foi
      R0_vec <- 1 - exp(-foi * age_groups)
      
      sim_out <- sirv_sim_coverageSwitch(
        T = T,
        A = 18,
        N = N,
        r = rep(0, 18),
        base_beta  = base_beta_draw,
        I0         = I0_draw,
        R0         = R0_vec,
        rho        = rho_draw,
        gamma      = gamma_draw,
        delay      = 2,   # Hard-coded delay
        VE_block   = VE_block,
        coverage_threshold = 1,
        target_age         = target,
        total_coverage     = total_coverage,
        weekly_delivery_speed = weekly_delivery_speed
      )
      
      draw_results[[d]] <- sim_out$age_stratified_cases  # A x T matrix
    }
    
    # Combine draws into a 3D array [18, T, n_draws]
    age_array <- array(unlist(draw_results), dim = c(18, T, n_draws))
    
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
    
    # If you still want YLL, keep this (otherwise remove):
    sim_df <- sim_df %>%
      dplyr::mutate(
        age_numeric = dplyr::case_when(
          stringr::str_detect(as.character(age_gr), "<5") ~ 0,
          TRUE ~ as.numeric(stringr::str_extract(as.character(age_gr), "^\\d+"))
        ),
        # Simple life-expectancy logic
        le_left = dplyr::case_when(
          age_numeric %in% c(0,5)   ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10,15) ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20,25) ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30,35) ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40,45) ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50,55) ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60,65) ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70,75) ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80,85) ~ quantile(le_sample$le_9, 0.5)
        ),
        yll     = fatal * le_left,
        yll_lo  = fatal_lo * le_left,
        yll_hi  = fatal_hi * le_left
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
    weekly_delivery_speed = NULL
    #base_beta,
    #I0,
    #gamma,
    #rho
) {
  # Number of weeks
  T <- nrow(observed)
  
  # Default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
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
      dplyr::group_by(AgeGroup) %>%  
      dplyr::summarise(
        Median          = sum(Cases, na.rm = TRUE),
        low95           = sum(post_low95, na.rm = TRUE),
        hi95            = sum(post_hi95, na.rm = TRUE),
        hospitalised    = sum(hospitalised, na.rm = TRUE),
        hospitalised_lo = sum(hospitalised_lo, na.rm = TRUE),
        hospitalised_hi = sum(hospitalised_hi, na.rm = TRUE),
        
        fatal      = sum(fatal, na.rm = TRUE),
        fatal_lo   = sum(fatal_lo, na.rm = TRUE),
        fatal_hi   = sum(fatal_hi, na.rm = TRUE),
        
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        # Attach pre-vacc values from pre_summary_cases_age
        pre_vacc  = pre_summary_cases_age$Median,
        pre_fatal = pre_summary_cases_age$fatal,
        
        # Differences for Cases
        diff   = pre_vacc - Median,
        impact = diff / pre_vacc * 100,
        
        # Differences for Fatal
        diff_fatal   = pre_fatal - fatal,
        impact_fatal = diff_fatal / pre_fatal * 100
      )
  })
  
  # Summarize by Week for entire population
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      dplyr::mutate(Scenario = paste0("Scenario_", i),
                    pre_vacc = pre_summary_cases$Median) %>% 
      dplyr::group_by(Week, Scenario) %>%
      dplyr::summarise(
        post_cases = sum(Cases, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        
        post_fatal = sum(fatal, na.rm = TRUE),
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        
        Scenario   = dplyr::first(Scenario),
        .groups    = "drop"
      )
  })
  
  summary_week_df <- do.call(rbind, summary_week)
  
  # Attach scenario-level weekly UI if it exists
  for (i in seq_len(n_scenarios)) {
    scenario_ui <- scenario_result[[i]]$weekly_df %>%
      dplyr::rename(
        post_weekly_median = weekly_median,
        post_weekly_low95  = weekly_low95,
        post_weekly_hi95   = weekly_hi95
      ) %>%
      dplyr::distinct(Week, .keep_all = TRUE)
    summary_week[[i]] <- dplyr::left_join(
      summary_week[[i]], scenario_ui, by = "Week", relationship = "many-to-many"
    )
  }
  summary_week_df <- do.call(rbind, summary_week)
  
  # Merge aggregated pre-vacc UI into the weekly summary
  if (all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      dplyr::distinct(Week, .keep_all = TRUE)
    summary_week_df <- dplyr::left_join(
      summary_week_df, pre_weekly_df, by = "Week", relationship = "many-to-many"
    )
  }
  
  summary_week_df$region <- region
  
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
      weekly_delivery_speed = weekly_delivery_speed,
      base_beta             = base_beta,
      gamma                 = gamma,
      I0                    = I0,
      rho                   = rho
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
    VE_block = 0.75
    
    #base_beta,
    #I0,
    #gamma,
    #rho
) {
  # ---- 0. Combine baseline + overrides ----
  params <- modifyList(baseline_params, param_overrides)
  
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
    target_age            = rep(0, 18),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    
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
    target_age            = rep(0, 18),
    coverage_threshold    = 0,
    total_coverage        = 0,
    weekly_delivery_speed = 0,
    
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
    
    total_coverage     = 0.5,
    
    # Overridden FOI, if you want post-vacc to also use a custom FOI:
    foi = params$foi,
    
    # Vaccine coverage parameters, possibly overridden
    VE_block           = params$VE_block,
    weekly_delivery_speed = params$weekly_delivery_speed
    
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
  
  pre_cases_sum <- sum(sim_res_pre$weekly_cases_median, na.rm = TRUE)
  
  
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
    total_coverage     = 0.5,
    foi                = baseline_params$foi,
    VE_block           = baseline_params$VE_block,
    weekly_delivery_speed = baseline_params$weekly_delivery_speed
    
    #base_beta        = baseline_params$base_beta,
    #I0               = baseline_params$I0,
    #gamma            = baseline_params$gamma,
    #rho              = baseline_params$rho
  )
  
  mid_week_sums <- lapply(mid_scenario_result, function(one_scenario) {
    sum(one_scenario$weekly_df$weekly_median, na.rm = TRUE)
  })
  mid_week_sums <- unlist(mid_week_sums)
  
  return(list(
    # Pre-vacc:
    pre_tot_cases    = pre_tot_cases,
    pre_tot_cases_owsa    = pre_tot_cases_owsa,
    # Scenarios:
    scenario_result  = scenario_result,
    scenario_week_sums = scenario_week_sums,
    
    #mid_scenario_result  = mid_scenario_result,
    mid_week_sums        = mid_week_sums,
    
    pre_cases_sum   = pre_cases_sum
    
  ))
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
sensitivity_results <- list()
for (param in names(sensitivity_ranges)) {
  param_values <- sensitivity_ranges[[param]]
  impacts <- sapply(param_values, function(val) {
    # Override the current parameter with the value 'val'
    run_simulation_with_params(
      target_age_list = target_age_list,
      param_overrides = setNames(list(val), param),
      baseline_params = baseline_params,
      posterior = posterior_bh,
      bra_foi_state_summ = bra_foi_state_summ,
      age_groups = age_groups,
      N = N_bahia$Bahia,
      region = region,      
      observed = observed_bh,
      age_gr_levels = age_gr_levels,
      lhs_sample_young = lhs_sample_young,
      lhs_old = lhs_old,
      le_sample = le_sample
    )
  })
  
  sensitivity_results[[param]] <- data.frame(
    Parameter = rep(param, length(param_values)),
    Value = param_values,
    Outcome = impacts
  )
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
