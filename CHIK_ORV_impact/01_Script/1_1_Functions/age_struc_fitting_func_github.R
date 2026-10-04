age_groups <- c(mean(0:1),
                mean(1:4),
                mean(5:9),
                mean(10:11),
                mean(12:17),
                mean(18:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)

age_gr_levels = c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
)

age_gr <- rep(c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
), 52)

age_gr_levels <- gsub("[–—]", "-", age_gr_levels)
age_gr        <- gsub("[–—]", "-", age_gr)

# function for post-process posterior prevacc data ----------------------------------------------------

create_summary_df <- function(
    fit_prevacc,       # 3D array: [iteration, week, age]
    bra_sum,
    weeks        = 1:nrow(bra_sum),   # Vector of week indices
    age_groups   = 1:20,   # Vector of age group indices
    scenario     = "Pre-Vaccination",
    region
) {
  # 1) Summarize the posterior samples across iterations
  posterior_array <- extract(fit_prevacc)
  age_cases <- posterior_array$age_stratified_cases
  
  summary_array <- apply(age_cases, c(2, 3), function(x) {
    c(
      median = median(x, na.rm = TRUE),
      lower  = quantile(x, 0.025, na.rm = TRUE),
      upper  = quantile(x, 0.975, na.rm = TRUE)
    )
  })
  
  # 2) Label the summary rows
  rownames(summary_array) <- c("median", "lower", "upper")
  
  # 3) Create a DataFrame for Week x AgeGroup
  df_out <- expand.grid(
    Week     = weeks,
    AgeGroup = age_groups
  )
  
  # 4) Fill the DataFrame from the summary array
  #    Notice the transpose-like call: aperm(..., c(2, 1))
  #    ensures we map [week, age] to the rows of df_out in the correct order
  df_out$Median <- as.vector(aperm(summary_array["median", , ], c(2, 1)))
  df_out$Lower  <- as.vector(aperm(summary_array["lower",  , ], c(2, 1)))
  df_out$Upper  <- as.vector(aperm(summary_array["upper",  , ], c(2, 1)))
  
  # 5) Add scenario name
  df_out$Scenario <- scenario
  
  df_summ <- df_out %>% group_by(Week) %>% summarise(
    Median = sum(Median),
    Lower =  sum(Lower),
    Upper = sum(Upper)
  )
  
  df_summ$Type <- "Predicted"
  df_summ$region <- region
  
  # Observed weekly cases
  observed <- data.frame(
    Week = 1:nrow(bra_sum),
    Observed = bra_sum$cases
  )
  
  observed$Type <- "Observed"
  
  
  # 6) Return the DataFrame
  return(list(df_out  = df_out,
              df_summ = df_summ,
              observed = observed
  ))
}

## overall fit 
overall_fit_gg <- function(
    observed,
    df_summ
) {
  
  g <-  ggplot() +
    # Observed points
    geom_point(data = observed, aes(x = Week, y = Observed, color = Type), size = 2) +
    
    # Predicted line
    geom_line(data = df_summ, aes(x = Week, y = Median, color = Type), size = 1) +
    
    # Prediction interval (ribbon)
    geom_ribbon(data = df_summ, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
                alpha = 0.2) +
    
    # Set colors using Brewer palette
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(labels = comma) + 
    # Labels & Theme
    labs(x = "Week", y = "Total cases") +
    theme_minimal()
  
  return(g)
}

## age strat
age_strat_gg <- function(
    df_age
){
  g <- ggplot(df_age, aes(x = Week, y = Median, color = Scenario)) +
    geom_line(size = 1) +
    facet_wrap(~ AgeGroup, ncol = 3)+
    theme_minimal()
  return(g)
}

## age struc fit
age_struc_fit <- function(
    fit_prevacc,
    model_age_bins,
    bra_expanded,
    bra_sum
){
  
  # 1) summary df
  posterior_array <- extract(fit_prevacc)
  age_cases <- posterior_array$age_stratified_cases
  
  summary_array <- apply(age_cases, c(2, 3), function(x) {
    c(
      median = median(x, na.rm = TRUE),
      lower  = quantile(x, 0.025, na.rm = TRUE),
      upper  = quantile(x, 0.975, na.rm = TRUE)
    )
  })
  
  # 2) Label the summary rows
  rownames(summary_array) <- c("median", "lower", "upper")
  
  
  # 3) wide format
  pred_df_wide <- as.data.frame(summary_array[1, , ])
  colnames(pred_df_wide) <- paste0("Week", 1:nrow(bra_sum))
  pred_df_wide$AgeGroup  <- model_age_bins
  
  pred_df_long <- pred_df_wide %>%
    pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Predicted") %>%
    mutate(Week = as.integer(gsub("Week", "", Week)))
  
  # Lower 95% interval
  lower_df_wide <- as.data.frame(summary_array[2, , ])
  colnames(lower_df_wide) <- paste0("Week", 1:nrow(bra_sum))
  lower_df_wide$AgeGroup <- model_age_bins
  lower_df_long <- lower_df_wide %>%
    pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Lower") %>%
    mutate(Week = as.integer(gsub("Week", "", Week)))
  
  # Upper 95% interval
  upper_df_wide <- as.data.frame(summary_array[3, , ])
  colnames(upper_df_wide) <- paste0("Week", 1:nrow(bra_sum))
  upper_df_wide$AgeGroup <- model_age_bins
  upper_df_long <- upper_df_wide %>%
    pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Upper") %>%
    mutate(Week = as.integer(gsub("Week", "", Week)))
  
  # Merge predicted, lower, and upper into one data frame:
  predicted_df <- pred_df_long %>%
    left_join(lower_df_long, by = c("AgeGroup", "Week")) %>%
    left_join(upper_df_long, by = c("AgeGroup", "Week"))
  
  observed_df <- bra_expanded[,c(4,6,5)]
  colnames(observed_df) <- c("Week", "Observed", "AgeGroup")
  
  combined_df <- predicted_df %>%
    left_join(observed_df, by = c("AgeGroup", "Week"))
  
  return(combined_df)
  
}

# plot age structured with 95% UI
agestrat_ui_gg <- function(combined_df){
  
  g <-  ggplot(combined_df, aes(x = Week)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.3) +
    geom_line(aes(y = Predicted), color = "blue", size = 1) +
    geom_point(aes(y = Observed), color = "red", size = 1) +
    facet_wrap(~ AgeGroup) +
    labs(title = "Observed vs Predicted Cases by Age Group",
         x = "Week",
         y = "Cases") +
    theme_minimal()
  
  return(g)
  
}

# extract stan params
extract_params <- function(fit_prevacc){
  
  posterior_prevacc <- rstan::extract(fit_prevacc)
  base_beta       <- apply(posterior_prevacc$base_beta, 2, median)       # length T
  I0              <- apply(posterior_prevacc$I0, 2, median)           
  gamma           <- median(posterior_prevacc$gamma)
  rho             <- median(posterior_prevacc$rho, 2, median)
  
  return(list(
    posterior_prevacc = posterior_prevacc,
    base_beta         = base_beta,
    I0                = I0,
    gamma             = gamma,
    rho               = rho
  ))
}

## prevacc simulation
simulate_pre_ui_age <- function(posterior, bra_foi_state_summ, age_groups, N, region, observed,
                                A = 20, 
                                delay = 53, VE_block = 0,
                                target_age = rep(0, A),
                                coverage_threshold = 0, total_coverage = 0,
                                total_supply = 0, weekly_delivery_speed = 0,
                                lhs_sample) {
  # Dynamically set T 
  T <- nrow(observed)
  
  # Determine the number of draws
  # 1. based on pre-defined posterior numbers
  posterior_idx <- posterior_idx_list[[region]]
  lhs_idx       <- lhs_idx_list[[region]]
  n_draws <- length(posterior_idx)
  # 2. based on posterior sample nums
  #n_draws <- length(posterior$gamma) 
  # 3. based on lhs_sapmle run 
  #n_draws <- nrow(lhs_sample)
  sim_results_list_rawinf <- vector("list", n_draws)
  sim_results_list_rawsymp <- vector("list", n_draws)
  
  set.seed(123)
  #foi_draws <- foi_draws_list[[region]]
  
  #I0_resampled        <- posterior$I0[posterior_idx, ]        
  #base_beta_resampled <- posterior$base_beta[posterior_idx, ] 
  
  # Run simulation for each draw and store the full age_stratified_cases matrix (dimensions: A x T)
  for (i in seq_len(n_draws)) {
    
    #idx <- posterior_idx[i] 
    idx_post <- posterior_idx[i]
    idx_lhs  <- lhs_idx[i]
    # 1. based on posterior draws
    #base_beta_draw <- posterior$base_beta[i, ]
    #I0_draw        <- posterior$I0[i, ]
    #gamma_draw     <- posterior$gamma[i]
    #rho_draw       <- posterior$rho[i]
    #sigma_draw     <- posterior$sigma[i]
    #FOI_rand_draw  <- foi_draws[i]
    
    # 2. based on lhs samples
    #base_beta_draw <- lhs_sample$base_beta[i]
    #I0_draw        <- lhs_sample$I0[i]
    #gamma_draw     <- lhs_sample$gamma[i]
    #rho_draw       <- lhs_sample$rho[i]
    #sigma_draw     <- lhs_sample$sigma[i]
    #FOI_rand_draw  <- lhs_sample$foi[i]
    
    # 3. based on pre-defined id
    base_beta_draw <- posterior$base_beta[idx_post, ]
    I0_draw        <- posterior$I0[idx_post, ]
    gamma_draw     <- lhs_sample$gamma[idx_lhs]
    rho_draw       <- lhs_sample$rho[idx_lhs]
    sigma_draw     <- lhs_sample$sigma[idx_lhs]
    FOI_rand_draw  <- lhs_sample$foi[idx_lhs]
    
    #R0_vec <- (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region] * age_groups))
    R0_vec <- (1 - exp(-FOI_rand_draw * age_groups))
    
    sim_out <- sirv_sim_coverageSwitch(
      T = T,
      A = 20,
      N = N,
      r = rep(0, 20),
      base_beta = base_beta_draw,
      I0_draw = I0_draw,
      R0 = R0_vec,
      rho = rho_draw,
      gamma = gamma_draw,
      sigma = sigma_draw, 
      delay = delay,
      VE_block = 0,
      VE_inf = 0,
      #coverage_threshold = 1,
      target_age = target_age,
      #total_coverage = 1e6 * 0.2 / sum(N),
      total_coverage = 0,
      #total_coverage = 180000 / sum(N), 
      #cap_supply =  1e6,
      weekly_delivery_speed = 0
    )
    
    # Save the age-stratified cases matrix (A rows x T columns) for this draw
    sim_results_list_rawinf[[i]] <- sim_out$age_stratified_cases_raw
    sim_results_list_rawsymp[[i]] <- sim_out$true_symptomatic
  }
  
  # Combine all draws into a 3D array with dimensions: [A, T, n_draws]
  age_array_rawinf <- array(unlist(sim_results_list_rawinf), dim = c(A, T, n_draws))
  age_array_rawsymp <- array(unlist(sim_results_list_rawsymp), dim = c(A, T, n_draws))
  
  # Initialize matrices to store the median and 95% UI for each age group and week
  median_by_age_rawinf <- matrix(0, nrow = A, ncol = T)
  low95_by_age_rawinf  <- matrix(0, nrow = A, ncol = T)
  hi95_by_age_rawinf   <- matrix(0, nrow = A, ncol = T)
  
  median_by_age_rawsymp <- matrix(0, nrow = A, ncol = T)
  low95_by_age_rawsymp  <- matrix(0, nrow = A, ncol = T)
  hi95_by_age_rawsymp   <- matrix(0, nrow = A, ncol = T)
  
  # For each age group and week, compute the median, low 95% and high 95% values across draws
  for (a in 1:A) {
    for (t in 1:T) {
      draws_rawinf <- age_array_rawinf[a, t, ]
      median_by_age_rawinf[a, t] <- median(draws_rawinf)
      low95_by_age_rawinf[a, t]  <- quantile(draws_rawinf, probs = 0.025)
      hi95_by_age_rawinf[a, t]   <- quantile(draws_rawinf, probs = 0.975)
      
      draws_rawsymp <- age_array_rawsymp[a, t, ]
      median_by_age_rawsymp[a, t] <- median(draws_rawsymp)
      low95_by_age_rawsymp[a, t]  <- quantile(draws_rawsymp, probs = 0.025)
      hi95_by_age_rawsymp[a, t]   <- quantile(draws_rawsymp, probs = 0.975)
      
    }
  }
  
  # Compute weekly totals by summing over the age groups for each simulation draw.
  # Summing over the first dimension (age groups) using margin = c(2, 3) returns a matrix of [T x n_draws].
  weekly_totals_rawinf  <- apply(age_array_rawinf, c(2, 3), sum)
  weekly_totals_rawsymp <- apply(age_array_rawsymp, c(2, 3), sum)
  
  # Now compute the median and 95% UI for the weekly totals (aggregated across ages) for each week
  weekly_cases_median_rawinf <- apply(weekly_totals_rawinf, 1, median)
  weekly_cases_low95_rawinf  <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.025)
  weekly_cases_hi95_rawinf   <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.975)
  
  weekly_cases_median_rawsymp <- apply(weekly_totals_rawsymp, 1, median)
  weekly_cases_low95_rawsymp  <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.025)
  weekly_cases_hi95_rawsymp   <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.975)
  
  dimnames(median_by_age_rawinf)  <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  dimnames(low95_by_age_rawinf)   <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  dimnames(hi95_by_age_rawinf)    <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  dimnames(median_by_age_rawsymp) <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  dimnames(low95_by_age_rawsymp)  <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  dimnames(hi95_by_age_rawsymp)   <- list(AgeGroup = as.character(seq_len(A)), Week = as.character(seq_len(T)))
  
  # Create a data frame summarizing the weekly totals with the UI
  df_rawinf <- data.frame(
    week   = 1:T,
    median = weekly_cases_median_rawinf,
    low95  = weekly_cases_low95_rawinf,
    hi95   = weekly_cases_hi95_rawinf
  )
  
  df_rawsymp <- data.frame(
    week   = 1:T,
    median = weekly_cases_median_rawsymp,
    low95  = weekly_cases_low95_rawsymp,
    hi95   = weekly_cases_hi95_rawsymp
  )
  
  # Return all desired outputs
  return(list(
    sim_out             = sim_out,
    sim_results_list_rawinf    = sim_results_list_rawinf,    # List of raw age-stratified matrices (each: A x T)
    sim_results_list_rawsymp   = sim_results_list_rawsymp,
    age_array_rawinf           = age_array_rawinf,           # 3D array of draws [A, T, n_draws]
    age_array_rawsymp          = age_array_rawsymp,
    median_by_age_rawinf       = median_by_age_rawinf,       # Age-stratified median estimates (A x T)
    median_by_age_rawsymp      = median_by_age_rawsymp,
    low95_by_age_rawinf        = low95_by_age_rawinf,        # Age-stratified lower 95% UI (A x T)
    low95_by_age_rawsymp       = low95_by_age_rawsymp,
    hi95_by_age_rawinf         = hi95_by_age_rawinf,         # Age-stratified upper 95% UI (A x T)
    hi95_by_age_rawsymp        = hi95_by_age_rawsymp,
    weekly_totals_rawinf       = weekly_totals_rawinf,       # Matrix of weekly totals for each draw [T x n_draws]
    weekly_cases_median_rawinf = weekly_cases_median_rawinf, # Weekly median across draws
    weekly_cases_low95_rawinf  = weekly_cases_low95_rawinf,  # Weekly lower 95% UI across draws
    weekly_cases_hi95_rawinf   = weekly_cases_hi95_rawinf,   # Weekly upper 95% UI across draws
    weekly_totals_rawsymp      = weekly_totals_rawsymp,
    weekly_cases_median_rawsymp = weekly_cases_median_rawsymp,
    weekly_cases_low95_rawsymp  = weekly_cases_low95_rawsymp,
    weekly_cases_hi95_rawsymp   = weekly_cases_hi95_rawsymp,
    df_rawinf                   = df_rawinf,                   # Data frame summarizing weekly totals with UI
    df_rawsymp                  = df_rawsymp
  ))
}

## infection added (final pre-vacc summarise function)--------------------------------------------------------------
summarise_presim_ui <- function(
    sim_result,        # output from simulate_pre_ui_age
    observed,
    age_gr_levels = c(
      "<1", "1-4", "5–9", "10-11", "12–17", "18–19", "20–24", "25–29",
      "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
      "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
    ),
    lhs_sample_young,  # data for younger ages (for YLD calculation)
    lhs_old,           # data for older ages (for YLD calculation)
    le_sample,         # life-expectancy sample
    hosp,              # hospitalization rate vector or single numeric
    fatal,             # fatality rate for hospitalized
    nh_fatal,           # fatality rate for non-hospitalized
    region
) {
  
  age_gr_levels <- gsub("[–—]", "-", age_gr_levels)
  
  # Set T dynamically
  T =  nrow(observed)
  A <- length(age_gr_levels)
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12–17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  # Create an age group vector for the age-by-week summary:
  age_gr <- rep(c(
    "<1", "1-4", "5–9", "10-11", "12–17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  ), 52)
  
  default_age_vector <- gsub("–", "-", default_age_vector)
  age_gr <- gsub("–", "-", age_gr)
  
  # --- Create Age-Stratified Summary Data Frame ---
  median_df <- as.data.frame.table(sim_result$median_by_age_rawsymp, responseName = "Median")
  colnames(median_df) <- c("AgeGroup", "Week", "Median")
  
  median_df <- median_df %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week = as.numeric(Week)
    )
  
  low95_df <- as.data.frame.table(sim_result$low95_by_age_rawsymp, responseName = "low95")
  hi95_df  <- as.data.frame.table(sim_result$hi95_by_age_rawsymp, responseName = "hi95")
  
  infection_df <- as.data.frame.table(sim_result$median_by_age_rawinf, responseName = "Infections")
  colnames(infection_df) <- c("AgeGroup", "Week", "Infections")
  
  infection_df <- infection_df %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week     = as.numeric(Week)
    )
  
  infection_lo_df <- as.data.frame.table(sim_result$low95_by_age_rawinf, responseName = "Infections_lo")
  infection_hi_df <- as.data.frame.table(sim_result$hi95_by_age_rawinf, responseName = "Infections_hi")
  
  infection_df <- infection_df %>%
    mutate(
      Infections_lo = infection_lo_df$Infections_lo,
      Infections_hi = infection_hi_df$Infections_hi
    )
  
  summary_cases_pre <- median_df %>%
    mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]  # assign provided age group labels
    )%>% left_join(infection_df, by = c("AgeGroup", "Week"))
  
  summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)
  
  # --- Weekly Summary ---
  weekly_df <- data.frame(
    Week = 1:length(sim_result$weekly_cases_median_rawsymp),
    weekly_median = sim_result$weekly_cases_median_rawsymp,
    weekly_low95  = sim_result$weekly_cases_low95_rawsymp,
    weekly_hi95   = sim_result$weekly_cases_hi95_rawsymp
  )
  
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(Median = sum(Median)) %>%
    mutate(Scenario = "Pre-vaccination") %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      lo95  = weekly_low95,
      hi95   = weekly_hi95
    )
  
  
  # --- Calculate Hospitalizations, Fatalities, and Cumulative Totals ---
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      hosp_rate = rep(hosp, T),
      hospitalised = Median * hosp_rate,
      hospitalised_lo = low95 * hosp_rate,
      hospitalised_hi = hi95 * hosp_rate,
      non_hospitalised = Median - hospitalised,
      non_hospitalised_lo = low95 - hospitalised_lo,  
      non_hospitalised_hi = hi95 - hospitalised_hi,
      fatality = rep(fatal, T),          # Ensure fatal is correctly replicated
      nh_fatality = rep(nh_fatal, T),
      fatal = (hospitalised * fatality +   
                 non_hospitalised * nh_fatality),
      fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
      fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality)
    )%>%
    arrange(Week, AgeGroup) %>%          # Order data by Week first, then AgeGroup
    group_by(Week, AgeGroup) %>%         # Group by Week and AgeGroup
    mutate(
      cum_fatal = cumsum(fatal),          # Calculate cumulative fatal cases
      cum_hosp  = cumsum(hospitalised)
    ) %>%
    ungroup()
  
  # --- Summarize by AgeGroup ---
  summary_cases_pre_age <- summary_cases_pre %>% group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Median),
      Infections   = sum(Infections),
      Infections_lo = sum(Infections_lo),
      Infections_hi  = sum(Infections_hi),
      hospitalised = sum(hospitalised),
      hospitalised_lo = sum(hospitalised_lo),
      hospitalised_hi = sum(hospitalised_hi),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    )
  
  summary_cases_pre_age$age_gr <- default_age_vector
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # --- Summarize by Week for Entire Population ---
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      Infections   = sum(Infections),
      Infections_lo = sum(Infections_lo),
      Infections_hi  = sum(Infections_hi),
      hospitalised = sum(hospitalised),
      fatality     = first(fatality),
      fatal        = sum(fatal),
      fatal_lo     = sum(fatal_lo),
      fatal_hi     = sum(fatal_hi)
    ) %>%
    mutate(
      Scenario     = "Pre-vaccination",
      cum_fatal    = cumsum(fatal),
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      cum_hosp     = cumsum(hospitalised)
    ) %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      pre_weekly_low95  = weekly_low95,
      pre_weekly_hi95   = weekly_hi95
    )
  
  # --- DALY Calculations ---
  summary_cases_pre <- summary_cases_pre %>% mutate(
    # YLD parameters
    dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
    dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
    dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
    dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
    dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
    dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
    dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5),
    subac_prop   = quantile(lhs_sample_young$subac, 0.5),
    chr_6m       = quantile(lhs_sample_young$chr6m, 0.5),
    chr_12m      = quantile(lhs_sample_young$chr12m, 0.5),
    chr_30m      = quantile(lhs_sample_young$chr30m, 0.5),
    chr_prop     = chr_6m + chr_12m + chr_30m 
  )
  summary_cases_pre <- summary_cases_pre %>%
    mutate(age_gr = gsub("[–—]", "-", age_gr)) %>%
    mutate(
      age_numeric = case_when(
        age_gr == "<1"    ~ 0,
        age_gr == "1-4"   ~ 2.5,
        age_gr == "5-9"   ~ 7,       
        age_gr == "10-11" ~ 10.5,
        age_gr == "12-17" ~ 14.5,     
        age_gr == "18-19" ~ 18.5,     
        age_gr == "20-24" ~ 22,
        age_gr == "25-29" ~ 27,
        age_gr == "30-34" ~ 32,
        age_gr == "35-39" ~ 37,
        age_gr == "40-44" ~ 42,
        age_gr == "45-49" ~ 47,
        age_gr == "50-54" ~ 52,
        age_gr == "55-59" ~ 57,
        age_gr == "60-64" ~ 62,
        age_gr == "65-69" ~ 67,
        age_gr == "70-74" ~ 72,
        age_gr == "75-79" ~ 77,
        age_gr == "80-84" ~ 82,
        age_gr == "85+"   ~ 87.5,
        TRUE              ~ NA_real_
      ),  # Extract numbers from age_gr
      subac_prop  = case_when(
        age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),   # Replace with the desired value for <39
        age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)   # Replace with the desired value for ≥40
      ),
      chr_6m  = case_when(
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
      chr_prop = chr_6m + chr_12m + chr_30m  # Recalculate chr_prop based on updated values
    ) %>% mutate(
      # YLD estimates
      yld_acute    = (hospitalised * dw_hosp * dur_acute) +
        (non_hospitalised * dw_nonhosp * dur_acute),
      yld_acute_lo = (hospitalised_lo * dw_hosp * dur_acute) +
        ((low95 - hospitalised_lo) * dw_nonhosp * dur_acute),
      yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
        ((hi95 - hospitalised_hi) * dw_nonhosp * dur_acute),
      
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
        ((low95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
        ((hi95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
      
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
        ((low95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
        ((hi95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
      
      yld_total    = yld_acute + yld_subacute + yld_chronic,
      yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
      yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
    ) %>% mutate(
      # le left 
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
    ) %>% mutate(
      # YLL  
      yll = fatal * le_left,
      yll_lo = fatal_lo * le_left,
      yll_hi = fatal_hi * le_left,
      # total DALY
      daly_tot = yld_total + yll,
      daly_tot_lo = yld_total_lo + yll_lo,
      daly_tot_hi = yld_total_hi + yll_hi,
      cum_daly = cumsum(daly_tot)
    )
  
  # --- Final Summaries ---
  summary_cases_pre_age <- summary_cases_pre %>% 
    group_by(AgeGroup) %>%
    summarise(
      Median           = sum(Median, na.rm = TRUE),
      Infections       = sum(Infections),
      Infections_lo    = sum(Infections_lo),
      Infections_hi    = sum(Infections_hi),
      low95            = sum(low95, na.rm = TRUE),
      hi95             = sum(hi95, na.rm = TRUE),
      hospitalised     = sum(hospitalised, na.rm = TRUE),
      hospitalised_lo  = sum(hospitalised_lo, na.rm = TRUE),
      hospitalised_hi  = sum(hospitalised_hi, na.rm = TRUE),
      yld_acute        = sum(yld_acute, na.rm = TRUE),
      yld_acute_lo     = sum(yld_acute_lo, na.rm = TRUE),
      yld_acute_hi     = sum(yld_acute_hi, na.rm = TRUE),
      yld_subacute     = sum(yld_subacute, na.rm = TRUE),
      yld_subacute_lo  = sum(yld_subacute_lo, na.rm = TRUE),
      yld_subacute_hi  = sum(yld_subacute_hi, na.rm = TRUE),
      yld_chronic      = sum(yld_chronic, na.rm = TRUE),
      yld_chronic_lo   = sum(yld_chronic_lo, na.rm = TRUE),
      yld_chronic_hi   = sum(yld_chronic_hi, na.rm = TRUE),
      yld_total        = sum(yld_total, na.rm = TRUE),
      yld_total_lo     = sum(yld_total_lo, na.rm = TRUE),
      yld_total_hi     = sum(yld_total_hi, na.rm = TRUE),
      yll              = sum(yll, na.rm = TRUE),
      yll_lo           = sum(yll_lo, na.rm = TRUE),
      yll_hi           = sum(yll_hi, na.rm = TRUE),
      daly_tot         = sum(daly_tot, na.rm = TRUE),
      daly_tot_lo      = sum(daly_tot_lo, na.rm = TRUE),
      daly_tot_hi      = sum(daly_tot_hi, na.rm = TRUE),
      fatal            = sum(fatal, na.rm = TRUE),
      fatal_lo         = sum(fatal_lo, na.rm = TRUE),
      fatal_hi         = sum(fatal_hi, na.rm = TRUE)
    ) %>% 
    mutate(
      Scenario = "Pre-vaccination"
    ) %>%
    ungroup()
  
  summary_cases_pre_age$age_gr <- age_gr[1:length(age_gr_levels)]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      Infections   = sum(Infections),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal),
      fatal_lo     = sum(fatal_lo),
      fatal_hi     = sum(fatal_hi)
    ) %>% mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal    = cumsum(fatal),
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      cum_hosp   = cumsum(hospitalised),
      cum_yld_acute = cumsum(yld_acute),
      cum_yld_subacute = cumsum(yld_subacute),
      cum_yld_chronic = cumsum(yld_chronic),
      cum_yld_total   = cumsum(yld_total),
      cum_yll         = cumsum(yll),
      cum_daly_tot    = cumsum(daly_tot)
    ) %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      lo95  = weekly_low95,
      hi95   = weekly_hi95
    )
  
  summary_cases_pre$region     <- region
  summary_cases_pre_all$region <- region
  summary_cases_pre_age$region <- region
  
  return(list(
    summary_cases_pre = summary_cases_pre,
    summary_cases_pre_all = summary_cases_pre_all,
    summary_cases_pre_age = summary_cases_pre_age
  ))
}

## post-vacc simulation---------------------------------------------------------
## v3. age group aligning to licensure
run_simulation_scenarios_ui_ixchiq <- function(target_age_list, 
                                               observed,
                                               N, bra_foi_state_summ, age_groups, region_name,
                                               hosp, fatal, nh_fatal,
                                               lhs_sample_young, lhs_old, le_sample,
                                               age_gr_levels,
                                               prevacc_ui = NULL,
                                               posterior,
                                               ve_inf = 0,
                                               total_coverage = 0,
                                               lhs_sample
) {
  # target_age_list: list of 0/1 vectors for each scenario
  # observed: observed data (used to set T)
  # N: population vector (length = number of age groups)
  # bra_foi_state_summ: data frame with avg_foi and NAME_1 for FOI calculation
  # age_groups: numeric vector of age group midpoints used in FOI calculation
  # hosp, fatal, nh_fatal: hospitalization, fatality rates (numeric or vector)
  # lhs_sample_young, lhs_old: samples for YLD parameters for young/old ages
  # le_sample: life-expectancy sample (for YLL calculation)
  # age_gr_levels: factor levels for age group labels
  # prevacc_ui: (optional) pre-vaccination UI data frame with columns: AgeGroup, Week, Median, low95, hi95
  # posterior: list with posterior draws (e.g., base_beta, I0, gamma, rho, etc.)
  
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  #n_draws <- length(posterior$gamma)
  #n_draws <- nrow(lhs_sample)
  posterior_idx <- posterior_idx_list[[region_name]]
  lhs_idx       <- lhs_idx_list[[region_name]]
  n_draws <- length(posterior_idx)
  
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  default_age_vector <- gsub("–", "-", default_age_vector)
  # Create an age group vector for the age-by-week summary:
  # Here, each default age label is repeated T times.
  age_gr <- rep(default_age_vector, T)
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    
    pct <- region_coverage[[region_name]]
    total_coverage_frac <- floor(pct) / 100
    
    #foi_row <- bra_foi_state_summ[bra_foi_state_summ$NAME_1 == region_name, ]
    #foi_mean <- foi_row$avg_foi
    foi_draws <- foi_draws_list[[region_name]]
    set.seed(123)
    #posterior_idx <- sample(seq_len(nrow(posterior$I0)), size = n_draws, replace = TRUE)
    #I0_resampled        <- posterior$I0[posterior_idx, ]        
    #base_beta_resampled <- posterior$base_beta[posterior_idx, ] 
    
    for (d in 1:n_draws) {
      
      idx_post <- posterior_idx[d] 
      idx_lhs  <- lhs_idx[d]
      # Extract d-th draw parameters
      #base_beta_draw <- posterior$base_beta[d, ]   
      #I0_draw        <- posterior$I0[d, ]          
      #gamma_draw     <- posterior$gamma[d]
      #rho_draw       <- posterior$rho[d]
      #sigma_draw     <- posterior$sigma[d]
      #ve_rand_draw   <- ve_draw[d]
      #wd_rand_draw   <- wd_draw[d]
      #FOI_rand_draw  <- foi_draws[d]
      #vc_rand_draw   <- vc_draw[d]
      #base_beta_draw <- base_beta_resampled[d, ]
      #I0_draw        <- I0_resampled[d, ]
      #gamma_draw     <- lhs_sample$gamma[d]
      #rho_draw       <- lhs_sample$rho[d]
      #sigma_draw     <- lhs_sample$sigma[d]
      #FOI_rand_draw  <- lhs_sample$foi[d]
      #ve_rand_draw <- lhs_sample$ve[d]
      #vc_rand_draw <- lhs_sample$vc[d]
      #wd_rand_draw <- lhs_sample$wd[d]
      
      # 2. predefined idx 
      base_beta_draw <- posterior$base_beta[idx_post, ]
      I0_draw        <- posterior$I0[idx_post, ]
      gamma_draw     <- lhs_sample$gamma[idx_lhs]
      rho_draw       <- lhs_sample$rho[idx_lhs]
      sigma_draw     <- lhs_sample$sigma[idx_lhs]
      FOI_rand_draw  <- lhs_sample$foi[idx_lhs]
      ve_rand_draw   <- lhs_sample$ve_ix[idx_lhs]
      #vc_rand_draw   <- lhs_sample$vc[idx_lhs]
      wd_rand_draw   <- lhs_sample$wd[idx_lhs]
      
      # sampling from lhs by defined VC levels 
      vc_rand_draw <- case_when(
        total_coverage == 0.10 ~ lhs_sample$vc10[idx_lhs],
        total_coverage == 0.50 ~ lhs_sample$vc50[idx_lhs],
        total_coverage == 0.90 ~ lhs_sample$vc90[idx_lhs],
        TRUE ~ NA_real_
      )
      
      # for ve x vc simulation
      VE_block_draw <- ve_rand_draw   ## always block disease                    
      
      VE_inf_draw <- if (is.na(ve_inf)) {                 
        ve_rand_draw                                      
      } else if (ve_inf == 0) {                           
        0                                                 
      } else {                                          
        ve_inf                                         
      }
      
      sim_out <- sirv_sim_coverageSwitch_ixchiq(
        T = T,
        A = length(age_gr_levels),
        N = N,
        r = rep(0, length(age_gr_levels)),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        #R0 = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
        R0 = 1 - exp(- FOI_rand_draw * age_groups),
        rho = rho_draw,
        gamma = gamma_draw,
        sigma = sigma_draw,
        delay = 2,
        VE_block = VE_block_draw,
        VE_inf = VE_inf_draw,
        #coverage_threshold = 1,
        target_age = target,
        #total_coverage = 2e6 * 0.2 / sum(N),
        #total_coverage = total_coverage_frac,
        total_coverage = vc_rand_draw,
        #total_coverage = 0.5,
        weekly_delivery_speed = wd_rand_draw
      )
      
      draw_results_raw_inf[[d]]  <- sim_out$age_stratified_cases_raw  # A x T matrix
      draw_results_raw_symp[[d]] <- sim_out$true_symptomatic  # A x T matrix
      
      if (anyNA(sim_out$age_stratified_cases_raw))
        stop("Draw ", d, ": sim_out produced NA")
    }
    
    # Combine draw results into a 3D array [18, T, n_draws]
    age_array_raw_inf  <- array(unlist(draw_results_raw_inf), dim = c(20, T, n_draws))
    age_array_raw_symp <- array(unlist(draw_results_raw_symp), dim = c(20, T, n_draws))
    
    # Compute quantiles for each age and week over draws:
    median_by_age_rawinf <- apply(age_array_raw_inf, c(1, 2), median)
    low95_by_age_rawinf  <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawinf   <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.975)
    
    median_by_age_rawsymp <- apply(age_array_raw_symp, c(1, 2), median)
    low95_by_age_rawsymp  <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawsymp   <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.975)
    
    # Weekly totals (sum over ages) for each draw:
    weekly_totals_rawinf       <- apply(age_array_raw_inf, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawinf <- apply(weekly_totals_rawinf, 1, median)
    weekly_cases_low95_rawinf  <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawinf   <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.975)
    
    weekly_totals_rawsymp       <- apply(age_array_raw_symp, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawsymp <- apply(weekly_totals_rawsymp, 1, median)
    weekly_cases_low95_rawsymp  <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawsymp   <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.975)
    
    # infection df
    infection_median_df <- as.data.frame.table(median_by_age_rawinf, responseName = "infection")
    infection_low_df    <- as.data.frame.table(low95_by_age_rawinf,  responseName = "infection_lo")
    infection_hi_df     <- as.data.frame.table(hi95_by_age_rawinf,   responseName = "infection_hi")
    
    colnames(infection_median_df) <- c("AgeGroup", "Week", "infection")
    colnames(infection_low_df)    <- c("AgeGroup", "Week", "infection_lo")
    colnames(infection_hi_df)     <- c("AgeGroup", "Week", "infection_hi")
    
    infection_df <- infection_median_df %>%
      mutate(
        infection_lo = infection_low_df$infection_lo,
        infection_hi = infection_hi_df$infection_hi,
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    # Create simulation data frame from median_by_age: (raw symptomatic cases)
    sim_df <- as.data.frame.table(median_by_age_rawsymp, responseName = "Cases")
    
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      mutate(
        Scenario = s,
        # Convert AgeGroup from factor to numeric via as.character
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    sim_df <- sim_df %>%
      mutate(region_full = region_name)
    
    # left join inf 
    sim_df <- left_join(sim_df, infection_df, by = c("AgeGroup", "Week"))
    
    # Attach UI for cases:
    low95_df <- as.data.frame.table(low95_by_age_rawsymp, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age_rawsymp, responseName = "post_hi95")
    
    sim_df <- sim_df %>%
      mutate(
        post_low95        = low95_df$post_low95,
        post_hi95         = hi95_df$post_hi95
      )
    
    # Weekly summary data frame:
    weekly_df <- data.frame(
      Week = 1:T,
      weekly_median = weekly_cases_median_rawsymp,
      weekly_low95 = weekly_cases_low95_rawsymp,
      weekly_hi95 = weekly_cases_hi95_rawsymp
    )
    
    # Attach age group labels using our default vector:
    sim_df <- sim_df %>%
      mutate(
        age_gr = rep(default_age_vector, T)
      )
    sim_df$age_gr <- factor(sim_df$age_gr, levels = age_gr_levels)
    
    # Calculate hospitalisations, fatalities, and uncertainty bounds:
    sim_df <- sim_df %>%
      mutate(
        hosp_rate = rep(hosp, T),
        hospitalised = Cases * hosp_rate,
        hospitalised_lo = post_low95 * hosp_rate,
        hospitalised_hi = post_hi95 * hosp_rate,
        non_hospitalised = Cases - hospitalised,
        non_hospitalised_lo = post_low95 - hospitalised_lo,   # conservative
        non_hospitalised_hi = post_hi95 - hospitalised_hi,
        fatality = rep(fatal, T),
        nh_fatality = rep(nh_fatal, T),
        fatal = (hospitalised * fatality + non_hospitalised * nh_fatality),
        fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
        fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality),
        
      ) %>%
      arrange(Week, AgeGroup) %>%
      group_by(Week, AgeGroup) %>%
      mutate(
        cum_fatal = cumsum(fatal),
        cum_hosp = cumsum(hospitalised)
      ) %>%
      ungroup()
    
    # Calculate age_numeric from age_gr for age-dependent parameters.
    sim_df <- sim_df %>%
      mutate(
        age_numeric =case_when(
          age_gr == "<1"    ~ 0,
          age_gr == "1-4"   ~ 2.5,
          age_gr == "5-9"   ~ 7,        
          age_gr == "10-11" ~ 10.5,
          age_gr == "12-17" ~ 14.5,     
          age_gr == "18-19" ~ 18.5,    
          age_gr == "20-24" ~ 22,
          age_gr == "25-29" ~ 27,
          age_gr == "30-34" ~ 32,
          age_gr == "35-39" ~ 37,
          age_gr == "40-44" ~ 42,
          age_gr == "45-49" ~ 47,
          age_gr == "50-54" ~ 52,
          age_gr == "55-59" ~ 57,
          age_gr == "60-64" ~ 62,
          age_gr == "65-69" ~ 67,
          age_gr == "70-74" ~ 72,
          age_gr == "75-79" ~ 77,
          age_gr == "80-84" ~ 82,
          age_gr == "85+"   ~ 87.5,
          TRUE              ~ NA_real_
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
          ((post_low95 - hospitalised_lo) * dw_nonhosp * dur_acute),
        yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
          ((post_hi95 - hospitalised_hi) * dw_nonhosp * dur_acute),
        
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
        
        yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
        
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
    
    if(!is.null(prevacc_ui)) {
      sim_df <- left_join(sim_df, 
                          dplyr::select(prevacc_ui, AgeGroup, Week, 
                                        pre_Median = Median, 
                                        pre_low95 = low95, 
                                        pre_hi95 = hi95), 
                          by = c("AgeGroup", "Week"))
    }
    
    scenario_result[[s]] <- list(sim_result = list(age_array_raw_symp           = age_array_raw_symp,
                                                   age_array_raw_inf            = age_array_raw_inf,
                                                   weekly_cases_median_rawsymp  = weekly_cases_median_rawsymp,
                                                   weekly_cases_median_rawinf   = weekly_cases_median_rawinf,
                                                   weekly_cases_low95_rawsymp   = weekly_cases_low95_rawsymp,
                                                   weekly_cases_low95_rawinf    = weekly_cases_low95_rawinf,
                                                   weekly_cases_hi95_rawsymp    = weekly_cases_hi95_rawsymp),
                                 sim_out                      = sim_out,
                                 sim_df                       = sim_df,
                                 weekly_df                    = weekly_df
    )
  }
  
  return(scenario_result)
}

### vimkunya
run_simulation_scenarios_ui_vimkun <- function(target_age_list, 
                                               observed,
                                               N, bra_foi_state_summ, age_groups, region_name,
                                               hosp, fatal, nh_fatal,
                                               lhs_sample_young, lhs_old, le_sample,
                                               age_gr_levels,
                                               prevacc_ui = NULL,
                                               posterior,
                                               ve_inf = 0,
                                               total_coverage = 0,
                                               lhs_sample
) {
  # target_age_list: list of 0/1 vectors for each scenario
  # observed: observed data (used to set T)
  # N: population vector (length = number of age groups)
  # bra_foi_state_summ: data frame with avg_foi and NAME_1 for FOI calculation
  # age_groups: numeric vector of age group midpoints used in FOI calculation
  # hosp, fatal, nh_fatal: hospitalization, fatality rates (numeric or vector)
  # lhs_sample_young, lhs_old: samples for YLD parameters for young/old ages
  # le_sample: life-expectancy sample (for YLL calculation)
  # age_gr_levels: factor levels for age group labels
  # prevacc_ui: (optional) pre-vaccination UI data frame with columns: AgeGroup, Week, Median, low95, hi95
  # posterior: list with posterior draws (e.g., base_beta, I0, gamma, rho, etc.)
  
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  #n_draws <- length(posterior$gamma)
  posterior_idx <- posterior_idx_list[[region_name]]
  lhs_idx       <- lhs_idx_list[[region_name]]
  n_draws <- length(posterior_idx)
  
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  default_age_vector <- gsub("–", "-", default_age_vector)
  
  # Create an age group vector for the age-by-week summary:
  # Here, each default age label is repeated T times.
  age_gr <- rep(default_age_vector, T)
  
  postwane_vec <- c(0.948, 0.948, 0.85, 0.841) 
  # test
  #postwane_vec <- c(0.1, 0.1, 0.1, 0.1) 
  VE_time_list <- list()
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    this_postwane <- postwane_vec[s]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    VE_time_draws <- matrix(NA, nrow = n_draws, ncol = 52) ## VE by target age group by draw
    
    pct <- region_coverage[[region_name]]
    total_coverage_frac <- floor(pct) / 100
    
    for (d in 1:n_draws) {
      
      #idx <- posterior_idx[d] 
      idx_post <- posterior_idx[d] 
      idx_lhs  <- lhs_idx[d]
      # Extract d-th draw parameters ------ calibration based 
      #base_beta_draw <- posterior$base_beta[d, ]   
      #I0_draw        <- posterior$I0[d, ]          
      #gamma_draw     <- posterior$gamma[d]
      #rho_draw       <- posterior$rho[d]
      #sigma_draw     <- posterior$sigma[d]
      
      # 2. predefined idx 
      base_beta_draw <- posterior$base_beta[idx_post, ]
      I0_draw        <- posterior$I0[idx_post, ]
      gamma_draw     <- lhs_sample$gamma[idx_lhs]
      rho_draw       <- lhs_sample$rho[idx_lhs]
      sigma_draw     <- lhs_sample$sigma[idx_lhs]
      FOI_rand_draw  <- lhs_sample$foi[idx_lhs]
      ve_rand_draw   <- lhs_sample$ve_vimkun[idx_lhs]
      vc_rand_draw   <- lhs_sample$vc[idx_lhs]
      wd_rand_draw   <- lhs_sample$wd[idx_lhs]
      
      # sampling from lhs by defined VC levels 
      #vc_rand_draw <- case_when(
      #  total_coverage == 0.10 ~ lhs_sample$vc10[idx_lhs],
      #  total_coverage == 0.50 ~ lhs_sample$vc50[idx_lhs],
      #  total_coverage == 0.90 ~ lhs_sample$vc90[idx_lhs],
      #  TRUE ~ NA_real_
      #)
      
      # waning curve helper function
      #ve_waning_curve <- function(T = 52, ve_week3, ve_week26) {
      #  slope <- (ve_week26 - ve_week3) / (26 - 3)
      
      #  VE_time <- sapply(1:T, function(t) {
      #    if (t <= 3) {
      #      ve_week3
      #    } else if (t <= 26) {
      #      ve_week3 + slope * (t - 3)
      #    } else {
      #      ve_week26 + slope * (t - 26)
      #    }
      #  })
      
      #  return(VE_time)
      #}
      
      ## exponential decay
      ve_waning_curve <- function(T = 52, ve_week3, ve_week26) {
        
        lambda <- log(ve_week3 / ve_week26) / (26 - 3)
        
        VE_time <- numeric(T)
        
        for (t in 1:T) {
          if (t <= 3) {
            VE_time[t] <- ve_week3
          } else {
            VE_time[t] <- ve_week3 * exp(-lambda * (t - 3))
          }
        }
        
        return(VE_time)
      }
      
      
      VE_time_base  <- ve_waning_curve(T         = 52, 
                                       ve_week3  = ve_rand_draw,
                                       ve_week26 = this_postwane)
      
      # 2) block always waning curve
      VE_block_vec <- VE_time_base
      
      # ve inf 
      VE_inf_vec <- case_when(
        is.na(ve_inf)             ~ VE_time_base,   
        ve_inf == 0               ~ rep(0, 52),     
        TRUE                      ~ rep(ve_inf, 52) 
      )
      
      VE_time_draws[d, ] <- VE_time_base
      
      sim_out <- sirv_sim_coverageSwitch_vimkun(
        T = T,
        A = length(age_gr_levels),
        N = N,
        r = rep(0, length(age_gr_levels)),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        #R0 = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
        R0 = 1 - exp(- FOI_rand_draw * age_groups),
        rho = rho_draw,
        gamma = gamma_draw,
        sigma = sigma_draw,
        delay = 2,
        VE_block = VE_time_base,
        VE_inf = VE_time_base,
        VE_inf_postwane   = this_postwane, 
        #coverage_threshold = 1,
        target_age = target,
        #total_coverage = 2e6 * 0.2 / sum(N),
        #total_coverage = total_coverage_frac,
        total_coverage = vc_rand_draw,
        #total_coverage = 0.01,
        weekly_delivery_speed = wd_rand_draw
      )
      
      draw_results_raw_inf[[d]]  <- sim_out$age_stratified_cases_raw  # A x T matrix
      draw_results_raw_symp[[d]] <- sim_out$true_symptomatic  # A x T matrix
      
    }
    
    VE_time_list[[s]] <- VE_time_draws
    
    # Combine draw results into a 3D array [18, T, n_draws]
    age_array_raw_inf  <- array(unlist(draw_results_raw_inf), dim = c(20, T, n_draws))
    age_array_raw_symp <- array(unlist(draw_results_raw_symp), dim = c(20, T, n_draws))
    
    # Compute quantiles for each age and week over draws:
    median_by_age_rawinf <- apply(age_array_raw_inf, c(1, 2), median)
    low95_by_age_rawinf  <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawinf   <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.975)
    
    median_by_age_rawsymp <- apply(age_array_raw_symp, c(1, 2), median)
    low95_by_age_rawsymp  <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawsymp   <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.975)
    
    # Weekly totals (sum over ages) for each draw:
    weekly_totals_rawinf       <- apply(age_array_raw_inf, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawinf <- apply(weekly_totals_rawinf, 1, median)
    weekly_cases_low95_rawinf  <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawinf   <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.975)
    
    weekly_totals_rawsymp       <- apply(age_array_raw_symp, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawsymp <- apply(weekly_totals_rawsymp, 1, median)
    weekly_cases_low95_rawsymp  <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawsymp   <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.975)
    
    
    # infection df
    infection_median_df <- as.data.frame.table(median_by_age_rawinf, responseName = "infection")
    infection_low_df    <- as.data.frame.table(low95_by_age_rawinf,  responseName = "infection_lo")
    infection_hi_df     <- as.data.frame.table(hi95_by_age_rawinf,   responseName = "infection_hi")
    
    colnames(infection_median_df) <- c("AgeGroup", "Week", "infection")
    colnames(infection_low_df)    <- c("AgeGroup", "Week", "infection_lo")
    colnames(infection_hi_df)     <- c("AgeGroup", "Week", "infection_hi")
    
    infection_df <- infection_median_df %>%
      mutate(
        infection_lo = infection_low_df$infection_lo,
        infection_hi = infection_hi_df$infection_hi,
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    # Create simulation data frame from median_by_age: (raw symptomatic cases)
    sim_df <- as.data.frame.table(median_by_age_rawsymp, responseName = "Cases")
    
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      mutate(
        Scenario = s,
        # Convert AgeGroup from factor to numeric via as.character
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    sim_df <- sim_df %>%
      mutate(region_full = region_name)
    
    # left join inf 
    sim_df <- left_join(sim_df, infection_df, by = c("AgeGroup", "Week"))
    
    # Attach UI for cases:
    low95_df <- as.data.frame.table(low95_by_age_rawsymp, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age_rawsymp, responseName = "post_hi95")
    
    sim_df <- sim_df %>%
      mutate(
        post_low95        = low95_df$post_low95,
        post_hi95         = hi95_df$post_hi95
      )
    
    # Weekly summary data frame:
    weekly_df <- data.frame(
      Week = 1:T,
      weekly_median = weekly_cases_median_rawsymp,
      weekly_low95 = weekly_cases_low95_rawsymp,
      weekly_hi95 = weekly_cases_hi95_rawsymp
    )
    
    # Attach age group labels using our default vector:
    sim_df <- sim_df %>%
      mutate(
        age_gr = rep(default_age_vector, T)
      )
    sim_df$age_gr <- factor(sim_df$age_gr, levels = age_gr_levels)
    
    # Calculate hospitalisations, fatalities, and uncertainty bounds:
    sim_df <- sim_df %>%
      mutate(
        hosp_rate = rep(hosp, T),
        hospitalised = Cases * hosp_rate,
        hospitalised_lo = post_low95 * hosp_rate,
        hospitalised_hi = post_hi95 * hosp_rate,
        non_hospitalised = Cases - hospitalised,
        non_hospitalised_lo = post_low95 - hospitalised_lo,   
        non_hospitalised_hi = post_hi95 - hospitalised_hi,
        fatality = rep(fatal, T),
        nh_fatality = rep(nh_fatal, T),
        fatal = (hospitalised * fatality + non_hospitalised * nh_fatality),
        fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
        fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality),
        
      ) %>%
      arrange(Week, AgeGroup) %>%
      group_by(Week, AgeGroup) %>%
      mutate(
        cum_fatal = cumsum(fatal),
        cum_hosp = cumsum(hospitalised)
      ) %>%
      ungroup()
    
    # Calculate age_numeric from age_gr for age-dependent parameters.
    sim_df <- sim_df %>%
      mutate(
        age_numeric = case_when(
          age_gr == "<1"    ~ 0,
          age_gr == "1-4"   ~ 2.5,
          age_gr == "5-9"   ~ 7,        
          age_gr == "10-11" ~ 10.5,
          age_gr == "12-17" ~ 14.5,     
          age_gr == "18-19" ~ 18.5,    
          age_gr == "20-24" ~ 22,
          age_gr == "25-29" ~ 27,
          age_gr == "30-34" ~ 32,
          age_gr == "35-39" ~ 37,
          age_gr == "40-44" ~ 42,
          age_gr == "45-49" ~ 47,
          age_gr == "50-54" ~ 52,
          age_gr == "55-59" ~ 57,
          age_gr == "60-64" ~ 62,
          age_gr == "65-69" ~ 67,
          age_gr == "70-74" ~ 72,
          age_gr == "75-79" ~ 77,
          age_gr == "80-84" ~ 82,
          age_gr == "85+"   ~ 87.5,
          TRUE              ~ NA_real_
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
          ((post_low95 - hospitalised_lo) * dw_nonhosp * dur_acute),
        yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
          ((post_hi95 - hospitalised_hi) * dw_nonhosp * dur_acute),
        
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
        
        yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
        
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
    
    if(!is.null(prevacc_ui)) {
      sim_df <- left_join(sim_df, 
                          dplyr::select(prevacc_ui, AgeGroup, Week, 
                                        pre_Median = Median, 
                                        pre_low95 = low95, 
                                        pre_hi95 = hi95), 
                          by = c("AgeGroup", "Week"))
    }
    
    scenario_result[[s]] <- list(sim_result = list(age_array_raw_symp           = age_array_raw_symp,
                                                   age_array_raw_inf            = age_array_raw_inf,
                                                   weekly_cases_median_rawsymp  = weekly_cases_median_rawsymp,
                                                   weekly_cases_median_rawinf   = weekly_cases_median_rawinf,
                                                   weekly_cases_low95_rawsymp   = weekly_cases_low95_rawsymp,
                                                   weekly_cases_low95_rawinf    = weekly_cases_low95_rawinf,
                                                   weekly_cases_hi95_rawsymp    = weekly_cases_hi95_rawsymp),
                                 sim_out                      = sim_out,
                                 sim_df                       = sim_df,
                                 weekly_df                    = weekly_df,
                                 VE_time_list                 = VE_time_list
    )
  }
  
  return(scenario_result)
}

# v4. age group aligning with licensure
postsim_all_ui <- function(scenario_result,   # list of scenario outputs: each with $sim_df and $weekly_df
                           observed,
                           age_gr_levels,
                           pre_summary_cases_age,   # aggregated pre-vacc by age; should include columns:
                           # Median, low95, hi95, yld_acute, yld_subacute, 
                           # yld_chronic, yld_total, yll, daly_tot, fatal
                           pre_summary_cases,       # pre-vacc data at age-week level (if needed)
                           pre_summary_cases_all,   # aggregated pre-vacc weekly UI data (with Week, Median, low95, hi95)
                           region) { 
  # -------- basic dims -------------------------------------------------
  T <- nrow(observed)
  # Define default age groups
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  default_age_vector <- gsub("–", "-", default_age_vector)
  # Age-week vector
  age_gr <- rep(default_age_vector, T)
  n_scenarios <- length(scenario_result)
  scenario_df  <- lapply(scenario_result, function(x) x$sim_df)
  
  # Build summary_list (age-week)
  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    # Attach aggregated pre-vaccination values (assumed to be at the same resolution)
    df$pre_infection      <- pre_summary_cases$Infections
    df$pre_infection_lo   <- pre_summary_cases$Infections_lo
    df$pre_infection_hi   <- pre_summary_cases$Infections_hi
    df$pre_vacc  <- pre_summary_cases$Median  
    df$pre_fatal <- pre_summary_cases$fatal
    df$pre_daly  <- pre_summary_cases$daly_tot 
    df <- df %>%
      mutate(
        diff   = pre_vacc - Cases,
        impact = diff / pre_vacc * 100
      )
    return(df)
  })
  
  summary_list_df <- do.call(rbind, summary_list)
  summary_list_df$age_gr <- rep(age_gr, length.out = nrow(summary_list_df))
  summary_list_df$age_gr <- factor(summary_list_df$age_gr, levels = age_gr_levels)
  
  # Merge aggregated pre-vacc UI (weekly) from pre_summary_cases_all into the age-week summary.
  if(all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_list_df <- left_join(summary_list_df, pre_weekly_df, by = "Week")
  }
  
  # Final age-group summary (sum → quantile)
  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%  
      summarise(
        infection        = sum(infection, na.rm = TRUE),
        infection_lo     = sum(infection_lo, na.rm = TRUE),
        infection_hi     = sum(infection_hi, na.rm = TRUE),
        Median           = sum(Cases, na.rm = TRUE),
        low95            = sum(post_low95, na.rm = TRUE),
        hi95             = sum(post_hi95, na.rm = TRUE),
        hospitalised     = sum(hospitalised, na.rm = TRUE),
        hospitalised_lo  = sum(hospitalised_lo, na.rm = TRUE),
        hospitalised_hi  = sum(hospitalised_hi, na.rm = TRUE),
        yld_acute        = sum(yld_acute, na.rm = TRUE),
        yld_acute_lo     = sum(yld_acute_lo, na.rm = TRUE),
        yld_acute_hi     = sum(yld_acute_hi, na.rm = TRUE),
        yld_subacute     = sum(yld_subacute, na.rm = TRUE),
        yld_subacute_lo  = sum(yld_subacute_lo, na.rm = TRUE),
        yld_subacute_hi  = sum(yld_subacute_hi, na.rm = TRUE),
        yld_chronic      = sum(yld_chronic, na.rm = TRUE),
        yld_chronic_lo   = sum(yld_chronic_lo, na.rm = TRUE),
        yld_chronic_hi   = sum(yld_chronic_hi, na.rm = TRUE),
        yld_total        = sum(yld_total, na.rm = TRUE),
        yld_total_lo     = sum(yld_total_lo, na.rm = TRUE),
        yld_total_hi     = sum(yld_total_hi, na.rm = TRUE),
        yll              = sum(yll, na.rm = TRUE),
        yll_lo           = sum(yll_lo, na.rm = TRUE),
        yll_hi           = sum(yll_hi, na.rm = TRUE),
        daly_tot         = sum(daly_tot, na.rm = TRUE),
        daly_tot_lo      = sum(daly_tot_lo, na.rm = TRUE),
        daly_tot_hi      = sum(daly_tot_hi, na.rm = TRUE),
        fatal            = sum(fatal, na.rm = TRUE),
        fatal_lo         = sum(fatal_lo, na.rm = TRUE),
        fatal_hi         = sum(fatal_hi, na.rm = TRUE)
      ) %>% 
      mutate(
        # Attach pre-vaccination values from pre_summary_cases_age
        pre_infection         = pre_summary_cases_age$Infections,
        pre_infection_lo      = pre_summary_cases_age$Infections_lo,
        pre_infection_hi      = pre_summary_cases_age$Infections_hi,
        
        pre_vacc              = pre_summary_cases_age$Median,
        pre_vacc_low95        = pre_summary_cases_age$low95,
        pre_vacc_hi95         = pre_summary_cases_age$hi95,
        
        pre_yld_acute         = pre_summary_cases_age$yld_acute,
        pre_yld_acute_low95   = pre_summary_cases_age$yld_acute_lo,
        pre_yld_acute_hi      = pre_summary_cases_age$yld_acute_hi,
        
        pre_yld_subacute      = pre_summary_cases_age$yld_subacute,
        pre_yld_subacute_low95= pre_summary_cases_age$yld_subacute_lo,
        pre_yld_subacute_hi   = pre_summary_cases_age$yld_subacute_hi,
        
        pre_yld_chronic       = pre_summary_cases_age$yld_chronic,
        pre_yld_chronic_low95 = pre_summary_cases_age$yld_chronic_lo,
        pre_yld_chronic_hi    = pre_summary_cases_age$yld_chronic_hi,
        
        pre_yld_total         = pre_summary_cases_age$yld_total,
        pre_yld_total_low95   = pre_summary_cases_age$yld_total_lo,
        pre_yld_total_hi      = pre_summary_cases_age$yld_total_hi,
        
        pre_yll               = pre_summary_cases_age$yll,
        pre_yll_low95         = pre_summary_cases_age$yll_lo,
        pre_yll_hi            = pre_summary_cases_age$yll_hi,
        
        pre_daly              = pre_summary_cases_age$daly_tot,
        pre_daly_low95        = pre_summary_cases_age$daly_tot_lo,
        pre_daly_hi           = pre_summary_cases_age$daly_tot_hi,
        
        pre_fatal             = pre_summary_cases_age$fatal,
        pre_fatal_low95       = pre_summary_cases_age$fatal_lo,
        pre_fatal_hi          = pre_summary_cases_age$fatal_hi,
        
        # Differences for infections:
        diff_inf    = pre_infection - infection,
        diff_inf_lo = pre_infection_lo - infection_lo,
        diff_inf_hi = pre_infection_hi - infection_hi,
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
        
        # Differences for YLD Acute:
        diff_yld_acute       = pre_yld_acute - yld_acute,
        diff_yld_acute_low   = pre_yld_acute_low95 - yld_acute_lo,
        diff_yld_acute_hi    = pre_yld_acute_hi - yld_acute_hi,
        impact_yld_acute     = diff_yld_acute / pre_yld_acute * 100,
        impact_yld_acute_low = diff_yld_acute_low / pre_yld_acute_low95 * 100,
        impact_yld_acute_hi  = diff_yld_acute_hi / pre_yld_acute_hi * 100,
        
        # Differences for YLD Subacute:
        diff_yld_subacute       = pre_yld_subacute - yld_subacute,
        diff_yld_subacute_low   = pre_yld_subacute_low95 - yld_subacute_lo,
        diff_yld_subacute_hi    = pre_yld_subacute_hi - yld_subacute_hi,
        impact_yld_subacute     = diff_yld_subacute / pre_yld_subacute * 100,
        impact_yld_subacute_low = diff_yld_subacute_low / pre_yld_subacute_low95 * 100,
        impact_yld_subacute_hi  = diff_yld_subacute_hi / pre_yld_subacute_hi * 100,
        
        # Differences for YLD Chronic:
        diff_yld_chronic       = pre_yld_chronic - yld_chronic,
        diff_yld_chronic_low   = pre_yld_chronic_low95 - yld_chronic_lo,
        diff_yld_chronic_hi    = pre_yld_chronic_hi - yld_chronic_hi,
        impact_yld_chronic     = diff_yld_chronic / pre_yld_chronic * 100,
        impact_yld_chronic_low = diff_yld_chronic_low / pre_yld_chronic_low95 * 100,
        impact_yld_chronic_hi  = diff_yld_chronic_hi / pre_yld_chronic_hi * 100,
        
        # Differences for YLD Total:
        diff_yld_total       = pre_yld_total - yld_total,
        diff_yld_total_low   = pre_yld_total_low95 - yld_total_lo,
        diff_yld_total_hi    = pre_yld_total_hi - yld_total_hi,
        impact_yld_total     = diff_yld_total / pre_yld_total * 100,
        impact_yld_total_low = diff_yld_total_low / pre_yld_total_low95 * 100,
        impact_yld_total_hi  = diff_yld_total_hi / pre_yld_total_hi * 100,
        
        # Differences for YLL:
        diff_yll       = pre_yll - yll,
        diff_yll_low   = pre_yll_low95 - yll_lo,
        diff_yll_hi    = pre_yll_hi - yll_hi,
        impact_yll     = diff_yll / pre_yll * 100,
        impact_yll_low = diff_yll_low / pre_yll_low95 * 100,
        impact_yll_hi  = diff_yll_hi / pre_yll_hi * 100,
        
        # Differences for DALY:
        diff_daly       = pre_daly - daly_tot,
        diff_daly_low   = pre_daly_low95 - daly_tot_lo,
        diff_daly_hi    = pre_daly_hi - daly_tot_hi,
        impact_daly     = diff_daly / pre_daly * 100,
        impact_daly_low = diff_daly_low / pre_daly_low95 * 100,
        impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100
      )
  })
  
  # Weekly summary
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median,
             pre_vacc_lo = pre_summary_cases$low95,
             pre_vacc_hi = pre_summary_cases$hi95
      ) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_infection = sum(infection),
        pre_infection = sum(pre_infection),
        post_cases = sum(Cases, na.rm = TRUE),
        post_cases_lo = sum(post_low95, na.rm = TRUE),
        post_cases_hi = sum(post_hi95, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        pre_cases_lo = sum(pre_vacc_lo, na.rm = TRUE),
        pre_cases_hi  = sum(pre_vacc_hi, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        post_fatal_lo = sum(fatal_lo, na.rm = TRUE),
        post_fatal_hi = sum(fatal_hi, na.rm = TRUE),
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
        post_daly_lo  = sum(daly_tot_lo, na.rm = TRUE),
        post_daly_hi  = sum(daly_tot_hi, na.rm = TRUE),
        pre_daly   = sum(pre_daly, na.rm = TRUE),
        diff_daly  = pre_daly - post_daly,
        Scenario   = first(Scenario),
        .groups = "drop"
      )
  })
  
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median,
             pre_vacc_lo = pre_summary_cases$low95,
             pre_vacc_hi = pre_summary_cases$hi95,
             pre_fatal   = pre_summary_cases$fatal,
             pre_fatal_lo   = pre_summary_cases$fatal_lo,
             pre_fatal_hi   = pre_summary_cases$fatal_hi,
             pre_daly       = pre_summary_cases$daly_tot,
             pre_daly_lo    = pre_summary_cases$daly_tot_lo,
             pre_daly_hi    = pre_summary_cases$daly_tot_hi,
      ) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_inf     = sum(infection),
        post_inf_lo  = sum(infection_lo),
        post_inf_hi  = sum(infection_hi),
        pre_inf      = sum(pre_infection),
        pre_inf_lo   = sum(pre_infection_lo),
        pre_inf_hi   = sum(pre_infection_hi),
        post_cases = sum(Cases, na.rm = TRUE),
        post_cases_lo = sum(post_low95, na.rm = TRUE),
        post_cases_hi = sum(post_hi95, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        pre_cases_lo = sum(pre_vacc_lo, na.rm = TRUE),
        pre_cases_hi  = sum(pre_vacc_hi, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        post_fatal_lo = sum(fatal_lo, na.rm = TRUE),
        post_fatal_hi = sum(fatal_hi, na.rm = TRUE),
        
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
        pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
        post_daly_lo  = sum(daly_tot_lo, na.rm = TRUE),
        post_daly_hi  = sum(daly_tot_hi, na.rm = TRUE),
        
        pre_daly   = sum(pre_daly, na.rm = TRUE),
        pre_daly_lo = sum(pre_daly_lo, na.rm = TRUE),
        pre_daly_hi = sum(pre_daly_hi, na.rm = TRUE),
        diff_daly  = pre_daly - post_daly,
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
        c(post_cases, post_cases_lo, post_cases_hi, pre_cases, pre_cases_lo, pre_cases_hi, diff, post_fatal, post_fatal_lo, post_fatal_hi,
          pre_fatal,  pre_fatal_lo, pre_fatal_hi, diff_fatal,
          post_daly,  post_daly_lo, post_daly_hi, pre_daly, pre_daly_lo, pre_daly_hi, diff_daly, post_weekly_median, post_weekly_low95,
          post_weekly_hi95, lo95, hi95),
        ~ .x / rho_p50
      ) 
    )
  
  # Global impact annotation
  global_impact_week <- summary_week_df %>%
    dplyr::group_by(Scenario) %>%
    dplyr::summarise(
      total_pre  = sum(pre_cases, na.rm=TRUE),
      total_post = sum(post_cases, na.rm=TRUE),
      total_pre_lo = sum(pre_cases_lo, na.rm = TRUE),
      total_pre_hi = sum(pre_cases_hi, na.rm = TRUE),
      total_post_lo = sum(post_cases_lo, na.rm = TRUE),
      total_post_hi = sum(post_cases_hi, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(impact = (total_pre - total_post)/total_pre *100,
                  impact_lo = (total_pre_lo - total_post_lo) / total_pre_lo * 100,
                  impact_hi = (total_pre_hi - total_post_hi) / total_pre_hi * 100
    )
  annotation_text <- paste0(global_impact_week$Scenario, ": ", round(global_impact_week$impact,1), "%", collapse="\n")
  global_impact_week$region <- region
  
  # Extract vaccination weeks
  vacc_start_week_s1 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s1   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  vacc_start_week_s2 <- scenario_result[[2]]$sim_out$vacc_start_week[5]
  vacc_end_week_s2   <- scenario_result[[2]]$sim_out$vacc_end_week[5]
  vacc_start_week_s3 <- scenario_result[[3]]$sim_out$vacc_start_week[13]
  vacc_end_week_s3   <- scenario_result[[3]]$sim_out$vacc_end_week[13]
  
  return(list(
    scenario_result    = scenario_result,
    scenario_data      = lapply(scenario_result, function(x) x$sim_result),
    scenario_df        = scenario_df,
    summary_list       = summary_list,
    summary_list_df    = summary_list_df,
    final_summ         = final_summ,
    summary_week       = summary_week,
    summary_week_df    = summary_week_df,
    global_impact_week = global_impact_week,
    annotation_text    = annotation_text,
    vacc_weeks = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
    )
  ))
}

# visualisation-------------------------------------------------------------------
## pre-post cases graph ----------------------------------------------------
epi_graph <- function(
    postsim_all,
    observed
){
  
  max_cases <- max(postsim_all$summary_week_df$hi95, na.rm = TRUE)
  
  # Extract vaccination week parameters from postsim_all
  vacc_start_week_s1 <- postsim_all$vacc_weeks$scenario1$start
  vacc_end_week_s1   <- postsim_all$vacc_weeks$scenario1$end
  vacc_start_week_s2 <- postsim_all$vacc_weeks$scenario2$start
  vacc_end_week_s2   <- postsim_all$vacc_weeks$scenario2$end
  vacc_start_week_s3 <- postsim_all$vacc_weeks$scenario3$start
  vacc_end_week_s3   <- postsim_all$vacc_weeks$scenario3$end
  
  g <- 
    ggplot(postsim_all$summary_week_df) +
    # Scenario shading ribbons (existing)
    geom_ribbon(data = postsim_all$summary_week_df %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
                aes(x = Week, ymin = 2/3 * max_cases, ymax = max_cases, fill = Scenario),
                alpha = 0.4) +
    geom_ribbon(data = postsim_all$summary_week_df %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
                aes(x = Week, ymin = 1/3 * max_cases, ymax = 2/3 * max_cases, fill = Scenario),
                alpha = 0.4) +
    geom_ribbon(data = postsim_all$summary_week_df %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
                aes(x = Week, ymin = 0, ymax = 1/3 * max_cases, fill = Scenario),
                alpha = 0.4) +
    # Post-vaccination UI ribbon (95% uncertainty interval)
    geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
                alpha = 0.2) +
    # Pre-vaccination UI ribbon (if available; omit if not)
    geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
                fill = "gray", alpha = 0.2) +
    # Post-vaccination case line (colored by scenario)
    geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases"), linewidth = 0.5) +
    # Pre-vaccination case line (dashed, in black)
    geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Pre Cases"), color = "black", linewidth = 0.5) +
    scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "1-20 years only", 
                                 "Scenario_2" = "20-59 years only", 
                                 "Scenario_3" = ">60 years only")) +
    scale_color_brewer(palette = "Set1",
                       labels = c("Scenario_1" = "1-20 years only", 
                                  "Scenario_2" = "20-59 years only", 
                                  "Scenario_3" = ">60 years only")) +
    scale_y_continuous(labels = scales::comma) +
    geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
    labs(color = "Vaccination strategy", linetype = "Type", fill = "Scenario",
         #title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
         x = "Week", y = "Predicted symptomatic reported cases") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      plot.margin = margin(5, 30, 5, 5)
    ) +
    annotate("text", x = 4, y = max(postsim_all$summary_week_df$hi95, na.rm = TRUE) * 1.1,
             label = "<----- Vaccine impact start (2 weeks after initiation)", 
             hjust = 0, vjust = 1, size = 3, color = "black") +
    annotate("text", x = Inf, y = Inf,
             hjust = 1, vjust = 1,
             label = postsim_all$annotation_text,
             color = "black", size = 3) +
    coord_cartesian(clip = "off") +
    geom_point(data = observed, aes(x = Week, y = Observed))
  theme(plot.margin = margin(5, 30, 5, 5)) 
  return(g)
}

epi_graph_nat <- function(
    postsim_all,
    observed
){
  
  max_cases <- max(postsim_all$summary_week_df$hi95, na.rm = TRUE)
  
  # Extract vaccination week parameters from postsim_all
  vacc_start_week_s1 <- postsim_all$vacc_start_week_s1$start
  vacc_end_week_s1   <- postsim_all$vacc_start_week_s1$end
  vacc_start_week_s2 <- postsim_all$vacc_start_week_s1$start
  vacc_end_week_s2   <- postsim_all$vacc_start_week_s1$end
  vacc_start_week_s3 <- postsim_all$vacc_start_week_s1$start
  vacc_end_week_s3   <- postsim_all$vacc_start_week_s1$end
  
  g <- 
    ggplot(postsim_all$summary_week_df) +
    # Scenario shading ribbons (existing)
    geom_ribbon(data = postsim_all$summary_week_df %>% filter(Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
                aes(x = Week, ymin = 0, ymax = Inf), fill = "grey70",
                alpha = 0.4) +
    #geom_ribbon(data = postsim_all$summary_week_df %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
    #            aes(x = Week, ymin = 1/3 * max_cases, ymax = 2/3 * max_cases, fill = Scenario),
    #            alpha = 0.4) +
    #geom_ribbon(data = postsim_all$summary_week_df %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
    #            aes(x = Week, ymin = 0, ymax = 1/3 * max_cases, fill = Scenario),
    #            alpha = 0.4) +
    # Post-vaccination UI ribbon (95% uncertainty interval)
    geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
                alpha = 0.2) +
    # Pre-vaccination UI ribbon (if available; omit if not)
    geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
                fill = "lightgrey", alpha = 0.2) +
    # Post-vaccination case line (colored by scenario)
    geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "With vaccination"), size = 0.5) +
    # Pre-vaccination case line (dashed, in black)
    geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Without vaccination"), color = "black", size = 0.5) +
    scale_linetype_manual(values = c("With vaccination" = "solid", "Without vaccination" = "dashed")) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "1-19 years only", 
                                 "Scenario_2" = "20-59 years only", 
                                 "Scenario_3" = "≥60 years only")) +
    scale_color_brewer(palette = "Set1",
                       labels = c("Scenario_1" = "1-19 years only", 
                                  "Scenario_2" = "20-59 years only", 
                                  "Scenario_3" = "≥60 years only")) +
    scale_y_continuous(labels = scales::comma) +
    geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
    labs(color = "Vaccination strategy", linetype = "Type", fill = "Vaccination strategy",
         x = "Week", y = "Predicted symptomatic cases") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      plot.margin = margin(5, 30, 5, 5)
    ) +
    annotate("text", x = 4, y = max(postsim_all$summary_week_df$hi95, na.rm = TRUE) * 1.1,
             label = "<----- Vaccine impact start (2 weeks after initiation)", 
             hjust = 0, vjust = 1, size = 3, color = "black") +
    annotate("text", x = Inf, y = Inf,
             hjust = 1, vjust = 1,
             label = postsim_all$annotation_text,
             color = "black", size = 3) +
    coord_cartesian(clip = "off") +
    geom_point(data = observed, aes(x = Week, y = Observed), size = 1)
  theme(plot.margin = margin(5, 30, 5, 5)) 
  
  return(g)
}

# nnv 
## v5: aligning with licensure age group
nnv_list <- function(vacc_allocation,
                     postsim_all_ui,
                     N,
                     region,
                     observed) {
  
  # Set T dynamically 
  T <- nrow(observed)
  n_scenarios <- length(target_age_list)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  default_age_vector <- gsub("[–—]", "-", default_age_vector)
  age_gr_levels <- gsub("[–—]", "-", age_gr_levels)
  
  # Create age-by-week lookup
  age_gr <- rep(default_age_vector, each = T)
  
  raw_allocation <- lapply(vacc_allocation$scenario_data, function(list){
    raw_allocation <- as.data.frame(list$raw_allocation_age)
    raw_allocation$tot_vacc <- rowSums(raw_allocation)
    list$raw_allocation <- raw_allocation
    tot_df <- setNames(as.data.frame(list$raw_allocation$tot_vacc), "tot_vacc")
  })
  
  
  raw_allocation_age <- lapply(seq_along(raw_allocation), function(id){
    df <- raw_allocation[[id]]
    tot_vacc <- sum(df$tot_vacc)
    return(tot_vacc)
  })
  
  # ② build final_summ_df
  final_summ_df <- do.call(rbind, lapply(seq_along(postsim_all_ui$final_summ), function(idx){
    
    # attach target col
    target <- c("<1 (novacc)", "1-11 years", "1-11 years", "1-11 years", "12-17 years",
                "18-59 years", "18-59 years", "18-59 years", "18-59 years",
                "18-59 years", "18-59 years", "18-59 years", "18-59 years", "18-59 years", "18-59 years",
                "60+ years", "60+ years", "60+ years", "60+ years", "60+ years")
    
    df <- postsim_all_ui$final_summ[[idx]]
    df <- df %>% mutate(
      tot_vacc = as.numeric(unlist(raw_allocation[[idx]])),
      scenario = paste0("Scenario_", idx),
      target   = target
    )
    
    #scenario_vacc <- sum(df$tot_vacc, na.rm = TRUE)
    #df <- df %>% mutate(scenario_vacc = scenario_vacc)
    
    summary_by_age <- df %>%
      group_by(scenario, AgeGroup) %>%
      summarise(
        # 1) age group당 투여된 총 백신량
        tot_vacc      = sum(tot_vacc,      na.rm = TRUE),
        pre_inf       = sum(pre_infection),
        pre_inf_lo    = sum(pre_infection_lo),
        pre_inf_hi    = sum(pre_infection_hi),
        post_inf      = sum(infection),
        post_inf_lo   = sum(infection_lo),
        post_inf_hi   = sum(infection_hi),
        pre_vacc      = sum(pre_vacc),
        pre_vacc_lo   = sum(pre_vacc_low95),
        pre_vacc_hi   = sum(pre_vacc_hi95),
        post_vacc     = sum(Median),
        post_vacc_lo  = sum(low95),
        post_vacc_hi  = sum(hi95),
        pre_fatal     = sum(pre_fatal),
        pre_fatal_lo  = sum(pre_fatal_low95),
        pre_fatal_hi  = sum(pre_fatal_hi),
        post_fatal    = sum(fatal),
        post_fatal_lo = sum(fatal_lo),
        post_fatal_hi = sum(fatal_hi),
        pre_daly      = sum(pre_daly),
        pre_daly_lo   = sum(pre_daly_low95),
        pre_daly_hi   = sum(pre_daly_hi),
        post_daly     = sum(daly_tot),
        post_daly_lo  = sum(daly_tot_lo),
        post_daly_hi  = sum(daly_tot_hi),
        # 2) age group당 averted diff 합계
        diff_inf      = sum(diff_inf),
        diff_inf_lo   = sum(diff_inf_lo),
        diff_inf_hi   = sum(diff_inf_hi),
        diff          = sum(diff,          na.rm = TRUE),
        diff_low      = sum(diff_low,      na.rm = TRUE),
        diff_hi       = sum(diff_hi,       na.rm = TRUE),
        # (fatal / daly 도 동일하게)
        diff_fatal    = sum(diff_fatal,    na.rm = TRUE),
        diff_fatal_low= sum(diff_fatal_low,na.rm = TRUE),
        diff_fatal_hi = sum(diff_fatal_hi, na.rm = TRUE),
        diff_daly     = sum(diff_daly,     na.rm = TRUE),
        diff_daly_low = sum(diff_daly_low, na.rm = TRUE),
        diff_daly_hi  = sum(diff_daly_hi,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # 2) 바로 스케일링
      mutate(region = region) %>%
      left_join(rho_df, by = "region") %>%
      mutate(
        across(
          .cols = c(pre_vacc, pre_vacc_lo, pre_vacc_hi, 
                    post_vacc, post_vacc_lo, post_vacc_hi,
                    pre_fatal, pre_fatal_lo, pre_fatal_hi,
                    post_fatal, post_fatal_lo, post_fatal_hi,
                    pre_daly, pre_daly_lo, pre_daly_hi,
                    post_daly, post_daly_lo, post_daly_hi,
                    diff, diff_low, diff_hi,
                    diff_fatal, diff_fatal_low, diff_fatal_hi,
                    diff_daly, diff_daly_low, diff_daly_hi),
          .fns  = ~ .x / rho_p50
        ), 
        pre_infection  = pre_vacc / 0.5242478,
        post_infection = post_vacc / 0.5242478
      ) %>%
      mutate(
        scenario_vacc = raw_allocation_age[[idx]]
      ) %>%
      mutate(
        nnv_inf     = scenario_vacc / diff_inf,
        nnv_inf_lo  = scenario_vacc / diff_inf_hi,
        nnv_inf_hi  = scenario_vacc / diff_inf_lo,
        
        nnv         = scenario_vacc / diff,
        nnv_lo      = scenario_vacc / diff_hi,
        nnv_hi      = scenario_vacc / diff_low,
        
        nnv_fatal       = scenario_vacc / diff_fatal,
        nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
        nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
        
        nnv_daly        = scenario_vacc / diff_daly,
        nnv_daly_lo     = scenario_vacc / diff_daly_hi,
        nnv_daly_hi     = scenario_vacc / diff_daly_low
        
      )
    summary_by_age <- summary_by_age %>% mutate(
      target = target
    )  %>%
      relocate(target, .after = scenario)  
    
  })
  )
  
  final_summ_df <- final_summ_df %>%
    mutate(
      tot_pop        = rep(N, n_scenarios),
      tot_vacc_prop  = tot_vacc / tot_pop,
      age_gr        = factor(default_age_vector[AgeGroup], 
                             levels = gsub("[–—]", "-", default_age_vector))
    )
  
  #final_summ_df$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
  final_summ_df$age_gr <- factor(final_summ_df$age_gr, levels = age_gr_levels)
  final_summ_df$region <- region
  
  # —— ② scenario × target별 per1M 요약 —— 
  per1M_summary <- final_summ_df %>%
    group_by(scenario, target) %>%
    summarise(
      scenario_vacc = first(scenario_vacc),
      tot_pop       = sum(tot_pop),
      vacc_prop     = scenario_vacc / tot_pop, 
      diff          = sum(diff,        na.rm = TRUE),
      diff_low      = sum(diff_low,    na.rm = TRUE),
      diff_hi       = sum(diff_hi,     na.rm = TRUE),
      diff_fatal    = sum(diff_fatal,  na.rm = TRUE),
      diff_fatal_low= sum(diff_fatal_low,na.rm = TRUE),
      diff_fatal_hi = sum(diff_fatal_hi, na.rm = TRUE),
      diff_daly     = sum(diff_daly,     na.rm = TRUE),
      diff_daly_low = sum(diff_daly_low, na.rm = TRUE),
      diff_daly_hi  = sum(diff_daly_hi,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      per1M_diff       = diff        / scenario_vacc * 1e6,
      per1M_diff_lo    = diff_low    / scenario_vacc * 1e6,
      per1M_diff_hi    = diff_hi     / scenario_vacc * 1e6,
      
      per1M_fatal      = diff_fatal      / scenario_vacc * 1e6,
      per1M_fatal_lo   = diff_fatal_low  / scenario_vacc * 1e6,
      per1M_fatal_hi   = diff_fatal_hi   / scenario_vacc * 1e6,
      
      per1M_daly       = diff_daly       / scenario_vacc * 1e6,
      per1M_daly_lo    = diff_daly_low   / scenario_vacc * 1e6,
      per1M_daly_hi    = diff_daly_hi    / scenario_vacc * 1e6,
      
      nnv         = scenario_vacc / diff,
      nnv_lo      = scenario_vacc / diff_hi,
      nnv_hi      = scenario_vacc / diff_low,
      
      nnv_fatal       = scenario_vacc / diff_fatal,
      nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
      nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
      
      nnv_daly        = scenario_vacc / diff_daly,
      nnv_daly_lo     = scenario_vacc / diff_daly_hi,
      nnv_daly_hi     = scenario_vacc / diff_daly_low
    )
  
  # ⑥ return
  list(
    raw_allocation     = raw_allocation,
    raw_allocation_age = raw_allocation_age,
    final_summ_df      = final_summ_df,
    per1M_summary      = per1M_summary
  )
}

## v2
nnv_gg <- function(final_nnv_df,
                   y_var = "nnv_fatal",  # name of y variable as a string
                   y_lab = "NNV to avert a single fatal case",
                   x_lab = "Age group",
                   title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2") {
  
  final_nnv_df <- final_nnv_df %>% filter(nnv !=0)
  
  # Convert the y_var string to a symbol for tidy evaluation
  y_sym <- sym(y_var)
  # Create symbols for the lower and upper bounds; assume they follow the naming pattern y_var_lo, y_var_hi
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))
  
  g <- ggplot(final_nnv_df, aes(x = age_gr, y = !!y_sym, fill = scenario)) +
    geom_bar(stat = "identity", alpha = 0.5) +
    geom_errorbar(aes(ymin = !!lower_sym, ymax = !!upper_sym),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    #facet_wrap(~scenario, nrow = 3, 
    #           labeller = as_labeller(c("Scenario_1" = "<20 years only", 
    #                                    "Scenario_2" = "20-59 years only", 
    #                                    "Scenario_3" = ">60 years only"))) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "<20 years only", 
                                 "Scenario_2" = "20-59 years only", 
                                 "Scenario_3" = ">60 years only")) +
    scale_y_continuous(labels = scales::comma) +
    scale_x_discrete(labels = function(x) gsub("years", "", x)) +
    theme_pubclean() +
    labs(
      y = y_lab,
      x = x_lab,
      title = title
    ) +
    theme(legend.position = "right", 
          title = element_text(size = 7, face = 'bold'),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(g)
}
##

nnv_nat_gg <- function(combined_df, 
                       y_var,
                       y_lab = "NNV to avert a single fatal case",
                       x_lab = "Age group",
                       title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2") {
  
  # Convert the y_var string to a symbol for tidy evaluation
  y_sym <- sym(y_var)
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))
  
  region_order <- c(
    "Ceará", "Paraíba", "Piauí", "Alagoas",               # High
    "Tocantins", "Pernambuco",                   # Moderate
    "Rio Grande do Norte", "Sergipe", "Bahia", "Minas Gerais", "Goiás"      # Low
  )
  
  setting_order <- c("High", "Moderate", "Low")  # case-sensitive
  
  combined_df <- combined_df %>%
    mutate(
      region = factor(region, levels = region_order),
      setting = factor(setting, levels = setting_order)
    )
  
  g <- 
    #ggplot(combined_df)+
    #geom_bar(aes(x = age_gr, y = !!y_sym, fill = scenario), stat = "identity", alpha = 0.7)+
    
    ggplot(combined_df, aes(x = age_gr, y = !!y_sym, fill = setting)) +
    geom_bar(stat = "identity", alpha = 0.5)+
    geom_errorbar(aes(ymin = !!lower_sym, ymax = !!upper_sym),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    facet_grid(scenario~region, 
               labeller = labeller(
                 scenario = c("Scenario_1" = "<20 years", 
                              "Scenario_2" = "20-59 years", 
                              "Scenario_3" = "≥60 years")
               ))+
    scale_x_discrete(labels = function(x) gsub("years", "", x)) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "<20 years", 
                                 "Scenario_2" = "20-59 years", 
                                 "Scenario_3" = "≥60 years"),
                      name = "Setting")+
    scale_y_continuous(labels = comma)+
    theme_pubclean()+
    labs(
      y = y_lab,
      x = x_lab,
      title = title
    )+
    theme(legend.position = "right", 
          title =element_text(size=7, face='bold'),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(g)
}


nnv_nat_setting <- function(combined_df, 
                            y_var,
                            y_lab = "NNV to avert a single symptomatic case",
                            x_lab = "Vaccination strategy",
                            title = "Coverage:50%, Delivery speed:10%, Deployment: Week 2") {
  
  y_sym     <- sym(y_var)
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))
  
  df <- combined_df %>%
    filter(!is.na(target)) %>%
    mutate(
      setting  = factor(setting,  levels = c("National","High","Moderate","Low")),
      # scenario와 target 레이블을 일치시킵니다
      scenario = factor(scenario,
                        levels = c("Scenario_1","Scenario_2","Scenario_3"),
                        labels = c("1–19 years","20–59 years","≥60 years")),
      target   = factor(target,
                        levels = c("<20 years","20-59 years",">60 years"),
                        labels = c("1–19 years","20–59 years","≥60 years"))
    ) %>%
    # 직접 vs 간접 효과: 레이블이 동일한 경우 Direct
    mutate(
      effect = if_else(scenario == target, "Direct", "Indirect")
    )
  
  sum_df <- df %>%
    group_by(setting, scenario, effect) %>%
    summarise(
      value = sum(!!y_sym,     na.rm = TRUE),
      lo    = sum(!!lower_sym, na.rm = TRUE),
      hi    = sum(!!upper_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4) 플롯
  ggplot(sum_df, aes(x = scenario, y = value, fill = effect)) +
    geom_col(
      position = position_dodge(width = 0.7),
      width    = 0.6,
      alpha    = 0.8
    ) +
    geom_errorbar(
      aes(ymin = lo, ymax = hi),
      position = position_dodge(width = 0.7),
      width    = 0.2,
      color    = "black"
    ) +
    facet_wrap(~ setting, nrow = 1) +
    scale_fill_manual(
      values = c("Direct" = "#66C2A5", "Indirect" = "#FC8D62"),
      labels = c("Direct"   = "NNV to avert a target group case",
                 "Indirect" = "NNV to avert a non-target group case"),
      name   = ""
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      title = title,
      x     = "Vaccination strategy",
      y     = y_lab
    ) +
    theme_pubclean(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}



nnv_nat_setting_summ <- function(combined_df, 
                                 y_var,
                                 y_lab = "NNV to avert a single fatal case",
                                 x_lab = "Scenario",
                                 title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2") {
  
  combined_df <- combined_df %>% filter(nnv !=0)
  
  # Convert the y_var string to a symbol for tidy evaluation
  y_sym <- sym(y_var)
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))
  
  setting_order <- c("National", "High", "Moderate", "Low")  # case-sensitive
  
  combined_df <- combined_df %>%
    mutate(
      setting = factor(setting, levels = setting_order)
    )
  
  combined_df <- combined_df  %>%
    mutate(
      setting  = factor(setting,  levels = c("National","High","Moderate","Low")),
      # scenario와 target 레이블을 일치시킵니다
      scenario = factor(scenario,
                        levels = c("Scenario_1","Scenario_2","Scenario_3", "Scenario_4"),
                        labels = c("1–11 years","12–17 years", "18-59 years", "60+ years"))
    )
  
  g <- 
    #ggplot(combined_df)+
    #geom_bar(aes(x = age_gr, y = !!y_sym, fill = scenario), stat = "identity", alpha = 0.7)+
    
    ggplot(combined_df, aes(x = scenario, y = !!y_sym, fill = scenario)) +
    geom_bar(stat = "identity", alpha = 0.5, position = position_dodge(width = 0.8), width = 0.8)+
    geom_errorbar(aes(ymin = !!lower_sym, ymax = !!upper_sym),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    facet_grid(~setting, 
               labeller = labeller(
                 scenario = c("Scenario_1" = "1-11 years", 
                              "Scenario_2" = "12-17 years", 
                              "Scenario_3" = "18-59 years",
                              "Scenario_4" = "60+ years")
               ))+
    #scale_x_discrete(labels = function(x) gsub("years", "", x)) +
    scale_x_discrete(
      labels = c(
        "1–11 years"  = "Strategy 1",
        "12-17 years" = "Strategy 2",
        "18-59 years"   = "Strategy 3",
        "60+ years"     = "Strategy 4"
      ))+
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "1-11 years", 
                                 "Scenario_2" = "12-17 years", 
                                 "Scenario_3" = "18-59 years",
                                 "Scenario_4" = "60+ years"),
                      name = "Vaccination strategy")+
    scale_y_continuous(labels = comma)+
    theme_pubclean()+
    labs(
      y = y_lab,
      x = x_lab,
      title = title
    )+
    theme(legend.position = "right", 
          axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
          title =element_text(size=7, face='bold'),
          strip.background  = element_rect(fill = "grey95", colour = NA)
    )
  
  return(g)
}





