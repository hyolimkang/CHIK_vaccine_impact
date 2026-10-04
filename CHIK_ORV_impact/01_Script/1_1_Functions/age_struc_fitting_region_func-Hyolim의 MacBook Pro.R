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


## pre-vacc sim (95%median point est) ----------------------------------------------------
summarize_simulation <- function(
    sim_result,        # e.g. output from sirv_sim_coverageSwitch()
    age_gr_levels = c("<5 years",
                      "5-9 years",
                      "10-14 years",
                      "15-19 years",
                      "20-24 years",
                      "25-29 years",
                      "30-34 years",
                      "35-39 years",
                      "40-44 years",
                      "45-49 years",
                      "50-54 years",
                      "55-59 years",
                      "60-64 years",
                      "65-69 years",
                      "70-74 years",
                      "75-79 years",
                      "80-84 years",
                      "85-89 years"),     # factor levels for age groups
    age_gr        = rep(c("<5 years",
                          "5-9 years",
                          "10-14 years",
                          "15-19 years",
                          "20-24 years",
                          "25-29 years",
                          "30-34 years",
                          "35-39 years",
                          "40-44 years",
                          "45-49 years",
                          "50-54 years",
                          "55-59 years",
                          "60-64 years",
                          "65-69 years",
                          "70-74 years",
                          "75-79 years",
                          "80-84 years",
                          "85-89 years"
    ), 52),            # repeated age group labels for each row in weekly data
    lhs_sample_young,  # data for younger ages (for YLD calculation)
    lhs_old,           # data for older ages (for YLD calculation)
    le_sample,         # life-expectancy sample
    hosp,              # hospitalization rate vector or single numeric
    fatal,             # fatality rate for hospitalized
    nh_fatal           # fatality rate for non-hospitalized
) {
  library(dplyr)
  
  # 1) Convert age_stratified_cases to a data frame
  summary_cases_pre <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Median")
  colnames(summary_cases_pre) <- c("AgeGroup", "Week", "Median")
  
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week = as.numeric(Week)
    )
  
  # Attach age group labels
  summary_cases_pre$age_gr <- age_gr
  summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)
  
  # 2) Summarize total cases by Week
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(Median = sum(Median)) %>%
    mutate(Scenario = "Pre-vaccination")
  
  # 3) Calculate hospitalizations, fatalities, YLD, YLL, etc.
  #    Here we replicate your stepwise approach
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      hosp_rate = rep(hosp, 52),  # ensure length matches rows
      hospitalised = Median * hosp_rate,
      non_hospitalised = Median - hospitalised,
      fatality = rep(hosp, 52),
      nh_fatality = rep(hosp, 52),
      fatal = (hospitalised * fatality) +
        (non_hospitalised * nh_fatality)
    ) %>%
    arrange(Week, AgeGroup) %>%
    group_by(Week, AgeGroup) %>%
    mutate(
      cum_fatal = cumsum(fatal),
      cum_hosp  = cumsum(hospitalised)
    ) %>%
    ungroup()
  
  # 4) Summarize by AgeGroup
  summary_cases_pre_age <- summary_cases_pre %>%
    group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Median),
      hospitalised = sum(hospitalised),
      fatality     = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>%
    mutate(Scenario = "Pre-vaccination")
  
  summary_cases_pre_age$age_gr <- age_gr[1:length(unique(summary_cases_pre$AgeGroup))]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # 5) Summarize by Week for entire population
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      hospitalised = sum(hospitalised),
      fatality     = first(fatality),
      fatal        = sum(fatal)
    ) %>%
    mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal  = cumsum(fatal),
      cum_hosp   = cumsum(hospitalised)
    )
  
  # 6) Example of DALY calculations (similar to your code)
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
    chr_prop     = chr_6m + chr_12m + chr_30m, 
  )
  
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      age_numeric = case_when(
        str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
        TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
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
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_total    = yld_acute + yld_subacute + yld_chronic
    ) %>% mutate(
      # le left 
      le_left = case_when(
        age_numeric %in% c(0, 5)   ~ quantile(le_sample$le_1, 0.5),
        age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
        age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
        age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
        age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
        age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
        age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
        age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
        age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
      ) 
    ) %>% mutate(
      # YLL  
      yll = fatal * le_left,
      
      # total DALY
      daly_tot = yld_total + yll,
      cum_daly = cumsum(daly_tot)
    )
  
  # Summaries by AgeGroup and by Week (DALY, YLD, etc.)
  summary_cases_pre_age <- summary_cases_pre %>% group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Median),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    )
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre$age_gr[1:18], levels = age_gr_levels)
  
  summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
    summarise(
      Median   = sum(Median),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination",
      cum_fatal    = cumsum(fatal),
      cum_hosp     = cumsum(hospitalised),
      cum_yld_acute    = cumsum(yld_acute),
      cum_yld_subacute = cumsum(yld_subacute),
      cum_yld_chronic  = cumsum(yld_chronic),
      cum_yld_total    = cumsum(yld_total),
      cum_yll          = cumsum(yll),
      cum_daly_tot     = cumsum(daly_tot)
    )
  
  # Return whatever data frames you need
  return(list(
    summary_cases_pre = summary_cases_pre,
    summary_cases_pre_all = summary_cases_pre_all,
    summary_cases_pre_age = summary_cases_pre_age
  ))
}

## postvacc sim  (95%median point est)----------------------------------------------------
run_simulation_scenarios <- function(target_age_list, supply, 
                                     N, param, bra_foi_state_summ, age_groups, region_name,
                                     hosp, fatal, nh_fatal,
                                     lhs_sample_young, lhs_old, le_sample,
                                     age_gr, age_gr_levels) {

  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    sim_result <- sirv_sim_coverageSwitch(
      T = 52,
      A = 18,
      N = N,
      r = rep(0, 18),
      base_beta = param$base_beta,
      I0_draw = param$I0,
      R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups)),
      rho = param$rho,
      gamma = param$gamma,
      delay = 53,                   # Vaccination starts from Week x
      VE_block = 0.989,
      coverage_threshold = 0.7,
      target_age = target,
      total_coverage = 1,
      total_supply = supply[[scenario_index]],
      weekly_delivery_speed = 0.1
    )
    
    sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
      # First, create the basic data frame
      mutate(
        Scenario  = scenario_index,
        AgeGroup  = as.numeric(Var1),
        Week      = as.numeric(Var2)
      ) %>%
      mutate(
        hosp_rate         = rep(hosp, 52),
        hospitalised      = Cases * hosp_rate,
        non_hospitalised  = Cases - hospitalised,
        fatality          = rep(fatal, 52),   # Ensure fatal is correctly replicated
        nh_fatality       = rep(nh_fatal, 52),
        fatal             = (hospitalised * fatality + non_hospitalised * nh_fatality)
      ) %>%
      mutate(
        age_numeric = case_when(
          str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
          TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
        )
      ) %>%
      # Add YLD parameters that are independent of age
      mutate(
        dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Add age-dependent parameters for subacute and chronic outcomes
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
      # Calculate YLD estimates
      mutate(
        yld_acute    = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_total    = yld_acute + yld_subacute + yld_chronic
      ) %>%
      # Assign life expectancy based on age_numeric and calculate YLL and DALY metrics
      mutate(
        le_left = case_when(
          age_numeric %in% c(0, 5)    ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10, 15)  ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20, 25)  ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30, 35)  ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40, 45)  ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50, 55)  ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60, 65)  ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70, 75)  ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80, 85)  ~ quantile(le_sample$le_9, 0.5)
        )
      ) %>%
      mutate(
        yll      = fatal * le_left,
        daly_tot = yld_total + yll,
        cum_daly = cumsum(daly_tot)
      )
    
    # Save simulation result and processed data frame for this scenario
    scenario_result[[scenario_index]] <- list(sim_result = sim_result, sim_df = sim_df)
  }
  
  return(scenario_result)
}

## postsim case aversion aggregate  (95%median point est) ----------------------------------------------------
postsim_all <- function(target_age_list, supply, N, param, 
                        bra_foi_state_summ, age_groups, region,
                        lhs_sample_young, 
                        lhs_old, 
                        le_sample,
                        hosp, 
                        fatal, 
                        nh_fatal,
                        age_gr, age_gr_levels,
                        pre_summary_cases_age,
                        pre_summary_cases) {
  # target_age_list: list of 0/1 vectors for each scenario
  # supply: list of total supply values for each scenario
  # N: population vector (e.g., N_pa)
  # param: a list with elements base_beta, I0, rho, gamma, etc.
  # bra_foi_state_summ: data frame with avg_foi and NAME_1 for FOI
  # age_groups: numeric vector of age group midpoints used in FOI calculation
  # lhs_sample_young, lhs_old: lists/data frames with YLD parameters for younger and older ages
  # le_sample: life expectancy sample (for YLL calculation)
  # hosp, fatal, nh_fatal: numeric values (or vectors) for hospitalization, fatality rates
  # age_gr: vector of age group labels repeated for each week (e.g., length = number of rows in sim_df)
  # age_gr_levels: factor levels for age_gr
  # pre_summary_cases_age: a data frame summarizing pre-vaccination by age (e.g., summary_cases_pre_age)
  
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    sim_result <- sirv_sim_coverageSwitch(
      T = 52,
      A = 18,
      N = N,
      r = rep(0, 18),
      base_beta = param$base_beta,
      I0_draw = param$I0,
      R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region] * age_groups),
      rho = param$rho,
      gamma = param$gamma,
      delay = 2,
      VE_block = 0.989,
      coverage_threshold = 1,
      target_age = target,
      total_coverage = 1,
      total_supply = supply[[scenario_index]],
      weekly_delivery_speed = 0.1
    )
    
    sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
      mutate(
        Scenario  = scenario_index,
        AgeGroup  = as.numeric(Var1),
        Week      = as.numeric(Var2)
      ) %>%
      mutate(
        hosp_rate         = rep(hosp, 52),
        hospitalised      = Cases * hosp_rate,
        non_hospitalised  = Cases - hospitalised,
        fatality          = rep(fatal, 52),
        nh_fatality       = rep(nh_fatal, 52),
        fatal             = (hospitalised * fatality + non_hospitalised * nh_fatality)
      ) %>%
      mutate(
        age_numeric = case_when(
          str_detect(age_gr, "<5") ~ 0,
          TRUE ~ as.numeric(str_extract(age_gr, "^\\d+"))
        )
      ) %>%
      # Add YLD parameters (assumed independent of age)
      mutate(
        dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Add age-dependent YLD parameters
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
      # Calculate YLD estimates
      mutate(
        yld_acute    = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_total    = yld_acute + yld_subacute + yld_chronic
      ) %>%
      # Assign life expectancy and calculate YLL and DALY metrics
      mutate(
        le_left = case_when(
          age_numeric %in% c(0, 5)    ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10, 15)  ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20, 25)  ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30, 35)  ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40, 45)  ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50, 55)  ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60, 65)  ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70, 75)  ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80, 85)  ~ quantile(le_sample$le_9, 0.5)
        )
      ) %>%
      mutate(
        yll      = fatal * le_left,
        daly_tot = yld_total + yll,
        cum_daly = cumsum(daly_tot)
      )
    
    scenario_result[[scenario_index]] <- list(sim_result = sim_result, sim_df = sim_df)
  }
  
  # Extract simulation data (raw sim_result) and processed data frame from each scenario
  scenario_data <- lapply(scenario_result, function(x) x$sim_result)
  scenario_df   <- lapply(scenario_result, function(x) x$sim_df)
  
  # Create a summary list that compares pre-vaccination (provided by pre_summary_cases) with simulated Cases
  summary_list <- lapply(seq_along(scenario_result), function(idx) {
    df <- scenario_result[[idx]]$sim_df
    df$pre_vacc <- pre_summary_cases$Median
    df$pre_fatal <- pre_summary_cases$fatal
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
  
  
  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%  # Include Scenario in grouping
      summarise(
        Median = sum(Cases, na.rm = TRUE),  # Summarize MedianCases
        yld_acute    = sum(yld_acute),
        yld_subacute = sum(yld_subacute),
        yld_chronic  = sum(yld_chronic),
        yld_total    = sum(yld_total),
        yll          = sum(yll),
        daly_tot     = sum(daly_tot),
        fatal  = sum(fatal),
        .groups = "drop"  # Drop grouping for clarity in the output
      ) %>% mutate(
        pre_vacc = pre_summary_cases_age$Median,
        pre_yld_acute = pre_summary_cases_age$yld_acute,
        pre_yld_subacute = pre_summary_cases_age$yld_subacute,
        pre_yld_chronic = pre_summary_cases_age$yld_chronic,
        pre_yld_total = pre_summary_cases_age$yld_total,
        pre_yll = pre_summary_cases_age$yll,
        pre_daly = pre_summary_cases_age$daly_tot,
        pre_fatal = pre_summary_cases_age$fatal,
        diff     = pre_vacc - Median,
        impact   = diff / pre_vacc * 100,
        diff_fatal = pre_fatal - fatal,
        impact_fatal = diff_fatal / pre_fatal * 100,
        # averted 
        diff_yld_acute    = pre_yld_acute - yld_acute,
        diff_yld_subacute = pre_yld_subacute - yld_subacute,
        diff_yld_chronic = pre_yld_chronic - yld_chronic,
        diff_yld_total   = pre_yld_total - yld_total,
        diff_yll         = pre_yll - yll,
        diff_daly        = pre_daly - daly_tot
      ) 
  })
  
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median) %>% group_by(Week, Scenario) %>%
      summarise(post_cases = sum(Cases),
                pre_cases  = sum(pre_vacc),
                diff       = pre_cases - post_cases,
                Scenario   = first(Scenario)
      ) 
    
  })
  
  
  summary_week_df <- do.call(rbind, summary_week)
  
  global_impact_week <- summary_week_df %>%
    group_by(Scenario) %>%
    summarise(
      total_post_cases = sum(post_cases),
      total_pre_case   = sum(pre_cases),
      .groups = "drop"
    ) %>%
    mutate(
      diff = total_pre_case - total_post_cases,
      impact = diff / total_pre_case  * 100
    )
  
  annotation_text <- paste0(
    global_impact_week$Scenario, ": ", 
    round(global_impact_week$impact, 1), "%",
    collapse = "\n"
  )
  
  # Extract vaccination start and end weeks for selected ages
  vacc_start_week_s1 <- scenario_data[[1]]$vacc_start_week[3]
  vacc_end_week_s1   <- scenario_data[[1]]$vacc_end_week[3]
  
  vacc_start_week_s2 <- scenario_data[[2]]$vacc_start_week[5]
  vacc_end_week_s2   <- scenario_data[[2]]$vacc_end_week[5]
  
  vacc_start_week_s3 <- scenario_data[[3]]$vacc_start_week[13]
  vacc_end_week_s3   <- scenario_data[[3]]$vacc_end_week[13]
  
  # Return all processed results as a list
  return(list(
    scenario_result = scenario_result,
    scenario_data   = scenario_data,
    scenario_df     = scenario_df,
    summary_list    = summary_list,
    summary_list_df = summary_list_df,
    final_summ      = final_summ,
    summary_week    = summary_week,
    summary_week_df = summary_week_df,
    global_impact_week = global_impact_week,
    annotation_text = annotation_text,
    vacc_weeks = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
      
    )
  ))
}

## prevacc posterior draws (95% UI) ----------------------------------------------------
make_region_foi_draws <- function(regions, bra_foi_states, n_draws = 1000) {
  draws_list <- list()
  set.seed(123)
  
  foi_cols <- paste0("foi", 1:100)  # foi1 ~ foi100
  
  for (reg in regions) {
    
    foi_row <- sf::st_drop_geometry(
      bra_foi_states[bra_foi_states$NAME_1 == reg, ]
    )
    
    foi_values <- unlist(foi_row[, foi_cols], use.names = FALSE)
    
    foi_values <- foi_values[!is.na(foi_values)]
    
    draws_list[[reg]] <- sample(foi_values, size = n_draws, replace = TRUE)
  }
  
  return(draws_list)
}

regions <- c(
  "Ceará",
  "Pernambuco",
  "Minas Gerais",
  "Bahia",
  "Paraíba",
  "Rio Grande do Norte",
  "Piauí",
  "Alagoas",
  "Tocantins",
  "Sergipe",
  "Goiás"
)

set.seed(123)
foi_draws_list <- make_region_foi_draws(regions, bra_foi_states)


######

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
  n_draws <- length(posterior_idx)
  # 2. based on posterior sample nums
  #n_draws <- length(posterior$gamma) 
  # 3. based on lhs_sapmle run 
  #n_draws <- nrow(lhs_sample)
  sim_results_list_rawinf <- vector("list", n_draws)
  sim_results_list_rawsymp <- vector("list", n_draws)
  
  set.seed(123)
  foi_draws <- foi_draws_list[[region]]
  
  #I0_resampled        <- posterior$I0[posterior_idx, ]        
  #base_beta_resampled <- posterior$base_beta[posterior_idx, ] 
  
  # Run simulation for each draw and store the full age_stratified_cases matrix (dimensions: A x T)
  for (i in seq_len(n_draws)) {
    
    idx <- posterior_idx[i] 
    
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
    base_beta_draw <- posterior$base_beta[idx, ]
    I0_draw        <- posterior$I0[idx, ]
    gamma_draw     <- lhs_sample$gamma[idx]
    rho_draw       <- lhs_sample$rho[idx]
    sigma_draw     <- lhs_sample$sigma[idx]
    FOI_rand_draw  <- lhs_sample$foi[idx]
    
    #R0_vec <- (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region] * age_groups))
    R0_vec <- (1 - exp(- FOI_rand_draw * age_groups))

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

## prevacc 95%UI summarise ----------------------------------------------------
summarise_presim_ui <- function(
    sim_result,        # output from simulate_pre_ui_age
    observed,
    age_gr_levels = c("<5 years",
                      "5-9 years",
                      "10-14 years",
                      "15-19 years",
                      "20-24 years",
                      "25-29 years",
                      "30-34 years",
                      "35-39 years",
                      "40-44 years",
                      "45-49 years",
                      "50-54 years",
                      "55-59 years",
                      "60-64 years",
                      "65-69 years",
                      "70-74 years",
                      "75-79 years",
                      "80-84 years",
                      "85-89 years"),
    lhs_sample_young,  # data for younger ages (for YLD calculation)
    lhs_old,           # data for older ages (for YLD calculation)
    le_sample,         # life-expectancy sample
    hosp,              # hospitalization rate vector or single numeric
    fatal,             # fatality rate for hospitalized
    nh_fatal,           # fatality rate for non-hospitalized
    region
) {
  
  # Set T dynamically
  T =  nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  age_gr <- rep(default_age_vector, each = T)
  
  # --- Create Age-Stratified Summary Data Frame ---
  median_df <- as.data.frame.table(sim_result$median_by_age, responseName = "Median")
  colnames(median_df) <- c("AgeGroup", "Week", "Median")
  
  median_df <- median_df %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week = as.numeric(Week)
    )
  
  low95_df <- as.data.frame.table(sim_result$low95_by_age, responseName = "low95")
  hi95_df  <- as.data.frame.table(sim_result$hi95_by_age, responseName = "hi95")
  
  summary_cases_pre <- median_df %>%
    mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]  # assign provided age group labels
    )
  summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)
  
  # --- Weekly Summary ---
  weekly_df <- data.frame(
    Week = 1:length(sim_result$weekly_cases_median),
    weekly_median = sim_result$weekly_cases_median,
    weekly_low95  = sim_result$weekly_cases_low95,
    weekly_hi95   = sim_result$weekly_cases_hi95
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
      non_hospitalised = Median - hospitalised,
      fatality = rep(fatal, T),          # Ensure fatal is correctly replicated
      nh_fatality = rep(nh_fatal, T),
      fatal = (hospitalised * fatality +   
                 non_hospitalised * nh_fatality)
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
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    )
  
  summary_cases_pre_age$age_gr <- age_gr[1:5]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # --- Summarize by Week for Entire Population ---
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      hospitalised = sum(hospitalised),
      fatality     = first(fatality),
      fatal        = sum(fatal)
    ) %>%
    mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal  = cumsum(fatal),
      cum_hosp   = cumsum(hospitalised)
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
    mutate(
      age_numeric = case_when(
        str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
        TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
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
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_total    = yld_acute + yld_subacute + yld_chronic
    ) %>% mutate(
      # le left 
      le_left = case_when(
        age_numeric %in% c(0, 5)   ~ quantile(le_sample$le_1, 0.5),
        age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
        age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
        age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
        age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
        age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
        age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
        age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
        age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
      ) 
    ) %>% mutate(
      # YLL  
      yll = fatal * le_left,
      
      # total DALY
      daly_tot = yld_total + yll,
      cum_daly = cumsum(daly_tot)
    )
  
  # --- Final Summaries ---
  summary_cases_pre_age <- summary_cases_pre %>% group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Median),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    ) %>%
    ungroup()
  
  summary_cases_pre_age$age_gr <- age_gr[1:5]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal  = cumsum(fatal),
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


## v2.------------------------
summarise_presim_ui <- function(
    sim_result,        # output from simulate_pre_ui_age
    observed,
    age_gr_levels = c(
      "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
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
  
  # Set T dynamically
  T =  nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  age_gr <- rep(default_age_vector, each = T)
  
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
  
  summary_cases_pre <- median_df %>%
    mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]  # assign provided age group labels
    )
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
      non_hospitalised_lo = low95 - hospitalised_hi,  # conservative: lower bound for non-hosp
      non_hospitalised_hi = hi95 - hospitalised_lo,
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
      hospitalised = sum(hospitalised),
      hospitalised_lo = sum(hospitalised_lo),
      hospitalised_hi = sum(hospitalised_hi),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    )
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # --- Summarize by Week for Entire Population ---
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(
      Median       = sum(Median),
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
    mutate(
      age_numeric = case_when(
        str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
        TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
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
        ((low95 - hospitalised_hi) * dw_nonhosp * dur_acute),
      yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
        ((hi95 - hospitalised_lo) * dw_nonhosp * dur_acute),
      
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
        ((low95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
      
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
        ((low95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
      
      yld_total    = yld_acute + yld_subacute + yld_chronic,
      yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
      yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
    ) %>% mutate(
      # le left 
      le_left = case_when(
        age_numeric %in% c(0, 5)   ~ quantile(le_sample$le_1, 0.5),
        age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
        age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
        age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
        age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
        age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
        age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
        age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
        age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
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
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
    summarise(
      Median       = sum(Median),
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

# v3. changed age groups 
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
  
  # Set T dynamically
  T =  nrow(observed)
  
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
  
  summary_cases_pre <- median_df %>%
    mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]  # assign provided age group labels
    )
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
      non_hospitalised_lo = low95 - hospitalised_hi,  # conservative: lower bound for non-hosp
      non_hospitalised_hi = hi95 - hospitalised_lo,
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
    mutate(
      age_numeric = case_when(
        age_gr == "<1" ~ 0,
        age_gr == "1-4" ~ 2.5,
        age_gr == "5–9" ~ 7,
        age_gr == "10-11" ~ 10.5,
        age_gr == "12–17" ~ 14.5,
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
        ((low95 - hospitalised_hi) * dw_nonhosp * dur_acute),
      yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
        ((hi95 - hospitalised_lo) * dw_nonhosp * dur_acute),
      
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
        ((low95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
      
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
        ((low95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
      
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


## postvacc sim (95%UI) ----------------------------------------------------
run_simulation_scenarios_ui <- function(target_age_list, 
                                        #supply, 
                                        observed,
                                        N, bra_foi_state_summ, age_groups, region_name,
                                        hosp, fatal, nh_fatal,
                                        lhs_sample_young, lhs_old, le_sample,
                                        age_gr_levels,
                                        prevacc_ui = NULL,
                                        posterior) {
  # target_age_list: list of 0/1 vectors for each scenario
  # supply: list of total supply values for each scenario
  # N: population vector (e.g., N_pa)
  # param: a list with elements base_beta, I0, gamma, rho, etc. (each with n_draws rows)
  # bra_foi_state_summ: data frame with avg_foi and NAME_1 for FOI
  # age_groups: numeric vector of age group midpoints used in FOI calculation
  # hosp, fatal, nh_fatal: numeric values (or vectors) for hospitalization, fatality rates
  # lhs_sample_young, lhs_old: YLD parameter samples for younger/older ages
  # le_sample: life expectancy sample (for YLL calculation)
  # age_gr: vector of age group labels repeated for each week (e.g., length = number of rows in sim_df)
  # age_gr_levels: factor levels for age_gr
  # prevacc_ui: a data frame with pre-vaccination summary UI (columns: AgeGroup, Week, Median, low95, hi95)
  # n_draws: number of posterior draws to run (default: number of rows in param$gamma)
  
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  n_draws <- length(posterior$gamma)
  
  # set T dynamically
  T = nrow(observed)
  
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    # Prepare a list to store results from each draw
    draw_results <- vector("list", n_draws)
    
    for (d in 1:n_draws) {
      # Extract the d-th draw parameters
      base_beta_draw <- posterior$base_beta[d, ]   # vector length = T
      I0_draw        <- posterior$I0[d, ]          # vector length = A
      gamma_draw     <- posterior$gamma[d]
      rho_draw       <- posterior$rho[d]
      
      sim_out <- sirv_sim_coverageSwitch(
        T = T,
        A = 18,
        N = N,
        r = rep(0, 18),
        base_beta = base_beta_draw,
        I0_draw   = I0_draw,
        R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
        rho = rho_draw,
        gamma = gamma_draw,
        delay = 2,
        VE_block = 0.989,
        coverage_threshold = 1,
        target_age = target,
        total_coverage = 0.5,
        #total_supply = supply[[s]],
        weekly_delivery_speed = 0.1
      )
      
      # Store the age-stratified cases matrix for this draw (dimensions: A x T)
      draw_results[[d]] <- sim_out$age_stratified_cases
    }
    
    # Combine draw results into a 3D array: dimensions [A, T, n_draws]
    age_array <- array(unlist(draw_results), dim = c(18, T, n_draws))
    
    # Compute quantiles for each age group and week over draws
    median_by_age <- apply(age_array, c(1, 2), median)
    low95_by_age  <- apply(age_array, c(1, 2), quantile, probs = 0.025)
    hi95_by_age   <- apply(age_array, c(1, 2), quantile, probs = 0.975)
    
    # Also compute weekly totals (by summing over ages) for each draw
    weekly_totals <- apply(age_array, c(2, 3), sum)  # dimensions: [T, n_draws]
    weekly_cases_median <- apply(weekly_totals, 1, median)
    weekly_cases_low95  <- apply(weekly_totals, 1, quantile, probs = 0.025)
    weekly_cases_hi95   <- apply(weekly_totals, 1, quantile, probs = 0.975)
    
    # Now create a simulation data frame from the median_by_age matrix
    sim_df <- as.data.frame.table(median_by_age, responseName = "Cases")
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      mutate(
        Scenario = s,
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    # Add the computed post-vaccination UI from the quantile matrices
    low95_df <- as.data.frame.table(low95_by_age, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age, responseName = "post_hi95")
    sim_df <- sim_df %>%
      mutate(
        post_low95 = low95_df$post_low95,
        post_hi95  = hi95_df$post_hi95
      )
    
    # Optionally, you might also add the computed weekly UI to a separate weekly summary data frame.
    weekly_df <- data.frame(
      Week = 1:T,
      weekly_median = weekly_cases_median,
      weekly_low95  = weekly_cases_low95,
      weekly_hi95   = weekly_cases_hi95
    )
    
    # Continue with additional processing as in your original function.
    sim_df <- sim_df %>%
      mutate(
        hosp_rate         = rep(hosp, T),
        hospitalised      = Cases * hosp_rate,
        non_hospitalised  = Cases - hospitalised,
        fatality          = rep(fatal, T),
        nh_fatality       = rep(nh_fatal, T),
        fatal             = (hospitalised * fatality + non_hospitalised * nh_fatality)
      ) %>%
      mutate(
        age_numeric = case_when(
          str_detect(age_gr, "<5") ~ 0,
          TRUE ~ as.numeric(str_extract(age_gr, "^\\d+"))
        )
      ) %>%
      # Add YLD parameters (age independent)
      mutate(
        dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Add age-dependent YLD parameters
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
      # Calculate YLD estimates
      mutate(
        yld_acute    = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_total    = yld_acute + yld_subacute + yld_chronic
      ) %>%
      # Assign life expectancy and calculate YLL and DALY metrics
      mutate(
        le_left = case_when(
          age_numeric %in% c(0, 5)    ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10, 15)  ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20, 25)  ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30, 35)  ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40, 45)  ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50, 55)  ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60, 65)  ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70, 75)  ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80, 85)  ~ quantile(le_sample$le_9, 0.5)
        )
      ) %>%
      mutate(
        yll      = fatal * le_left,
        daly_tot = yld_total + yll,
        cum_daly = cumsum(daly_tot)
      )
    
    # If a pre-vaccination UI data frame is provided, merge its selected columns.
    if(!is.null(prevacc_ui)) {
      sim_df <- left_join(sim_df, 
                          dplyr::select(prevacc_ui, AgeGroup, Week, 
                                        pre_Median = Median, 
                                        pre_low95 = low95, 
                                        pre_hi95 = hi95), 
                          by = c("AgeGroup", "Week"))
    }
    
    scenario_result[[s]] <- list(sim_result = list(age_array = age_array,
                                                   weekly_cases_median = weekly_cases_median,
                                                   weekly_cases_low95 = weekly_cases_low95,
                                                   weekly_cases_hi95 = weekly_cases_hi95),
                                 sim_out    = sim_out,
                                 sim_df = sim_df,
                                 weekly_df = weekly_df)
  }
  
  return(scenario_result)
}


### v2.
run_simulation_scenarios_ui <- function(target_age_list, 
                                        observed,
                                        N, bra_foi_state_summ, age_groups, region_name,
                                        hosp, fatal, nh_fatal,
                                        lhs_sample_young, lhs_old, le_sample,
                                        age_gr_levels,
                                        prevacc_ui = NULL,
                                        posterior,
                                        ve_inf = 0
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
  n_draws <- length(posterior$gamma)
  
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Here, each default age label is repeated T times.
  age_gr <- rep(default_age_vector, T)
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    
    pct <- region_coverage[[region_name]]
    total_coverage_frac <- floor(pct) / 100
    
    for (d in 1:n_draws) {
      # Extract d-th draw parameters
      base_beta_draw <- posterior$base_beta[d, ]   # vector length = T
      I0_draw        <- posterior$I0[d, ]          # vector length = 18 (A)
      gamma_draw     <- posterior$gamma[d]
      rho_draw       <- posterior$rho[d]
      sigma_draw     <- posterior$sigma[d]
      
      sim_out <- sirv_sim_coverageSwitch(
        T = T,
        A = 18,
        N = N,
        r = rep(0, 18),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        R0 = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
        rho = rho_draw,
        gamma = gamma_draw,
        delay = 1,
        VE_block = 0.989,
        VE_inf = VE_inf,
        #coverage_threshold = 1,
        target_age = target,
        #total_coverage = 1e6 * 0.2 / sum(N),
        total_coverage = total_coverage_frac,
        #total_coverage = 180000 / sum(N), 
        #cap_supply =  1e6,
        weekly_delivery_speed = 0.1,
        VE_inf_assumed = ve_inf
      )
      
      draw_results_raw_inf[[d]]  <- sim_out$age_stratified_cases_raw  # A x T matrix
      draw_results_raw_symp[[d]] <- sim_out$true_symptomatic  # A x T matrix
      
    }
    
    # Combine draw results into a 3D array [18, T, n_draws]
    age_array_raw_inf  <- array(unlist(draw_results_raw_inf), dim = c(18, T, n_draws))
    age_array_raw_symp <- array(unlist(draw_results_raw_symp), dim = c(18, T, n_draws))
    
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
        non_hospitalised_lo = post_low95 - hospitalised_hi,   # conservative
        non_hospitalised_hi = post_hi95 - hospitalised_lo,
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
          str_detect(as.character(age_gr), "<5") ~ 0,
          TRUE ~ as.numeric(str_extract(as.character(age_gr), "^\\d+"))
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
          age_numeric %in% c(0,5) ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10,15) ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20,25) ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30,35) ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40,45) ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50,55) ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60,65) ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70,75) ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80,85) ~ quantile(le_sample$le_9, 0.5)
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



## v3. age group aligning to licensure
run_simulation_scenarios_ui <- function(target_age_list, 
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
      
      idx <- posterior_idx[d] 
      
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
      base_beta_draw <- posterior$base_beta[idx, ]
      I0_draw        <- posterior$I0[idx, ]
      gamma_draw     <- lhs_sample$gamma[idx]
      rho_draw       <- lhs_sample$rho[idx]
      sigma_draw     <- lhs_sample$sigma[idx]
      FOI_rand_draw  <- lhs_sample$foi[idx]
      ve_rand_draw   <- lhs_sample$ve[idx]
      vc_rand_draw   <- lhs_sample$vc[idx]
      wd_rand_draw   <- lhs_sample$wd[idx]
      
      # for ve x vc simulation
      if (is.na(ve_inf) || ve_inf != 0) {  
        VE_block_draw <- ve_rand_draw
        VE_inf_draw   <- ve_rand_draw
      } else {                            
        VE_block_draw <- 0
        VE_inf_draw   <- 0
      }
      
      sim_out <- sirv_sim_coverageSwitch(
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
        #total_coverage = 0.01,
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
run_simulation_scenarios_ui <- function(target_age_list, 
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
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    this_postwane <- postwane_vec[s]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    
    pct <- region_coverage[[region_name]]
    total_coverage_frac <- floor(pct) / 100
    
    for (d in 1:n_draws) {
      
      idx <- posterior_idx[d] 
      
      # Extract d-th draw parameters ------ calibration based 
      #base_beta_draw <- posterior$base_beta[d, ]   
      #I0_draw        <- posterior$I0[d, ]          
      #gamma_draw     <- posterior$gamma[d]
      #rho_draw       <- posterior$rho[d]
      #sigma_draw     <- posterior$sigma[d]
      
      # 2. predefined idx 
      base_beta_draw <- posterior$base_beta[idx, ]
      I0_draw        <- posterior$I0[idx, ]
      gamma_draw     <- lhs_sample$gamma[idx]
      rho_draw       <- lhs_sample$rho[idx]
      sigma_draw     <- lhs_sample$sigma[idx]
      FOI_rand_draw  <- lhs_sample$foi[idx]
      ve_rand_draw   <- lhs_sample$ve[idx]
      vc_rand_draw   <- lhs_sample$vc[idx]
      wd_rand_draw   <- lhs_sample$wd[idx]
      
      # for ve x vc simulation
      if (is.na(ve_inf) || ve_inf != 0) {  
        VE_block_draw <- ve_rand_draw
        VE_inf_draw   <- ve_rand_draw
      } else {                            
        VE_block_draw <- 0
        VE_inf_draw   <- 0
      }
      
      ve_waning_curve <- function(T = 52, ve_week3, ve_week26) {
        slope <- (ve_week26 - ve_week3) / (26 - 3)
        
        VE_time <- sapply(1:T, function(t) {
          if (t <= 3) {
            ve_week3
          } else if (t <= 26) {
            ve_week3 + slope * (t - 3)
          } else {
            ve_week26 + slope * (t - 26)
          }
        })
        
        return(VE_time)
      }
      
      VE_time <- ve_waning_curve(T = 52,
                                 ve_week3  = ve_rand_draw,
                                 ve_week26 = this_postwane)
      
      sim_out <- sirv_sim_coverageSwitch(
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
        VE_block = VE_time,
        VE_inf = VE_time,
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
                                 weekly_df                    = weekly_df
    )
  }
  
  return(scenario_result)
}



###
postsim_all_ui <- function(scenario_result,   # list of scenario outputs: each with $sim_df and $weekly_df
                           observed,
                           age_gr_levels,
                           pre_summary_cases_age,   # aggregated pre-vacc by age
                           pre_summary_cases,       # pre-vacc data at age-week level (if needed)
                           pre_summary_cases_all,   # aggregated pre-vacc weekly UI data
                           region) { 
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # Number of scenarios
  n_scenarios <- length(scenario_result)
  
  # Extract the age-week simulation data for each scenario
  scenario_df <- lapply(scenario_result, function(x) x$sim_df)
  
  # Build a summary list comparing pre-vaccination vs. post-vaccination cases (by age-week)
  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    # Use the aggregated pre-vacc cases from pre_summary_cases (assumed at the same resolution)
    df$pre_vacc <- pre_summary_cases$Median  
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
  
  # Merge aggregated pre-vacc UI from pre_summary_cases_all into the age-week summary.
  # pre_summary_cases_all is expected to have columns: Week, Median, low95, hi95.
  if(all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_list_df <- left_join(summary_list_df, pre_weekly_df, by = "Week")
  }
  
  # Summarize by AgeGroup (aggregating across weeks)
  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%  
      summarise(
        Median       = sum(Cases, na.rm = TRUE),
        yld_acute    = sum(yld_acute),
        yld_subacute = sum(yld_subacute),
        yld_chronic  = sum(yld_chronic),
        yld_total    = sum(yld_total),
        yll          = sum(yll),
        daly_tot     = sum(daly_tot),
        fatal        = sum(fatal),
        .groups = "drop"
      ) %>% mutate(
        pre_vacc         = pre_summary_cases_age$Median,
        pre_yld_acute    = pre_summary_cases_age$yld_acute,
        pre_yld_subacute = pre_summary_cases_age$yld_subacute,
        pre_yld_chronic  = pre_summary_cases_age$yld_chronic,
        pre_yld_total    = pre_summary_cases_age$yld_total,
        pre_yll          = pre_summary_cases_age$yll,
        pre_daly         = pre_summary_cases_age$daly_tot,
        pre_fatal        = pre_summary_cases_age$fatal,
        diff             = pre_vacc - Median,
        impact           = diff / pre_vacc * 100,
        diff_fatal       = pre_fatal - fatal,
        impact_fatal     = diff_fatal / pre_fatal * 100,
        diff_yld_acute    = pre_yld_acute - yld_acute,
        diff_yld_subacute = pre_yld_subacute - yld_subacute,
        diff_yld_chronic  = pre_yld_chronic - yld_chronic,
        diff_yld_total    = pre_yld_total - yld_total,
        diff_yll          = pre_yll - yll,
        diff_daly         = pre_daly - daly_tot
      )
  })
  
  # Summarize by Week for the entire population
  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>% 
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_cases = sum(Cases),
        pre_cases  = sum(pre_vacc),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal),
        pre_fatal  = sum(pre_fatal),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot),
        pre_daly   = sum(pre_daly),
        diff_daly  = pre_daly - post_daly,
        Scenario   = first(Scenario),
        .groups = "drop"
      )
  })
  
  summary_week_df <- do.call(rbind, summary_week)
  
  # Incorporate the weekly UI from each scenario.
  # Each scenario's weekly_df is stored in scenario_result[[i]]$weekly_df.
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
  
  # Also, merge the aggregated pre-vacc UI from pre_summary_cases_all into the weekly summary.
  if(all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_week_df <- left_join(summary_week_df, pre_weekly_df, by = "Week", relationship = "many-to-many")
  }
  
  summary_week_df$region <- region
  
  # Compute global weekly impact
  global_impact_week <- summary_week_df %>%
    group_by(Scenario) %>%
    summarise(
      total_post_cases = sum(post_cases),
      total_pre_case   = sum(pre_cases),
      total_post_fatal = sum(post_fatal),
      total_pre_fatal  = sum(pre_fatal),
      total_post_daly  = sum(post_daly),
      total_pre_daly   = sum(pre_daly),
      .groups = "drop"
    ) %>%
    mutate(
      diff         = total_pre_case - total_post_cases,
      impact       = diff / total_pre_case * 100,
      diff_fatal   = total_pre_fatal - total_post_fatal,
      impact_fatal = diff_fatal / total_pre_fatal * 100,
      diff_daly    = total_pre_daly - total_post_daly,
      impact_daly  = diff_daly / total_pre_daly * 100 
    )
  
  annotation_text <- paste0(
    global_impact_week$Scenario, ": ", 
    round(global_impact_week$impact, 1), "%",
    collapse = "\n"
  )
  
  global_impact_week$region <- region
  
  # Extract vaccination start and end weeks (if available in sim_result; adjust indices as needed)
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
    vacc_weeks         = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
    )
  ))
}

### v2 -------------------------------------------------------------------------
postsim_all_ui <- function(scenario_result,   # list of scenario outputs: each with $sim_df and $weekly_df
                           observed,
                           age_gr_levels,
                           pre_summary_cases_age,   # aggregated pre-vacc by age; should include columns:
                           # Median, low95, hi95, yld_acute, yld_subacute, 
                           # yld_chronic, yld_total, yll, daly_tot, fatal
                           pre_summary_cases,       # pre-vacc data at age-week level (if needed)
                           pre_summary_cases_all,   # aggregated pre-vacc weekly UI data (with Week, Median, low95, hi95)
                           region) { 
  # Set T dynamically 
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # Number of scenarios
  n_scenarios <- length(scenario_result)
  
  # Extract the age-week simulation data for each scenario
  scenario_df <- lapply(scenario_result, function(x) x$sim_df)
  
  # Build a summary list comparing pre-vaccination vs. post-vaccination cases (by age-week)
  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    # Attach aggregated pre-vaccination values (assumed to be at the same resolution)
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
  
  # Final summarization by AgeGroup: aggregate (sum) over weeks for each scenario.
  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%  
      summarise(
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
  

  # Summarize by Week for the entire population for each scenario.
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
        post_daly  = sum(daly_tot, na.rm = TRUE),
        pre_daly   = sum(pre_daly, na.rm = TRUE),
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
  
  # Compute global weekly impact
  global_impact_week <- summary_week_df %>%
    group_by(Scenario) %>%
    summarise(
      total_post_cases = sum(post_cases, na.rm = TRUE),
      total_pre_case   = sum(pre_cases, na.rm = TRUE),
      total_post_fatal = sum(post_fatal, na.rm = TRUE),
      total_pre_fatal  = sum(pre_fatal, na.rm = TRUE),
      total_post_daly  = sum(post_daly, na.rm = TRUE),
      total_pre_daly   = sum(pre_daly, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      diff         = total_pre_case - total_post_cases,
      impact       = diff / total_pre_case * 100,
      diff_fatal   = total_pre_fatal - total_post_fatal,
      impact_fatal = diff_fatal / total_pre_fatal * 100,
      diff_daly    = total_pre_daly - total_post_daly,
      impact_daly  = diff_daly / total_pre_daly * 100 
    )
  
  annotation_text <- paste0(
    global_impact_week$Scenario, ": ", 
    round(global_impact_week$impact, 1), "%",
    collapse = "\n"
  )
  
  global_impact_week$region <- region
  
  vacc_start_week_s1 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s1   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  
  vacc_start_week_s2 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s2   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  
  vacc_start_week_s3 <- scenario_result[[3]]$sim_out$vacc_start_week[13]
  vacc_end_week_s3   <- scenario_result[[3]]$sim_out$vacc_end_week[13]
  
  return(list(
    scenario_result = scenario_result,
    scenario_data = lapply(scenario_result, function(x) x$sim_result),
    scenario_df = scenario_df,
    summary_list = summary_list,
    summary_list_df = summary_list_df,
    final_summ = final_summ,
    summary_week = summary_week,
    summary_week_df = summary_week_df,
    global_impact_week = global_impact_week,
    annotation_text = annotation_text,
    vacc_weeks = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
    )
  ))
}


### v3 (truesymptomatic -> post-scaling: /rho_p50: original: I[t] * symprate * rho)
postsim_all_ui <- function(scenario_result,   # list of scenario outputs: each with $sim_df and $weekly_df
                           observed,
                           age_gr_levels,
                           pre_summary_cases_age,   # aggregated pre-vacc by age; should include columns:
                           # Median, low95, hi95, yld_acute, yld_subacute, 
                           # yld_chronic, yld_total, yll, daly_tot, fatal
                           pre_summary_cases,       # pre-vacc data at age-week level (if needed)
                           pre_summary_cases_all,   # aggregated pre-vacc weekly UI data (with Week, Median, low95, hi95)
                           region) { 
  # Set T dynamically 
  T <- nrow(observed)
  # Define default age groups
  default_age_vector <- c("<5 years","5-9 years","10-14 years","15-19 years",
                          "20-24 years","25-29 years","30-34 years","35-39 years",
                          "40-44 years","45-49 years","50-54 years","55-59 years",
                          "60-64 years","65-69 years","70-74 years","75-79 years",
                          "80-84 years","85-89 years")
  # Age-week vector
  age_gr <- rep(default_age_vector, T)
  n_scenarios <- length(scenario_result)
  scenario_df  <- lapply(scenario_result, function(x) x$sim_df)
  
  # Build summary_list (age-week)
  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    # Attach aggregated pre-vaccination values (assumed to be at the same resolution)
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
             pre_vacc = pre_summary_cases$Median) %>% 
      group_by(Week, Scenario) %>%
      summarise(
        post_cases = sum(Cases, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
        pre_daly   = sum(pre_daly, na.rm = TRUE),
        diff_daly  = pre_daly - post_daly,
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
        post_daly  = sum(daly_tot, na.rm = TRUE),
        pre_daly   = sum(pre_daly, na.rm = TRUE),
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
        c(post_cases, pre_cases, diff, post_fatal, pre_fatal, diff_fatal,
          post_daly, pre_daly, diff_daly, post_weekly_median, post_weekly_low95,
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
      .groups    = "drop"
    ) %>%
    dplyr::mutate(impact = (total_pre - total_post)/total_pre *100)
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
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
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
             pre_vacc_hi = pre_summary_cases$hi95
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
        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
        pre_daly   = sum(pre_daly, na.rm = TRUE),
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
        c(post_cases, post_cases_lo, post_cases_hi, pre_cases, pre_cases_lo, pre_cases_hi, diff, post_fatal, pre_fatal, diff_fatal,
          post_daly, pre_daly, diff_daly, post_weekly_median, post_weekly_low95,
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


## target week
target_week <- function(postsim_all,
                        prevacc_ui,
                        region){
  # 1. 시나리오별 타겟 연령 정의
  age_targets <- list(
    Scenario_1 = 1:4,
    Scenario_2 = 5:12,
    Scenario_3 = 13:18
  )
  
  # 2. 시나리오별 + 연령별 + 주별 집계
  summary_age_week <- lapply(seq_along(postsim_all$summary_list), function(i) {
    df <- postsim_all$summary_list[[i]]
    
    df <- df %>%
      mutate(
        Scenario        = paste0("Scenario_", i),
        pre_vacc_cases  = prevacc_ui$Median,
        pre_vacc_lo     = prevacc_ui$low95,
        pre_vacc_hi     = prevacc_ui$hi95,
      )
    
    df_sum <- df %>%
      group_by(Week, AgeGroup, Scenario) %>%
      summarise(
        post_cases         = sum(Cases, na.rm = TRUE),
        post_weekly_low95  = sum(post_low95, na.rm = TRUE),
        post_weekly_hi95   = sum(post_hi95,  na.rm = TRUE),
        pre_cases          = sum(pre_vacc_cases, na.rm = TRUE),
        pre_cases_lo       = sum(pre_vacc_lo, na.rm = TRUE),   # optional
        pre_cases_hi       = sum(pre_vacc_hi, na.rm = TRUE),   # optional
        .groups = "drop"
      )
    
    return(df_sum)
  }) %>% bind_rows()
  
  # 3. 타겟 연령 플래그 추가 및 타겟 연령군만 필터
  summary_age_week <- summary_age_week %>%
    mutate(target = map2_lgl(Scenario, AgeGroup, ~ .y %in% age_targets[[.x]])) %>%
    filter(target)
  
  # 4. 시나리오별 주간 총합 계산
  summary_target_week <- summary_age_week %>%
    group_by(Week, Scenario) %>%
    summarise(
      post_cases = sum(post_cases, na.rm = TRUE),
      post_cases_lo = sum(post_weekly_low95),
      post_cases_hi = sum(post_weekly_hi95),
      pre_cases  = sum(pre_cases,  na.rm = TRUE),
      pre_cases_lo = sum(pre_cases_lo),
      pre_cases_hi = sum(pre_cases_hi),
      .groups    = "drop",
    )
  
  summary_target_week$region <- region

  summary_target_week <- summary_target_week %>%
    left_join(rho_df, by = "region") %>%
    mutate(
      across(
        c(post_cases, post_cases_lo, post_cases_hi, pre_cases, pre_cases_lo, pre_cases_hi),
        ~ .x / rho_p50
      )
    )
  return(summary_target_week)
}
## only direct group
epi_graph_direct <- function(postsim_all, prevacc_ui, observed) {
  
  # 1. 시나리오별 타겟 연령 정의
  age_targets <- list(
    Scenario_1 = 1:4,
    Scenario_2 = 5:12,
    Scenario_3 = 13:18
  )
  
  # 2. 시나리오별 + 연령별 + 주별 집계
  summary_age_week <- lapply(seq_along(postsim_all_ce_ui$summary_list), function(i) {
    df <- postsim_all_ce_ui$summary_list[[i]]
    
    df <- df %>%
      mutate(
        Scenario        = paste0("Scenario_", i),
        pre_vacc_cases  = prevacc_ui_ce$Median,
        pre_vacc_lo     = prevacc_ui_ce$low95,
        pre_vacc_hi     = prevacc_ui_ce$hi95,
      )
    
    df_sum <- df %>%
      group_by(Week, AgeGroup, Scenario) %>%
      summarise(
        post_cases         = sum(Cases, na.rm = TRUE),
        post_weekly_low95  = sum(post_low95, na.rm = TRUE),
        post_weekly_hi95   = sum(post_hi95,  na.rm = TRUE),
        pre_cases          = sum(pre_vacc_cases, na.rm = TRUE),
        pre_cases_lo       = sum(pre_vacc_lo, na.rm = TRUE),   # optional
        pre_cases_hi       = sum(pre_vacc_hi, na.rm = TRUE),   # optional
        .groups = "drop"
      )
    
    return(df_sum)
  }) %>% bind_rows()
  
  # 3. 타겟 연령 플래그 추가 및 타겟 연령군만 필터
  summary_age_week <- summary_age_week %>%
    mutate(target = map2_lgl(Scenario, AgeGroup, ~ .y %in% age_targets[[.x]])) %>%
    filter(target)
  
  # 4. 시나리오별 주간 총합 계산
  summary_target_week <- summary_age_week %>%
    group_by(Week, Scenario) %>%
    summarise(
      post_cases = sum(post_cases, na.rm = TRUE),
      post_cases_lo = sum(post_weekly_low95),
      post_cases_hi = sum(post_weekly_hi95),
      pre_cases  = sum(pre_cases,  na.rm = TRUE),
      pre_cases_lo = sum(pre_cases_lo),
      pre_cases_hi = sum(pre_cases_hi),
      .groups    = "drop",
    )
  
  # 5. 백신 주차 정보 + 최대 케이스 계산
  vacc_start_week <- postsim_all$vacc_weeks$scenario1$start
  vacc_end_week   <- postsim_all$vacc_weeks$scenario1$end
  
  scenario_max_df <- summary_target_week %>%
    group_by(Scenario) %>%
    summarise(max_cases = max(post_cases, na.rm = TRUE), .groups = "drop")
  
  # summary_target_week에 join
  summary_target_week <- summary_target_week %>%
    left_join(scenario_max_df, by = "Scenario") %>%
    mutate(
      ribbon_ymin = case_when(
        Scenario == "Scenario_1" ~ 2/3 * max_cases,
        Scenario == "Scenario_2" ~ 1/3 * max_cases,
        Scenario == "Scenario_3" ~ 0
      ),
      ribbon_ymax = case_when(
        Scenario == "Scenario_1" ~ max_cases,
        Scenario == "Scenario_2" ~ 2/3 * max_cases,
        Scenario == "Scenario_3" ~ 1/3 * max_cases
      ),
      annotation_y = 0.9 * max_cases
    )
  
  summary_target_week <- summary_target_week %>%
    left_join(rho_df, by = "region") %>%
    mutate(
      across(
        c(post_cases, post_cases_lo, post_cases_hi, pre_cases, pre_cases_lo, pre_cases_hi),
        ~ .x / rho_p50
      )
    )
  
  global_impact_all <- summary_target_week %>%
    group_by(Scenario) %>%
    summarise(
      total_post_cases = sum(post_cases),
      total_pre_case   = sum(pre_cases),
      .groups = "drop"
    ) %>%
    mutate(
      diff   = total_pre_case - total_post_cases,
      impact = diff / total_pre_case * 100
    )
  
  annotation_df <- global_impact_all %>%
    mutate(
      Strategy = dplyr::recode(Scenario,
                               "Scenario_1" = "Strategy 1",
                               "Scenario_2" = "Strategy 2",
                               "Scenario_3" = "Strategy 3"),
      ann = paste0(Strategy, ": ", round(impact, 1), "%")
    ) %>%
    select(Scenario, ann)
  
  # Create an annotation data frame that has one row per region
  vaccine_ann_df <- global_impact_all  %>%
    mutate(x = 4, 
           y = overall_max * 0.95, 
           label = "<----- Vaccine impact start (+ 2 wks after initiation)")
  
  # 6. 시각화
  g <- 
    ggplot(summary_target_week) +
    geom_ribbon(data = summary_target_week %>% filter(Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
                aes(x = Week, ymin = 0, ymax = Inf), fill = "grey70",
                alpha = 0.3)  +
    geom_ribbon(aes(x = Week, ymin = post_cases_lo, ymax = post_cases_hi, fill = Scenario),
                alpha = 0.3) +
    geom_line(aes(x = Week, y = pre_cases, linetype = "Without vaccination"), color = "black", linewidth = 0.4) +
    geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "With vaccination"), linewidth = 0.4) + 
    geom_ribbon(aes(x = Week, ymin = pre_cases_lo, ymax = pre_cases_hi),
                fill = "lightgrey", alpha = 0.6) +
    geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
    geom_text(
      data = summary_target_week %>%
        filter(Week == 4) %>%
        distinct(Scenario, annotation_y),
      mapping = aes(x = 4, y = annotation_y, label = "<----- Vaccine impact start (+2 wks)"),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 3
    ) +
    facet_wrap(~Scenario) +
    geom_text(data = annotation_df, 
              aes(x = Inf, y = Inf, label = ann), 
              inherit.aes = FALSE,
              hjust = 1.1, 
              vjust = 1.1, 
              size = 3) +
    theme(plot.margin = margin(5, 30, 5, 5)) +
    #geom_point(data = observed_all, aes(x = Week, y = Observed),size = 0.05) + 
    facet_wrap(~Scenario, scales = "free_y")+
    theme_pubclean()+
    theme(legend.position = "right")
  
  
  return(g)
}


## vaccine allocations ----------------------------------------------
vacc_allocation <- function(postsim_all_ui, observed, region) {
  
  T = nrow(observed)
    
  scenario_data <- lapply(seq_along(postsim_all_ui$scenario_result), function(idx){
    list <- postsim_all_ui$scenario_result[[idx]]
    data <- list$sim_out
    return(data)
  })
  
  raw_allocation_week <- lapply(scenario_data, function(scenario){
    scenario$raw_allocation_age
  })
  scenario_names <- c("Scenario1", "Scenario2", "Scenario3", "Scenario4")
  
  weekly_allocation_df <- do.call(rbind, lapply(seq_along(raw_allocation_week), function(idx) {
    df <- as.data.frame(raw_allocation_week[[idx]])  # Convert matrix to dataframe
    df$Scenario <- scenario_names[idx]               # Assign scenario label
    df$AgeGroup <- 1:nrow(df)                        # Add age group identifier
    df                                               # Return dataframe
  }))
  
  colnames(weekly_allocation_df)[1:(ncol(weekly_allocation_df)-2)] <- as.character(1:(ncol(weekly_allocation_df)-2))
  
  weekly_allocation_df$tot_sum <- rowSums(weekly_allocation_df[,1:T])
  
  weekly_allocation_long <- weekly_allocation_df %>%
    pivot_longer(
      cols = all_of(as.character(1:T)),  # Use the character names of the week columns
      names_to = "Week",
      values_to = "Vaccinated"
    ) %>%
    mutate(
      Week = as.numeric(Week)  # Convert week names to numeric for plotting
    )%>%   # Ensure ordering is correct
    mutate(age_gr = rep(age_gr[1:length(age_gr_levels)], each = T, times = n_scenarios)) %>%
    mutate(age_gr = factor(age_gr, levels = unique(age_gr)))
  
  weekly_allocation_long$region <- region
  
  paired <- brewer.pal(12, "Paired")
  extra <- brewer.pal(8, "Dark2")[1:6]
  full_palette <- c(paired, extra)
  
  total_allocations <- weekly_allocation_long %>% group_by(Scenario) %>% summarise(
    tot_sum = sum(tot_sum)
  )
  
  return(list(
    scenario_data      = scenario_data,
    raw_allocation_week = raw_allocation_week,
    weekly_allocation_df = weekly_allocation_df,
    weekly_allocation_long = weekly_allocation_long,
    total_allocations = total_allocations
  ))
  
}


## nnv  ----------------------------------------------------
nnv_list <- function(vacc_allocation,
                     postsim_all_ui,
                     N,
                     region,
                     observed){
  
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
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
  
  
  final_summ_df <- do.call(rbind, lapply(seq_along(postsim_all_ui$final_summ), function(idx){
    df <- postsim_all_ui$final_summ[[idx]]
    df <- df %>% mutate(
      tot_vacc = as.numeric(unlist(raw_allocation[[idx]])) 
    )
    df$scenario <- paste0("Scenario_", idx)
    df <- df %>% mutate(
      nnv = tot_vacc / diff,
      nnv_fatal = tot_vacc / diff_fatal,
      nnv_yld_acute    =  tot_vacc / diff_yld_acute,
      nnv_yld_subacute = tot_vacc / diff_yld_subacute,
      nnv_yld_chronic = tot_vacc / diff_yld_chronic,
      nnv_yld_total =tot_vacc / diff_yld_total,
      nnv_yll = tot_vacc / diff_yll,
      nnv_daly = tot_vacc / diff_daly
    )
    colnames(df)[colnames(df) == "tot_vacc$tot_vacc"] <- "tot_vacc"
    colnames(df)[colnames(df) == "nnv$tot_vacc"] <- "nnv"
    return(df)
  }))
  
  final_summ_df <- final_summ_df %>% mutate(
    tot_pop = rep(N, n_scenarios),
    pre_vacc_10k = pre_vacc / tot_vacc * 100000,
    averted_per10k = diff / tot_vacc * 100000,
    tot_vacc_per10k = tot_vacc / tot_pop,
    nnv_per = tot_vacc_per10k / averted_per10k
  )
  
  final_summ_df$age_gr <- rep(age_gr[1:18],n_scenarios)
  final_summ_df$age_gr <- factor(final_summ_df$age_gr, levels = age_gr_levels)
  final_summ_df$region <- region
  
  return(list(
    raw_allocation     = raw_allocation,
    raw_allocation_age = raw_allocation_age,
    final_summ_df      = final_summ_df
  ))

}


### v3
nnv_list <- function(vacc_allocation,
                     postsim_all_ui,
                     N,
                     region,
                     observed){
  
  # Set T dynamically 
  T = nrow(observed)
  
  n_scenarios <- length(target_age_list)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
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
  
  
  final_summ_df <- do.call(rbind, lapply(seq_along(postsim_all_ui$final_summ), function(idx){
    
    # attach target col
    target <- c("<20 years", "<20 years", "<20 years", "<20 years",
                "20-59 years", "20-59 years", "20-59 years", "20-59 years",
                "20-59 years", "20-59 years", "20-59 years", "20-59 years", 
                ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years")

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
        # 2) age group당 averted diff 합계
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
      mutate(
        scenario_vacc = raw_allocation_age[[idx]]
      ) %>%
      mutate(
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
    )

  })
  )
  
  # —— ② scenario × target별 per1M 요약 —— 
  per1M_summary <- final_summ_df %>%
    group_by(scenario, target) %>%
    summarise(
      scenario_vacc = first(scenario_vacc),
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
      per1M_daly_hi    = diff_daly_hi    / scenario_vacc * 1e6
    )
  
  final_summ_df <- final_summ_df %>%
    mutate(
      tot_pop        = rep(N, n_scenarios),
      tot_vacc_prop  = tot_vacc / tot_pop
    )
  
  final_summ_df$age_gr <- rep(age_gr[1:18],n_scenarios)
  final_summ_df$age_gr <- factor(final_summ_df$age_gr, levels = age_gr_levels)
  final_summ_df$region <- region
  
  return(list(
    raw_allocation     = raw_allocation,
    raw_allocation_age = raw_allocation_age,
    final_summ_df      = final_summ_df,
    per1M_summary      = per1M_summary
  ))
  
}

### v4 (add rho df)

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
    "<5 years","5-9 years","10-14 years","15-19 years",
    "20-24 years","25-29 years","30-34 years","35-39 years",
    "40-44 years","45-49 years","50-54 years","55-59 years",
    "60-64 years","65-69 years","70-74 years","75-79 years",
    "80-84 years","85-89 years"
  )
  
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
    target <- c("<20 years", "<20 years", "<20 years", "<20 years",
                "20-59 years", "20-59 years", "20-59 years", "20-59 years",
                "20-59 years", "20-59 years", "20-59 years", "20-59 years", 
                ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years")
    
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
        ) 
      ) %>%
      mutate(
        scenario_vacc = raw_allocation_age[[idx]]
      ) %>%
      mutate(
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
    )
    
  })
  )
  
  final_summ_df <- final_summ_df %>%
    mutate(
      tot_pop        = rep(N, n_scenarios),
      tot_vacc_prop  = tot_vacc / tot_pop
    )
  
  final_summ_df$age_gr <- rep(age_gr[1:18],n_scenarios)
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

## nnv graph ----------------------------------------------------


nnv_gg <- function(final_nnv_df,
                   y_var = "nnv_fatal",  # name of y variable as a string
                   y_lab = "NNV to avert a single fatal case",
                   x_lab = "Age group",
                   title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2") {
  
  # Convert the y_var string to a symbol for tidy evaluation
  y_sym <- sym(y_var)
  
  g <- ggplot(final_nnv_df) +
    geom_bar(aes(x = age_gr, y = !!y_sym, fill = scenario),
             stat = "identity", alpha = 0.7) +
    facet_wrap(~scenario, nrow = 3, 
               labeller = as_labeller(c("Scenario_1" = "<20 years only", 
                                        "Scenario_2" = "20-59 years only", 
                                        "Scenario_3" = ">60 years only"))) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "<20 years only", 
                                 "Scenario_2" = "20-59 years only", 
                                 "Scenario_3" = ">60 years only")) +
    scale_y_continuous(labels = comma) +
    theme_minimal() +
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
                       y_lab = "NNV to avert a single fatal case",
                       x_lab = "Age group",
                       title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2") {
  
  combined_df <- combined_df %>% filter(nnv !=0)
  
  # Convert the y_var string to a symbol for tidy evaluation
  y_sym <- sym(y_var)
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))

  setting_order <- c("National", "High", "Moderate", "Low")  # case-sensitive
  
  combined_df <- combined_df %>%
    mutate(
      setting = factor(setting, levels = setting_order),
    )
  
  g <- 
    #ggplot(combined_df)+
    #geom_bar(aes(x = age_gr, y = !!y_sym, fill = scenario), stat = "identity", alpha = 0.7)+
    
    ggplot(combined_df, aes(x = age_gr, y = !!y_sym, fill = scenario)) +
    geom_bar(stat = "identity", alpha = 0.5)+
    geom_errorbar(aes(ymin = !!lower_sym, ymax = !!upper_sym),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    facet_grid(~setting, 
               labeller = labeller(
                 scenario = c("Scenario_1" = "1-19 years", 
                              "Scenario_2" = "20-59 years", 
                              "Scenario_3" = "≥60 years")
               ))+
    scale_x_discrete(labels = function(x) gsub("years", "", x)) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "1-19 years", 
                                 "Scenario_2" = "20-59 years", 
                                 "Scenario_3" = "≥60 years"),
                      name = "Vaccination strategy")+
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
## v2
nnv_nat_setting <- function(combined_df, 
                            y_var,
                            y_lab = "NNV to avert a single symptomatic case",
                            x_lab = "Vaccination strategy",
                            title = "Coverage:50%, Delivery speed:10%, Deployment: Week 2") {
  
  # tidy eval 준비
  y_sym     <- sym(y_var)
  lower_sym <- sym(paste0(y_var, "_lo"))
  upper_sym <- sym(paste0(y_var, "_hi"))
  
  combined_df <- combined_df %>%
    filter(!is.na(target)) %>%
    mutate(
      setting  = factor(setting,  levels = c("National","High","Moderate","Low")),
      scenario = factor(scenario, levels = c("Scenario_1","Scenario_2","Scenario_3"),
                        labels = c("1–19 years","20–59 years","≥60 years"))
    )
  
  ggplot(combined_df, aes(x = scenario, y = !!y_sym, fill = target)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = !!lower_sym, ymax = !!upper_sym),
                  position = position_dodge(width = 0.8),
                  width = 0.2, color = "black") +
    facet_wrap(~ setting, nrow = 1) +
    scale_fill_brewer("Age target", palette = "Set2") +
    scale_y_continuous(labels = comma) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )
}
# v3. 
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

make_per1m_plot <- function(df,
                            y_var,      # ex: "per1M_diff"
                            y_lo,       # ex: "per1M_diff_lo"
                            y_hi,       # ex: "per1M_diff_hi"
                            y_label,    # ex: "Cases averted per 1 M doses"
                            x_label     = "Vaccination strategy",
                            legend_title = "Type"
) {
  # tidy-eval 준비
  y_sym    <- sym(y_var)
  ylo_sym  <- sym(y_lo)
  yhi_sym  <- sym(y_hi)
  
  # 공통 dodge 설정
  pd <- position_dodge(width = 0.8)
  
  custom_labeller <- function(x) {
    # 1,000 이상인 경우만 k 붙이고 반올림
    labs <- ifelse(
      abs(x) >= 1000,
      paste0(round(x / 1000), "k"),
      as.character(x)
    )
    labs
  }
  g1 <- ggplot(df, aes(x = scenario, y = !!y_sym, fill = effect_group)) +
    geom_col(position = pd, width = 0.8, alpha = 0.8) +
    geom_errorbar(aes(ymin = !!ylo_sym, ymax = !!yhi_sym),
                  position = pd, width = 0.2) +
    facet_wrap(~ setting, nrow = 1, strip.position = "top") +
    scale_fill_brewer(palette = "Set2", name = legend_title) +
    scale_y_continuous(labels = custom_labeller) +
    labs(x = x_label, y = y_label) +
    theme_pubclean() +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95", colour = NA),
      axis.text.x      = element_text(angle = 90, hjust = 1)
    )
  
  g1 <- g1 + shared_theme + theme(
    legend.position = "right",
    plot.title = element_blank(),
    legend.text     = element_text(size = 8),  
    legend.title    = element_text(size = 9),
    strip.text      = element_blank(),
    strip.background  = element_blank(),
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
  ) 
  
  return(g1)
  
}

## v2.1 tracking indirect effect
make_per1m_plot <- function(df, y_var, y_lo, y_hi,
                            y_label,
                            x_label = "Vaccination strategy",
                            legend_title = "Effect type") {
  
  y_sym   <- rlang::sym(y_var)
  ylo_sym <- rlang::sym(y_lo)
  yhi_sym <- rlang::sym(y_hi)
  
  level_order <- c(
    "Total effect in target group",
    "Indirect — <1 (novacc)",
    "Indirect — 1-11 years",
    "Indirect — 12-17 years",
    "Indirect — 18-59 years",
    "Indirect — 60+ years"
  )
  
  df <- df %>%
    mutate(effect_group = factor(effect_group, levels = level_order))
  
  plot_summary <- df %>%
    group_by(setting, scenario, effect_group) %>%
    summarise(
      mid_val = sum(!!y_sym,   na.rm = TRUE),
      lo_val  = sum(!!ylo_sym, na.rm = TRUE),
      hi_val  = sum(!!yhi_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  error_df <- plot_summary %>%
    group_by(setting, scenario) %>%             # 스택 꼭대기 찾기
    summarise(
      y    = sum(mid_val),
      ymin = sum(lo_val),
      ymax = sum(hi_val),
      .groups = "drop"
    )
  
  ggplot(plot_summary,
         aes(x = scenario, y = mid_val, fill = effect_group)) +
    geom_col(position = "stack", width = 0.8) +
    geom_errorbar(
      data        = error_df,
      inherit.aes = FALSE,
      aes(x = scenario, y = y, ymin = ymin, ymax = ymax),
      width  = 0.25, colour = "black") +
    facet_wrap(~ setting, nrow = 1, strip.position = "top") +
    scale_fill_brewer(
      palette = "Set2",
      name    = legend_title,
      breaks  = level_order,
      labels  = level_order
    ) +
    scale_y_continuous(
      labels = function(x) ifelse(abs(x) >= 1000,
                                  paste0(round(x/1000), "k"), x)) +
    labs(x = x_label, y = y_label) +
    theme_pubclean() +
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}



## v2.
make_per1m_plot <- function(df,
                            y_var,         # 예: "per1M_diff"
                            y_lo,          # 예: "per1M_diff_lo"
                            y_hi,          # 예: "per1M_diff_hi"
                            y_label,       # 예: "Cases averted per 1M doses"
                            x_label      = "Vaccination strategy",
                            legend_title = "Type"
) {
  # 1) tidy evaluation을 위한 심볼 준비
  y_sym   <- sym(y_var)
  ylo_sym <- sym(y_lo)
  yhi_sym <- sym(y_hi)
  
  custom_labeller <- function(x) {
    # 1,000 이상인 경우만 k 붙이고 반올림
    labs <- ifelse(
      abs(x) >= 1000,
      paste0(round(x / 1000), "k"),
      as.character(x)
    )
    labs
  }
  
  level_order <- c(
    "Indirect — >60 years",
    "Indirect — 20-59 years",
    "Indirect — <20 years",
    "Total effect in target group"
  )
  
  df <- df %>% 
    mutate(effect_group = factor(effect_group, levels = rev(level_order)))%>%
    # keep only the total‐effect slice (all others are zero)
    filter(effect_group == "Total effect in target group")
  
  # 3) 세그먼트별 중앙값(mid_val) 계산
  plot_summary <- df %>%
    group_by(setting, scenario) %>%
    summarise(
      mid_val = sum(!!y_sym,   na.rm = TRUE),
      lo_val  = sum(!!ylo_sym, na.rm = TRUE),
      hi_val  = sum(!!yhi_sym, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # 4) 전체(총효과) 에러바용 테이블: “Total effect”만 필터링
  error_df <- plot_summary %>%
    #filter(effect_group == "Total effect in target group") %>%
    group_by(setting, scenario) %>%
    summarise(
      y       = cumsum(mid_val),       # 스택 전체 높이
      ymin    = cumsum(lo_val),
      ymax    = cumsum(hi_val),
      .groups = "drop"
    )
  
  plot_summary <- df %>%
    group_by(setting, scenario, effect_group) %>%
    summarise(
      mid_val = sum(!!y_sym,   na.rm = TRUE),
      lo_val  = sum(!!ylo_sym, na.rm = TRUE),
      hi_val  = sum(!!yhi_sym, na.rm = TRUE),
      .groups  = "drop"
    )
  # 5) 그리기
  ggplot(plot_summary, 
         aes(x = scenario, y = mid_val, fill = scenario)) +
    # (a) 스택형 막대
    geom_col(position = "stack", width = 0.8, alpha = 0.5) +
    # (b) 총효과 에러바: 스택 꼭대기(vjust=1)에 고정
    geom_errorbar(
      data        = error_df,
      inherit.aes = F,
      aes(x = scenario,
          y = y,
          ymin = ymin,
          ymax = ymax),
      #position = position_stack(vjust = 1),
      width    = 0.2,
      colour   = "black"
    ) +
    # (c) 페싯 & 범례 & 테마
    facet_wrap(~ setting, nrow = 1, strip.position = "top") +
    scale_fill_brewer(palette = "Set1", name = legend_title, 
                      labels = c("Strategy 1 (1-11 yrs)", "Strategy 2 (12-17 yrs)","Strategy 3 (18-59 yrs)", "Strategy 4 (60+ yrs)"
    )) +
    scale_y_continuous(labels = custom_labeller) +
    labs(x = x_label, y = y_label) +
    theme_pubclean() +
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

combine_nnv_per1m <-function(nnv_plot, per1m_plot, tag = "A"){
  
  custom_labeller <- function(x) {
    v <- abs(x)
    ifelse(v >= 1e6,
           paste0(round(x / 1e6, 1), "M"),
           ifelse(v >= 1e3,
                  paste0(round(x / 1e3, 1), "k"),
                  as.character(x)
           )
    )
  }
  
  shared_theme  <- theme_pubclean() +
    theme(
      axis.text.y     = element_text(size = 10),
      axis.title.y    = element_text(size = 9),
      legend.text    = element_text(size = 8),
      legend.title   = element_text(size = 9),
      plot.title      = element_text(size = 9),
      plot.tag       = element_text(size = 10, face = "bold")
    )
  
  p1 <- nnv_plot + shared_theme + scale_y_continuous(labels = custom_labeller) +
    shared_theme +
    labs(tag = tag) +
    theme(
      axis.title.x  = element_blank(),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      plot.title    = element_blank(),
      legend.position = "right"
    )
  
  g1 <- per1m_plot + shared_theme + scale_y_continuous(labels = custom_labeller) +
    shared_theme +
    #labs(tag = "B") +
    theme(
      strip.text        = element_blank(),
      strip.background  = element_blank(),
      # x-axis는 세로 회전 유지 (필요 시)
      axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
      plot.title        = element_blank(),
      legend.position   = "right"
    )+
    guides(fill = "none", colour = "none") 
    
  combined <-  (p1) / 
    (g1) +
    plot_layout(guides="collect") &
   theme(legend.position="right")
  
  return(combined)
}

# risk plot
compute_nnv_summary <- function(df_risk,
                                y_scale_factor = 1e6) {
  df_long <- df_risk %>%
    mutate(
      cases_averted_direct = if_else(is_target_group, cases_averted_direct, 0),
      cases_averted_direct_lo = if_else(is_target_group, cases_averted_direct_lo, 0),
      cases_averted_direct_hi = if_else(is_target_group, cases_averted_direct_hi, 0)
    ) %>%
    select(scenario, target, is_target_group,
           cases_averted_direct, cases_averted_direct_lo, cases_averted_direct_hi,
           cases_averted_indirect, cases_averted_indirect_lo, cases_averted_indirect_hi,
           tot_vacc
    ) %>%
    pivot_longer(
      cols = starts_with("cases_averted_"),
      names_to = "effect_full",
      values_to = "value"
    ) %>%
    # 4) effect_type, suffix 분리
    mutate(
      effect_type = case_when(
        grepl("indirect", effect_full)          ~ if_else(is_target_group,
                                                          "Indirect — target age",
                                                          paste0("Indirect — ", target)),
        grepl("cases_averted_direct", effect_full) ~ "Direct effect",
        TRUE                                     ~ NA_character_
      ),
      suffix = case_when(
        grepl("_lo$", effect_full) ~ "lo",
        grepl("_hi$", effect_full) ~ "hi",
        TRUE                       ~ "med"
      )
    ) %>%
    # 5) 이제 effect_full은 필요 없으니 제거
    select(-effect_full)
  
  df_summary <- df_long %>%
    # 1) 시나리오별 전체 tot_vacc 계산
    group_by(scenario) %>%
    mutate(tot_vacc_scenario = max(tot_vacc, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # 2) scenario × effect_type × suffix 단위로 value 합산
    group_by(scenario, effect_type, suffix, tot_vacc_scenario) %>%
    summarise(sum_value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    
    # 3) 1M FVP당 값으로 변환
    mutate(per1M = sum_value / tot_vacc_scenario * 1e6) %>%
    
    # 4) suffix (med/lo/hi) 기준으로 와이드 포맷
    pivot_wider(
      names_from  = suffix,
      values_from = per1M,
      values_fill = list(per1M = 0)
    )
  
  df_plot <- df_summary %>%
    # 만약 중복된 행이 있다면 summarize 로 한 행으로 묶어주세요.
    # (suffix 단계에서 이미 묶지 않으셨다면 여기서 한 번 더)
    group_by(scenario, effect_type) %>%
    summarise(
      med = sum(med),
      lo  = sum(lo),
      hi  = sum(hi),
      .groups="drop"
    ) %>%
    mutate(
      effect_type = factor(
        effect_type,
        levels = c(
          "Direct effect",
          "Indirect — target age",
          "Indirect — <20 years",
          "Indirect — 20-59 years",
          "Indirect — >60 years"
        )
      )
    ) 
  
  return(df_plot)
}

## vaccine allocation graph ----------------------------------------

vacc_alloc_graph <- function(vacc_alloc_list){
  
  paired <- brewer.pal(12, "Paired")
  extra <- brewer.pal(8, "Dark2")
  full_palette <- c(paired, extra)
  
  total_allocations <- vacc_alloc_list$weekly_allocation_long %>% group_by(Scenario) %>% summarise(
    tot_sum = sum(tot_sum)
  )
  
  p <- 
    ggplot(vacc_alloc_list$weekly_allocation_long, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
    geom_bar(stat = "identity") +   # Bars are automatically stacked
    labs(
      x = "Week",
      y = "Total doses of vaccines distributed per age group",
      fill = "Age Group"
    ) +
    scale_fill_manual(values = full_palette) + 
    scale_y_continuous(label=comma)+
    facet_wrap(~ Scenario, 
               labeller = labeller(Scenario = c("Scenario1" = "1-11 years only", 
                                                "Scenario2" = "12-17 years only", 
                                                "Scenario3" = "18-59 years only",
                                                "Scenario4" = "60+ years only"
                                                )))  +         # Facet by Scenario if you want separate plots for each
    theme_minimal()
  
  return(p)
}


## case region map -------------------------------------------------
region_prepost_map <- function(post_case, title){
  
  overall_range <- range(c(br_states$case_22, br_states$post_case), na.rm = TRUE)
 
  g1 <-  
   ggplot(br_states) +
    geom_sf(aes(fill = case_22)) +
    scale_fill_distiller(
      palette = "Reds", 
      direction = 1, 
      name = "Symptomatic cases (2022)",
      limits = overall_range  # sets the same scale for all maps
    ) +
    labs(title = "Pre-vaccination symptomatic cases (2022)",
         size  = 5) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_sf_text(aes(label = region), size = 2, color = "black")
  
  g2 <- 
    ggplot(br_states) +
    geom_sf(aes(fill = .data[[post_case]])) +
    scale_fill_distiller(
      palette = "Reds", 
      direction = 1, 
      name = "Symptomatic cases (2022)",
      limits = overall_range  # sets the same scale for all maps
    ) +
    labs(title = title,
         size  = 5) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_sf_text(aes(label = region), size = 2, color = "black")
  
  combined_plot <- ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
  
  return(combined_plot)
}


## R0 estimation -------------------------------------------------------
r0_estim <- function(posterior){
  
  base_beta_draw <- posterior$base_beta
  beta_ui <- apply(posterior$base_beta, 2, quantile, probs = c(0.025, 0.5, 0.975))
  
  beta_df <- as.data.frame(t(beta_ui))
  colnames(beta_df) <- c("low", "median", "hi")
  
  # Create a time variable (assuming time points 1:T)
  beta_df$time <- 1:nrow(beta_df)
  
  gamma_draw <- posterior$gamma
  gamma_mat <- matrix(gamma_draw, nrow = nrow(base_beta_draw), ncol = ncol(base_beta_draw), byrow = FALSE)
  
  r0 <- base_beta_draw / gamma_mat
  r0_ui <- apply(r0, 2, quantile, probs = c(0.025, 0.5, 0.975))
  r0_df <- as.data.frame(t(r0_ui))
  colnames(r0_df) <- c("low", "median", "hi")
  r0_df$time <- 1:nrow(r0_df)
  
  return(list(
    beta_ui = beta_ui,
    beta_df = beta_df,
    r0_ui   = r0_ui,
    r0_df   = r0_df
  ))
}



## heatmap for delay -------------------------------------------------
heatmap_delay <- function(target_age_list,
                          observed,
                          posterior, 
                          N, 
                          bra_foi_state_summ, 
                          age_groups, 
                          region_name,
                          hosp, 
                          fatal, 
                          nh_fatal,
                          lhs_sample_young, 
                          lhs_old, 
                          le_sample, 
                          pre_summary_cases) {

  
  # no vacc scenario 
  pre_ui <- simulate_pre_ui_age(
    posterior, bra_foi_state_summ, age_groups, N, region_name, observed,
    delay     = Inf,       
    VE_block  = 0,        
    target_age = rep(0, length(N)),
    coverage_threshold = 0,
    total_coverage     = 0,
    weekly_delivery_speed = 0
  )
  
  pre_summary_cases <- summarise_presim_ui(
    sim_result       = pre_ui,
    observed         = observed,
    age_gr_levels    = c("<5 years","5-9 years","10-14 years","15-19 years",
                         "20-24 years","25-29 years","30-34 years","35-39 years",
                         "40-44 years","45-49 years","50-54 years","55-59 years",
                         "60-64 years","65-69 years","70-74 years","75-79 years",
                         "80-84 years","85-89 years"),
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = hosp,
    fatal            = fatal,
    nh_fatal         = nh_fatal,
    region           = region_name
  )$summary_cases_pre 
  
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # supply and delay
  supply <- seq(0, 0.1, by = 0.005)
  delay_steps  <- seq(1, T, by = 1)
  
  # Initialize an empty list to store impact results for each scenario, supply, delay combination.
  impact_data <- list()
  
  # Number of scenarios (each target_age_list element defines a scenario)
  n_scenarios <- length(target_age_list)
  
  # extract median params
  base_beta       <- apply(posterior$base_beta, 2, median)      
  I0              <- apply(posterior$I0, 2, median)           
  gamma           <- median(posterior$gamma)
  rho             <- median(posterior$rho)
  
  # Loop over scenarios, supply rates, and delay steps
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    for (supply_rate in supply) { 
      for (delay_step in delay_steps) {
        
        # Run simulation with given parameters
        sim_result <- sirv_sim_coverageSwitch(
          T = T, 
          A = 18,
          N = N,
          r = rep(0, 18),
          base_beta = base_beta,   
          I0 = I0,            
          R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
          rho = rho,               
          gamma = gamma,           
          delay = delay_step,
          VE_block = 0.989,
          coverage_threshold = 1,
          target_age = target,
          total_coverage = supply_rate,
          weekly_delivery_speed = 0.1
        )
        
        # Create a data frame from the simulated age-stratified cases
        sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
          mutate(
            Scenario = scenario_index,
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          ) %>%
          mutate(
            hosp_rate = rep(hosp, T),
            hospitalised = Cases * hosp_rate,
            non_hospitalised = Cases - hospitalised,
            fatality = rep(fatal, T),
            nh_fatality = rep(nh_fatal, T),
            fatal = (hospitalised * fatality + non_hospitalised * nh_fatality)
          ) %>%
          mutate(
            age_numeric = case_when(
              str_detect(age_gr, "<5") ~ 0,
              TRUE ~ as.numeric(str_extract(age_gr, "^\\d+"))
            )
          ) %>%
          # (Additional YLD calculations omitted for brevity – include your full chain here)
          mutate(
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
            dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
            dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
          ) %>%
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
          mutate(
            yld_acute = (hospitalised * dw_hosp * dur_acute) +
              (non_hospitalised * dw_nonhosp * dur_acute),
            yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
              (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
            yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
              (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
            yld_total = yld_acute + yld_subacute + yld_chronic
          ) %>%
          mutate(
            le_left = case_when(
              age_numeric %in% c(0, 5) ~ quantile(le_sample$le_1, 0.5),
              age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
              age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
              age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
              age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
              age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
              age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
              age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
              age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
            )
          ) %>%
          mutate(
            yll = fatal * le_left,
            daly_tot = yld_total + yll,
            cum_daly = cumsum(daly_tot)
          )
       
         # Summarize pre-vaccination cases from the simulation output (using a separate summary data frame)
        pre_vacc_cases <- pre_summary_cases %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 1–4
            tot_pre_cases_agegr1  = sum(Median[AgeGroup %in% 1:4],  na.rm = TRUE),
            tot_pre_fatal_agegr1  = sum(fatal[AgeGroup   %in% 1:4],  na.rm = TRUE),
            tot_pre_daly_agegr1   = sum(daly_tot[AgeGroup   %in% 1:4],  na.rm = TRUE),
            # AgeGroup 5–12
            tot_pre_cases_agegr2  = sum(Median[AgeGroup %in% 5:12], na.rm = TRUE),
            tot_pre_fatal_agegr2  = sum(fatal[AgeGroup   %in% 5:12], na.rm = TRUE),
            tot_pre_daly_agegr2   = sum(daly_tot[AgeGroup   %in% 5:12], na.rm = TRUE),
            # AgeGroup 13–18
            tot_pre_cases_agegr3  = sum(Median[AgeGroup %in% 13:18], na.rm = TRUE),
            tot_pre_fatal_agegr3  = sum(fatal[AgeGroup   %in% 13:18], na.rm = TRUE),
            tot_pre_daly_agegr3   = sum(daly_tot[AgeGroup   %in% 13:18], na.rm = TRUE),
            .groups = "drop"
          )
        
        # Summarize post-vaccination cases from sim_df
        post_vacc_cases <- sim_df %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 1–4
            tot_post_cases_agegr1  = sum(Cases[AgeGroup %in% 1:4],      na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(fatal[AgeGroup %in% 1:4],      na.rm = TRUE),
            tot_post_daly_agegr1   = sum(daly_tot[AgeGroup %in% 1:4],   na.rm = TRUE),
            # AgeGroup 5–12
            tot_post_cases_agegr2  = sum(Cases[AgeGroup %in% 5:12],     na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(fatal[AgeGroup %in% 5:12],     na.rm = TRUE),
            tot_post_daly_agegr2   = sum(daly_tot[AgeGroup %in% 5:12],  na.rm = TRUE),
            # AgeGroup 13–18
            tot_post_cases_agegr3  = sum(Cases[AgeGroup %in% 13:18],    na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(fatal[AgeGroup %in% 13:18],    na.rm = TRUE),
            tot_post_daly_agegr3   = sum(daly_tot[AgeGroup %in% 13:18], na.rm = TRUE),
            .groups = "drop"
          )
        
        # Join the two and compute overall impact
        impact_data_summ <- pre_vacc_cases %>%
          left_join(post_vacc_cases, by = "Week") %>%
          summarise(
            # — Age group 1 (1–4)
            tot_pre_cases_agegr1   = sum(tot_pre_cases_agegr1,   na.rm = TRUE),
            tot_post_cases_agegr1  = sum(tot_post_cases_agegr1,  na.rm = TRUE),
            diff_cases_agegr1      = tot_pre_cases_agegr1   - tot_post_cases_agegr1,
            tot_pre_fatal_agegr1   = sum(tot_pre_fatal_agegr1,   na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(tot_post_fatal_agegr1,  na.rm = TRUE),
            diff_fatal_agegr1      = tot_pre_fatal_agegr1   - tot_post_fatal_agegr1,
            tot_pre_daly_agegr1    = sum(tot_pre_daly_agegr1,    na.rm = TRUE),
            tot_post_daly_agegr1   = sum(tot_post_daly_agegr1,   na.rm = TRUE),
            diff_daly_agegr1       = tot_pre_daly_agegr1    - tot_post_daly_agegr1,

            # — Age group 2 (5–12)
            tot_pre_cases_agegr2   = sum(tot_pre_cases_agegr2,   na.rm = TRUE),
            tot_post_cases_agegr2  = sum(tot_post_cases_agegr2,  na.rm = TRUE),
            diff_cases_agegr2      = tot_pre_cases_agegr2   - tot_post_cases_agegr2,
            tot_pre_fatal_agegr2   = sum(tot_pre_fatal_agegr2,   na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(tot_post_fatal_agegr2,  na.rm = TRUE),
            diff_fatal_agegr2      = tot_pre_fatal_agegr2   - tot_post_fatal_agegr2,
            tot_pre_daly_agegr2    = sum(tot_pre_daly_agegr2,    na.rm = TRUE),
            tot_post_daly_agegr2   = sum(tot_post_daly_agegr2,   na.rm = TRUE),
            diff_daly_agegr2       = tot_pre_daly_agegr2    - tot_post_daly_agegr2,

            # — Age group 3 (13–18)
            tot_pre_cases_agegr3   = sum(tot_pre_cases_agegr3,   na.rm = TRUE),
            tot_post_cases_agegr3  = sum(tot_post_cases_agegr3,  na.rm = TRUE),
            diff_cases_agegr3      = tot_pre_cases_agegr3   - tot_post_cases_agegr3,
            tot_pre_fatal_agegr3   = sum(tot_pre_fatal_agegr3,   na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(tot_post_fatal_agegr3,  na.rm = TRUE),
            diff_fatal_agegr3      = tot_pre_fatal_agegr3   - tot_post_fatal_agegr3,
            tot_pre_daly_agegr3    = sum(tot_pre_daly_agegr3,    na.rm = TRUE),
            tot_post_daly_agegr3   = sum(tot_post_daly_agegr3,   na.rm = TRUE),
            diff_daly_agegr3       = tot_pre_daly_agegr3    - tot_post_daly_agegr3,

            .groups = "drop"
          )
        
        
        # Save results with a unique name in impact_data
        impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_Delay_", delay_step)]] <- list(
          pre_vacc_cases = pre_vacc_cases,
          post_vacc_cases = post_vacc_cases,
          impact_data_summ = impact_data_summ
        )
      }
    }
  }
  
  # Combine impact summaries into one data frame
  impact_summ <- do.call(rbind, lapply(names(impact_data), function(idx) {
    list <- impact_data[[idx]]  # Access the element using its name
    df <- list$impact_data_summ           
    df$scenario_id <- idx          
    return(df)                        
  }))
  
  
  impact_summ <- impact_summ %>%
    mutate(
      Scenario = sub("_Supply.*", "", scenario_id),
      Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
      Delay = as.numeric(sub("Delay_", "", stringr::str_extract(scenario_id, "Delay_\\d+(\\.\\d+)?")))
    )
  
  impact_summ$region <- region_name
  
  return(impact_summ)
}


### sim pre med
simulate_pre_ui_age_med <- function(posterior, bra_foi_state_summ, age_groups, N, region, observed,
                                    A = 20,
                                    delay = 53, VE_block = 0, VE_inf = 0,
                                    target_age = rep(0, A),
                                    coverage_threshold = 0, total_coverage = 0,
                                    weekly_delivery_speed = 0) {
  # 1) 시간 길이 설정
  T <- nrow(observed)
  
  # 2) posterior에서 median 파라미터 추출
  base_beta_med <- apply(posterior$base_beta, 2, median)
  I0_med        <- apply(posterior$I0,      2, median)
  gamma_med     <- median(posterior$gamma)
  rho_med       <- median(posterior$rho)
  sigma_med     <- median(posterior$sigma)
  
  # 3) 단일 median 기반 시뮬레이션
  sim_out <- sirv_sim_coverageSwitch(
    T                    = T,
    A                    = A,
    N                    = N,
    r                    = rep(0, A),
    base_beta            = base_beta_med,
    I0_draw              = I0_med,
    R0                   = 1 - exp(- bra_foi_state_summ$avg_foi[
      bra_foi_state_summ$NAME_1 == region
    ] * age_groups),
    rho                  = rho_med,
    gamma                = gamma_med,
    sigma                = sigma_med,
    delay                = delay,
    VE_block             = VE_block,
    VE_inf               = VE_inf,
    target_age           = target_age,
    coverage_threshold   = coverage_threshold,
    total_coverage       = total_coverage,
    weekly_delivery_speed= weekly_delivery_speed
  )
  
  # 4) age-stratified symptomatic median (A×T)
  median_by_age_rawsymp <- sim_out$true_symptomatic
  # no uncertainty at median-run
  low95_by_age_rawsymp  <- median_by_age_rawsymp
  hi95_by_age_rawsymp   <- median_by_age_rawsymp
  
  # 5) 주별 합계
  weekly_cases_rawsymp <- colSums(median_by_age_rawsymp)
  
  # 6) weekly UI = median only
  weekly_low95_rawsymp <- weekly_cases_rawsymp
  weekly_hi95_rawsymp  <- weekly_cases_rawsymp
  
  # 7) Data frames
  df_rawsymp <- data.frame(
    week   = 1:T,
    median = weekly_cases_rawsymp,
    low95  = weekly_low95_rawsymp,
    hi95   = weekly_hi95_rawsymp
  )
  
  # 반환
  list(
    sim_out                     = sim_out,
    median_by_age_rawsymp       = median_by_age_rawsymp,
    low95_by_age_rawsymp        = low95_by_age_rawsymp,
    hi95_by_age_rawsymp         = hi95_by_age_rawsymp,
    weekly_cases_median_rawsymp = df_rawsymp$median,
    weekly_cases_low95_rawsymp  = df_rawsymp$low95,
    weekly_cases_hi95_rawsymp   = df_rawsymp$hi95,
    df_rawsymp                  = df_rawsymp
  )
}
## v2 heatmap (using true symp)
heatmap_delay <- function(target_age_list,
                          observed,
                          posterior, 
                          N, 
                          bra_foi_state_summ, 
                          age_groups, 
                          region_name,
                          hosp, 
                          fatal, 
                          nh_fatal,
                          lhs_sample_young, 
                          lhs_old, 
                          le_sample, 
                          pre_summary_cases,
                          VE_inf = 0 ) {
  
  # no vacc scenario 
  pre_ui <- simulate_pre_ui_age_med(
    posterior, bra_foi_state_summ, age_groups, N, region_name, observed,
    delay     = Inf,       
    VE_block  = 0,   
    VE_inf    = 0,
    target_age = rep(0, length(N)),
    coverage_threshold = 0,
    total_coverage     = 0,
    weekly_delivery_speed = 0
  )
  
  pre_summary_cases <- summarise_presim_ui(
    sim_result       = pre_ui,
    observed         = observed,
    age_gr_levels    =  c(
      "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
      "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
      "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
    ),
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = hosp,
    fatal            = fatal,
    nh_fatal         = nh_fatal,
    region           = region_name
  )$summary_cases_pre 
  
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <-  c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # supply and delay

  pct <- region_coverage[region_name]
  #max_coverage <- floor(pct) / 100
  
  supply <- seq(0, 1, by = 0.1)
  
  #supply <- seq(0, max_coverage, by = 0.005)
  delay_steps  <- seq(1, T, by = 1)
  
  # Initialize an empty list to store impact results for each scenario, supply, delay combination.
  impact_data <- list()
  
  # Number of scenarios (each target_age_list element defines a scenario)
  n_scenarios <- length(target_age_list)
  
  # extract median params
  base_beta       <- apply(posterior$base_beta, 2, median)      
  I0              <- apply(posterior$I0, 2, median)           
  gamma           <- median(posterior$gamma)
  rho             <- median(posterior$rho)
  sigma           <- median(posterior$sigma)
  
  # Loop over scenarios, supply rates, and delay steps
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    for (supply_rate in supply) { 
      for (delay_step in delay_steps) {
        
        # Run simulation with given parameters
        sim_result <- sirv_sim_coverageSwitch(
          T = T, 
          A = length(age_gr_levels),
          N = N,
          r = rep(0, 20),
          base_beta = base_beta,   
          I0 = I0,            
          R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
          rho = rho,               
          gamma = gamma,        
          sigma = sigma,
          delay = delay_step,
          VE_block = 0.989,
          VE_inf = 0.989,
          coverage_threshold = 1,
          target_age = target,
          total_coverage = supply_rate,
          weekly_delivery_speed = 0.1
        )
        
        # Create a data frame from the simulated age-stratified cases
        sim_df <- as.data.frame.table(sim_result$true_symptomatic, responseName = "Cases") %>%
          mutate(
            Scenario = scenario_index,
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          ) %>%
          mutate(
            hosp_rate = rep(hosp, T),
            hospitalised = Cases * hosp_rate,
            non_hospitalised = Cases - hospitalised,
            fatality = rep(fatal, T),
            nh_fatality = rep(nh_fatal, T),
            fatal = (hospitalised * fatality + non_hospitalised * nh_fatality)
          ) %>%
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
          ) %>%
          # (Additional YLD calculations omitted for brevity – include your full chain here)
          mutate(
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
            dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
            dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
          ) %>%
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
          mutate(
            yld_acute = (hospitalised * dw_hosp * dur_acute) +
              (non_hospitalised * dw_nonhosp * dur_acute),
            yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
              (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
            yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
              (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
            yld_total = yld_acute + yld_subacute + yld_chronic
          ) %>%
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
            daly_tot = yld_total + yll,
            cum_daly = cumsum(daly_tot)
          )
        
        
        # Summarize pre-vaccination cases from the simulation output (using a separate summary data frame)
        pre_vacc_cases <- pre_summary_cases %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 2-4 (strategy1: 1-11)
            tot_pre_cases_agegr1  = sum(Median[AgeGroup %in% 2:4],  na.rm = TRUE),
            tot_pre_fatal_agegr1  = sum(fatal[AgeGroup   %in% 2:4],  na.rm = TRUE),
            tot_pre_daly_agegr1   = sum(daly_tot[AgeGroup   %in% 2:4],  na.rm = TRUE),
            # AgeGroup 5 (strategy2: 12-17)
            tot_pre_cases_agegr2  = sum(Median[AgeGroup %in% 5], na.rm = TRUE),
            tot_pre_fatal_agegr2  = sum(fatal[AgeGroup   %in% 5], na.rm = TRUE),
            tot_pre_daly_agegr2   = sum(daly_tot[AgeGroup   %in% 5], na.rm = TRUE),
            # AgeGroup 6-15 (strategy3: 18-59)
            tot_pre_cases_agegr3  = sum(Median[AgeGroup %in% 6:15], na.rm = TRUE),
            tot_pre_fatal_agegr3  = sum(fatal[AgeGroup   %in% 6:15], na.rm = TRUE),
            tot_pre_daly_agegr3   = sum(daly_tot[AgeGroup   %in% 6:15], na.rm = TRUE),
            # AgeGroup 16-20 (strategy3: 60+)
            tot_pre_cases_agegr4  = sum(Median[AgeGroup %in% 16:20], na.rm = TRUE),
            tot_pre_fatal_agegr4  = sum(fatal[AgeGroup   %in% 16:20], na.rm = TRUE),
            tot_pre_daly_agegr4   = sum(daly_tot[AgeGroup   %in% 16:20], na.rm = TRUE),

            .groups = "drop"
          )
        
        # scaling
        pre_vacc_cases <- pre_vacc_cases %>% 
          mutate(region = region_name) %>%
          left_join(rho_df, by = "region") %>%
          mutate(
            across(
              c(tot_pre_cases_agegr1, tot_pre_fatal_agegr1, tot_pre_daly_agegr1,
                tot_pre_cases_agegr2, tot_pre_fatal_agegr2, tot_pre_daly_agegr2,
                tot_pre_cases_agegr3, tot_pre_fatal_agegr3, tot_pre_daly_agegr3,
                tot_pre_cases_agegr4, tot_pre_fatal_agegr4, tot_pre_daly_agegr4
              ),
              ~ .x / rho_p50
            )
          )
        
        # Summarize post-vaccination cases from sim_df
        post_vacc_cases <- sim_df %>%
          mutate(region = region_name) %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 2-4 (strategy1: 1-11)
            tot_post_cases_agegr1  = sum(Cases[AgeGroup %in% 2:4],  na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(fatal[AgeGroup   %in% 2:4],  na.rm = TRUE),
            tot_post_daly_agegr1   = sum(daly_tot[AgeGroup   %in% 2:4],  na.rm = TRUE),
            # AgeGroup 5 (strategy2: 12-17)
            tot_post_cases_agegr2  = sum(Cases[AgeGroup %in% 5], na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(fatal[AgeGroup   %in% 5], na.rm = TRUE),
            tot_post_daly_agegr2   = sum(daly_tot[AgeGroup   %in% 5], na.rm = TRUE),
            # AgeGroup 6-15 (strategy3: 18-59)
            tot_post_cases_agegr3  = sum(Cases[AgeGroup %in% 6:15], na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(fatal[AgeGroup   %in% 6:15], na.rm = TRUE),
            tot_post_daly_agegr3   = sum(daly_tot[AgeGroup   %in% 6:15], na.rm = TRUE),
            # AgeGroup 16-20 (strategy3: 60+)
            tot_post_cases_agegr4  = sum(Cases[AgeGroup %in% 16:20], na.rm = TRUE),
            tot_post_fatal_agegr4  = sum(fatal[AgeGroup   %in% 16:20], na.rm = TRUE),
            tot_post_daly_agegr4   = sum(daly_tot[AgeGroup   %in% 16:20], na.rm = TRUE),

            .groups = "drop"
          )
        
        post_vacc_cases <- post_vacc_cases %>%
          mutate(region = region_name) %>%
          left_join(rho_df, by = "region") %>%
          mutate(
            across(
              c(tot_post_cases_agegr1, tot_post_fatal_agegr1, tot_post_daly_agegr1,
                tot_post_cases_agegr2, tot_post_fatal_agegr2, tot_post_daly_agegr2,
                tot_post_cases_agegr3, tot_post_fatal_agegr3, tot_post_daly_agegr3,
                tot_post_cases_agegr4, tot_post_fatal_agegr4, tot_post_daly_agegr4),
              ~ .x / rho_p50
            )
          )
        
        # Join the two and compute overall impact
        impact_data_summ <- pre_vacc_cases %>%
          left_join(post_vacc_cases, by = "Week") %>%
          summarise(
            # — Age group 1 
            tot_pre_cases_agegr1   = sum(tot_pre_cases_agegr1,   na.rm = TRUE),
            tot_post_cases_agegr1  = sum(tot_post_cases_agegr1,  na.rm = TRUE),
            diff_cases_agegr1      = tot_pre_cases_agegr1   - tot_post_cases_agegr1,
            tot_pre_fatal_agegr1   = sum(tot_pre_fatal_agegr1,   na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(tot_post_fatal_agegr1,  na.rm = TRUE),
            diff_fatal_agegr1      = tot_pre_fatal_agegr1   - tot_post_fatal_agegr1,
            tot_pre_daly_agegr1    = sum(tot_pre_daly_agegr1,    na.rm = TRUE),
            tot_post_daly_agegr1   = sum(tot_post_daly_agegr1,   na.rm = TRUE),
            diff_daly_agegr1       = tot_pre_daly_agegr1    - tot_post_daly_agegr1,
            
            # — Age group 2 
            tot_pre_cases_agegr2   = sum(tot_pre_cases_agegr2,   na.rm = TRUE),
            tot_post_cases_agegr2  = sum(tot_post_cases_agegr2,  na.rm = TRUE),
            diff_cases_agegr2      = tot_pre_cases_agegr2   - tot_post_cases_agegr2,
            tot_pre_fatal_agegr2   = sum(tot_pre_fatal_agegr2,   na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(tot_post_fatal_agegr2,  na.rm = TRUE),
            diff_fatal_agegr2      = tot_pre_fatal_agegr2   - tot_post_fatal_agegr2,
            tot_pre_daly_agegr2    = sum(tot_pre_daly_agegr2,    na.rm = TRUE),
            tot_post_daly_agegr2   = sum(tot_post_daly_agegr2,   na.rm = TRUE),
            diff_daly_agegr2       = tot_pre_daly_agegr2    - tot_post_daly_agegr2,
            
            # — Age group 3 
            tot_pre_cases_agegr3   = sum(tot_pre_cases_agegr3,   na.rm = TRUE),
            tot_post_cases_agegr3  = sum(tot_post_cases_agegr3,  na.rm = TRUE),
            diff_cases_agegr3      = tot_pre_cases_agegr3   - tot_post_cases_agegr3,
            tot_pre_fatal_agegr3   = sum(tot_pre_fatal_agegr3,   na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(tot_post_fatal_agegr3,  na.rm = TRUE),
            diff_fatal_agegr3      = tot_pre_fatal_agegr3   - tot_post_fatal_agegr3,
            tot_pre_daly_agegr3    = sum(tot_pre_daly_agegr3,    na.rm = TRUE),
            tot_post_daly_agegr3   = sum(tot_post_daly_agegr3,   na.rm = TRUE),
            diff_daly_agegr3       = tot_pre_daly_agegr3    - tot_post_daly_agegr3,
            
            # — Age group 4 
            tot_pre_cases_agegr4   = sum(tot_pre_cases_agegr4,   na.rm = TRUE),
            tot_post_cases_agegr4  = sum(tot_post_cases_agegr4,  na.rm = TRUE),
            diff_cases_agegr4      = tot_pre_cases_agegr4   - tot_post_cases_agegr4,
            tot_pre_fatal_agegr4   = sum(tot_pre_fatal_agegr4,   na.rm = TRUE),
            tot_post_fatal_agegr4  = sum(tot_post_fatal_agegr4,  na.rm = TRUE),
            diff_fatal_agegr4      = tot_pre_fatal_agegr4   - tot_post_fatal_agegr4,
            tot_pre_daly_agegr4    = sum(tot_pre_daly_agegr4,    na.rm = TRUE),
            tot_post_daly_agegr4   = sum(tot_post_daly_agegr4,   na.rm = TRUE),
            diff_daly_agegr4       = tot_pre_daly_agegr4    - tot_post_daly_agegr4,
            
            .groups = "drop"
          )
        
        total_vacc <- sum(sim_result$raw_allocation_age, na.rm = TRUE)
        
        # Save results with a unique name in impact_data
        impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_Delay_", delay_step)]] <- list(
          pre_vacc_cases = pre_vacc_cases,
          post_vacc_cases = post_vacc_cases,
          impact_data_summ = impact_data_summ,
          total_vacc = total_vacc
        )
      }
    }
  }
  
  # Combine impact summaries into one data frame
  impact_summ <- do.call(rbind, lapply(names(impact_data), function(idx) {
    list <- impact_data[[idx]]  # Access the element using its name
    df <- list$impact_data_summ           
    df$scenario_id <- idx   
    df$total_vacc <- list$total_vacc
    return(df)                        
  }))
  
  
  impact_summ <- impact_summ %>%
    mutate(
      Scenario = sub("_Supply.*", "", scenario_id),
      Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
      Delay = as.numeric(sub("Delay_", "", stringr::str_extract(scenario_id, "Delay_\\d+(\\.\\d+)?")))
    )
  
  impact_summ$region <- region_name
  
  return(impact_summ)
}


prepare_heatmap_summary <- function(heatmap_df) {
  heatmap_df %>%
    mutate(
      scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
    )
}

plot_heatmap <- function(data, prefix, legend_title) {
  
  custom_labeller <- function(x) {
    sapply(x, function(v) {
      if (abs(v) >= 1e6) {
        paste0(round(v / 1e6, 1), "M")
      } else if (abs(v) >= 1e3) {
        formatC(v, format = "d", big.mark = ",")
      } else {
        as.character(v)
      }
    })
  }
  
  data %>%
    select(Supply, Delay, Scenario, scenario_num, starts_with(prefix)) %>%
    pivot_longer(
      cols         = starts_with(prefix),  # prefix로 시작하는 모든 agegr 열들
      names_to     = "agegr",              # 변수명에서 숫자 부분만 뽑아 agegr에
      names_prefix = prefix,               # prefix(예: "impact_cases_agegr") 제거
      values_to    = "impact"              # 값은 impact라는 단일 열로
    ) %>%
    mutate(agegr = as.integer(agegr)) %>% # "1","2","3" → 정수형으로
    filter(agegr == scenario_num) %>%     # 시나리오 번호와 같은 agegr만 남김
    ggplot(aes(x = Supply, y = Delay, fill = impact)) +
    geom_tile() +
    scale_fill_distiller(
      palette   = "Spectral",
      direction = -1,
      name      = legend_title,
      labels = custom_labeller
    ) + 
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),  
      labels = function(x) sprintf("%.0f", x * 100)
    ) +
    scale_y_continuous(labels = custom_labeller) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_wrap(~ Scenario, labeller = as_labeller(c(
      "Scenario_1" = "1-19 years only",
      "Scenario_2" = "20-59 years only",
      "Scenario_3" = "≥60 years only"
    ))) +
    #geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", size = 0.2) +
    #geom_hline(yintercept = 2,   linetype = "dashed", color = "black", size = 0.2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(),
      panel.grid   = element_blank()
    ) +
    labs(x = "Vaccine supply (% of total population)")
}


# contour map 
ve_contour <- function(target_age_list,
                          observed,
                          posterior, 
                          N, 
                          bra_foi_state_summ, 
                          age_groups, 
                          region_name,
                          hosp, 
                          fatal, 
                          nh_fatal,
                          lhs_sample_young, 
                          lhs_old, 
                          le_sample, 
                          pre_summary_cases) {
  
  # no vacc scenario 
  pre_ui <- simulate_pre_ui_age_med(
    posterior, bra_foi_state_summ, age_groups, N, region_name, observed,
    delay     = Inf,       
    VE_block  = 0,        
    target_age = rep(0, length(N)),
    coverage_threshold = 0,
    total_coverage     = 0,
    weekly_delivery_speed = 0
  )
  
  pre_summary_cases <- summarise_presim_ui(
    sim_result       = pre_ui,
    observed         = observed,
    age_gr_levels    = c("<5 years","5-9 years","10-14 years","15-19 years",
                         "20-24 years","25-29 years","30-34 years","35-39 years",
                         "40-44 years","45-49 years","50-54 years","55-59 years",
                         "60-64 years","65-69 years","70-74 years","75-79 years",
                         "80-84 years","85-89 years"),
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = hosp,
    fatal            = fatal,
    nh_fatal         = nh_fatal,
    region           = region_name
  )$summary_cases_pre 
  
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # VE_i and VE_b
  VE_block <- seq(0.7, 0.989, by = 0.05)
  VE_inf   <- seq(0, 0.989, by = 0.05)
  
  # Initialize an empty list to store impact results for each scenario, supply, delay combination.
  impact_data <- list()
  
  # Number of scenarios (each target_age_list element defines a scenario)
  n_scenarios <- length(target_age_list)
  
  # extract median params
  base_beta       <- apply(posterior$base_beta, 2, median)      
  I0              <- apply(posterior$I0, 2, median)           
  gamma           <- median(posterior$gamma)
  rho             <- median(posterior$rho)
  
  # Loop over scenarios, supply rates, and delay steps
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    for (VE_inf_rate in VE_inf) { 
      for (VE_block_rate in VE_block) {
        
        # Run simulation with given parameters
        sim_result <- sirv_sim_coverageSwitch(
          T = T, 
          A = 18,
          N = N,
          r = rep(0, 18),
          base_beta = base_beta,   
          I0 = I0,            
          R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
          rho = rho,               
          gamma = gamma,           
          delay = 2,
          VE_block = VE_block_rate,
          VE_inf = VE_inf_rate,
          coverage_threshold = 1,
          target_age = target,
          total_coverage = 0.12,
          weekly_delivery_speed = 0.1,
          VE_inf_assumed = 0.989
        )
        
        # Create a data frame from the simulated age-stratified cases
        sim_df <- as.data.frame.table(sim_result$true_symptomatic, responseName = "Cases") %>%
          mutate(
            Scenario = scenario_index,
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          ) %>%
          mutate(
            hosp_rate = rep(hosp, T),
            hospitalised = Cases * hosp_rate,
            non_hospitalised = Cases - hospitalised,
            fatality = rep(fatal, T),
            nh_fatality = rep(nh_fatal, T),
            fatal = (hospitalised * fatality + non_hospitalised * nh_fatality)
          ) %>%
          mutate(
            age_numeric = case_when(
              str_detect(age_gr, "<5") ~ 0,
              TRUE ~ as.numeric(str_extract(age_gr, "^\\d+"))
            )
          ) %>%
          # (Additional YLD calculations omitted for brevity – include your full chain here)
          mutate(
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
            dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
            dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
          ) %>%
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
          mutate(
            yld_acute = (hospitalised * dw_hosp * dur_acute) +
              (non_hospitalised * dw_nonhosp * dur_acute),
            yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
              (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
            yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
              (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
            yld_total = yld_acute + yld_subacute + yld_chronic
          ) %>%
          mutate(
            le_left = case_when(
              age_numeric %in% c(0, 5) ~ quantile(le_sample$le_1, 0.5),
              age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
              age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
              age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
              age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
              age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
              age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
              age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
              age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
            )
          ) %>%
          mutate(
            yll = fatal * le_left,
            daly_tot = yld_total + yll,
            cum_daly = cumsum(daly_tot)
          )
        
        
        # Summarize pre-vaccination cases from the simulation output (using a separate summary data frame)
        pre_vacc_cases <- pre_summary_cases %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 1–4
            tot_pre_cases_agegr1  = sum(Median[AgeGroup %in% 1:4],  na.rm = TRUE),
            tot_pre_fatal_agegr1  = sum(fatal[AgeGroup   %in% 1:4],  na.rm = TRUE),
            tot_pre_daly_agegr1   = sum(daly_tot[AgeGroup   %in% 1:4],  na.rm = TRUE),
            # AgeGroup 5–12
            tot_pre_cases_agegr2  = sum(Median[AgeGroup %in% 5:12], na.rm = TRUE),
            tot_pre_fatal_agegr2  = sum(fatal[AgeGroup   %in% 5:12], na.rm = TRUE),
            tot_pre_daly_agegr2   = sum(daly_tot[AgeGroup   %in% 5:12], na.rm = TRUE),
            # AgeGroup 13–18
            tot_pre_cases_agegr3  = sum(Median[AgeGroup %in% 13:18], na.rm = TRUE),
            tot_pre_fatal_agegr3  = sum(fatal[AgeGroup   %in% 13:18], na.rm = TRUE),
            tot_pre_daly_agegr3   = sum(daly_tot[AgeGroup   %in% 13:18], na.rm = TRUE),
            .groups = "drop"
          )
        
        # scaling
        pre_vacc_cases <- pre_vacc_cases %>% 
          mutate(region = region_name) %>%
          left_join(rho_df, by = "region") %>%
          mutate(
            across(
              c(tot_pre_cases_agegr1, tot_pre_fatal_agegr1, tot_pre_daly_agegr1,
                tot_pre_cases_agegr2, tot_pre_fatal_agegr2, tot_pre_daly_agegr2,
                tot_pre_cases_agegr3, tot_pre_fatal_agegr3, tot_pre_daly_agegr3),
              ~ .x / rho_p50
            )
          )
        
        # Summarize post-vaccination cases from sim_df
        post_vacc_cases <- sim_df %>%
          group_by(Week) %>%
          summarise(
            # AgeGroup 1–4
            tot_post_cases_agegr1  = sum(Cases[AgeGroup %in% 1:4],      na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(fatal[AgeGroup %in% 1:4],      na.rm = TRUE),
            tot_post_daly_agegr1   = sum(daly_tot[AgeGroup %in% 1:4],   na.rm = TRUE),
            # AgeGroup 5–12
            tot_post_cases_agegr2  = sum(Cases[AgeGroup %in% 5:12],     na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(fatal[AgeGroup %in% 5:12],     na.rm = TRUE),
            tot_post_daly_agegr2   = sum(daly_tot[AgeGroup %in% 5:12],  na.rm = TRUE),
            # AgeGroup 13–18
            tot_post_cases_agegr3  = sum(Cases[AgeGroup %in% 13:18],    na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(fatal[AgeGroup %in% 13:18],    na.rm = TRUE),
            tot_post_daly_agegr3   = sum(daly_tot[AgeGroup %in% 13:18], na.rm = TRUE),
            .groups = "drop"
          )
        
        post_vacc_cases <- post_vacc_cases %>%
          mutate(region = region_name) %>%
          left_join(rho_df, by = "region") %>%
          mutate(
            across(
              c(tot_post_cases_agegr1, tot_post_fatal_agegr1, tot_post_daly_agegr1,
                tot_post_cases_agegr2, tot_post_fatal_agegr2, tot_post_daly_agegr2,
                tot_post_cases_agegr3, tot_post_fatal_agegr3, tot_post_daly_agegr3),
              ~ .x / rho_p50
            )
          )
        
        # Join the two and compute overall impact
        impact_data_summ <- pre_vacc_cases %>%
          left_join(post_vacc_cases, by = "Week") %>%
          summarise(
            # — Age group 1 (1–4)
            tot_pre_cases_agegr1   = sum(tot_pre_cases_agegr1,   na.rm = TRUE),
            tot_post_cases_agegr1  = sum(tot_post_cases_agegr1,  na.rm = TRUE),
            diff_cases_agegr1      = tot_pre_cases_agegr1   - tot_post_cases_agegr1,
            tot_pre_fatal_agegr1   = sum(tot_pre_fatal_agegr1,   na.rm = TRUE),
            tot_post_fatal_agegr1  = sum(tot_post_fatal_agegr1,  na.rm = TRUE),
            diff_fatal_agegr1      = tot_pre_fatal_agegr1   - tot_post_fatal_agegr1,
            tot_pre_daly_agegr1    = sum(tot_pre_daly_agegr1,    na.rm = TRUE),
            tot_post_daly_agegr1   = sum(tot_post_daly_agegr1,   na.rm = TRUE),
            diff_daly_agegr1       = tot_pre_daly_agegr1    - tot_post_daly_agegr1,
            
            # — Age group 2 (5–12)
            tot_pre_cases_agegr2   = sum(tot_pre_cases_agegr2,   na.rm = TRUE),
            tot_post_cases_agegr2  = sum(tot_post_cases_agegr2,  na.rm = TRUE),
            diff_cases_agegr2      = tot_pre_cases_agegr2   - tot_post_cases_agegr2,
            tot_pre_fatal_agegr2   = sum(tot_pre_fatal_agegr2,   na.rm = TRUE),
            tot_post_fatal_agegr2  = sum(tot_post_fatal_agegr2,  na.rm = TRUE),
            diff_fatal_agegr2      = tot_pre_fatal_agegr2   - tot_post_fatal_agegr2,
            tot_pre_daly_agegr2    = sum(tot_pre_daly_agegr2,    na.rm = TRUE),
            tot_post_daly_agegr2   = sum(tot_post_daly_agegr2,   na.rm = TRUE),
            diff_daly_agegr2       = tot_pre_daly_agegr2    - tot_post_daly_agegr2,
            
            # — Age group 3 (13–18)
            tot_pre_cases_agegr3   = sum(tot_pre_cases_agegr3,   na.rm = TRUE),
            tot_post_cases_agegr3  = sum(tot_post_cases_agegr3,  na.rm = TRUE),
            diff_cases_agegr3      = tot_pre_cases_agegr3   - tot_post_cases_agegr3,
            tot_pre_fatal_agegr3   = sum(tot_pre_fatal_agegr3,   na.rm = TRUE),
            tot_post_fatal_agegr3  = sum(tot_post_fatal_agegr3,  na.rm = TRUE),
            diff_fatal_agegr3      = tot_pre_fatal_agegr3   - tot_post_fatal_agegr3,
            tot_pre_daly_agegr3    = sum(tot_pre_daly_agegr3,    na.rm = TRUE),
            tot_post_daly_agegr3   = sum(tot_post_daly_agegr3,   na.rm = TRUE),
            diff_daly_agegr3       = tot_pre_daly_agegr3    - tot_post_daly_agegr3,
            
            .groups = "drop"
          )
        
        
        # Save results with a unique name in impact_data
        impact_data[[paste0("Scenario_", scenario_index, "_VEinf_", VE_inf_rate, "_VEblock_", VE_block_rate)]] <- list(
          pre_vacc_cases = pre_vacc_cases,
          post_vacc_cases = post_vacc_cases,
          impact_data_summ = impact_data_summ
        )
      }
    }
  }
  
  # Combine impact summaries into one data frame
  impact_summ <- do.call(rbind, lapply(names(impact_data), function(idx) {
    list <- impact_data[[idx]]  # Access the element using its name
    df <- list$impact_data_summ           
    df$scenario_id <- idx          
    return(df)                        
  }))
  
  
  impact_summ <- impact_summ %>%
    mutate(
      Scenario   = sub("_VEinf.*", "", scenario_id),
      VE_inf     = as.numeric(str_extract(scenario_id, "(?<=VEinf_)\\d*\\.?\\d+")),
      VE_block   = as.numeric(str_extract(scenario_id, "(?<=VEblock_)\\d*\\.?\\d+"))
    )
  
  impact_summ$region <- region_name
  
  return(impact_summ)
}


# all scenario result summarise to compare pre-post cases --------------
attach_target_column <- function(df_nnv, observed, region) {

  T = nrow(observed)
  
  weeks <- T
  
  target_list <- list(
    c("<20 years", "<20 years", "<20 years", "<20 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years"),
    
    c("<20 years", "<20 years", "<20 years", "<20 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years"),
    
    c("<20 years", "<20 years", "<20 years", "<20 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years")
  )
  
  target1 <- rep(target_list[[1]])
  target2 <- rep(target_list[[2]])
  target3 <- rep(target_list[[3]])
  target  <- c(target1, target2, target3)
  
  df <- df_nnv
  
  df$target <- target
  
  df_summ  <- df %>% 
    group_by(target, scenario) %>% 
    summarise(
      pre_cases    = sum(pre_vacc, na.rm = TRUE),
      post_cases   = sum(Median, na.rm = TRUE),
      post_low95   = sum(low95, na.rm = TRUE),
      post_hi95    = sum(hi95, na.rm = TRUE),
      pre_low95    = sum(pre_vacc_low95, na.rm = TRUE),
      pre_hi95     = sum(pre_vacc_hi95, na.rm = TRUE),
      
      post_fatal   = sum(fatal, na.rm = TRUE),
      post_fatal_lo = sum(fatal_lo, na.rm = TRUE),
      post_fatal_hi = sum(fatal_hi, na.rm = TRUE),
      pre_fatal    = sum(pre_fatal, na.rm = TRUE),
      pre_fatal_lo = sum(pre_fatal_low95, na.rm = TRUE),
      pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
      
      post_daly    = sum(daly_tot, na.rm = TRUE),
      post_daly_lo = sum(daly_tot_lo, na.rm = TRUE),
      post_daly_hi = sum(daly_tot_hi, na.rm = TRUE),
      pre_daly     = sum(pre_daly, na.rm = TRUE),
      pre_daly_lo  = sum(pre_daly_low95, na.rm = TRUE),
      pre_daly_hi  = sum(pre_daly_hi, na.rm = TRUE),
      
      tot_vacc     = sum(tot_vacc, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(
      ## For Cases:
      diff_case_mid   = pre_cases - post_cases,
      diff_case_low   = pre_low95 - post_low95,
      diff_case_hi    = pre_hi95 - post_hi95,
      impact_case     = diff_case_mid / pre_cases * 100,
      impact_case_low = diff_case_low / pre_low95 * 100,
      impact_case_hi  = diff_case_hi / pre_hi95 * 100,
      
      ## For Fatal outcomes:
      diff_fatal      = pre_fatal - post_fatal,
      diff_fatal_low  = pre_fatal_lo - post_fatal_lo,
      diff_fatal_hi   = pre_fatal_hi - post_fatal_hi,
      impact_fatal    = diff_fatal / pre_fatal * 100,
      impact_fatal_low = diff_fatal_low / pre_fatal_lo * 100,
      impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,
      
      ## For DALY:
      diff_daly       = pre_daly - post_daly,
      diff_daly_low   = pre_daly_lo - post_daly_lo,
      diff_daly_hi    = pre_daly_hi - post_daly_hi,
      impact_daly     = diff_daly / pre_daly * 100,
      impact_daly_low = diff_daly_low / pre_daly_lo * 100,
      impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100,
      
      ## Optionally, you can compute NNV values for Cases and Fatal:
      nnv             = tot_vacc / diff_case_mid,
      nnv_lo          = tot_vacc / diff_case_low,
      nnv_hi          = tot_vacc / diff_case_hi,
      
      nnv_fatal       = tot_vacc / diff_fatal,
      nnv_fatal_lo    = tot_vacc / diff_fatal_low,
      nnv_fatal_hi    = tot_vacc / diff_fatal_hi,
      
      nnv_daly        = tot_vacc / diff_daly,
      nnv_daly_lo     = tot_vacc / diff_daly_low,
      nnv_daly_hi     = tot_vacc / diff_daly_hi
    )
  
  
  df_summ$region <- region

  return(df_summ)
}

make_tot_pop_summary <- function(N_pop) {
  
  target_list <- list(
    c("<20 years", "<20 years", "<20 years", "<20 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years",
      "20-59 years", "20-59 years", "20-59 years", "20-59 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years", ">60 years"))
    
  target1 <- rep(target_list[[1]])
  
  N_pop$age_gr <- age_gr[1:18]
  N_pop$target <- target1
  
  N_pop_sum <- N_pop
  colnames(N_pop_sum) <- c("pop", "age_gr", "target")
  
  N_pop_sum <- N_pop_sum %>% group_by(target) %>% summarise(
    tot_pop = sum(pop)
  )
  
  return(N_pop_sum)
}

# cumulative fatal case graph 
fatal_summ_func <- function(postsim_all_ui, N, observed, region) {
  library(dplyr)
  # vacc doses
  totvacc_s1 <- 10056615
  totvacc_s2 <- 22029741
  totvacc_s3 <- 5800257
  
  # 1) 시뮬레이션 결과 불러오기
  scenario_results <- postsim_all_ui$scenario_result
  
  mutated_list <- lapply(scenario_results, function(scn) {
    scn$sim_df %>%
      mutate(
        post_cases  = Cases,
        post_cases_lo = post_low95,
        post_cases_hi = post_hi95,
        # post-vax 사망
        post_fatal     = fatal,
        post_fatal_lo   = fatal_lo,
        post_fatal_hi  = fatal_hi,
        # post vacc daly
        post_daly = daly_tot,
        post_daly_lo = daly_tot_lo,
        post_daly_hi = daly_tot_hi,
        # pre-vax 사망
        pre_fatal      = pre_Median * hosp_rate * fatality +
          (pre_Median * (1 - hosp_rate)) * nh_fatality,
        pre_fatal_lo      = pre_low95 * hosp_rate * fatality +
          (pre_low95 * (1 - hosp_rate)) * nh_fatality,
        pre_fatal_hi      = pre_hi95 * hosp_rate * fatality +
          (pre_hi95 * (1 - hosp_rate)) * nh_fatality,
        # 누적 사망
        cum_post_fatal = cumsum(post_fatal),
        cum_post_fatal_lo = cumsum(post_fatal_lo),
        cum_post_fatal_hi = cumsum(post_fatal_hi),
        cum_pre_fatal  = cumsum(pre_fatal),
        cum_pre_fatal_lo  = cumsum(pre_fatal_lo),
        cum_pre_fatal_hi  = cumsum(pre_fatal_hi)
      )
  })
  
  # 3) 주차·연령그룹별 누적 사망 집계
  fatal_week_df <- bind_rows(lapply(seq_along(mutated_list), function(i) {
    mutated_list[[i]] %>%
      group_by(Week, AgeGroup) %>%
      summarise(
        cum_post_fatal = last(cum_post_fatal),
        cum_post_fatal_lo = last(cum_post_fatal_lo),
        cum_post_fatal_hi = last(cum_post_fatal_hi),
        cum_pre_fatal  = last(cum_pre_fatal),
        cum_pre_fatal_lo  = last(cum_pre_fatal_lo),
        cum_pre_fatal_hi  = last(cum_pre_fatal_hi),
        .groups = "drop"
      ) %>%
      mutate(Scenario = paste0("Scenario_", i))
  })) %>%
    ungroup() %>%
    mutate(region = region)
  
  # 4) 목표 연령그룹(1-4, 5-12, 13-18)별 최종 누적 사망 & diff per million 계산
  summary_df <- fatal_week_df %>%
    mutate(
      agegrp3 = case_when(
        AgeGroup %in% 1:4   ~ "1-4",
        AgeGroup %in% 5:12  ~ "5-12",
        AgeGroup %in% 13:18 ~ "13-18",
        TRUE                ~ NA_character_
      )
    ) %>%
    group_by(Scenario, agegrp3) %>%
    summarise(
      final_cum_post_fatal = last(cum_post_fatal),
      final_cum_post_fatal_lo = last(cum_post_fatal_lo),
      final_cum_post_fatal_hi = last(cum_post_fatal_hi),
      final_cum_pre_fatal  = last(cum_pre_fatal),
      final_cum_pre_fatal_lo  = last(cum_pre_fatal_lo),
      final_cum_pre_fatal_hi  = last(cum_pre_fatal_hi),
      .groups = "drop"
    ) %>%
    # 시나리오별 목표 그룹만 남기기
    filter(
      (Scenario == "Scenario_1" & agegrp3 == "1-4") |
        (Scenario == "Scenario_2" & agegrp3 == "5-12") |
        (Scenario == "Scenario_3" & agegrp3 == "13-18")
    ) %>%
    mutate(
      totvacc = case_when(
        Scenario == "Scenario_1" ~ totvacc_s1,
        Scenario == "Scenario_2" ~ totvacc_s2,
        Scenario == "Scenario_3" ~ totvacc_s3
      ),
      # (사망자수 백신 이전 – 사망자수 백신 이후) ÷ 투여된 백신수 × 1e6
      diff_per_million = (final_cum_pre_fatal - final_cum_post_fatal) / totvacc * 1e6,
      diff_per_million_lo = (final_cum_pre_fatal_lo - final_cum_post_fatal_lo) / totvacc * 1e6,
      diff_per_million_hi = (final_cum_pre_fatal_hi - final_cum_post_fatal_hi) / totvacc * 1e6
    )%>%
    mutate(
      region = region
    )
  
  return(list(
    fatal_week_df     = fatal_week_df,
    fatal_agegrp3_sum = summary_df
  ))
}

fatal_summ_age <- function(postsim_all_ui, N, observed, region) {
  # sim_result_list: a list of scenario results (each should have a sim_df element)
  
  T = 18
  
  scenario_results  <- postsim_all_ui$scenario_result
  
  mutated_list <- lapply(scenario_results, function(scenario) {
    
    df <- scenario$sim_df
    df <- df %>% mutate(
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      pre_hosp     = pre_Median * hosp_rate,
      pre_hosp_lo  = pre_low95 * hosp_rate,
      pre_hosp_hi  = pre_hi95 * hosp_rate,
      pre_nonhosp  = pre_Median - pre_hosp,
      pre_nonhosp_lo = pre_low95 - pre_hosp_lo,
      pre_nonhosp_hi = pre_hi95 - pre_hosp_hi,
      pre_fatal = (pre_hosp * fatality + pre_nonhosp * nh_fatality),
      pre_fatal_lo = (pre_hosp_lo * fatality + pre_nonhosp_lo * nh_fatality),
      pre_fatal_hi = (pre_hosp_hi * fatality + pre_nonhosp_hi * nh_fatality),
      cum_fatal_pre = cumsum(pre_fatal),
      cum_fatal_pre_lo = cumsum(pre_fatal_lo),
      cum_fatal_pre_hi = cumsum(pre_fatal_hi)
    )
    return(df)
  })
  
  fatal_week_list <- lapply(seq_along(mutated_list), function(i) {
    mutated_list[[i]] %>% 
      group_by(age_gr) %>% 
      summarise(
        fatal        = sum(fatal, na.rm = TRUE),
        fatal_lo     = sum(fatal_lo, na.rm = TRUE),
        fatal_hi     = sum(fatal_hi, na.rm = TRUE),
        pre_fatal    = sum(pre_fatal, na.rm = TRUE),
        pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
        pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      mutate(
        Scenario        = paste0("Scenario_", i),
        cum_fatal       = cumsum(fatal),
        cum_fatal_lo    = cumsum(fatal_lo),
        cum_fatal_hi    = cumsum(fatal_hi),
        cum_fatal_pre   = cumsum(pre_fatal),
        cum_fatal_pre_lo= cumsum(pre_fatal_lo),
        cum_fatal_pre_hi= cumsum(pre_fatal_hi)
      )
  })
  
  # Bind all scenario summaries together
  fatal_week_df <- do.call(rbind, fatal_week_list)
  
  fatal_week_df <- fatal_week_df %>% 
    mutate(
      tot_pop   = rep(N, 3),
      post_fatal_prop = (cum_fatal / tot_pop) * 100,
      post_fatal_prop_lo = (cum_fatal_lo / tot_pop) * 100,
      post_fatal_prop_hi = (cum_fatal_hi / tot_pop) * 100,
      pre_fatal_prop  = (cum_fatal_pre / tot_pop) * 100,
      pre_fatal_prop_lo = (cum_fatal_pre_lo / tot_pop) * 100, 
      pre_fatal_prop_hi = (cum_fatal_pre_hi / tot_pop) * 100
    )
  
  fatal_impact_week <- fatal_week_df %>%
    group_by(Scenario) %>%
    summarise(
      total_post_fatal = sum(cum_fatal),  
      total_post_fatal_lo = sum(cum_fatal_lo),
      total_post_fatal_hi = sum(cum_fatal_hi),
      total_pre_fatal   = sum(cum_fatal_pre),
      total_pre_fatal_lo = sum(cum_fatal_pre_lo),
      total_pre_fatal_hi = sum(cum_fatal_pre_hi),
      .groups = "drop"
    ) %>%
    mutate(
      diff_fatal   = `total_pre_fatal` - `total_post_fatal`,
      impact_fatal = diff_fatal / `total_pre_fatal`  * 100,
      diff_fatal_lo = `total_pre_fatal_lo` - `total_post_fatal_lo`,
      impact_fatal_lo = diff_fatal_lo / `total_pre_fatal_lo`  * 100,
      diff_fatal_hi = `total_pre_fatal_hi` - `total_post_fatal_hi`,
      impact_fatal_hi = diff_fatal_hi / `total_pre_fatal_hi`  * 100
    )
  
  
  annotation_text_fatal <- paste0(
    fatal_impact_week$Scenario, ": ", 
    round(fatal_impact_week$impact_fatal, 1), "%",
    collapse = "\n"
  )
  
  fatal_week_df$region <- region
  
  return(list(mutated_list = mutated_list,
              fatal_week_list   = fatal_week_list,
              fatal_week_df  = fatal_week_df,
              fatal_impact_week = fatal_impact_week,
              annotation_text_fatal = annotation_text_fatal)
  )
}

## cum fatal graph
cum_fatal_plot <-  function(cum_fatal, 
                            postsim_all
                            ) {
  
  fatal_week_df <- cum_fatal$fatal_week_df
  
  annotation_text_fatal <- cum_fatal$annotation_text_fatal
    
  # Compute maximum pre-vaccination fatal proportion for scaling the shading ribbons
  max_fatal <- max(fatal_week_df$pre_fatal_prop, na.rm = TRUE)

  vacc_start_week_s1 <- postsim_all$vacc_weeks$scenario1$start
  vacc_end_week_s1   <- postsim_all$vacc_weeks$scenario1$end
  vacc_start_week_s2 <- postsim_all$vacc_weeks$scenario2$start
  vacc_end_week_s2   <- postsim_all$vacc_weeks$scenario2$end
  vacc_start_week_s3 <- postsim_all$vacc_weeks$scenario3$start
  vacc_end_week_s3   <- postsim_all$vacc_weeks$scenario3$end
  
  p <- ggplot(fatal_week_df) +
    # Scenario 1 shading (Top third: from 2/3*max_fatal to max_fatal)
    geom_ribbon(data = fatal_week_df %>% 
                  filter(Week >= vacc_start_week_s1, Week <= vacc_end_week_s1),
                aes(x = Week, ymin = 0, ymax = Inf), fill = "grey90",
                alpha = 0.3) +
    
    # Scenario 2 shading (Middle third: from 1/3*max_fatal to 2/3*max_fatal)
   # geom_ribbon(data = fatal_week_df %>% 
     #             filter(Scenario == "Scenario_2",
    #                     Week >= vacc_start_week_s2, Week <= vacc_end_week_s2),
      #          aes(x = Week, ymin = (1/3)*max_fatal, ymax = (2/3)*max_fatal, fill = Scenario),
       #         alpha = 0.3) +
    
    # Scenario 3 shading (Bottom third: from 0 to 1/3*max_fatal)
    #geom_ribbon(data = fatal_week_df %>% 
     #             filter(Scenario == "Scenario_3",
      #                   Week >= vacc_start_week_s3, Week <= vacc_end_week_s3),
       #         aes(x = Week, ymin = 0, ymax = (1/3)*max_fatal, fill = Scenario),
        #        alpha = 0.3) +
    
    # 95% UI ribbon for post-fatal proportion (assuming columns: post_fatal_prop_low & post_fatal_prop_hi)
    geom_ribbon(aes(x = Week, ymin = post_fatal_prop_lo, ymax = post_fatal_prop_hi, fill = Scenario),
                alpha = 0.2) +
    
    # Post-vaccination fatal proportion line
    geom_line(aes(x = Week, y = post_fatal_prop, color = Scenario, linetype = "With vaccination"), size = 0.5) +
    
    # Pre-vaccination fatal proportion line (dashed, in black)
    geom_line(aes(x = Week, y = pre_fatal_prop, group = Scenario, linetype = "Without vaccination"), 
              color = "black", size = 0.5) +
    # 95% UI ribbon for pre-vaccination fatal proportion (in gray)
    geom_ribbon(aes(x = Week, ymin = pre_fatal_prop_lo, ymax = pre_fatal_prop_hi),
                fill = "#D8BFD8", alpha = 0.2) +
    
    scale_linetype_manual(values = c("With vaccination" = "solid", "Without vaccination" = "dashed")) +
    
    scale_fill_brewer(palette = "Set1",
                      labels = c("Scenario_1" = "<20 years only", 
                                 "Scenario_2" = "20-59 years only", 
                                 "Scenario_3" = ">60 years only")) +
    
    scale_color_brewer(palette = "Set1",
                       labels = c("Scenario_1" = "<20 years only", 
                                  "Scenario_2" = "20-59 years only", 
                                  "Scenario_3" = ">60 years only")) +
    
    # Vertical dashed line indicating vaccine impact start (adjust as needed)
    geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
    
    # Format the y-axis as percentages
    scale_y_continuous(labels = scales::percent) +
    
    labs(
      title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
      x = "Week",
      y = "Cumulative Fatal Cases (%)",
      color = "Scenario",
      linetype = "Type",
      fill = "Scenario"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      plot.margin = margin(5, 0, 1, 10)
    ) +
    
    # Annotation at the top-right
    annotate("text",
             x = Inf, y = Inf,
             hjust = 1, vjust = 1,
             label = annotation_text_fatal,
             color = "black", size = 3) +
    
    coord_cartesian(clip = "off") +
    
    # Annotation for vaccine impact start near x=4
    annotate("text",
             x = 4, y = max_fatal * 1.1,
             label = "<----- Vaccine impact start (+ 2 wks after initiation)",
             hjust = 0, vjust = 1,
             size = 2, color = "black")
  
  return(p)
}



## vaccine efficacy and contour map --------------------------------------------
contour_ve <- function(target_age_list,
                          observed,
                          posterior, 
                          N, 
                          bra_foi_state_summ, 
                          age_groups, 
                          region_name,
                          hosp, 
                          fatal, 
                          nh_fatal,
                          lhs_sample_young, 
                          lhs_old, 
                          le_sample, 
                          pre_summary_cases) {
  # Set T dynamically 
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # supply and VE
  supply <- seq(0.1, 1, by = 0.05)
  ve  <- seq(0.5, 1, by = 0.05)
  
  # Initialize an empty list to store impact results for each scenario, supply, delay combination.
  ve_impact_data <- list()
  
  # Number of scenarios (each target_age_list element defines a scenario)
  n_scenarios <- length(target_age_list)
  
  # extract median params
  base_beta       <- apply(posterior$base_beta, 2, median)      
  I0              <- apply(posterior$I0, 2, median)           
  gamma           <- median(posterior$gamma)
  rho             <- median(posterior$rho)
  
  # Loop over scenarios, supply rates, and delay steps
  for (scenario_index in seq_along(target_age_list)) {
    target <- target_age_list[[scenario_index]]
    
    for (supply_rate in supply) { 
      for (ve_rate in ve) {
        
        # Run simulation with given parameters
        sim_result <- sirv_sim_coverageSwitch(
          T = T, 
          A = 18,
          N = N,
          r = rep(0, 18),
          base_beta = base_beta,   
          I0 = I0,            
          R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
          rho = rho,               
          gamma = gamma,           
          delay = 2,
          VE_block = ve_rate,
          coverage_threshold = 1,
          target_age = target,
          total_coverage = supply_rate,
          weekly_delivery_speed = 0.1
        )
        
        # Create a data frame from the simulated age-stratified cases
        sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
          mutate(
            Scenario = scenario_index,
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          ) %>%
          mutate(
            hosp_rate = rep(hosp, T),
            hospitalised = Cases * hosp_rate,
            non_hospitalised = Cases - hospitalised,
            fatality = rep(fatal, T),
            nh_fatality = rep(nh_fatal, T),
            fatal = (hospitalised * fatality + non_hospitalised * nh_fatality)
          ) %>%
          mutate(
            age_numeric = case_when(
              str_detect(age_gr, "<5") ~ 0,
              TRUE ~ as.numeric(str_extract(age_gr, "^\\d+"))
            )
          ) %>%
          # (Additional YLD calculations omitted for brevity – include your full chain here)
          mutate(
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
            dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
            dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
          ) %>%
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
          mutate(
            yld_acute = (hospitalised * dw_hosp * dur_acute) +
              (non_hospitalised * dw_nonhosp * dur_acute),
            yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
              (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
            yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
              (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
            yld_total = yld_acute + yld_subacute + yld_chronic
          ) %>%
          mutate(
            le_left = case_when(
              age_numeric %in% c(0, 5) ~ quantile(le_sample$le_1, 0.5),
              age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
              age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
              age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
              age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
              age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
              age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
              age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
              age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
            )
          ) %>%
          mutate(
            yll = fatal * le_left,
            daly_tot = yld_total + yll,
            cum_daly = cumsum(daly_tot)
          )
        
        # Summarize pre-vaccination cases from the simulation output (using a separate summary data frame)
        pre_vacc_cases <- pre_summary_cases %>%
          group_by(Week) %>%
          summarise(
            tot_pre_cases = sum(Median, na.rm = TRUE),
            tot_pre_fatal = sum(fatal, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Summarize post-vaccination cases from sim_df
        post_vacc_cases <- sim_df %>%
          group_by(Week) %>%
          summarise(
            tot_post_cases = sum(Cases, na.rm = TRUE),
            tot_post_fatal = sum(fatal, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Join the two and compute overall impact
        impact_data_summ <- pre_vacc_cases %>%
          left_join(post_vacc_cases, by = "Week") %>%
          summarise(
            tot_pre_cases = sum(tot_pre_cases, na.rm = TRUE),
            tot_post_cases = sum(tot_post_cases, na.rm = TRUE),
            tot_diff = tot_pre_cases - tot_post_cases,
            impact = tot_diff / tot_pre_cases * 100,
            tot_pre_fatal = sum(tot_pre_fatal, na.rm = TRUE),
            tot_post_fatal = sum(tot_post_fatal, na.rm = TRUE),
            diff_fatal = tot_pre_fatal - tot_post_fatal,
            impact_fatal = diff_fatal / tot_pre_fatal * 100,
            .groups = "drop"
          )
        
        # Save results with a unique name in impact_data
        ve_impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_VE_", ve_rate)]] <- list(
          pre_vacc_cases = pre_vacc_cases,
          post_vacc_cases = post_vacc_cases,
          impact_data_summ = impact_data_summ
        )
      }
    }
  }
  
  # Combine impact summaries into one data frame
  ve_impact_data_summ <- do.call(rbind, lapply(names(ve_impact_data), function(idx) {
    list <- ve_impact_data[[idx]]  # Access the element using its name
    df <- list$impact_data_summ           
    df$scenario_id <- idx          
    return(df)                        
  }))
  
  
  # Parse out scenario, supply, and VE
  ve_impact_data_summ <- ve_impact_data_summ %>%
    dplyr::mutate(
      Scenario = sub("_Supply.*", "", scenario_id),
      Supply   = as.numeric(sub("Supply_", "",
                                stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
      VE       = as.numeric(sub("VE_", "",
                                stringr::str_extract(scenario_id, "VE_\\d+(\\.\\d+)?")))
    )
  
  ve_impact_data_summ$region <- region_name
  
  return(ve_impact_data_summ)
}

contour_graph <- function(contour_ve) { 
  
  ve_impact_summ <- contour_ve
  
  p <- 
  ggplot(data = ve_impact_summ, aes(x = Supply, y = VE, z = impact)) +
    geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
    scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
    theme_minimal()+
    geom_vline(xintercept = 0.75, linetype = "dashed", color = "black", size = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.2) +
    facet_wrap(~Scenario) + 
    facet_wrap(~ Scenario, 
               labeller = as_labeller(c("Scenario_1" = "<20 years only", 
                                        "Scenario_2" = "20-59 years only", 
                                        "Scenario_3" = ">60 years only")))+
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))
  
  return(p)
  
  }


rho_estim <- function(posterior) {
  
  rho_mat <- as.matrix(posterior$rho)
  rho_mat <- apply(rho_mat, 2, quantile, probs = c(0.5, 0.025, 0.975))
  
  return(rho_mat)
}


beta_estim <- function(posterior) {
  
  beta_mat <- as.matrix(posterior$base_beta)
  beta_mat <- apply(beta_mat, 2, quantile, probs = c(0.5, 0.025, 0.975))
  
  return(beta_mat)
}

gamma_estim <- function(posterior) {
  
  gamma_mat <- as.matrix(posterior$gamma)
  gamma_mat <- quantile(posterior$gamma, c(0.5, 0.025, 0.975))
  
  return(gamma_mat)
}

# posterior check fns 
prep_ppc_interval_data <- function(obs_list, yrep_list, quantiles = c(0.1, 0.9)) {
  do.call(rbind, lapply(names(obs_list), function(region) {
    y <- as.numeric(obs_list[[region]])
    yrep <- yrep_list[[region]][1:100, ]  # limit to 100 draws
    
    lower <- apply(yrep, 2, quantile, probs = quantiles[1])
    upper <- apply(yrep, 2, quantile, probs = quantiles[2])
    median_pred <- apply(yrep, 2, median)
    
    data.frame(
      index = seq_along(y),
      observed = y,
      lower = lower,
      upper = upper,
      predicted = median_pred,
      region = region
    )
  }))
}

plot_ppc_interval_facet <- function(df) {
  ggplot(df, aes(x = index)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.4) +
    geom_line(aes(y = predicted), color = "blue", linetype = "dashed") +
    geom_point(aes(y = observed), color = "black") +
    facet_wrap(~ region, scales = "free_y") +
    labs(
      x = "Index",
      y = "Cases",
      title = "Posterior Predictive Intervals by Region"
    ) +
    theme_minimal()
}

prep_obs_pred_data <- function(obs_list, yrep_list) {
  do.call(rbind, lapply(names(obs_list), function(region) {
    y_obs <- as.numeric(obs_list[[region]])
    y_pred_median <- apply(yrep_list[[region]], 2, median)
    
    data.frame(
      observed = y_obs,
      predicted = y_pred_median,
      region = region
    )
  }))
}

plot_obs_pred_facet <- function(df) {
  ggplot(df, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ region, scales = "free") +
    labs(
      x = "Observed cases",
      y = "Predicted cases (95% median)"
    ) +
    theme_minimal()
}

compute_all_r2 <- function(obs_list, yrep_list) {
  do.call(rbind, lapply(names(obs_list), function(region) {
    y <- as.numeric(obs_list[[region]])
    yrep <- yrep_list[[region]]
    
    bayes_R2 <- function(y, yrep_draws) {
      var_fit <- apply(yrep_draws, 1, var)
      var_resid <- apply(yrep_draws, 1, function(yp) var(y - yp))
      var_fit / (var_fit + var_resid)
    }
    
    r2_samples <- bayes_R2(y, yrep)
    
    data.frame(
      region = region,
      r2_mean = mean(r2_samples),
      r2_lower = quantile(r2_samples, 0.025),
      r2_upper = quantile(r2_samples, 0.975)
    )
  }))
}

plot_ppc_dens <- function(y, yrep, region_name = "") {
  bayesplot::ppc_dens_overlay(y = as.numeric(y), yrep = yrep[1:100, ]) +
    ggtitle(paste("Posterior Predictive Density:", region_name))
}

plot_ppc_intervals <- function(y, yrep, region_name = "") {
  bayesplot::ppc_intervals(y = as.numeric(y), yrep = yrep[1:100, ]) +
    ggtitle(paste("Predictive Intervals:", region_name))
}

plot_obs_pred <- function(y, yrep, region_name = "") {
  y_pred_median <- apply(yrep, 2, median)
  ggplot(data.frame(observed = as.numeric(y), predicted = y_pred_median), aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Observed Cases",
      y = "Predicted Cases (Posterior Median)",
      title = paste("Observed vs. Predicted:", region_name)
    ) +
    theme_minimal()
}

compute_bayes_R2 <- function(y, yrep) {
  bayes_R2 <- function(y, yrep_draws) {
    var_fit <- apply(yrep_draws, 1, var)
    var_resid <- apply(yrep_draws, 1, function(yp) var(y - yp))
    var_fit / (var_fit + var_resid)
  }
  y <- as.numeric(y)
  r2_samples <- bayes_R2(y, yrep)
  list(
    mean_r2 = mean(r2_samples),
    ci_r2 = quantile(r2_samples, c(0.025, 0.975)),
    all_samples = r2_samples
  )
}
