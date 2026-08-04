N <- as.matrix(brazil_pop_transposed)
N_2023 <- brazil_pop_transposed[, 13]
N_tot <- brazil_pop_transposed$tot_pop

bra_foi <- allfoi %>% filter(country == "Brazil")
lon_target <- -43.7908
lat_target <- -17.9302
tol <- 0.05

# total time step
T = 365

# Filter the dataframe using subset
bra_mg <- subset(bra_foi, abs(x - lon_target) <= tol & abs(y - lat_target) <= tol)

foi_bra_mg <- mean(bra_mg$foi_mid)

# age groups
age_groups <- c(mean(0:4),
                mean(5:9),
                mean(10:14),
                mean(15:19),
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

sim_result_novacc <- sirv_sim_coverageSwitch(
  T = 365, 
  A = 18,
  N = N_2023,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = beta_time,
  I0_draw = rep(1, 18),
  R0  = 1 - exp(-foi_bra_mg * age_groups),
  rho = 0.30,
  gamma = 1 / 7,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  daily_delivery_speed = 0
)

age_gr <- rep(c("<5 years",
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
), T)
age_gr_levels <- c("<5 years",
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
                   "85-89 years")

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Median")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )
summary_cases_pre$age_gr <- age_gr
summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)

summary_cases_pre$Scenario <- "Pre-Vaccination"

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_R_pre <- as.data.frame.table(sim_result_novacc$R, responseName = "R")
colnames(summary_R_pre) <- c("AgeGroup", "Week","Median")

summary_R_pre <- summary_R_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week),
    Median_rho = Median * rho
  )

summary_R_pre$Scenario <- "Pre-Vaccination"

summary_R_pre_all <- summary_R_pre %>% group_by(Week) %>% 
  summarise(
    Median_raw = sum(Median),
    Median   = Median_raw * rho
  )

p <- 
  ggplot()+
  geom_line(data = summary_cases_pre, aes(x = Week, y = Median, color = Scenario))+
  facet_wrap(~age_gr)+
  theme_bw()+
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/pre_vacc_age.jpg", p, width = 7, height = 4, dpi = 1200)

### loop through all possible beta_time
N_list <- as.list(brazil_pop_transposed[,1:27])
  
sim_results_list <- vector("list", length(beta_time_list))

for (i in seq_along(beta_time_list)) {
  cat(sprintf("Running simulation for beta_time_list element %d\n", i))
  
  sim_results_list[[i]] <- sirv_sim_coverageSwitch(
    T = 365,                          # Number of weeks (here, actually days if running daily)
    A = 18,                           # Number of age groups
    N = N_list[[i]],                       # Population by age
    r = rep(0, 18),                   # Aging rate (set to 0 for all age groups)
    base_beta = beta_time_list[[i]],  # Use the ith beta_time vector
    I0_draw = rep(1, 18),             # Initial infected in each age group
    R0 = 1 - exp(-foi_state_list[[i]] * age_groups),  # Underlying immunity from catalytic model
    rho = 0.30,                       # Detection probability ~ 0.3
    gamma = 1 / 7,                    # Recovery rate (weekly or daily? In your current weekly model, this should be consistent)
    delay = 53,                       # Vaccination starts at week 53 (or day 53 if daily)
    VE_block = 0,                     # Vaccine efficacy (here 0, meaning no vaccine effect)
    target_age = rep(0, 18),          # Targeting specific age groups (all 0 here)
    coverage_threshold = 0,
    total_coverage = 0,
    daily_delivery_speed = 0
  )
}

summary_list <- lapply(seq_along(sim_results_list), function(i) {
  sim_res <- sim_results_list[[i]]
  # Convert the age-stratified cases table to a data frame
  summary_cases <- as.data.frame.table(sim_res$age_stratified_cases, responseName = "Median")
  colnames(summary_cases) <- c("AgeGroup", "Day", "Median")
  
  # Convert AgeGroup and Week to numeric
  summary_cases <- summary_cases %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Day = as.numeric(Day),
      State = paste0("State", i),
      age_gr = rep(age_gr, length.out = nrow(summary_cases))
    )
  
  return(summary_cases)
})
for(i in seq_along(summary_list)) {
  summary_list[[i]]$FOI <- foi_list[i]
}
summary_cases_all <- bind_rows(summary_list)

# Optionally, if you also want a summary aggregated over age groups per week:
summary_cases_all_agg <- summary_cases_all %>%
  group_by(State, Day) %>%
  summarise(Median = sum(Median), .groups = "drop")

ggplot(summary_cases_all, aes(x = Day, y = Median, color = State))+
  geom_line() +
  facet_wrap(~age_gr) +
  theme_bw() +
  theme_minimal()
# post-process epidemics --------
# symptomatic cases: assumed that rho multiplied values are already sypmtomatic
# hospitalised cases
# fatal cases 

# can we define fatal cases/hospitalisation occur at the same week of incidence?
# actually, there should be a lag between initial incidence - hospitalisation - fatal
# not exactly 4% of cases will be hospitalised at that week
# the 
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

summary_cases_pre_age$age_gr <- age_gr[1:18]
summary_cases_pre_age$age_gr <- factor(summary_cases_pre$age_gr[1:18], levels = age_gr_levels)

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median),
    hospitalised = sum(hospitalised),
    fatalilty    = first(fatality),
    fatal        = sum(fatal)
  ) %>% mutate(
    Scenario = "Pre-vaccination",
    cum_fatal    = cumsum(fatal),
    cum_hosp     = cumsum(hospitalised)
  )


# daly estiamtion

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

# loop through strategy ------------------------------------------------------
target_age_list <- list(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1), # >60 yrs old
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) # <20 yrs old

scenario_result <- list()

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  sim_result <- sirv_sim_coverageSwitch(
    T = T,
    A = 18,
    N = N_2023,
    r = rep(0, 18),
    base_beta = beta_time,
    I0_draw = rep(1, 18),
    R0  = 1 - exp(-foi_bra_mg * age_groups),
    rho = 0.3,
    gamma = 1/7,
    delay = 100,
    VE_block = 0.75,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = 0.6,
    daily_delivery_speed = 0.01
  )
  
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
    # First, create the basic data frame
    mutate(
      Scenario  = scenario_index,
      AgeGroup  = as.numeric(Var1),
      Week      = as.numeric(Var2)
    ) %>%
    mutate(
      hosp_rate         = rep(hosp, T),
      hospitalised      = Cases * hosp_rate,
      non_hospitalised  = Cases - hospitalised,
      fatality          = rep(fatal, T),   # Ensure fatal is correctly replicated
      nh_fatality       = rep(nh_fatal, T),
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
  
  scenario_result[[scenario_index]] <- list(sim_result = sim_result, 
                                            sim_df     = sim_df)
}

scenario_data <- lapply(seq_along(scenario_result), function(idx){
  list <- scenario_result[[idx]]
  data <- list$sim_result
  return(data)
})

scenario_df <- lapply(seq_along(scenario_result), function(idx){
  list <- scenario_result[[idx]]
  data <- list$sim_df
  return(data)
})


summary_list <- lapply(seq_along(scenario_result), function(idx){
  list <- scenario_result[[idx]]
  df <- list$sim_df
  df$pre_vacc <- summary_cases_pre$Median
  df$pre_fatal <- summary_cases_pre$fatal
  df <- df %>% mutate(
    diff     = pre_vacc - Cases,
    impact   = diff / pre_vacc * 100
  )
  return(df)
})

summary_list_df <- do.call(rbind, summary_list)
summary_list_df$age_gr <- rep(age_gr, 3)
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
      pre_vacc = summary_cases_pre_age$Median,
      pre_yld_acute = summary_cases_pre_age$yld_acute,
      pre_yld_subacute = summary_cases_pre_age$yld_subacute,
      pre_yld_chronic = summary_cases_pre_age$yld_chronic,
      pre_yld_total = summary_cases_pre_age$yld_total,
      pre_yll = summary_cases_pre_age$yll,
      pre_daly = summary_cases_pre_age$daly_tot,
      pre_fatal = summary_cases_pre_age$fatal,
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
           pre_vacc = summary_cases_pre$Median) %>% group_by(Week, Scenario) %>%
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
    total_post_cases = sum(post_cases),  # Sum cumulative fatal cases for each scenario and type
    total_pre_case   = sum(pre_cases),
    .groups = "drop"
  ) %>%
  mutate(
    diff = `total_pre_case` - `total_post_cases`,
    impact = diff / `total_pre_case`  * 100
  )

annotation_text <- paste0(
  global_impact_week$Scenario, ": ", 
  round(global_impact_week$impact, 1), "%",
  collapse = "\n"
)

vacc_start_day <- scenario_data[[1]]$vacc_start_day[13] 
vacc_end_day <- scenario_data[[1]]$vacc_end_day[13]

p <- 
  ggplot(summary_week_df) +
  geom_ribbon(data = summary_week_df[summary_week_df$Week >= vacc_start_day & 
                                       summary_week_df$Week <= vacc_end_day, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_cases, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 100 + 15, linetype = "dashed", color = "black", size = 0.4)+
  #scale_color_brewer(type = "qual", palette = 1)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, delivery speed: 10%, deployment: week 6",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Day", 
    y = "Cases"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10)
  ) + 
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text,  # from the paste step
           color = "black", size = 3)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5))+
  annotate(
    "text", x = 100 + 15, y = max(summary_week_df$pre_cases, na.rm = TRUE) * 0.9, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )

