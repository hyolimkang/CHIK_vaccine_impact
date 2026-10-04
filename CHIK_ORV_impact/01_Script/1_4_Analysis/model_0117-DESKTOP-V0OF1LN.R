source("01_Script/1_1_Functions/sir_functions.R")
source("01_Script/1_1_Functions/sim_functions.R")
source("01_Script/1_1_Functions/library.R")
load("00_Data/0_2_Processed/fit_prevacc.RData")
load("00_Data/0_2_Processed/fit_prevacc_bra23.RData")

# prevacc data------------------------------------------------------------------
posterior_prevacc <- rstan::extract(fit_prevacc)

base_beta       <- apply(posterior_prevacc$base_beta, 2, median)       # length T
I0              <- apply(posterior_prevacc$I0, 2, median)           
gamma           <- median(posterior_prevacc$gamma)
#rho            <- apply(posterior_prevacc$rho, 2, median)
rho             <- median(posterior_prevacc$rho, 2, median)

epi_report_plisa$week <- 1:nrow(epi_report_plisa)

epi_report_plisa <- epi_report_plisa %>%
  mutate(year_week = paste(year, epi_week, sep = "-"))

epi_report_plisa <- epi_report_plisa %>%
  mutate(year_week = factor(year_week, levels = unique(year_week), ordered = TRUE))

brazil_2019 <- epi_report_plisa %>% filter(year == 2019)
brazil_2020 <- epi_report_plisa %>% filter(year == 2020)
brazil_2021 <- epi_report_plisa %>% filter(year == 2021)
brazil_2022 <- epi_report_plisa %>% filter(year == 2022)
brazil_2023 <- epi_report_plisa %>% filter(year == 2023)
brazil_2024 <- epi_report_plisa %>% filter(year == 2024)


N <- as.matrix(brazil_pop_transposed)
N_2023 <- brazil_pop_transposed[, 13]

# cruzeiro de sul
lambda <- 0.0076 # 0.00279100 0.005134
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
  T = 52, 
  A = 18,
  N = N_2023,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  R0  = 1 - exp(-lambda * age_groups),
  #R0 = rep(0, 18),
  rho = rho,
  gamma = gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  weekly_delivery_speed = 0
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
                ), 52)
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
  ylab("Predicted reported symptomatic cases (95% median)")+
  theme_minimal()

  ggsave(filename = "02_Outputs/2_1_Figures/pre_vacc_age.jpg", p, width = 7, height = 4, dpi = 1200)


# post-process epidemics 
# symptomatic cases: assumed that rho multiplied values are already sypmtomatic
# hospitalised cases
# fatal cases 

# can we define fatal cases/hospitalisation occur at the same week of incidence?
# actually, there should be a lag between initial incidence - hospitalisation - fatal
# not exactly 4% of cases will be hospitalised at that week
# the 
summary_cases_pre <- summary_cases_pre %>%
    mutate(
      hosp_rate = rep(hosp, 52),
      hospitalised = Median * hosp_rate,
      non_hospitalised = Median - hospitalised,
      fatality = rep(fatal, 52),          # Ensure fatal is correctly replicated
      nh_fatality = rep(nh_fatal, 52),
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

target_age_list <- list(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), # >60 yrs old
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) # <20 yrs old

target_age_list <- list(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1), # >60 yrs old
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)) # <20 yrs old

n_scenarios <- length(target_age_list)

scenario_result <- list()

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]

  sim_result <- sirv_sim_coverageSwitch(
    T = 52,
    A = 18,
    N = N_2023,
    r = rep(0, 18),
    base_beta = base_beta,
    I0_draw = I0,
    R0  = 1 - exp(-lambda * age_groups),
    #R0  = rep(0, 18),
    rho = rho,
    gamma = gamma,
    delay = 1,
    VE_block = 0.75,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = 0.6,
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
summary_list_df$age_gr <- rep(age_gr, n_scenarios)
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

vacc_start_week_s1 <- scenario_data[[1]]$vacc_start_week[13] # originally 13
vacc_end_week_s1 <- scenario_data[[1]]$vacc_end_week[13]

vacc_start_week_s2 <- scenario_data[[2]]$vacc_start_week[5] # originally 13
vacc_end_week_s2 <- scenario_data[[2]]$vacc_end_week[5]

vacc_start_week_s3 <- scenario_data[[3]]$vacc_start_week[3] # originally 13
vacc_end_week_s3 <- scenario_data[[3]]$vacc_end_week[3]

vacc_start_week_s4 <- scenario_data[[4]]$vacc_start_week[3] # originally 13
vacc_end_week_s4 <- scenario_data[[4]]$vacc_end_week[3]


p <- 
  ggplot(summary_week_df) +
  geom_ribbon(data = summary_week_df[summary_week_df$Week >= vacc_start_week & 
                                    summary_week_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_cases, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "black", size = 0.4)+
  #scale_color_brewer(type = "qual", palette = 1)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, delivery speed: 10%, deployment: week 11",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
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
    "text", x = 13, y = max(summary_week_df$pre_cases, na.rm = TRUE) * 0.9, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig1.jpg", p, width = 7, height = 4, dpi = 1200)

max_cases <- max(summary_week_df$pre_cases, na.rm = TRUE)

p <- 
ggplot(summary_week_df) +
  # Scenario 1 shading
  geom_ribbon(data = summary_week_df %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 3/4 * max_cases, ymax = max_cases, fill = Scenario),
              alpha = 0.4) +
  
  # Scenario 2 shading
  geom_ribbon(data = summary_week_df %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
              aes(x = Week, ymin = 2/4 * max_cases, ymax = 3/4 * max_cases, fill = Scenario),
              alpha = 0.4) +
  
  # Scenario 3 shading
  geom_ribbon(data = summary_week_df %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
              aes(x = Week, ymin = 1/4, ymax = 2/4 * max_cases, fill = Scenario),
              alpha = 0.4) +
  
  # Scenario 3 shading
  #geom_ribbon(data = summary_week_df %>% filter(Scenario == "Scenario_4", Week >= vacc_start_week_s4 & Week <= vacc_end_week_s4),
  #            aes(x = Week, ymin = 0, ymax = 1/4 * max_cases, fill = Scenario),
  #            alpha = 0.4) +
  
  # Post-vaccination case lines with scenario colors
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases"), size = 0.5) +
  
  # Pre-vaccination case lines in gray
  geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Pre Cases"), color = "black", size = 0.5) +
  
  # Maintain correct linetypes
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  
  # Automatically assign colors based on Scenario
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  # Vertical line to indicate vaccine impact start
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4) +
  
  # Labels and titles
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       title = "Coverage: 60%, Delivery Speed: 10%, Deployment: Week 1",
       x = "Week", y = "Predicted reported symptomatic cases (95% median)") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    plot.margin = margin(5, 30, 5, 5)
  ) + 
  
  # Annotation for vaccine impact
  annotate("text", x = 3, y = max(summary_week_df$pre_cases, na.rm = TRUE) * 0.9,
           label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
           hjust = 0, vjust = 1, size = 3, color = "black") + 
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text,  # from the paste step
           color = "black", size = 3)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5))

ggsave(filename = "02_Outputs/2_1_Figures/fig1_upd.jpg", p, width = 8, height = 6, dpi = 1200)

###### age group specific graph #####
global_impact_week <- summary_list_df %>%
  group_by(Scenario, age_gr) %>%
  summarise(
    total_post_cases = sum(Cases),  # Sum cumulative fatal cases for each scenario and type
    total_pre_case   = sum(pre_vacc),
    .groups = "drop"
  ) %>%
  mutate(
    diff = `total_pre_case` - `total_post_cases`,
    impact = diff / `total_pre_case`  * 100
  )

annotation_df <- global_impact_week %>%
  filter(Scenario == 3) %>%
  mutate(label = paste0("Impact: ", round(impact, 1), "%"))

p <- 
  ggplot(summary_list_df[summary_list_df$Scenario == 3,]) +
  geom_ribbon(data = summary_list_df[summary_list_df$Week >= vacc_start_week_s2 & 
                                     summary_list_df$Week <= vacc_end_week_s2, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = Cases, color = factor(Scenario), linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_vacc, color = factor(Scenario), linetype = "Pre Cases"), color = "black")+
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  scale_color_manual(
    name = "Scenario 3",  
    values = c("3" = "#990000")  # Darker shade of blue from Brewer #084594 #990000  #006633
  ) + 
  #scale_color_brewer(palette = "Dark2")+ # Set2, Dark 2, RdYlBu
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4)+
  #scale_color_brewer(type = "qual", palette = 1)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Scenario3: Coverage:60%, Delivery speed: 10%, Deployment: Week 1",  # Add title
    color = "Scenario3", 
    linetype = "Type", 
    x = "Week", 
    y = "Cases"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10),
    legend.position="bottom"
  ) + 
  geom_text(data = annotation_df, aes(x = 40, y = max(summary_list_df[summary_list_df$Scenario == 3,]$Cases, na.rm = TRUE), 
                                      label = label), inherit.aes = FALSE, size = 3, hjust = 1)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5))+
  #annotate(
  #  "text", x = 8, y = max(summary_week_df$pre_cases, na.rm = TRUE), # Adjust x/y for placement
  #  label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
  #  hjust = 0, vjust = 1, 
  #  size = 2, color = "black"
  #)+
  facet_wrap(~age_gr, scales = "free_y")

ggsave(filename = "02_Outputs/2_1_Figures/fig_s3_age_inf.jpg", p, width = 12, height = 6, dpi = 1200)



# cumulative deaths ------------------------------------------------------------

fatal_list <- list()

fatal_summ <- lapply(seq_along(scenario_data), function(index) {
  
  list <- scenario_data[[index]]
  
  age_strat_cases <- as.data.frame.table(list$age_stratified_cases, responseName = "Cases")%>%
    mutate(
      Scenario = paste0("Scenario_", index),
      AgeGroup = as.numeric(Var1),
      Week = as.numeric(Var2),
      hosp_rate    = rep(hosp, 52),
      hospitalised = Cases * hosp_rate,
      non_hospitalised = Cases - hospitalised,
      nh_fatality  = rep(nh_fatal, 52),
      fatality     = rep(fatal, 52)
    )%>% mutate(
      fatal = (hospitalised * fatality + 
               non_hospitalised * nh_fatality),
      pre_fatal = summary_cases_pre$fatal,
      pre_infection = summary_cases_pre$Median,
      pre_hosp      = summary_cases_pre$hospitalised
    )

  fatal_age <- age_strat_cases %>% group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Cases),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = paste0("Scenario_", index) 
    )
  
  fatal_age$age_gr <- age_gr[1:18]
  fatal_age$age_gr <- factor(fatal_age$age_gr[1:18], levels = age_gr_levels)
  
  fatal_all <- age_strat_cases %>% group_by(Week) %>%
    summarise(
      Median   = sum(Cases),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario     = paste0("Scenario_", index),
      cum_fatal    = cumsum(fatal),
      cum_hosp     = cumsum(hospitalised)
    )
  list(age_strat_cases = age_strat_cases, fatal_age = fatal_age, fatal_all = fatal_all)
})

fatal_week_df <- do.call(rbind, lapply(seq_along(fatal_summ), function(index){
  
  df <- fatal_summ[[index]]$fatal_all
  
  return(df)

}))

fatal_week_df <- fatal_week_df %>% mutate(
  # fatal data
  pre_fatal = rep(summary_cases_pre_all$cum_fatal, length(scenario_data)),
  tot_pop   = rep(sum(N_2023), 52*length(scenario_data)),
  post_fatal_prop = cum_fatal / tot_pop,
  pre_fatal_prop  = pre_fatal / tot_pop,
  # hospitalised data
  pre_hosp  = rep(summary_cases_pre_all$cum_hosp, length(scenario_data)),
  post_hosp_prop = cum_hosp / tot_pop,
  pre_hosp_prop  = pre_hosp / tot_pop,
)

fatal_impact_week <- fatal_week_df %>%
  group_by(Scenario) %>%
  summarise(
    total_post_fatal = sum(cum_fatal),  # Sum cumulative fatal cases for each scenario and type
    total_pre_fatal   = sum(pre_fatal),
    total_post_hosp  = sum(cum_hosp),
    total_pre_hosp  = sum(pre_hosp),
    .groups = "drop"
  ) %>%
  mutate(
    diff_fatal   = `total_pre_fatal` - `total_post_fatal`,
    impact_fatal = diff_fatal / `total_pre_fatal`  * 100,
    diff_hosp    = `total_pre_hosp` - `total_post_hosp`,
    impact_hosp  = diff_hosp / `total_pre_hosp`  * 100
  )

annotation_text_fatal <- paste0(
  fatal_impact_week$Scenario, ": ", 
  round(fatal_impact_week$impact_fatal, 1), "%",
  collapse = "\n"
)

annotation_text_hosp <- paste0(
  fatal_impact_week$Scenario, ": ", 
  round(fatal_impact_week$impact_hosp, 1), "%",
  collapse = "\n"
)

## cum fatal graph 
p <- 
  ggplot(fatal_week_df) +
  geom_ribbon(data = fatal_week_df[fatal_week_df$Week >= vacc_start_week & 
                                   fatal_week_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_fatal_prop, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_fatal_prop, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4)+
  scale_y_continuous(labels = percent)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, delivery speed: 10%, Deployment: Week1",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
    y = "Cumulative fatal cases (%)"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10)
  )+
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text_fatal,  # from the paste step
           color = "black", size = 2)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 0, 1, 10))+
  annotate(
    "text", x = 3, y = max(fatal_week_df$pre_fatal_prop, na.rm = TRUE) * 1.1, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig5.jpg", p, width = 6, height = 4, dpi = 1200)

# Calculate the global maximum for fatal proportions from the entire dataset
max_fatal <- max(fatal_week_df$pre_fatal_prop, na.rm = TRUE)

p <- 
ggplot(fatal_week_df) +
  # Scenario 1 shading (Top third: from 2/3*max_fatal to max_fatal)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
    aes(x = Week, ymin = (2/3) * max_fatal, ymax = max_fatal, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Scenario 2 shading (Middle third: from 1/3*max_fatal to 2/3*max_fatal)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
    aes(x = Week, ymin = (1/3) * max_fatal, ymax = (2/3) * max_fatal, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Scenario 3 shading (Bottom third: from 0 to 1/3*max_fatal)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
    aes(x = Week, ymin = 0, ymax = (1/3) * max_fatal, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Plot post-fatal proportion lines with scenario colors
  geom_line(aes(x = Week, y = post_fatal_prop, color = Scenario, linetype = "Post Cases"), size = 0.5) +
  
  # Plot pre-fatal proportion lines in black (grouped by Scenario so lines are separate)
  geom_line(aes(x = Week, y = pre_fatal_prop, group = Scenario, linetype = "Pre Cases"), color = "black", size = 0.5) +
  
  # Specify linetype mapping for the two line types
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  
  # Let ggplot assign colors automatically based on Scenario (using a Brewer palette)
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  # Vertical line indicating vaccine impact start (adjust x-intercept as needed)
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4) +
  
  # Format the y-axis as percentages
  scale_y_continuous(labels = percent) +
  
  # Labels and title
  labs(
    title = "Coverage: 60%, Delivery Speed: 10%, Deployment: Week 1",
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
  
  # Annotation (adjust placement as needed)
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text_fatal,  # your annotation text
           color = "black", size = 2) +
  coord_cartesian(clip = "off") +
  annotate("text",
           x = 3, y = max(fatal_week_df$pre_fatal_prop, na.rm = TRUE) * 1.1,
           label = "<----- Vaccine impact start (+ 2 wks after initiation)",
           hjust = 0, vjust = 1,
           size = 2, color = "black")

ggsave(filename = "02_Outputs/2_1_Figures/fig6.jpg", p, width = 6, height = 4, dpi = 1200)


## cum hosp graph 
p <- 
  ggplot(fatal_week_df) +
  geom_ribbon(data = fatal_week_df[fatal_week_df$Week >= vacc_start_week & 
                                   fatal_week_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_hosp_prop, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_hosp_prop, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4)+
  scale_y_continuous(labels = percent)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, Delivery speed: 10%, Deployment: Week 1",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
    y = "Cumulative hospitalised cases (%)"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10)
  )+
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text_hosp,  # from the paste step
           color = "black", size = 2)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 0, 1, 10))+
  annotate(
    "text", x = 3, y = max(fatal_week_df$pre_hosp_prop, na.rm = TRUE) * 1.1, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )


# Calculate the global maximum for post_hosp_prop across the entire fatal_week_df
max_hosp <- max(fatal_week_df$pre_hosp_prop, na.rm = TRUE)

p <- 
ggplot(fatal_week_df) +
  # Scenario 1 shading (Top third: from 2/3*max_hosp to max_hosp)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
    aes(x = Week, ymin = (2/3) * max_hosp, ymax = max_hosp, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Scenario 2 shading (Middle third: from 1/3*max_hosp to 2/3*max_hosp)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
    aes(x = Week, ymin = (1/3) * max_hosp, ymax = (2/3) * max_hosp, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Scenario 3 shading (Bottom third: from 0 to 1/3*max_hosp)
  geom_ribbon(
    data = fatal_week_df %>% 
      filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
    aes(x = Week, ymin = 0, ymax = (1/3) * max_hosp, fill = Scenario),
    alpha = 0.3
  ) +
  
  # Post-hospitalisation proportion lines, colored by Scenario
  geom_line(aes(x = Week, y = post_hosp_prop, color = Scenario, linetype = "Post Cases"), size = 0.5) +
  
  # Pre-hospitalisation proportion lines in black (grouped by Scenario)
  geom_line(aes(x = Week, y = pre_hosp_prop, group = Scenario, linetype = "Pre Cases"), color = "black", size = 0.5) +
  
  # Set line type mapping
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  
  # Use a Brewer palette for fill and color (this automatically assigns colors by Scenario)
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  # Add vertical dashed line indicating vaccine impact start (Week 8)
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 0.4) +
  
  # Set y-axis to show percentages
  scale_y_continuous(labels = scales::percent_format()) +
  
  # Labels and title
  labs(
    title = "Coverage: 60%, Delivery Speed: 10%, Delay: Week 1",
    x = "Week", 
    y = "Cumulative Hospitalised Cases (%)",
    color = "Scenario",
    linetype = "Type",
    fill = "Scenario"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    plot.margin = margin(5, 0, 1, 10)
  ) +
  
  # Add annotation text for vaccine impact start
  annotate("text", x = 3, y = max(fatal_week_df$pre_hosp_prop, na.rm = TRUE) * 1.1,
           label = "<----- Vaccine impact start (+ 2 wks after initiation)",
           hjust = 0, vjust = 1, size = 2, color = "black") +
  
  coord_cartesian(clip = "off") +
  
  # Optionally, add another annotation (if needed)
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text_hosp,
           color = "black", size = 2)


ggsave(filename = "02_Outputs/2_1_Figures/fig_cumhosp.jpg", p, width = 6, height = 4, dpi = 1200)


# sub-graphs: age stratified cumulative fatal cases by scenario-----------------------------------

fatal_age_list <- lapply(seq_along(fatal_summ), function(idx) {
  
  list <- fatal_summ[[idx]]
  
  age_strat_cases <- list$age_strat_cases
  
  age_strat_cases <-age_strat_cases %>%
    group_by(AgeGroup, Scenario) %>%  # Group by AgeGroup and Scenario
    arrange(Week) %>%                 # Ensure ordered by Week
    mutate(
      cum_fatal = cumsum(fatal),       # Calculate cumulative fatal cases
      pre_fatal = cumsum(pre_fatal)
    ) %>%
    ungroup()
  
   return(age_strat_cases)
})

fatal_age_df <- do.call(rbind, fatal_age_list)

target_df <- as.data.frame(do.call(cbind, target_age_list))

target <- c(rep(target_df[,1], 52),
            rep(target_df[,2], 52),
            rep(target_df[,3], 52))

fatal_age_df <- fatal_age_df %>% mutate(
  target = target
)

# Plot cumulative fatalities by AgeGroup
pre_vacc_target_df <- fatal_age_df %>%
  filter(target == 1)  # Retain rows where 'target' is 1 for pre-vaccination

post_vacc_target_df <- fatal_age_df %>%
  filter(target == 1)  # Retain rows where 'target' is 1 for post-vaccination

combined_df <- bind_rows(
  pre_vacc_target_df %>% mutate(Type = "Pre-Vaccination"),
  post_vacc_target_df %>% mutate(Type = "Post-Vaccination")
)

global_impact_fatal <- combined_df %>%
  group_by(Scenario) %>%
  summarise(
    total_post_cases = sum(cum_fatal),  # Sum cumulative fatal cases for each scenario and type
    total_pre_case   = sum(pre_fatal),
    .groups = "drop"
  ) %>%
  mutate(
    diff = `total_pre_case` - `total_post_cases`,
    impact = diff / `total_pre_case`  * 100
  )

combined_df <- combined_df %>%
  left_join(global_impact_fatal %>% select(Scenario, impact), by = "Scenario")

facet_y_max <- pre_vacc_target_df %>%
  group_by(Scenario) %>%
  summarize(y_max = max(pre_fatal, cum_fatal, na.rm = TRUE), .groups = "drop")

global_impact_fatal$annotation_text <- paste0(
  "Impact: ", round(global_impact_fatal$impact, 1), "%"
) 

global_impact_fatal <- global_impact_fatal %>%
  left_join(facet_y_max, by = "Scenario")

# Create the plot
p <- 
  ggplot() +
  # Plot pre-vaccination cumulative fatal cases with dashed linetype
  geom_line(
    data = pre_vacc_target_df,
    aes(
      x = Week, y = pre_fatal, group = AgeGroup, linetype = "Pre-Vaccination"
    ),
    color = "grey", size = 1, alpha = 0.5
  ) +
  # Plot post-vaccination cumulative fatal cases with solid linetype
  geom_line(
    data = post_vacc_target_df,
    aes(
      x = Week, y = cum_fatal, group = AgeGroup, color = as.factor(AgeGroup), linetype = "Post-Vaccination"
    ),
    size = 1
  ) +
  geom_ribbon(data = pre_vacc_target_df[pre_vacc_target_df$Week >= vacc_start_week & 
                                        pre_vacc_target_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) + 
  facet_wrap(~ Scenario, scales = "free_y", ncol = 3) +  # Facet by Scenario, 3 columns
  #scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentage
  scale_color_discrete(name = "Target Age Group") +  # Legend for Age Groups
  scale_linetype_manual(
    values = c("Pre-Vaccination" = "dashed", "Post-Vaccination" = "solid"),  # Set linetypes
    name = "Type"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Facet title styling
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center title
    legend.position = "bottom"  # Place legend at the bottom
  )+
  geom_text(
    data = global_impact_fatal,
    aes(x = max(pre_vacc_target_df$Week), y = y_max * 1.1, 
        label = annotation_text),
    inherit.aes = FALSE,  # Avoid inheriting other aesthetics
    size = 4, color = "black", hjust = 1, vjust = 1
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 0, 1, 10)) +
  labs(
    title = "Coverage: 60%, weekly delivery speed: 10%",
    y = "Fatal cases"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig7.jpg", plot = p,
       width = 12, height = 6)

# sub-graphs: age stratified infections by scenario-----------------------------------
infection_age_list <- lapply(seq_along(fatal_summ), function(idx) {
  
  list <- fatal_summ[[idx]]
  
  age_strat_cases <- list$age_strat_cases
  
  age_strat_cases <-age_strat_cases %>%
    group_by(AgeGroup, Scenario) %>%  # Group by AgeGroup and Scenario
    arrange(Week) %>%                 # Ensure ordered by Week
    mutate(
      post_infection = Cases,       # Calculate cumulative fatal cases
      pre_infection  = pre_infection
    ) %>%
    ungroup()
  
  return(age_strat_cases)
})

infection_age_df <- do.call(rbind, infection_age_list)

target <- as.data.frame(do.call(cbind, target_age_list))

target <- c(rep(target[,1], 52),
            rep(target[,2], 52),
            rep(target[,3], 52))

infection_age_df <- infection_age_df %>% mutate(
  target = target
)

# Plot cumulative fatalities by AgeGroup
pre_vacc_inf_df <- infection_age_df %>%
  filter(target == 1)  # Retain rows where 'target' is 1 for pre-vaccination

post_vacc_inf_df <- infection_age_df %>%
  filter(target == 1)  # Retain rows where 'target' is 1 for post-vaccination

combined_inf <- bind_rows(
  pre_vacc_inf_df %>% mutate(Type = "Pre-Vaccination"),
  post_vacc_inf_df %>% mutate(Type = "Post-Vaccination")
)

combined_inf <- combined_inf %>% mutate(
  diff = pre_infection - post_infection
)

global_impact_inf <- combined_inf %>%
  filter(target == 1) %>%  # Filter to include only target age groups
  group_by(Scenario) %>%  # Group by Scenario
  summarise(
    total_diff = sum(diff),  # Sum all the weekly differences (pre - post) across age groups and weeks
    total_pre_cases = sum(pre_infection),  # Sum pre-infection cases for the denominator
    global_impact = total_diff / total_pre_cases * 100,  # Calculate the global impact as a percentage
    .groups = "drop"  # Drop grouping for a clean output
  )

combined_inf <- combined_inf %>%
  left_join(global_impact_inf %>% select(Scenario, global_impact), by = "Scenario")

annotation_inf <- combined_inf %>%
  group_by(Scenario) %>%
  summarise(
    max_x = max(Week),               # Max x for annotation placement
    max_y = max(pre_infection, na.rm = TRUE) * 1.05, # Max y for annotation placement
    global_impact = unique(global_impact) # Global impact percentage
  ) %>%
  mutate(label = paste0("Impact: ", round(global_impact, 2), "%"))

reference_x <- annotation_inf$max_x[annotation_inf$Scenario == "Scenario_1"]
reference_y <- annotation_inf$max_y[annotation_inf$Scenario == "Scenario_1"]

annotation_inf <- annotation_inf %>%
  mutate(
    max_x = reference_x,  # Use the same x position for all facets
    max_y = reference_y   # Use the same y position for all facets
  )

total_duration <- scenario_data[[1]]$total_vaccine_duration_age[1]

# Create the plot
custom_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(18)

p <- 
  ggplot() +
  # Plot pre-vaccination cumulative fatal cases with dashed linetype
  geom_line(
    data = pre_vacc_inf_df,
    aes(
      x = Week, y = pre_infection, group = AgeGroup, linetype = "Pre-Vaccination"
    ),
    color = "grey", size = 1, alpha = 0.5
  ) +
  # Plot post-vaccination cumulative fatal cases with solid linetype
  geom_line(
    data = post_vacc_inf_df,
    aes(
      x = Week, y = post_infection, group = AgeGroup, color = as.factor(AgeGroup), linetype = "Post-Vaccination"
    ),
    size = 1
  ) +
  geom_ribbon(data = pre_vacc_inf_df[pre_vacc_inf_df$Week >= vacc_start_week & 
                                     pre_vacc_inf_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) + 
  facet_wrap(~ Scenario,  ncol = 3) +  # Facet by Scenario, 3 columns
  scale_color_manual(
    values = custom_colors,  # Use the Dark2 palette with 18 colors
    name = "Target Age Group"
  ) +
  scale_linetype_manual(
    values = c("Pre-Vaccination" = "dashed", "Post-Vaccination" = "solid"),  # Set linetypes
    name = "Type"
  ) +
  geom_text(
    data = annotation_inf,
    aes(
      x = max_x,  # Place annotation at the max Week
      y = max_y * 1.05,  # Place annotation slightly above the max Value
      label = label
    ),
    inherit.aes = FALSE,
    size = 4, hjust = 1, vjust = 1  # Styling for annotation text
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Facet title styling
    plot.title = element_text(size = 11),  # Center title
    legend.position = "bottom"  # Place legend at the bottom
  )+
  labs(
    title = "Coverage: 60%, weekly delivery speed: 10%",
    y = "Infections"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig8.jpg", plot = p,
       width = 14, height = 7)

pre_vacc_inf_s1 <- pre_vacc_inf_df %>% filter(Scenario == "Scenario_1")
pre_vacc_inf_s2 <- pre_vacc_inf_df %>% filter(Scenario == "Scenario_2")
pre_vacc_inf_s3 <- pre_vacc_inf_df %>% filter(Scenario == "Scenario_3")


## colored lines for targeted and grey lines for untargeted age groups within each scenario
p <- ggplot() +
  # Plot non-targeted (target == 0) post-vaccination cumulative cases with dashed linetype
  geom_line(
    data = infection_age_df %>% filter(target == 0),
    aes(
      x = Week, y = post_infection, group = AgeGroup, linetype = "Untargeted"
    ),
    color = "grey", size = 1, alpha = 0.5
  ) +
  # Plot targeted (target == 1) post-vaccination cumulative cases with solid linetype
  geom_line(
    data = infection_age_df %>% filter(target == 1),
    aes(
      x = Week, y = post_infection, group = AgeGroup, color = as.factor(AgeGroup), linetype = "Targeted"
    ),
    size = 1
  ) +
  # Highlight the intervention period with a ribbon for targeted groups
  geom_ribbon(
    data = infection_age_df %>% filter(target == 1, Week <= vacc_end_week & Week >= vacc_start_week),
    aes(x = Week, ymin = -Inf, ymax = Inf),
    fill = "grey", alpha = 0.3
  ) +
  # Facet by Scenario
  facet_wrap(~ Scenario, ncol = 3) +
  scale_color_manual(
    values = custom_colors,  # Use your Dark2 palette with 18 colors
    name = "Target Age Group"
  ) +
  scale_linetype_manual(
    values = c("Untargeted" = "dashed", "Targeted" = "solid"),  # Set linetypes
    name = "Group"
  ) +
  geom_text(
    data = annotation_inf,
    aes(
      x = max_x,  # Place annotation at the max Week
      y = max_y * 1.05,  # Place annotation slightly above the max Value
      label = label
    ),
    inherit.aes = FALSE,
    size = 4, hjust = 1, vjust = 1  # Styling for annotation text
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Facet title styling
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center title
    legend.position = "bottom"  # Place legend at the bottom
  )


ggsave(filename = "02_Outputs/2_1_Figures/fig8.jpg", plot = p,
       width = 14, height = 7)


# DALY -------------------------------------------------------------------------

daly_summ <- lapply(seq_along(scenario_df), function(index) {
  
  # Avoid using 'list' as a variable name; using 'dat' instead.
  dat <- scenario_df[[index]]
  
  # Summarise by AgeGroup
  daly_age <- dat %>% 
    group_by(AgeGroup) %>%
    summarise(
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot)
    ) %>% 
    mutate(
      cum_yld_acute = cumsum(yld_acute),
      cum_yld_subacute = cumsum(yld_subacute),
      cum_yld_chronic = cumsum(yld_chronic),
      cum_yld_total = cumsum(yld_total),
      cum_yll = cumsum(yll),
      cum_daly = cumsum(daly_tot),
      Scenario = paste0("Scenario_", index)
    )
  
  # Assign age_gr from an external object (ensure these exist)
  daly_age$age_gr <- age_gr[1:18]
  daly_age$age_gr <- factor(daly_age$age_gr[1:18], levels = age_gr_levels)
  
  # Summarise by Week.  
  # Note: If you want to use 'fatal' and 'hospitalised' later, consider summarising them as well.
  daly_all <- dat %>% 
    group_by(Week) %>%
    summarise(
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot)
    ) %>% 
    arrange(Week) %>% 
    mutate(
      cum_yld_acute = cumsum(yld_acute),
      cum_yld_subacute = cumsum(yld_subacute),
      cum_yld_chronic = cumsum(yld_chronic),
      cum_yld_total = cumsum(yld_total),
      cum_yll = cumsum(yll),
      cum_daly = cumsum(daly_tot),
      Scenario = paste0("Scenario_", index)
    )
  
  list(daly_age = daly_age, daly_all = daly_all)
})

daly_week_df <- do.call(rbind, lapply(seq_along(daly_summ), function(index){
  
  df <- daly_summ[[index]]$daly_all
  
  return(df)
  
}))

daly_week_df <- daly_week_df %>% mutate(
  # fatal data
  total_pre_yld_acute    = rep(summary_cases_pre_all$cum_yld_acute, length(scenario_data)),
  total_pre_yld_subacute = rep(summary_cases_pre_all$cum_yld_subacute, length(scenario_data)),
  total_pre_yld_chronic  = rep(summary_cases_pre_all$cum_yld_chronic, length(scenario_data)),
  total_pre_yld_total    = rep(summary_cases_pre_all$cum_yld_total, length(scenario_data)),
  total_pre_yll          = rep(summary_cases_pre_all$cum_yll, length(scenario_data)),
  total_pre_daly_tot     = rep(summary_cases_pre_all$cum_daly_tot, length(scenario_data))
)

daly_impact_week <- daly_week_df %>%
  group_by(Scenario) %>%
  summarise(
    total_pre_yld_total    = sum(total_pre_yld_total),
    total_pre_daly_tot     = sum(total_pre_daly_tot),
    
    total_post_yld_acute    = sum(cum_yld_acute),
    total_post_yld_subacute = sum(cum_yld_subacute),
    total_post_yld_chronic  = sum(cum_yld_chronic),
    total_post_yld_total    = sum(cum_yld_total),
    total_post_yll          = sum(cum_yll),
    total_post_daly_tot     = sum(cum_daly),
    .groups = "drop"
  ) %>%
  mutate(
    diff_yld_total    = total_pre_yld_total - total_post_yld_total,
    impact_yld_total  = diff_yld_total / total_pre_yld_total * 100,
    diff_daly         = total_pre_daly_tot - total_post_daly_tot,
    impact_daly       = diff_daly / total_pre_daly_tot * 100
  )

annotation_text_daly <- paste0(
  daly_impact_week$Scenario, ": ", 
  round(daly_impact_week$impact_daly, 1), "%",
  collapse = "\n"
)

vacc_start_week <- scenario_data[[1]]$vacc_start_week[13]
vacc_end_week <- scenario_data[[1]]$vacc_end_week[13]

## cum fatal graph 
p <- 
  ggplot(daly_week_df) +
  geom_ribbon(data = daly_week_df[daly_week_df$Week >= vacc_start_week & 
                                    daly_week_df$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = cum_daly, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = total_pre_daly_tot, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 8, linetype = "dashed", color = "black", size = 0.4)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, delivery speed: 10%",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
    y = "Cumulative DALY"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10)
  )+
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text_daly,  # from the paste step
           color = "black", size = 2)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 0, 1, 10))+
  annotate(
    "text", x = 8, y = max(daly_week_df$total_pre_daly_tot, na.rm = TRUE) * 1.1, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )

# reduction by supply rate -----------------------------------------------------

reduction_result <- list()

supply <- seq(0, 1, by = 0.05)

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  for (supply_rate in supply) { 
  sim_result <- sirv_sim_coverageSwitch(
    T = 52,
    A = 18,
    N = N_2023,
    r = rep(0, 18),
    base_beta = base_beta,
    I0_draw = I0,
    R0  = 1 - exp(-lambda * age_groups),
    rho = rho,
    gamma = gamma,
    delay = 6,
    VE_block = 0.75,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = supply_rate,
    weekly_delivery_speed = 0.1
  )
  
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
    mutate(
      Scenario            = scenario_index,
      AgeGroup            = as.numeric(Var1),
      Week                = as.numeric(Var2)
    )
  
  reduction_result[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate)]] <- sim_result
  }
}

summ_reduction <- list()

for (scenario_index in seq_along(reduction_result)) {
  
  # Extract the data for the current scenario index using extract_single_scenario_data
  scenario_df <- extract_single_scenario_data(scenario_index, reduction_result)
  
  # Compute median across iterations
  median_summary <- compute_median(scenario_df)
  
  # Get the corresponding scenario name from reduction_result
  scenario_name <- names(reduction_result)[scenario_index]
  
  # Parse the scenario number and supply rate from the scenario name
  parsed_name <- strsplit(scenario_name, "_")[[1]]
  scenario_num <- parsed_name[2]  # Extract "1" from "Scenario_1"
  supply_rate  <- as.numeric(parsed_name[4])  # Extract "0.1" from "Supply_0.1"
  
  # Add scenario and supply rate information to the median summary
  median_summary <- median_summary %>% 
    mutate(
      Scenario = paste0("Scenario_", scenario_num),
      Supply   = supply_rate
    )
  
  # Store the summarized data in summ_reduction
  summ_reduction[[scenario_name]] <- median_summary
}

combined_reduction <- do.call(rbind, lapply(names(summ_reduction), function(name) {
  df <- summ_reduction[[name]]
  df$ScenarioSupply <- name  # Tag with Scenario-Supply identifier
  return(df)
}))

target_agegr <- c(rep(target_df[, 1], 52 * length(supply)),
                  rep(target_df[, 2], 52 * length(supply)),
                  rep(target_df[, 3], 52 * length(supply)))

replicated_pre <- summary_cases_pre %>%
  slice(rep(row_number(), times = length(supply) * length(unique(combined_reduction$Scenario))))

combined_reduction <- combined_reduction %>% mutate(
  target = target_agegr
) %>% mutate(
  pre_vacc = replicated_pre$Median
)

global_impact_summary <- combined_reduction %>%
  filter(target == 1) %>%  # Include only target age groups
  group_by(Scenario, Supply) %>%  # Group by Scenario and Supply
  summarise(
    total_post_cases = sum(MedianCases, na.rm = TRUE),  # Sum of post-vaccination cases
    total_pre_cases = sum(pre_vacc, na.rm = TRUE),      # Sum of pre-vaccination cases
    total_diff = total_pre_cases - total_post_cases,    # Difference (impact)
    global_impact = total_diff / total_pre_cases * 100, # Percentage impact
    .groups = "drop"                                    # Ungroup after summarization
  )

p <- ggplot(global_impact_summary)+
  geom_line(aes(x = Supply, y = global_impact, color = Scenario))+
  theme_minimal()+
  #geom_vline(xintercept = 0.1, linetype = "dashed", color = "black")+
  labs(
    y = "Total impact (%)",
    x = "Coverage",
    title = "Coverage:60%, delivery speed: 10%, deployment: week 6"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'))

ggsave(filename = "02_Outputs/2_1_Figures/fig9.jpg", plot = p,
       width = 5, height = 4)


# by coverage and delay (heatmap generation) -----------------------------------------------------


supply <- seq(0.1, 1, by = 0.05)
delay_steps  <- seq(1, 52, by = 1)

debug_impact_data <- list()
# Iterate over scenarios, supply rates, and delays
for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  for (supply_rate in supply) { 
    for (delay_step in delay_steps) {
      # Run the simulation
      sim_result <- sirv_sim_coverageSwitch(
        T = 52,
        A = 18,
        N = N_2023,
        r = rep(0, 18),
        base_beta = base_beta,
        I0_draw = I0,
        R0  = 1 - exp(-lambda * age_groups),
        rho = rho,
        gamma = gamma,
        delay = delay_step,
        VE_block = 0.75,
        coverage_threshold = 1,
        target_age = target,
        total_coverage = supply_rate,
        weekly_delivery_speed = 0.1
      )
      
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
        mutate(
          Scenario = scenario_index,
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )%>%
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
      
      pre_vacc_cases <- summary_cases_pre %>%
        group_by(Week) %>%
        summarise(tot_pre_cases = sum(Median, na.rm = TRUE),
                  tot_pre_fatal = sum(fatal, na.rm = TRUE), .groups = "drop")
      
      post_vacc_cases <- sim_df %>%
        group_by(Week) %>%
        summarise(tot_post_cases = sum(Cases, na.rm = TRUE),
                  tot_post_fatal = sum(fatal, na.rm = TRUE), .groups = "drop")
      
      impact_data <- pre_vacc_cases %>%
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
      
      debug_impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_Delay_", delay_step)]] <- list(
        pre_vacc_cases = pre_vacc_cases,
        post_vacc_cases = post_vacc_cases,
        impact_data = impact_data
      )
      
    }
  }
}

impact_summ <- do.call(rbind, lapply(names(debug_impact_data), function(idx) {
  list <- debug_impact_data[[idx]]  # Access the element using its name
  df <- list$impact_data           
  df$scenario_id <- idx          
  return(df)                        
}))

impact_summ <- impact_summ %>%
  mutate(
    Scenario = sub("_Supply.*", "", scenario_id),  # Extract only the Scenario part
    Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
    Delay = as.numeric(sub("Delay_", "", stringr::str_extract(scenario_id, "Delay_\\d+(\\.\\d+)?")))
  )

p <- 
  ggplot(impact_summ, aes(x = Supply, y = Delay, fill = impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% Cases averted") +  # Add color legend breaks
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "black", size = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~Scenario)+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig11.jpg", plot = p,
       width = 8, height = 3)

p <- 
  ggplot(impact_summ, aes(x = Supply, y = Delay, fill = impact_fatal)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% Fatal cases averted") +  # Add color legend breaks
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", size = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~Scenario)+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig11_fatal.jpg", plot = p,
       width = 8, height = 3)


supply_60 <- impact_summ %>% filter(Supply == 0.6)

p <- 
ggplot(supply_60)+
  geom_line(aes(x = Delay, y = impact, color = factor(Scenario)))+
  theme_minimal()+
  labs(
    x = "Week",
    y = "Median case reduction (%)",
    color = "Scenario"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig18.jpg", plot = p,
       width = 5, height = 4)


# summarise by scenario ---------------------------------------------------------
summary_scenario_df <- summary_week_df %>% group_by(Scenario) %>% 
                          summarise(
                            post_cases = sum(post_cases),
                            pre_cases  = sum(pre_cases),
                            diff       = pre_cases - post_cases,
                            impact     = (diff / pre_cases) * 100
                          )

# by ve and supply (contour map generation) -----------------------------------------------------

supply <- seq(0.1, 1, by = 0.05)
ve  <- seq(0.1, 1, by = 0.05)

ve_impact_data <- list()
# Iterate over scenarios, supply rates, and delays
for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  for (supply_rate in supply) { 
    for (ve_rate in ve) {
      # Run the simulation
      sim_result <- sirv_sim_coverageSwitch(
        T = 52,
        A = 18,
        N = N_2023,
        r = rep(0, 18),
        base_beta = base_beta,
        I0_draw = I0,
        R0  = 1 - exp(-lambda * age_groups),
        rho = rho,
        gamma = gamma,
        delay = 10,
        VE_block = ve_rate,
        coverage_threshold = 1,
        target_age = target,
        total_coverage = supply_rate,
        weekly_delivery_speed = 0.05
      )
      
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
        mutate(
          Scenario = scenario_index,
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      pre_vacc_cases <- summary_cases_pre %>%
        group_by(Week) %>%
        summarise(tot_pre_cases = sum(Median, na.rm = TRUE), .groups = "drop")
      
      post_vacc_cases <- sim_df %>%
        group_by(Week) %>%
        summarise(tot_post_cases = sum(Cases, na.rm = TRUE), .groups = "drop")
      
      impact_data <- pre_vacc_cases %>%
        left_join(post_vacc_cases, by = "Week") %>%
        summarise(
          tot_pre_cases = sum(tot_pre_cases, na.rm = TRUE),
          tot_post_cases = sum(tot_post_cases, na.rm = TRUE),
          tot_diff = tot_pre_cases - tot_post_cases,
          impact = tot_diff / tot_pre_cases * 100,
          .groups = "drop"
        )
      
      ve_impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_VE_", ve_rate)]] <- list(
        pre_vacc_cases = pre_vacc_cases,
        post_vacc_cases = post_vacc_cases,
        impact_data = impact_data
      )
      
    }
  }
}

ve_impact_summ <- do.call(rbind, lapply(names(ve_impact_data), function(idx) {
  list <- ve_impact_data[[idx]]  # Access the element using its name
  df <- list$impact_data           
  df$scenario_id <- idx          
  return(df)                        
}))

ve_impact_summ <- ve_impact_summ %>%
  mutate(
    Scenario = sub("_Supply.*", "", scenario_id),  # Extract only the Scenario part
    Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
    VE = as.numeric(sub("VE_", "", stringr::str_extract(scenario_id, "VE_\\d+(\\.\\d+)?")))
  )

p <- 
  ggplot(data = ve_impact_summ, aes(x = VE, y = Supply, z = impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 0.2) +
  facet_wrap(~Scenario)

ggsave(filename = "02_Outputs/2_1_Figures/fig_ve.jpg", plot = p,
       width = 8, height = 3)


# raw allocations --------------------------------------------------------------
raw_allocation <- lapply(scenario_data, function(list){
  raw_allocation <- as.data.frame(list$raw_allocation_age)
  raw_allocation$tot_vacc <- rowSums(raw_allocation)
  list$raw_allocation <- raw_allocation
  tot_df <- setNames(as.data.frame(list$raw_allocation$tot_vacc), "tot_vacc")
})

tot_vacc_df <- do.call(rbind, raw_allocation)
tot_vacc_df <- tot_vacc_df %>% filter(tot_vacc > 0)
tot_vacc_df <- tot_vacc_df %>% mutate(
  tot_alloc = c(rep(sum(tot_vacc[1:6]),6),
                rep(sum(tot_vacc[7:14]),8),
                rep(sum(tot_vacc[15:18]),4)
                ),
  vacc_prop = tot_vacc / tot_alloc
)
scenario <- c(rep("Scenario1", 6),
              rep("Scenario2", 8),
              rep("Scenario3", 4)
)
tot_vacc_df$Scenario <- scenario
tot_vacc_df$age_gr <- c(age_gr[13:18],age_gr[5:12],age_gr[1:4])
levels <- c(age_gr[1:4], age_gr[5:12], age_gr[13:18])
tot_vacc_df$age_gr <- factor(tot_vacc_df$age_gr[1:18], levels = levels)


p <- 
ggplot() +
  geom_bar(data = tot_vacc_df, aes(x = age_gr, y = vacc_prop, fill = Scenario), stat = "identity", position = "dodge", alpha = 0.8) +
  facet_wrap(~Scenario, ncol = 1, scales = "free_y") +
  theme_bw()+
  scale_x_discrete(limits = levels)+
  scale_y_continuous(labels = scales::percent)+
  ylab("Weekly vaccine distribution per age group (%)")+
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(size = 10), # Adjust facet label text size
    axis.text.x = element_text(angle = 90, hjust = 1) 
  )+
  labs(x = "Age group")+
  scale_fill_brewer(palette = "Set1")

ggsave(filename = "02_Outputs/2_1_Figures/fig12.jpg", plot = p,
       width = 7, height = 4)

raw_allocation_week <- lapply(scenario_data, function(scenario){
  scenario$raw_allocation_age
})
scenario_names <- c("Scenario1", "Scenario2", "Scenario3")

weekly_allocation_df <- do.call(rbind, lapply(seq_along(raw_allocation_week), function(idx) {
  df <- as.data.frame(raw_allocation_week[[idx]])  # Convert matrix to dataframe
  df$Scenario <- scenario_names[idx]               # Assign scenario label
  df$AgeGroup <- 1:nrow(df)                        # Add age group identifier
  df                                             # Return dataframe
}))

colnames(weekly_allocation_df)[1:(ncol(weekly_allocation_df)-2)] <- as.character(1:(ncol(weekly_allocation_df)-2))

weekly_allocation_df$tot_sum <- rowSums(weekly_allocation_df[,1:52])

weekly_allocation_long <- weekly_allocation_df %>%
  pivot_longer(
    cols = all_of(as.character(1:52)),  # Use the character names of the week columns
    names_to = "Week",
    values_to = "Vaccinated"
  ) %>%
  mutate(
    Week = as.numeric(Week)  # Convert week names to numeric for plotting
  )%>%   # Ensure ordering is correct
  mutate(age_gr = rep(age_gr[1:18], each = 52, times = 3)) %>%
  mutate(age_gr = factor(age_gr, levels = unique(age_gr)))

paired <- brewer.pal(12, "Paired")
extra <- brewer.pal(8, "Dark2")[1:6]
full_palette <- c(paired, extra)

p <- 
ggplot(weekly_allocation_long, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=comma)+
  facet_wrap(~ Scenario) +         # Facet by Scenario if you want separate plots for each
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/fig17.jpg", plot = p,
       width = 8, height = 6)

# foi tracking by scenario -----------------------------------------------------

foi_df <- do.call(rbind, lapply(seq_along(scenario_data), function(idx){
  list <- scenario_data[[idx]]
  foi_vec <- list$phi
  data.frame(
    scenario = paste0("Scenario_",idx),
    week     = seq_along(foi_vec),
    foi      = foi_vec
  )
}))

foi_pre <- data.frame(
  scenario = "Pre-vacc",
  week = seq_along(sim_result_novacc$phi),
  foi = sim_result_novacc$phi
)

foi_comb <- rbind(foi_df, foi_pre)

foi_impact <- foi_df %>% mutate(
  diff = foi_pre$foi - foi_df$foi,
  impact = diff / foi_pre$foi
)

foi_impact_avg <- foi_impact %>% group_by(scenario)%>%
  summarise(
  avg = mean(impact, na.rm = TRUE)
)

p <- 
  ggplot(foi_comb) +
  geom_line(aes(x = week, y = foi, color = scenario, linetype = ifelse(scenario == "Pre-vacc", "Pre-vacc", "Post-vacc")), size = 0.7)+
  geom_vline(xintercept = 12, linetype = "dashed", color = "black", size = 0.4)+
  theme_minimal()+
  scale_linetype_manual(
    name = "Vaccination scenario",
    values = c("Pre-vacc" = "dashed",
               "Post-vacc" = "solid")
  )+
  annotate(
    "text", x = 6, y = max(foi_comb$foi, na.rm = TRUE) * 0.9, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 3, color = "black"
  )+
  labs(y = "FOI")

ggsave(filename = "02_Outputs/2_1_Figures/fig14.jpg", plot = p,
       width = 10, height = 7)


# nnv --------------------------------------------------------------------------

# attack rate nominator: total cumulative new infection = recovered + fatal cases (by age group)
summary_R_pre_age <- summary_R_pre %>%
  filter(Week == 52 | Week == 1) %>%
  group_by(AgeGroup) %>%
  summarize(Difference = Median_rho[Week == 52] - Median_rho[Week == 1])

summary_fatal_pre_age <- fatal_age_df %>%
  filter(Week == 52 | Week == 1)%>%
  group_by(Scenario, AgeGroup) %>%
  summarize(Difference = pre_fatal[Week == 52] - pre_fatal[Week == 1])

raw_allocation_age <- lapply(seq_along(raw_allocation), function(id){
  df <- raw_allocation[[id]]
  tot_vacc <- sum(df$tot_vacc)
  return(tot_vacc)
})


final_summ_df <- do.call(rbind, lapply(seq_along(final_summ), function(idx){
  df <- final_summ[[idx]]
  df <- df %>% mutate(
    tot_vacc = as.numeric(unlist(raw_allocation[[idx]])) 
  )
  df$scenario <- paste0("Scenario_", idx)
  df <- df %>% mutate(
    nnv = tot_vacc / diff,
    nnv_fatal = tot_vacc / diff_fatal,
    nnv_yld_acute    =  tot_vacc / diff_yld_acute,
    nnv_yld_subacutt = tot_vacc / diff_yld_subacute,
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
  tot_pop = rep(N_2023, n_scenarios),
  pre_vacc_10k = pre_vacc / tot_pop * 100000,
  averted_per10k = diff / tot_pop * 100000,
  tot_vacc_per10k = tot_vacc / tot_pop,
  nnv_per = tot_vacc_per10k / averted_per10k
)

final_summ_df$age_gr <- rep(age_gr[1:18],n_scenarios)
final_summ_df$age_gr <- factor(final_summ_df$age_gr, levels = age_gr_levels)


filtered_df <- final_summ_df %>% filter(scenario == "Scenario_1")
filtered_df <- filtered_df %>% mutate(
  age_gr = age_gr[1:18]
)
filtered_df$age_gr <- factor(filtered_df$age_gr, levels = age_gr_levels)

p <- ggplot(filtered_df) +
  geom_bar(aes(x = age_gr, y = pre_vacc_10k, fill = scenario), stat = "identity")+
  labs(
    y = "Pre-vaccination cases per 100k population"
  )+
  scale_y_continuous(labels = comma) + 
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig15.jpg", plot = p,
       width = 10, height = 7)


ggplot(final_summ_df) +
  geom_bar(aes(x = AgeGroup, y = tot_vacc), stat = "identity")+
  facet_wrap(~scenario)

ggplot(final_summ_df) +
  geom_bar(aes(x = AgeGroup, y = tot_vacc_per10k), stat = "identity")+
  facet_wrap(~scenario)

ggplot(final_summ_df) +
  geom_bar(aes(x = AgeGroup, y = diff), stat = "identity")+
  facet_wrap(~scenario)

p <- 
  ggplot(final_summ_df) +
  geom_bar(aes(x = age_gr, y = averted_per10k, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~scenario)+
  scale_y_continuous(labels = comma) +
  theme_minimal()+
  labs(
    y = "Cases averted per 100k population"
  )+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig_avert_per.jpg", plot = p,
       width = 10, height = 6)


p <- 
  ggplot(final_summ_df)+
  geom_bar(aes(x = age_gr, y = nnv_daly, fill = scenario), stat = "identity", alpha = 0.7)+
  facet_wrap(~scenario, nrow = 3)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "NNV to avert DALY",
    x = "Age group",
    title = "Coverage:60%, Delivery speed: 10%, Deployment: Week 6"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(filename = "02_Outputs/2_1_Figures/fig_nnv_4scenarios.jpg", plot = p,
       width = 7, height = 4)

# population dynamics
brazil_2023_pop <- as.data.frame(matrix(NA, nrow = 18, ncol = 2))
names(brazil_2023_pop) <- c("age_gr", "pop")
brazil_2023_pop[, 1] <- 1:18
brazil_2023_pop[, 2] <- brazil_pop_transposed$`Minas Gerais`

ggplot(data = brazil_2023_pop) +
  geom_bar(aes(x = age_gr, y = pop), stat = "identity")

brazil_2023_pop$scenario <- NA
brazil_2023_pop$scenario[1:6] <- 1
brazil_2023_pop$scenario[7:10] <- 2
brazil_2023_pop$scenario[11:18] <- 3

g1 <- filter(brazil_2023_pop, scenario == 1)
g2 <- filter(brazil_2023_pop, scenario == 2)
g3 <- filter(brazil_2023_pop, scenario == 3)

final_summ_nonzero <- final_summ_df %>% filter(nnv != 0)
final_summ_nonzero_s4 <- final_summ_df %>% filter(nnv != 0 & scenario == "Scenario_4")

p1 <- p <- 
ggplot(final_summ_nonzero)+
  geom_bar(aes(x = age_gr, y = nnv_daly, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "NNV to avert a single DALY",
    x = "Age group",
    title = "Coverage:60%, Delivery speed: 10%, Deployment: Week 1"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "02_Outputs/2_1_Figures/nnv_daly.jpg", plot = p,
       width = 5, height = 4)


p2 <- 
  ggplot(final_summ_nonzero)+
  geom_bar(aes(x = age_gr, y = nnv_yld_chronic, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "NNV (YLD chronic)",
    x = "Age group",
    title = "Coverage:60%, weekly speed: 10%"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))

p3 <- 
  ggplot(final_summ_nonzero)+
  geom_bar(aes(x = age_gr, y = nnv_yll, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "NNV (YLL)",
    x = "Age group",
    title = "Coverage:60%, weekly speed: 10%"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))

p4 <- 
  ggplot(final_summ_nonzero)+
  geom_bar(aes(x = age_gr, y = nnv, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "NNV to avert a single symptomatic case",
    x = "Age group",
    title = "Coverage:60%, weekly speed: 10%"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "02_Outputs/2_1_Figures/nnv_symp.jpg", p4,
       width = 8, height = 4)


p1 <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p2 <- p2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Remove legends from all plots except one (keep it in p4)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

# Arrange the plots in a 2x2 grid
final_plot <- (p1 + p2) / (p3 + p4 + plot_layout(guides = "collect")) + 
  plot_annotation(tag_levels = 'A')

ggsave(filename = "02_Outputs/2_1_Figures/final_plot.jpg", plot = final_plot,
       width = 8, height = 4)

p <- 
ggplot(final_summ_nonzero)+
  geom_bar(aes(x = age_gr, y = tot_vacc_per10k, fill = scenario), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "fraction of pop vaccinated",
    x = "Age group",
    title = "Coverage:60%, weekly speed: 10%"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "02_Outputs/2_1_Figures/fig_vaccprop.jpg", plot = p,
       width = 5, height = 4)

group_totals <- brazil_2023_pop %>%
  group_by(scenario) %>%
  summarize(
    total_pop = sum(pop),
    min_x = min(age_gr),   # Leftmost age grouphttp://127.0.0.1:20937/graphics/4b21c84a-f8f4-4bfa-b4b5-04ac74fee866.png
    max_x = max(age_gr),   # Rightmost age group
    label_y = max(brazil_2023_pop$pop) * 1.2,  # Position text above bars
    .groups = "drop"
  )

p <- 
ggplot(brazil_2023_pop)+
  geom_bar(aes(x = age_gr, y = pop, fill = factor(scenario)), stat = "identity", alpha = 0.7)+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = comma)+
  theme_minimal()+
  labs(
    y = "total population",
    x = "Age group"
  )+
  theme(legend.position = "right", 
        title =element_text(size=7, face='bold'),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  #geom_curve(data = group_totals, 
  #           aes(x = min_x, xend = max_x, y = label_y, yend = label_y ),
  #           curvature = 0, size = 0.6, color = "black")
  geom_text(data = group_totals, aes(x = (min_x + max_x) / 2, y = label_y + 500000, 
                                     label = comma(total_pop)), 
            size = 3, vjust = 0)
ggsave(filename = "02_Outputs/2_1_Figures/fig_demog.jpg", plot = p,
       width = 5, height = 4)


## scenario 4
s4_long <- final_summ_nonzero_s4 %>%
  pivot_longer(cols = c(nnv_fatal, nnv_daly, nnv),  # Select columns to gather
               names_to = "NNV_Type",  # New column for facetting
               values_to = "NNV_Value") 

legend_labels <- c("nnv_fatal" = "NNV to avert a fatal case", 
                   "nnv_daly" = "NNV to avert a DALY",
                   "nnv" = "NNV to avert a case")
custom_labels <- c(
  "nnv" = "NNV to avert a single symptomatic case",
  "nnv_daly" = "NNV to avert a single DALY",
  "nnv_fatal" = "NNV to avert a single fatal case"
)

p <- 
ggplot(s4_long, aes(x = age_gr, y = NNV_Value, fill = NNV_Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +  # Change "dodge" to "stack" if needed
  facet_wrap(~NNV_Type, scales = "free_y", labeller = labeller(NNV_Type = custom_labels)) +  # Separate facets for each NNV type
  scale_fill_brewer(palette = "Set1", labels = legend_labels) +  
  scale_y_continuous(labels = comma) +
  labs(title = "Coverage:60%, Delivery speed: 10%, Deployment: Week 6",
       x = "Age Group",
       y = "NNV",
       fill = "NNV Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 9)) 

ggsave(filename = "02_Outputs/2_1_Figures/fig_nnv_allvacc.jpg", plot = p,
       width = 12, height = 4)


# susceptible over time by age group to check diminishing returns---------------

s_pop <- lapply(seq_along(scenario_data), function(idx){
  
  list <- scenario_data[[idx]]$S
  
  s_df <- as.data.frame.table(list, responseName = "susceptible") %>%
    mutate(
      Scenario            = paste0("Scenario_", idx),
      AgeGroup            = as.numeric(Var1),
      Week                = as.numeric(Var2)
    )
  
  return(s_df)

})

s_1_data <- s_pop[[1]]  # Select the first scenario
s_1_data <- s_1_data %>% mutate(
  target = rep(c(target_df[,1]), 52)
)

s_2_data <- s_pop[[2]]  # Select the first scenario
s_2_data <- s_2_data %>% mutate(
  target = rep(c(target_df[,2]), 52)
)


ggplot(s_1_data, aes(x = Week, y = susceptible, group = AgeGroup, color = ifelse(target == 1, as.factor(AgeGroup), "Non-Target"))) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    x = "Week",
    y = "Susceptible Population",
    color = "Age Group"
  ) +
  scale_color_manual(
    values = c(
      setNames(scales::hue_pal()(length(unique(s_1_data$AgeGroup[s_1_data$target == 1]))), 
               unique(as.character(s_1_data$AgeGroup[s_1_data$target == 1]))), # Automatic colors for target groups
      "Non-Target" = "grey" # Gray for non-target groups
    )
  )  +
  theme(legend.position = "bottom")


# nnv heatmap by different coverag level and weekly speed -------------------------------

delay <- seq(1, 52, by = 1) 
supply  <- seq(0.1, 1, by = 0.1)

nnv_list <- list()

# Iterate over scenarios, supply rates, and delays
for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  for (supply_rate in supply) { 
    for (delay_step in delay) {
      
      scenario_name <- paste0("Scenario_", scenario_index, 
                              "_Supply_", round(supply_rate, 2), 
                              "_Delay_", round(delay_step, 2))
      
      # Run the simulation
      sim_result <- sirv_sim_coverageSwitch(
        T = 52,
        A = 18,
        N = N_region1,
        r = rep(0, 18),
        base_beta = base_beta,
        I0_draw = I0,
        R0  = 1 - exp(-lambda * age_groups),
        rho = rho,
        gamma = gamma,
        delay = delay_step,
        VE_block = 1,
        coverage_threshold = 1,
        target_age = target,
        total_coverage = supply_rate,
        weekly_delivery_speed = 0.05
      )
      
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
        mutate(
          Scenario = scenario_index,
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      tot_supply <- data.frame(
        AgeGroup = seq_along(sim_result$total_supply_age),  # Assign indices as AgeGroup
        Supply = sim_result$total_supply_age               # Supply values
      )
      
      pre_vacc_cases <- summary_cases_pre %>%
        group_by(AgeGroup) %>%
        summarise(tot_pre_cases = sum(Median, na.rm = TRUE), .groups = "drop")
      
      post_vacc_cases <- sim_df %>%
        group_by(AgeGroup) %>%
        summarise(tot_post_cases = sum(Cases, na.rm = TRUE), .groups = "drop")
      
      nnv_data <- pre_vacc_cases %>%
        left_join(post_vacc_cases, by = "AgeGroup") %>%
        mutate(tot_vacc = tot_supply$Supply,
               diff = tot_pre_cases - tot_post_cases,
               nnv  = tot_vacc / diff,
               scenario = scenario_name) 
      
      nnv_list[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_Delay_", delay_step)]] <- list(
        pre_vacc_cases = pre_vacc_cases,
        post_vacc_cases = post_vacc_cases,
        nnv_data = nnv_data
      )
      
    }
  }
}

# extract only nnv data
nnv_summ <- do.call(rbind, lapply(seq_along(nnv_list), function(idx){
  list <- nnv_list[[idx]]$nnv_data
  df <- list %>% summarise(
    tot_pre_cases  = sum(tot_pre_cases),
    tot_post_cases = sum(tot_post_cases),
    diff  = tot_pre_cases - tot_post_cases,
    tot_vacc = sum(tot_vacc),
    nnv = tot_vacc / diff,
    scenario = first(scenario)
  ) %>% mutate(
    Scenario = sub("_Supply.*", "", scenario),  # Extract only the Scenario part
    Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario, "Supply_\\d+(\\.\\d+)?"))),
    Delay = as.numeric(sub("Delay_", "", stringr::str_extract(scenario, "Delay_\\d+(\\.\\d+)?")))
  )
  
  return(df)
}))
nnv_summ <- nnv_summ %>% filter(!is.nan(nnv))

p <- ggplot(nnv_summ, aes(x = Supply, y = Delay, fill = nnv)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "NNV") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~Scenario)+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig13.jpg", plot = p,
       width = 8, height = 3)


# run multiple years -----------------------------------------------------------
base_beta_2yr <- rep(base_beta, 2)  # length 156

## no vacc
sim_result_novacc_2yrs <- sirv_sim_coverageSwitch(
  T = 52*2, 
  A = 18,
  N = N_2023,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta_2yr,
  I0_draw = I0,
  #R0  = 1 - exp(-lambda * age_groups),
  R0 = rep(0, 18),
  rho = rho,
  gamma = gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,             # Vaccine efficacy
  target_age = c(rep(0,18)),  # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  weekly_delivery_speed = 0
)

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc_2yrs$age_stratified_cases, responseName = "Median")

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
##

# 2) Set total simulation time to 2 years * 52 weeks
n_years <- 2
T_2yr <- 52 * 2

target_age_list <- list(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1), # >60 yrs old
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) # <20 yrs old

scenario_result_2yrs <- list()

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  sim_result <- sirv_sim_coverageSwitch(
    T = T_2yr,
    A = 18,
    N = N_2023,
    r = rep(0, 18),
    base_beta = base_beta_2yr,
    I0_draw = I0,
    #R0  = 1 - exp(-lambda * age_groups),
    R0 = rep(0, 18),
    rho = rho,
    gamma = gamma,
    delay = 40,
    VE_block = 0.75,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = 0.6,
    weekly_delivery_speed = 0.05
  )
  
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
    mutate(
      Scenario            = scenario_index,
      AgeGroup            = as.numeric(Var1),
      Week                = as.numeric(Var2)
    ) %>% mutate(
      hosp_rate = rep(hosp, 52*n_years),
      hospitalised = Cases * hosp_rate,
      non_hospitalised = Cases - hospitalised,
      fatality = rep(fatal, 52*n_years),          # Ensure fatal is correctly replicated
      nh_fatality = rep(nh_fatal, 52*n_years),
      fatal = (hospitalised * fatality +   
                 non_hospitalised * nh_fatality)
    )
  
  scenario_result_2yrs[[scenario_index]] <- list(sim_result = sim_result, 
                                            sim_df     = sim_df)
}

scenario_data_2yrs <- lapply(seq_along(scenario_result_2yrs), function(idx){
  list <- scenario_result_2yrs[[idx]]
  data <- list$sim_result
  return(data)
})

summary_list_2yrs <- lapply(seq_along(scenario_result_2yrs), function(idx){
  list <- scenario_result_2yrs[[idx]]
  df <- list$sim_df
  df$pre_vacc <- summary_cases_pre$Median
  df$pre_fatal <- summary_cases_pre$fatal
  df <- df %>% mutate(
    diff     = pre_vacc - Cases,
    impact   = diff / pre_vacc * 100
  )
  return(df)
})

summary_list_df_2yrs <- do.call(rbind, summary_list_2yrs)
summary_list_df_2yrs$age_gr <- rep(age_gr, 3)
summary_list_df_2yrs$age_gr <- factor(summary_list_df_2yrs$age_gr, levels = age_gr_levels)

final_summ_2yrs <- lapply(summary_list_2yrs, function(df) {
  df %>%
    group_by(AgeGroup) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(Cases, na.rm = TRUE),  # Summarize MedianCases
      fatal  = sum(fatal),
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_cases_pre_age$Median,
      pre_fatal = summary_cases_pre_age$fatal,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100,
      diff_fatal = pre_fatal - fatal,
      impact_fatal = diff_fatal / pre_fatal * 100
    ) 
})

summary_week_2yrs <- lapply(seq_along(summary_list_2yrs), function(i) {
  summary_list_2yrs[[i]] %>% 
    mutate(Scenario = paste0("Scenario_", i),
           pre_vacc = summary_cases_pre$Median) %>% group_by(Week, Scenario) %>%
    summarise(post_cases = sum(Cases),
              pre_cases  = sum(pre_vacc),
              diff       = pre_cases - post_cases,
              Scenario   = first(Scenario)
    ) 
  
})

summary_week_df_2yrs <- do.call(rbind, summary_week_2yrs)

global_impact_week_2yrs <- summary_week_df_2yrs %>%
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
  global_impact_week_2yrs$Scenario, ": ", 
  round(global_impact_week_2yrs$impact, 1), "%",
  collapse = "\n"
)

vacc_start_week <- scenario_data_2yrs[[1]]$vacc_start_week[13]
vacc_end_week <- scenario_data_2yrs[[1]]$vacc_end_week[13]

p <- 
  ggplot(summary_week_df_2yrs) +
  geom_ribbon(data = summary_week_df_2yrs[summary_week_df_2yrs$Week >= vacc_start_week & 
                                       summary_week_df_2yrs$Week <= vacc_end_week, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_cases, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "black", size = 0.4)+
  #scale_color_brewer(type = "qual", palette = 1)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:60%, delivery speed: 10%, deployment: week 6",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
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
    "text", x = 42, y = max(summary_week_df$pre_cases, na.rm = TRUE) * 0.9, # Adjust x/y for placement
    label = "<----- Vaccine impact start (+ 2 wks after initiation)", 
    hjust = 0, vjust = 1, 
    size = 2, color = "black"
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig_2yr_outbreaks_delay40.jpg", plot = p,
       width = 8, height = 3)
