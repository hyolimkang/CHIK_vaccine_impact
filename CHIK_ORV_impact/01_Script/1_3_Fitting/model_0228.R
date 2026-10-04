posterior_prevacc <- rstan::extract(fit_prevacc_22)

base_beta_21       <- apply(posterior_prevacc$base_beta, 2, median)       # length T
I0_21              <- apply(posterior_prevacc$I0, 2, median)           
gamma_21           <- median(posterior_prevacc$gamma)
rho_21             <- median(posterior_prevacc$rho, 2, median)

base_beta_22       <- apply(posterior_prevacc$base_beta, 2, median)       # length T
I0_22              <- apply(posterior_prevacc$I0, 2, median)           
gamma_22           <- median(posterior_prevacc$gamma)
rho_22             <- median(posterior_prevacc$rho, 2, median)

base_beta <- c(base_beta_21, base_beta_22)
I0 <- c(I0_21, I0_22)
rho <- mean(rho_21, rho_22)
gamma <- mean(gamma_21, gamma_22)

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
  T = 95, 
  A = 18,
  N = N_2023,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  R0  = 1 - exp(-lambda * age_groups),
 # R0 = rep(0, 18),
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
), 95)
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

ggsave(filename = "02_Outputs/2_1_Figures/fig_prevacc_age.jpg", p, width = 7, height = 4, dpi = 1200)




summary_cases_pre <- summary_cases_pre %>%
  mutate(
    hosp_rate = rep(hosp, 44),
    hospitalised = Median * hosp_rate,
    non_hospitalised = Median - hospitalised,
    fatality = rep(fatal, 44),          # Ensure fatal is correctly replicated
    nh_fatality = rep(nh_fatal, 44),
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

n_scenarios <- length(target_age_list)

scenario_result <- list()

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  sim_result <- sirv_sim_coverageSwitch(
    T = 44,
    A = 18,
    N = N_2023,
    r = rep(0, 18),
    base_beta = base_beta,
    I0_draw = I0,
    R0  = 1 - exp(-lambda * age_groups),
    #R0  = rep(0, 18),
    rho = rho,
    gamma = gamma,
    delay = 2,
    VE_block = 0.75,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = 0.5,
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
      hosp_rate         = rep(hosp, 44),
      hospitalised      = Cases * hosp_rate,
      non_hospitalised  = Cases - hospitalised,
      fatality          = rep(fatal, 44),   # Ensure fatal is correctly replicated
      nh_fatality       = rep(nh_fatal, 44),
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
