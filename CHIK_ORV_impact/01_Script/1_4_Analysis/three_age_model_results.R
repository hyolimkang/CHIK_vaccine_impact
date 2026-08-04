# cruzeiro de sul
lambda <- 0.003843 # 0.00279100 0.005134
# 3 age groups
age_groups <- c(mean(0:20),
                mean(21:59),
                mean(60:89)
)

N_agegr <- c(sum(N_region1[1:4]),
             sum(N_region1[5:12]),
             sum(N_region1[13:18])
)

sim_result_novacc <- sirv_sim_coverageSwitch_Fatal(
  T = 52, 
  A = 3,
  N = N_agegr,
  r = rep(0, 3),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  R0  = 1 - exp(-lambda * age_groups),
  rho = rho,
  gamma = gamma,
  
  delay = 52,                   # Vaccination starts from Week x
  VE_block = 0,             # Vaccine efficacy
  target_age = c(rep(0,3)),  # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  weekly_delivery_speed = 0,
  VE_p = 0, 
  FR_infection = fatal
)

age_gr <- rep(c("< 20 years old",
            "20-59 years old",
            "> 60 years old"), 52)
age_gr_levels <- c("< 20 years old", "20-59 years old", "> 60 years old")

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Median")
summary_Ifatal_pre <- as.data.frame.table(sim_result_novacc$IFatal_detect, responseName = "IFatal")
summary_IVfatal_pre <- as.data.frame.table(sim_result_novacc$IVFatal_detect, responseName = "IVFatal")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_Ifatal_pre) <- c("AgeGroup", "Week","Median")
colnames(summary_IVfatal_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

# Convert the age_gr column to a factor with the specified order
summary_cases_pre$age_gr <- age_gr
summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)

summary_cases_pre$Scenario <- "Pre-Vaccination"

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_Ifatal_pre <- summary_Ifatal_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_Ifatal_pre$Scenario <- "Pre-Vaccination"

summary_Ifatal_pre_all <- summary_Ifatal_pre %>% group_by(Week) %>%
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

ggplot()+
  geom_line(data = summary_cases_pre, aes(x = Week, y = Median, color = Scenario))+
  facet_wrap(~age_gr)+
  theme_bw()

