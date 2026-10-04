
#-------------------------------------------------------------------------------
# Data 
#-------------------------------------------------------------------------------
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

p <- ggplot(data = epi_report_plisa) +
  geom_point(aes(x = year_week, y = incidence, color = factor(year), group = year)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(breaks = levels(epi_report_plisa$year_week)[seq(1, nlevels(epi_report_plisa$year_week), by = 10)])+
  labs(color = "Year")

ggsave(filename = "02_Outputs/2_1_Figures/bra_201924.jpg", p, width = 11, height = 6, dpi = 1200)

observed_cases <- brazil_2023$incidence
T <- length(observed_cases)  

N_tot <- as.matrix(brazil_pop_transposed$tot_pop)
N_2023 <- as.matrix(brazil_pop_transposed[,13]) # minas gerais

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_bra_v3.stan")
indexP_single <- indexP_matrix[,1]
indexP_singlemat <- matrix(indexP_single, nrow=52, ncol=1)

stan_data_vacc <- list(
  T = length(observed_cases),             # Number of weeks
  A = 18,                                 # Number of age groups
  N = as.vector(N_region1),                                  # Population size (age x regions)
  observed_cases = observed_cases,        # Weekly observed cases
  #r = rep(1 / (5 * 52), 18),             # Ageing rates for each group
  r = rep(0, 18),
  indexP = indexP_matrix[,1],           # Vector capacity indices for regions
  delay = 1,                              # Vaccination starts in week T
  VE_block = 1,                        # Vaccine efficacy
  indexP = as.vector(indexP_matrix[,1]),
  eta = 0.4,
  vaccination_rate = 1 - (0.6)^(1/6),  # Weekly vaccination rate for 40% over 6 weeks
  vaccine_target_age = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0)  # Targeting age group 4
)

stan_data_prevacc <- list(
  T = length(observed_cases),
  A = 18,
  N = as.vector(N_2023), 
  observed_cases = observed_cases,
  r = rep(1 / (5 * 52), 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0, 52),
  eta = 0
)


fit_prevacc <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc,
  iter = 1000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 100,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

save(fit_prevacc, file = "00_Data/0_2_Processed/fit_prevacc_bra23.RData")

fit_vacc <- sampling(
  object = stan_model_age,
  data = stan_data_vacc,
  iter = 1000,            # Reduced from 4000
  chains = 2,             # Reduced from 4
  warmup = 500,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

posterior_cases_prevacc <- extract(fit_prevacc)
age_cases_prevacc <- posterior_cases_prevacc$age_stratified_cases

summary_cases_prevacc <- apply(age_cases_prevacc, c(2, 3), function(x) {
  c(
    meain = median(x, na.rm = TRUE),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975)
  )
})
rownames(summary_cases_prevacc) <- c("median", "lower", "upper")

age_groups <- 1:18  # Age group labels
weeks <- 1:52       # Weeks

df_prevacc <- expand.grid(
  Week = weeks,
  AgeGroup = age_groups
)

# Extract mean, lower, and upper from the summary array
df_prevacc$Median <- as.vector(aperm(summary_cases_prevacc["median", , ], c(2, 1)))  # Weeks x Age
df_prevacc$Lower <- as.vector(aperm(summary_cases_prevacc["lower", , ], c(2, 1)))   # Weeks x Age
df_prevacc$Upper <- as.vector(aperm(summary_cases_prevacc["upper", , ], c(2, 1)))   # Weeks x Age
df_prevacc$Scenario <- "Pre-Vaccination"

ggplot(df_prevacc, aes(x = Week, y = Median, color = Scenario)) +
  geom_line(size = 1) +
  ylab("reported symptomatic cases (95%mid)")+
  facet_wrap(~ AgeGroup, ncol = 3)+
  theme_bw()

#-------------------------------------------------------------------------------
posterior_vacc <- extract(fit_vacc)
age_cases_vacc <- posterior_vacc$age_stratified_cases

summary_cases_vacc <- apply(age_cases_vacc, c(2, 3), function(x) {
  c(
    meain = median(x, na.rm = TRUE),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975)
  )
})
rownames(summary_cases_vacc) <- c("median", "lower", "upper")

age_groups <- 1:18  # Age group labels
weeks <- 1:52       # Weeks

df_vacc <- expand.grid(
  Week = weeks,
  AgeGroup = age_groups
)

# Extract mean, lower, and upper from the summary array
df_vacc$Median <- as.vector(aperm(summary_cases_vacc["median", , ], c(2, 1)))  # Weeks x Age
df_vacc$Lower <- as.vector(aperm(summary_cases_vacc["lower", , ], c(2, 1)))   # Weeks x Age
df_vacc$Upper <- as.vector(aperm(summary_cases_vacc["upper", , ], c(2, 1)))   # Weeks x Age
df_vacc$Scenario <- "Post-Vaccination"

ggplot(df_vacc, aes(x = Week, y = Median, color = Scenario)) +
  geom_line(size = 1) +
  facet_wrap(~ AgeGroup, ncol = 3)+
  theme_bw()

combined_df <- rbind(df_prevacc, df_vacc)

ggplot(combined_df, aes(x = Week, y = Median, color = Scenario)) +
  geom_line(size = 0.8) +
  #geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Scenario), alpha = 0.2) +
  facet_wrap(~ AgeGroup, ncol = 3)+
  theme_bw()

ggplot()+
  geom_line(data = df_vacc, aes(x = Week, y = Median, color = Scenario))+
  facet_wrap(~ AgeGroup, ncol = 3)+
  geom_line(data = df_prevacc, aes(x = Week, y = Median, color = Scenario))+
  theme_bw()


## validate
iv_pred <- as.data.frame(posterior_vacc$IV_pred[1,1,])
i_pred <- as.data.frame(posterior_vacc$I_pred[1,1,])
infection_all <- cbind(iv_pred, i_pred)
t <- c(1:52)
infection_all$t <- t 
colnames(infection_all) <- c("iv_pred","i_pred","time")

ggplot(infection_all)+
  geom_line(aes(x = time, y = iv_pred, color = "blue"))+
  geom_line(aes(x = time, y = i_pred, color = "red"))



predicted_cases_samples <- posterior_cases_prevacc$pred_cases 
pre_vacc_cases <- posterior_cases_prevacc$age_stratified_cases

predicted_cases_summary <- apply(predicted_cases_samples, 2, function(x) {
  c(median = median(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
})
rownames(predicted_cases_summary) <- c("median", "lower", "upper")


predicted_cases_df <- data.frame(
  Week = 1:ncol(predicted_cases_samples),
  Median = predicted_cases_summary["median", ],
  Lower = predicted_cases_summary["lower", ],
  Upper = predicted_cases_summary["upper", ]
)
# Observed weekly cases
observed_df <- data.frame(
  Week = 1:length(stan_data_prevacc$observed_cases),
  Observed = stan_data_prevacc$observed_cases
)

observed_df$Type <- "Observed"
predicted_cases_df$Type <- "Predicted"

p <- 
ggplot() +
  # Observed points
  geom_point(data = observed_df, aes(x = Week, y = Observed, color = Type), size = 2) +
  
  # Predicted line
  geom_line(data = predicted_cases_df, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = predicted_cases_df, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2) +
  
  # Set colors using Brewer palette
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = comma) + 
  # Labels & Theme
  labs(x = "Week", y = "Weekly reported total symptomatic cases") +
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/pre_vacc_total.jpg", p, width = 7, height = 4, dpi = 1200)


ggsave(filename = "02_Outputs/2_1_Figures/fitted_curve_bra23.jpg", p, width = 7, height = 4, dpi = 1200)


library(bayesplot)

ppc_dens_overlay(stan_data_prevacc$observed_cases, predicted_cases_samples[1:100, ]) +
  ggtitle("Posterior Predictive Check: Observed vs. Predicted Cases")

ppc_intervals(stan_data_prevacc$observed_cases,predicted_cases_samples[1:100, ]) +
  ggtitle("Posterior Predictive Check: Predictive Intervals of Cases")

ppc_stat(stan_data_prevacc$observed_cases, predicted_cases_samples, stat = "mean")
ppc_stat(stan_data_prevacc$observed_cases, predicted_cases_samples, stat = "sd")

print(fit_prevacc, pars = c("base_beta[3]", "base_beta[5]"))  # Check Rhat and ESS
mcmc_trace(as.array(fit_prevacc), pars = "rho")