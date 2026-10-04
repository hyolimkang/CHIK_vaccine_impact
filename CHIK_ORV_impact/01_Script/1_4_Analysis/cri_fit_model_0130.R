# monthly data fitting: all LAC countries / descriptive analysis first --------------------------------------

# all countries
lac_report_plisa <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                               sheet = "Sheet3")

total_cases <- lac_report_plisa %>%
  group_by(country) %>%
  summarize(total_cases = sum(cases, na.rm = TRUE),
            attack_rate = total_cases / tot_pop)

p <- ggplot(lac_report_plisa)+
  geom_line(aes(x = date, y = case_per_100k, color = country))+
  theme_minimal()+
  facet_wrap(~country, scales = "free_y")+
  theme(legend.position = "none") 

ggsave(filename = "02_Outputs/2_1_Figures/outbreaks_lac.jpg", p, width = 15, height = 7, dpi = 1200)

classified_countries <- lac_report_plisa %>%
  group_by(country) %>%
  summarize(max_cases = max(cases, na.rm = TRUE)) %>%
  mutate(group = case_when(
    max_cases > 10000 ~ "High",
    max_cases > 1000 & max_cases <= 10000 ~ "Medium",
    TRUE ~ "Low"
  )) %>%
  arrange(desc(max_cases))

classified_data <- lac_report_plisa %>%
  left_join(classified_countries, by = "country")

classified_data <- classified_data %>%
  mutate(group = factor(group, levels = c("High", "Medium", "Low")))

# remove imperfect data
imperfect_countries <- c("Colombia", "Jamaica", "Belize", 
                         "Bermuda", "Bonaire Saint Eustatius and Saba",
                         "British Virgin Islands", "Chile", "Cuba",
                         "United States of America", "Saint Vincent and the Grenadines") 
classified_data <- classified_data %>%
  filter(!country %in% imperfect_countries)


p <- ggplot(classified_data, aes(x = date, y = cases, color = group, group = country)) +
  geom_line(size = 1) +
  theme_minimal() +
  facet_wrap(~ group + country, scales = "free_y")  +  # Adjusts y-axis scale per country
  theme(legend.position = "none")

high_countries <- classified_data %>% filter(group == "High")

# sample 
costa_rica <- high_countries %>% filter(country == "Costa Rica")

ggplot(costa_rica[1:12,], aes(x = date, y = cases, color = group, group = country)) +
  geom_line(size = 1) +
  theme_minimal() +
  facet_wrap(~ group + country, scales = "free_y")  +  # Adjusts y-axis scale per country
  theme(legend.position = "none")

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_scenario_cri.stan")


# 
N_cri <- as.matrix(cri_pop_transposed[,1])

stan_data_prevacc <- list(
  T = length(costa_rica[1:12,]$cases),
  A = 18,
  N = as.vector(N_cri), 
  observed_cases = costa_rica[1:12,]$cases,
  r = rep(0,18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,12),
  eta = 0
)

fit_prevacc <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc,
  iter = 1000,            # Reduced from 4000
  chains = 2,             # Reduced from 4
  warmup = 100,           # Specify warmup period
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
month <- 1:12       # Weeks

df_prevacc <- expand.grid(
  Month = month,
  AgeGroup = age_groups
)

# Extract mean, lower, and upper from the summary array
df_prevacc$Median <- as.vector(aperm(summary_cases_prevacc["median", , ], c(2, 1)))  # Weeks x Age
df_prevacc$Lower <- as.vector(aperm(summary_cases_prevacc["lower", , ], c(2, 1)))   # Weeks x Age
df_prevacc$Upper <- as.vector(aperm(summary_cases_prevacc["upper", , ], c(2, 1)))   # Weeks x Age
df_prevacc$Scenario <- "Pre-Vaccination"

ggplot(df_prevacc, aes(x = Month, y = Median, color = Scenario)) +
  geom_line(size = 1) +
  facet_wrap(~ AgeGroup, ncol = 3)+
  theme_bw()

observed_df <- data.frame(
  Month = 1:length(stan_data_prevacc$observed_cases),
  Observed = stan_data_prevacc$observed_cases
)

posterior_prevacc <- extract(fit_prevacc)
predicted_cases_samples <- posterior_prevacc$pred_cases 
pre_vacc_cases <- posterior_prevacc$age_stratified_cases

predicted_cases_summary <- apply(predicted_cases_samples, 2, function(x) {
  c(median = median(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
})
rownames(predicted_cases_summary) <- c("median", "lower", "upper")
# Convert to data frame
predicted_cases_df <- data.frame(
  Month = 1:ncol(predicted_cases_samples),
  Median = predicted_cases_summary["median", ],
  Lower = predicted_cases_summary["lower", ],
  Upper = predicted_cases_summary["upper", ]
)


p <- ggplot() +
  geom_point(data = observed_df, aes(x = Month, y = Observed), color = "red", size = 2) +
  geom_line(data = predicted_cases_df, aes(x = Month, y = Median), color = "blue", size = 1) +
  geom_ribbon(data = predicted_cases_df, aes(x = Month, ymin = Lower, ymax = Upper), 
              fill = "blue", alpha = 0.2) +
  #labs(title = "aggregated model predictions vs. observed cases",
  #     x = "Week", y = "Number of infectious individuals") +
  labs(title = "Epi-curve fitted to monthly CRI data (2014-15)",
       x = "Month", y = "Number of infectious individuals") +
  scale_x_continuous(breaks = unique(predicted_cases_df$Month))+
  theme_minimal() 

ggsave(filename = "02_Outputs/2_1_Figures/cri_fit.jpg", p, width = 9, height = 6, dpi = 1200)

# post check
ppc_dens_overlay(stan_data_prevacc$observed_cases, predicted_cases_samples[1:100, ]) +
  ggtitle("Posterior Predictive Check: Observed vs. Predicted Cases")

ppc_intervals(stan_data_prevacc$observed_cases,predicted_cases_samples[1:100, ]) +
  ggtitle("Posterior Predictive Check: Predictive Intervals of Cases")

beta_df <- data.frame(Month = 1:12, Beta = colMeans(posterior_prevacc$beta))

ggplot(beta_df, aes(x = Month, y = Beta)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = apply(posterior_prevacc$beta, 2, quantile, probs = 0.025),
                  ymax = apply(posterior_prevacc$beta, 2, quantile, probs = 0.975)),
              fill = "blue", alpha = 0.2) +
  labs(title = "Estimated Transmission Rate (Beta)", x = "Week", y = "Beta") +
  theme_minimal()

traceplot(fit_prevacc, pars = c("beta[1]", "beta[10]", "beta[12]", "rho"))

posterior_array <- as.array(fit_prevacc) 
posterior_matrix <- posterior_array[, , c("beta[1]", "beta[2]", "beta[3]",
                                          "beta[4]", "beta[5]", "beta[6]",
                                          "beta[7]", "beta[8]", "beta[9]",
                                          "beta[10]", "beta[11]", "beta[12]")]

mcmc_areas(
  posterior_matrix,
  pars = c("beta[1]", "beta[2]", "beta[3]")
)
mcmc_areas(
  posterior_matrix,
  pars = c("beta[4]", "beta[5]", "beta[6]")
)
mcmc_areas(
  posterior_matrix,
  pars = c("beta[7]", "beta[8]", "beta[9]")
)
mcmc_areas(
  posterior_matrix,
  pars = c("beta[10]", "beta[11]", "beta[12]")
)