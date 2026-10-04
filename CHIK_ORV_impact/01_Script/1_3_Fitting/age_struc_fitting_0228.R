
observed_cases_age_mat_21 <- as.matrix(bra_expanded_21[,4])
observed_cases_age_mat_21 <- matrix(observed_cases_age_mat_21, nrow = 18, ncol = 51, byrow = FALSE)

observed_cases_age_mat_22 <- as.matrix(bra_expanded_22[,4])
observed_cases_age_mat_22 <- matrix(observed_cases_age_mat_22, nrow = 18, ncol = 44, byrow = FALSE)

N_tot <- as.matrix(brazil_pop_transposed$tot_pop)
N_2022 <- as.matrix(brazil_pop_transposed[,13]) # minas gerais

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_bra_v4.1.stan")

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

stan_data_prevacc <- list(
  T = 44,
  A = 18,
  N = as.vector(N_tot), 
  observed_cases_by_age = round(observed_cases_age_mat_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0, 44),
  eta = 0,
  prior_I0 = round(observed_cases_age_mat_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(-lambda * age_groups)
)

#prior_rho = rep(brazil_symp$detect_prop, each = 2),
#prior_sd_rho = 0.3

fit_prevacc_21 <- sampling(
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

fit_prevacc_22 <- sampling(
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

save(fit_prevacc_21, file = "00_Data/0_2_Processed/fit_prevacc_bra21_agestrat.RData")
save(fit_prevacc_22, file = "00_Data/0_2_Processed/fit_prevacc_bra22_agestrat.RData")

posterior_cases_prevacc <- extract(fit_prevacc_22)
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
weeks <- 1:44       # Weeks

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
  facet_wrap(~ AgeGroup, ncol = 3)+
  theme_minimal()

df_prevacc <- df_prevacc %>% group_by(Week) %>% summarise(
  Median = sum(Median),
  Lower =  sum(Lower),
  Upper = sum(Upper)
)
# Observed weekly cases
observed_df <- data.frame(
  Week = 1:44,
  Observed = bra_all_sum_22$tot_cases
)

observed_df$Type <- "Observed"
df_prevacc$Type <- "Predicted"

p <- 
ggplot() +
  # Observed points
  geom_point(data = observed_df, aes(x = Week, y = Observed, color = Type), size = 2) +
  
  # Predicted line
  geom_line(data = df_prevacc, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = df_prevacc, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2) +
  
  # Set colors using Brewer palette
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = comma) + 
  # Labels & Theme
  labs(x = "Week", y = "Total cases") +
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/fig_prevacc_fit.jpg", plot = p,
       width = 5, height = 4)


# median
pred_df_wide <- as.data.frame(summary_cases_prevacc[1, , ])
colnames(pred_df_wide) <- paste0("Week", 1:44)
pred_df_wide$AgeGroup <- model_age_bins
pred_df_long <- pred_df_wide %>%
  pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Predicted") %>%
  mutate(Week = as.integer(gsub("Week", "", Week)))

# Lower 95% interval
lower_df_wide <- as.data.frame(summary_cases_prevacc[2, , ])
colnames(lower_df_wide) <- paste0("Week", 1:44)
lower_df_wide$AgeGroup <- model_age_bins
lower_df_long <- lower_df_wide %>%
  pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Lower") %>%
  mutate(Week = as.integer(gsub("Week", "", Week)))

# Upper 95% interval
upper_df_wide <- as.data.frame(summary_cases_prevacc[3, , ])
colnames(upper_df_wide) <- paste0("Week", 1:44)
upper_df_wide$AgeGroup <- model_age_bins
upper_df_long <- upper_df_wide %>%
  pivot_longer(cols = starts_with("Week"), names_to = "Week", values_to = "Upper") %>%
  mutate(Week = as.integer(gsub("Week", "", Week)))

# Merge predicted, lower, and upper into one data frame:
predicted_df <- pred_df_long %>%
  left_join(lower_df_long, by = c("AgeGroup", "Week")) %>%
  left_join(upper_df_long, by = c("AgeGroup", "Week"))

observed_df <- bra_expanded[,c(1,4,5)]
colnames(observed_df) <- c("Week", "Observed", "AgeGroup")

combined_df <- predicted_df %>%
  left_join(observed_df, by = c("AgeGroup", "Week"))

# Now plot with facets by AgeGroup.
ggplot(combined_df, aes(x = Week)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = Predicted), color = "blue", size = 1) +
  geom_point(aes(y = Observed), color = "red", size = 1) +
  facet_wrap(~ AgeGroup) +
  labs(title = "Observed vs Predicted Cases by Age Group",
       x = "Week",
       y = "Cases") +
  theme_minimal()

# brazil age specific symp
brazil_symp <- all_countries %>% filter(country == "Brazil" & type == "symptomatic")

bra_summ <- bra_expanded %>% group_by(age_bin_label) %>% summarise(
  cases = sum(cases)
)

model_age_bins <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                    "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                    "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                    "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                    "80-84 years", "85-89 years")
bra_indexed <- bra_summ %>%
  mutate(age_index = match(age_bin_label, model_age_bins))

bra_indexed <- bra_indexed %>%
  mutate(pairGroup = ceiling(age_index / 2))

# 3) Now group by pairGroup and Week (plus any other grouping vars you need, e.g. region)
#    and sum the Cases:
bra_paired <- bra_indexed %>%
  group_by(pairGroup) %>%
  summarise(cases = sum(cases), .groups = "drop")

pair_labels <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89")
bra_paired <- bra_paired %>%
  mutate(pairGroupLabel = pair_labels[pairGroup])


brazil_symp$obs_symp <- bra_paired$cases

brazil_symp <- brazil_symp %>% mutate(
  detect_prop = obs_symp / tot_med
)



