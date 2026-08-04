load("00_Data/0_2_Processed/bra_cases_2022_cleaned_otherstates.RData")
load("00_Data/0_2_Processed/bra_pop_2022_cleaned_otherstates.RData")

observed_cases_go_22 <- as.matrix(bra_go_22[,6])
observed_cases_go_22 <- matrix(observed_cases_go_22, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_se_22 <- as.matrix(bra_se_22[,6])
observed_cases_se_22 <- matrix(observed_cases_se_22, nrow = 18, ncol = 52, byrow = FALSE)

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_bra_v4.1.stan")

# state specific foi
bra_foi <- allfoi %>% filter(country == "Brazil")
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))
bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)

bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = mean(foi_mid, na.rm = TRUE))%>%
  filter(!is.na(NAME_1))

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

## prevacc data 
stan_data_prevacc_se_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_se$Sergipe), 
  observed_cases_by_age = round(observed_cases_se_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_se_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Sergipe"] * age_groups)
)

stan_data_prevacc_go_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_go$Goiás), 
  observed_cases_by_age = round(observed_cases_go_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_go_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Goiás"] * age_groups)
)

### fitting 
fit_prevacc_se_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_se_22,
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

fit_prevacc_go_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_go_22,
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


save(fit_prevacc_se_22, file = "00_Data/0_2_Processed/fit_prevacc_se_22.RData")
save(fit_prevacc_go_22, file = "00_Data/0_2_Processed/fit_prevacc_go_22.RData")

## post processing 
list_se_22 <- create_summary_df(fit_prevacc_se_22,
                                bra_sum_se,
                                region = "Sergipe")


list_go_22 <- create_summary_df(fit_prevacc_go_22,
                                bra_sum_go,
                                region = "Goiás")


## Sergipe
df_se_22 <- list_se_22$df_out

df_se_summ_22 <- list_se_22$df_summ

observed_se_22 <- list_se_22$observed

age_strat_gg(df_se_22)

i <- overall_fit_gg(observed_se_22, df_se_summ_22)


## goias
df_go_22 <- list_go_22$df_out

df_go_summ_22 <- list_go_22$df_summ

observed_go_22 <- list_go_22$observed

age_strat_gg(df_go_22)

f <- overall_fit_gg(observed_go_22, df_go_summ_22)

combined_go <- age_struc_fit(fit_prevacc_go_22,
                             model_age_bins,
                             bra_go_22,
                             bra_sum_go_22)
agestrat_ui_gg(combined_go)




## overall graph 
observed_se_22$region <- "Sergipe" 
observed_go_22$region <- "Goiás"

observed_all_22_med <- bind_rows(observed_se_22, 
                                 observed_go_22)

pred_all <- bind_rows(
   df_se_summ_22,
   df_go_summ_22
)

p <- 
  
  ggplot()+
  geom_point(data = observed_all_22_med, aes(x = Week, y = Observed), size = 1)+
  facet_wrap(~region, scales = "free_y")+
  # Predicted line
  geom_line(data = pred_all, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2)+
  theme_light()+
  scale_y_continuous(label = comma)+
  ylab("Observed and predicted symptomatic cases")

ggsave(filename = "02_Outputs/2_1_Figures/fig_2022_med_fit.jpg", p, width = 10, height = 6, dpi = 1200)

save("observed_sg_22", "observed_sp_22",
     "observed_go_22", "observed_ms_22", "observed_mh_22", "observed_all_22_med", "bra_all_sum_22_med",
     "N_go", "N_sp", "N_sg", "N_ms", "N_mh",
     file = "00_Data/0_2_Processed/observed_2022_med.RData")


