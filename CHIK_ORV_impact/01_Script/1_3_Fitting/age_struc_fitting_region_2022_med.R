
observed_cases_go_22 <- as.matrix(bra_go_22[,4])
observed_cases_go_22 <- matrix(observed_cases_go_22, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_sg_22 <- as.matrix(bra_sg_22[,4])
observed_cases_sg_22 <- matrix(observed_cases_sg_22, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_sp_22 <- as.matrix(bra_sp_22[,4])
observed_cases_sp_22 <- matrix(observed_cases_sp_22, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_ms_22 <- as.matrix(bra_ms_22[,4])
observed_cases_ms_22 <- matrix(observed_cases_ms_22, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_mh_22 <- as.matrix(bra_mh_22[,4])
observed_cases_mh_22 <- matrix(observed_cases_mh_22, nrow = 18, ncol = 52, byrow = FALSE)

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
stan_data_prevacc_sg_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_sg), 
  observed_cases_by_age = round(observed_cases_sg_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_sg_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Sergipe"] * age_groups)
)

stan_data_prevacc_go_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_go), 
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

stan_data_prevacc_sp_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_sp), 
  observed_cases_by_age = round(observed_cases_sp_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_sp_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "São Paulo"] * age_groups)
)

stan_data_prevacc_ms_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_ms), 
  observed_cases_by_age = round(observed_cases_ms_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ms_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Mato Grosso do Sul"] * age_groups)
)

stan_data_prevacc_mh_22 <- list(
  T = 52,
  A = 18,
  N = as.vector(N_mh), 
  observed_cases_by_age = round(observed_cases_mh_22),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_mh_22[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Maranhão"] * age_groups)
)

### fitting 
fit_prevacc_sg_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_sg_22,
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

fit_prevacc_sp_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_sp_22,
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

fit_prevacc_ms_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ms_22,
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

fit_prevacc_mh_22 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mh_22,
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

save(fit_prevacc_sg_22, file = "00_Data/0_2_Processed/fit_prevacc_sg_22.RData")
save(fit_prevacc_sp_22, file = "00_Data/0_2_Processed/fit_prevacc_sp_22.RData")
save(fit_prevacc_go_22, file = "00_Data/0_2_Processed/fit_prevacc_go_22.RData")
save(fit_prevacc_ms_22, file = "00_Data/0_2_Processed/fit_prevacc_ms_22.RData")
save(fit_prevacc_mh_22, file = "00_Data/0_2_Processed/fit_prevacc_mh_22.RData")

## post processing 
list_sg_22 <- create_summary_df(fit_prevacc_se_22,
                                bra_sum_se,
                                region = "Sergipe")

list_sp_22 <- create_summary_df(fit_prevacc_sp_22,
                                bra_sum_sp_22,
                                region = "São Paulo")

list_go_22 <- create_summary_df(fit_prevacc_go_22,
                                bra_sum_go,
                                region = "Goiás")

list_ms_22 <- create_summary_df(fit_prevacc_ms_22,
                                bra_sum_ms_22,
                                region = "Mato Grosso do Sul")

list_mh_22 <- create_summary_df(fit_prevacc_mh_22,
                                bra_sum_mh_22,
                                region = "Maranhão")

## Sergipe
df_sg_22 <- list_sg_22$df_out

df_sg_summ_22 <- list_sg_22$df_summ

observed_sg_22 <- list_sg_22$observed

age_strat_gg(df_sg_22)

i <- overall_fit_gg(observed_sg_22, df_sg_summ_22)

## São Paulo
df_sp_22 <- list_sp_22$df_out

df_sp_summ_22 <- list_sp_22$df_summ

observed_sp_22 <- list_sp_22$observed

age_strat_gg(df_sp_22)

f <- overall_fit_gg(observed_sp_22, df_sp_summ_22)

combined_sp <- age_struc_fit(fit_prevacc_sp_22,
                             model_age_bins,
                             bra_sp_22,
                             bra_sum_sp_22)
agestrat_ui_gg(combined_sp)

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



## Mato Grosso do Sul
df_ms_22 <- list_ms_22$df_out

df_ms_summ_22 <- list_ms_22$df_summ

observed_ms_22 <- list_ms_22$observed

age_strat_gg(df_ms_22)

g <- overall_fit_gg(observed_ms_22, df_ms_summ_22)

combined_ms <- age_struc_fit(fit_prevacc_ms_22,
                             model_age_bins,
                             bra_ms_22,
                             bra_sum_ms_22)
agestrat_ui_gg(combined_ms)


## Maranhão
df_mh_22 <- list_mh_22$df_out

df_mh_summ_22 <- list_mh_22$df_summ

observed_mh_22 <- list_mh_22$observed

age_strat_gg(df_mh_22)

g <- overall_fit_gg(observed_mh_22, df_mh_summ_22)

combined_mh <- age_struc_fit(fit_prevacc_mh_22,
                             model_age_bins,
                             bra_mh_22,
                             bra_sum_mh_22)
agestrat_ui_gg(combined_mh)



## overall graph 
observed_sg_22$region <- "Sergipe" 
observed_sp_22$region <- "São Paulo" 
observed_go_22$region <- "Goiás"
observed_ms_22$region <- "Mato Grosso do Sul"
observed_mh_22$region <- "Maranhão"

observed_all_22_med <- bind_rows(observed_sg_22, 
                                   observed_sp_22,
                                   observed_go_22,
                                   observed_ms_22,
                                   observed_mh_22)

pred_all <- bind_rows(
   df_sg_summ_22,
   df_sp_summ_22,
   df_go_summ_22,
   df_ms_summ_22,
   df_mh_summ_22
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

save("observed_se_22",
     "observed_go_22",
     "N_go", "N_se",
     file = "00_Data/0_2_Processed/observed_2022_med.RData")


