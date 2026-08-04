
observed_cases_bh_23 <- as.matrix(bra_bh_23[,4])
observed_cases_bh_23 <- matrix(observed_cases_bh_23, nrow = 18, ncol = 41, byrow = FALSE)

observed_cases_ce_23 <- as.matrix(bra_ce_23[,4])
observed_cases_ce_23 <- matrix(observed_cases_ce_23, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_mg_23 <- as.matrix(bra_mg_23[,4])
observed_cases_mg_23 <- matrix(observed_cases_mg_23, nrow = 18, ncol = 36, byrow = FALSE)

observed_cases_ms_23 <- as.matrix(bra_ms_23[,4])
observed_cases_ms_23 <- matrix(observed_cases_ms_23, nrow = 18, ncol = 41, byrow = FALSE)

observed_cases_sp_23 <- as.matrix(bra_sp_23[,4])
observed_cases_sp_23 <- matrix(observed_cases_sp_23, nrow = 18, ncol = 31, byrow = FALSE)

observed_cases_tc_23 <- as.matrix(bra_tc_23[,4])
observed_cases_tc_23 <- matrix(observed_cases_tc_23, nrow = 18, ncol = 30, byrow = FALSE)


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
stan_data_prevacc_bh <- list(
  T = 41,
  A = 18,
  N = as.vector(N_bahia), 
  observed_cases_by_age = round(observed_cases_bh_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,41),
  eta = 0,
  prior_I0 = round(observed_cases_bh_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Bahia"] * age_groups)
)

stan_data_prevacc_mg <- list(
  T = 36,
  A = 18,
  N = as.vector(N_mg), 
  observed_cases_by_age = round(observed_cases_mg_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,36),
  eta = 0,
  prior_I0 = round(observed_cases_mg_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Minas Gerais"] * age_groups)
)

stan_data_prevacc_ms <- list(
  T = 41,
  A = 18,
  N = as.vector(N_ms), 
  observed_cases_by_age = round(observed_cases_ms_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,41),
  eta = 0,
  prior_I0 = round(observed_cases_ms_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Mato Grosso do Sul"] * age_groups)
)

stan_data_prevacc_sp <- list(
  T = 31,
  A = 18,
  N = as.vector(N_sp), 
  observed_cases_by_age = round(observed_cases_sp_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,31),
  eta = 0,
  prior_I0 = round(observed_cases_sp_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "São Paulo"] * age_groups)
)


stan_data_prevacc_ce <- list(
  T = 52,
  A = 18,
  N = as.vector(N_ceara), 
  observed_cases_by_age = round(observed_cases_ce_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ce_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)
)

stan_data_prevacc_tc <- list(
  T = 30,
  A = 18,
  N = as.vector(N_tc), 
  observed_cases_by_age = round(observed_cases_tc_23),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,30),
  eta = 0,
  prior_I0 = round(observed_cases_tc_23[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Tocantins"] * age_groups)
)

### fitting 
fit_prevacc_bh_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_bh,
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

fit_prevacc_ms_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ms,
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

fit_prevacc_sp_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_sp,
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

fit_prevacc_mg_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mg,
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

fit_prevacc_tc_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_tc,
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


fit_prevacc_ce_23 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ce,
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

save(fit_prevacc_bh_23, file = "00_Data/0_2_Processed/fit_prevacc_bh_23.RData")
save(fit_prevacc_ce_23, file = "00_Data/0_2_Processed/fit_prevacc_ce_23.RData")
save(fit_prevacc_mg_23, file = "00_Data/0_2_Processed/fit_prevacc_mg_23.RData")
save(fit_prevacc_tc_23, file = "00_Data/0_2_Processed/fit_prevacc_tc_23.RData")
save(fit_prevacc_ms_23, file = "00_Data/0_2_Processed/fit_prevacc_ms_23.RData")
save(fit_prevacc_sp_23, file = "00_Data/0_2_Processed/fit_prevacc_sp_23.RData")

## post processing 
list_bh_23 <- create_summary_df(fit_prevacc_bh_23,
                                bra_sum_bh_23,
                                region = "Bahia")

list_mg_23 <- create_summary_df(fit_prevacc_mg_23,
                                bra_sum_mg_23,
                                region = "Minas Gerais")

list_ms_23 <- create_summary_df(fit_prevacc_ms_23,
                                bra_sum_ms_23,
                                region = "Mato Grosso do Sul")

list_sp_23 <- create_summary_df(fit_prevacc_sp_23,
                                bra_sum_sp_23,
                                region = "São Paulo")

list_ce_23 <- create_summary_df(fit_prevacc_ce_23,
                                bra_sum_ce_23,
                                region = "Ceará")

list_tc_23 <- create_summary_df(fit_prevacc_tc_23,
                                bra_sum_tc_23,
                                region = "Tocantins")


## bahia
df_bh_23 <- list_bh_23$df_out

df_bh_summ_23 <- list_bh_23$df_summ

observed_bh_23 <- list_bh_23$observed

age_strat_gg(df_bh_23)

i <- overall_fit_gg(observed_bh_23, df_bh_summ_23)

## ceara
df_ce_23 <- list_ce_23$df_out

df_ce_summ_23 <- list_ce_23$df_summ

observed_ce_23 <- list_ce_23$observed

age_strat_gg(df_ce_23)

a <- overall_fit_gg(observed_ce_23, df_ce_summ_23)

combined_ce <- age_struc_fit(fit_prevacc_ce_23,
                             model_age_bins,
                             bra_ce,
                             bra_sum_ce_23)
agestrat_ui_gg(combined_ce)

## minas gerais
df_mg_23 <- list_mg_23$df_out

df_mg_summ_23 <- list_mg_23$df_summ

observed_mg_23 <- list_mg_23$observed

age_strat_gg(df_mg_23)

b <- overall_fit_gg(observed_mg_23, df_mg_summ_23)

combined_mg <- age_struc_fit(fit_prevacc_mg_23,
                             model_age_bins,
                             bra_mg,
                             bra_sum_mg_23)
agestrat_ui_gg(combined_mg)


## Mato Grosso do Sul
df_ms_23 <- list_ms_23$df_out

df_ms_summ_23 <- list_ms_23$df_summ

observed_ms_23 <- list_ms_23$observed

age_strat_gg(df_ms_23)

d <- overall_fit_gg(observed_ms_23, df_ms_summ_23)

combined_ms <- age_struc_fit(fit_prevacc_ms_23,
                             model_age_bins,
                             bra_ms,
                             bra_sum_ms_23)
agestrat_ui_gg(combined_ms)


## São Paulo
df_sp_23 <- list_sp_23$df_out

df_sp_summ_23 <- list_sp_23$df_summ

observed_sp_23 <- list_sp_23$observed

age_strat_gg(df_sp_23)

f <- overall_fit_gg(observed_sp_23, df_sp_summ_23)

combined_sp <- age_struc_fit(fit_prevacc_sp_23,
                             model_age_bins,
                             bra_sp,
                             bra_sum_sp_23)
agestrat_ui_gg(combined_sp)


## Tocantins
df_tc_23 <- list_tc_23$df_out

df_tc_summ_23 <- list_tc_23$df_summ

observed_tc_23 <- list_tc_23$observed

age_strat_gg(df_tc_23)

g <- overall_fit_gg(observed_tc_23, df_tc_summ_23)

combined_tc <- age_struc_fit(fit_prevacc_tc_23,
                             model_age_bins,
                             bra_tc,
                             bra_sum_tc_23)
agestrat_ui_gg(combined_tc)


## overall graph 
observed_ce_23$region <- "Ceará" 
observed_bh_23$region <- "Bahia" 
observed_mg_23$region <- "Minas Gerais" 
observed_ms_23$region <- "Mato Grosso do Sul"
observed_sp_23$region <- "São Paulo" 
observed_tc_23$region <- "Tocantins"

observed_all <- bind_rows(observed_bh_23, observed_ce_23, observed_mg_23,
                          observed_ms_23, observed_tc_23, observed_sp_23)

pred_all <- bind_rows(
  df_bh_summ_23, df_mg_summ_23, df_ms_summ_23,
  df_sp_summ_23, df_ce_summ_23, df_tc_summ_23
)

p <- 
  
  ggplot()+
  geom_point(data = observed_all, aes(x = Week, y = Observed), size = 1)+
  facet_wrap(~region, scales = "free_y")+
  # Predicted line
  geom_line(data = pred_all, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2)+
  theme_light()+
  scale_y_continuous(label = comma)+
  ylab("Observed and predicted symptomatic cases")

ggsave(filename = "02_Outputs/2_1_Figures/fig_2023_fit.jpg", p, width = 10, height = 8, dpi = 1200)


