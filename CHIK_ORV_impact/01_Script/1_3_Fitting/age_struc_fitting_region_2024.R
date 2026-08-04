
observed_cases_bh_24 <- as.matrix(bra_bh_24[,4])
observed_cases_bh_24 <- matrix(observed_cases_bh_24, nrow = 18, ncol = 26, byrow = FALSE)

observed_cases_ce_24 <- as.matrix(bra_ce_24[,4])
observed_cases_ce_24 <- matrix(observed_cases_ce_24, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_mg_24 <- as.matrix(bra_mg_24[,4])
observed_cases_mg_24 <- matrix(observed_cases_mg_24, nrow = 18, ncol = 30, byrow = FALSE)

observed_cases_pn_24 <- as.matrix(bra_pn_24[,4])
observed_cases_pn_24 <- matrix(observed_cases_pn_24, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_go_24 <- as.matrix(bra_go_24[,4])
observed_cases_go_24 <- matrix(observed_cases_go_24, nrow = 18, ncol = 52, byrow = FALSE)

observed_cases_mt_24 <- as.matrix(bra_mt_24[,4])
observed_cases_mt_24 <- matrix(observed_cases_mt_24, nrow = 18, ncol = 40, byrow = FALSE)

observed_cases_ms_24 <- as.matrix(bra_ms_24[,4])
observed_cases_ms_24 <- matrix(observed_cases_ms_24, nrow = 18, ncol = 25, byrow = FALSE)

observed_cases_sp_24 <- as.matrix(bra_sp_24[,4])
observed_cases_sp_24 <- matrix(observed_cases_sp_24, nrow = 18, ncol = 35, byrow = FALSE)


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
  T = 26,
  A = 18,
  N = as.vector(N_bahia), 
  observed_cases_by_age = round(observed_cases_bh_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,26),
  eta = 0,
  prior_I0 = round(observed_cases_bh_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Bahia"] * age_groups)
)

stan_data_prevacc_pn <- list(
  T = 52,
  A = 18,
  N = as.vector(N_pemam), 
  observed_cases_by_age = round(observed_cases_pn_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pn_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Pernambuco"] * age_groups)
)

stan_data_prevacc_mg <- list(
  T = 30,
  A = 18,
  N = as.vector(N_mg), 
  observed_cases_by_age = round(observed_cases_mg_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,30),
  eta = 0,
  prior_I0 = round(observed_cases_mg_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Minas Gerais"] * age_groups)
)

stan_data_prevacc_ms <- list(
  T = 25,
  A = 18,
  N = as.vector(N_ms), 
  observed_cases_by_age = round(observed_cases_ms_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,25),
  eta = 0,
  prior_I0 = round(observed_cases_ms_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Mato Grosso do Sul"] * age_groups)
)

stan_data_prevacc_sp <- list(
  T = 35,
  A = 18,
  N = as.vector(N_sp), 
  observed_cases_by_age = round(observed_cases_sp_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,35),
  eta = 0,
  prior_I0 = round(observed_cases_sp_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "São Paulo"] * age_groups)
)

stan_data_prevacc_mt <- list(
  T = 40,
  A = 18,
  N = as.vector(N_mt), 
  observed_cases_by_age = round(observed_cases_mt_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,40),
  eta = 0,
  prior_I0 = round(observed_cases_mt_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Mato Grosso"] * age_groups)
)

stan_data_prevacc_ce <- list(
  T = 52,
  A = 18,
  N = as.vector(N_ceara), 
  observed_cases_by_age = round(observed_cases_ce_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ce_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)
)

stan_data_prevacc_go <- list(
  T = 52,
  A = 18,
  N = as.vector(N_go), 
  observed_cases_by_age = round(observed_cases_go_24),
  r = rep(0, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_go_24[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Goiás"] * age_groups)
)

### fitting 
fit_prevacc_bh_24 <- sampling(
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

fit_prevacc_ms_24 <- sampling(
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

fit_prevacc_mg_24 <- sampling(
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

fit_prevacc_pn_24 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pn,
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

fit_prevacc_mt_24 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mt,
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

fit_prevacc_sp_24 <- sampling(
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

fit_prevacc_ce_24 <- sampling(
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

fit_prevacc_go_24 <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_go,
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

save(fit_prevacc_bh_24, file = "00_Data/0_2_Processed/fit_prevacc_bh_24.RData")
save(fit_prevacc_ce_24, file = "00_Data/0_2_Processed/fit_prevacc_ce_24.RData")
save(fit_prevacc_mg_24, file = "00_Data/0_2_Processed/fit_prevacc_mg_24.RData")
save(fit_prevacc_pn_24, file = "00_Data/0_2_Processed/fit_prevacc_pn_24.RData")
save(fit_prevacc_ms_24, file = "00_Data/0_2_Processed/fit_prevacc_ms_24.RData")
save(fit_prevacc_mt_24, file = "00_Data/0_2_Processed/fit_prevacc_mt_24.RData")
save(fit_prevacc_sp_24, file = "00_Data/0_2_Processed/fit_prevacc_sp_24.RData")
save(fit_prevacc_go_24, file = "00_Data/0_2_Processed/fit_prevacc_go_24.RData")

## post processing 
list_bh_24 <- create_summary_df(fit_prevacc_bh_24,
                                bra_sum_bh_24,
                                region = "Bahia")

list_mt_24 <- create_summary_df(fit_prevacc_mt_24,
                                bra_sum_mt_24,
                                region = "Mato Grosso")

list_mg_24 <- create_summary_df(fit_prevacc_mg_24,
                                bra_sum_mg_24,
                                region = "Minas Gerais")

list_pn_24 <- create_summary_df(fit_prevacc_pn_24,
                                bra_sum_pn_24,
                                region = "Pernambuco")

list_ms_24 <- create_summary_df(fit_prevacc_ms_24,
                                bra_sum_ms_24,
                                region = "Mato Grosso do Sul")

list_sp_24 <- create_summary_df(fit_prevacc_sp_24,
                                bra_sum_sp_24,
                                region = "São Paulo")

list_ce_24 <- create_summary_df(fit_prevacc_ce_24,
                                bra_sum_ce_24,
                                region = "Ceará")

list_go_24 <- create_summary_df(fit_prevacc_go_24,
                                bra_sum_go_24,
                                region = "Goiás")


## bahia
df_bh_24 <- list_bh_24$df_out

df_bh_summ_24 <- list_bh_24$df_summ

observed_bh_24 <- list_bh_24$observed

age_strat_gg(df_bh_24)

i <- overall_fit_gg(observed_bh_24, df_bh_summ_24)

## ceara
df_ce_24 <- list_ce_24$df_out

df_ce_summ_24 <- list_ce_24$df_summ

observed_ce_24 <- list_ce_24$observed

age_strat_gg(df_ce_24)

a <- overall_fit_gg(observed_ce_24, df_ce_summ_24)

combined_ce <- age_struc_fit(fit_prevacc_ce_24,
                             model_age_bins,
                             bra_ce_24,
                             bra_sum_ce_24)
agestrat_ui_gg(combined_ce)

## minas gerais
df_mg_24 <- list_mg_24$df_out

df_mg_summ_24 <- list_mg_24$df_summ

observed_mg_24 <- list_mg_24$observed

age_strat_gg(df_mg_24)

b <- overall_fit_gg(observed_mg_24, df_mg_summ_24)

combined_mg <- age_struc_fit(fit_prevacc_mg_24,
                             model_age_bins,
                             bra_mg_24,
                             bra_sum_mg_24)
agestrat_ui_gg(combined_mg)

## pernambuco
df_pn_24 <- list_pn_24$df_out

df_pn_summ_24 <- list_pn_24$df_summ

observed_pn_24 <- list_pn_24$observed

age_strat_gg(df_pn)

c <- overall_fit_gg(observed_pn_24, df_pn_summ_24)

combined_pn <- age_struc_fit(fit_prevacc_pn_24,
                             model_age_bins,
                             bra_pn_24,
                             bra_sum_pn_24)
agestrat_ui_gg(combined_pn)


## Mato Grosso do Sul
df_ms_24 <- list_ms_24$df_out

df_ms_summ_24 <- list_ms_24$df_summ

observed_ms_24 <- list_ms_24$observed

age_strat_gg(df_ms_24)

d <- overall_fit_gg(observed_ms_24, df_ms_summ_24)

combined_ms <- age_struc_fit(fit_prevacc_ms_24,
                             model_age_bins,
                             bra_ms_24,
                             bra_sum_ms_24)
agestrat_ui_gg(combined_ms)


## Mato Grosso
df_mt_24 <- list_mt_24$df_out

df_mt_summ_24 <- list_mt_24$df_summ

observed_mt_24 <- list_mt_24$observed

age_strat_gg(df_mt_24)

e <- overall_fit_gg(observed_mt_24, df_mt_summ_24)

combined_mt <- age_struc_fit(fit_prevacc_mt_24,
                             model_age_bins,
                             bra_mt_24,
                             bra_sum_mt_24)
agestrat_ui_gg(combined_mt)

## São Paulo
df_sp_24 <- list_sp_24$df_out

df_sp_summ_24 <- list_sp_24$df_summ

observed_sp_24 <- list_sp_24$observed

age_strat_gg(df_sp_24)

f <- overall_fit_gg(observed_sp_24, df_sp_summ_24)

combined_sp <- age_struc_fit(fit_prevacc_sp_24,
                             model_age_bins,
                             bra_sp_24,
                             bra_sum_sp_24)
agestrat_ui_gg(combined_sp)

## Goiás
df_go_24 <- list_go_24$df_out

df_go_summ_24 <- list_go_24$df_summ

observed_go_24 <- list_go_24$observed

age_strat_gg(df_go_24)

g <- overall_fit_gg(observed_go_24, df_go_summ_24)

combined_go <- age_struc_fit(fit_prevacc_go_24,
                             model_age_bins,
                             bra_go_24,
                             bra_sum_go_24)
agestrat_ui_gg(combined_go)


## overall graph 
observed_ce_24$region <- "Ceará" 
observed_bh_24$region <- "Bahia" 
observed_mt_24$region <- "Mato Grosso" 
observed_mg_24$region <- "Minas Gerais" 
observed_ms_24$region <- "Mato Grosso do Sul"
observed_sp_24$region <- "São Paulo" 
observed_pn_24$region <- "Pernambuco"
observed_go_24$region <- "Goiás"

observed_all <- bind_rows(observed_bh_24, observed_pn_24, observed_mg_24,
                          observed_ms_24, observed_mt_24, observed_sp_24,
                          observed_ce_24, observed_go_24)

pred_all <- bind_rows(
  df_bh_summ_24, df_pn_summ_24, df_mg_summ_24,
  df_ms_summ_24, df_mt_summ_24, df_sp_summ_24,
  df_ce_summ_24, df_go_summ_24
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

ggsave(filename = "02_Outputs/2_1_Figures/fig_2024_fit.jpg", p, width = 10, height = 8, dpi = 1200)

save("observed_ce_24", "observed_bh_24", "observed_mt_24", 
     "observed_mg_24", "observed_sp_24", "observed_ms_24",
     "observed_go_24", "observed_pn_24",
     file = "00_Data/0_2_Processed/observed_2024.RData")


