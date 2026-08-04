load("00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")

observed_cases_bahia <- as.matrix(bra_bh_22[,6])
observed_cases_bahia <- matrix(observed_cases_bahia, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_ce <- as.matrix(bra_ce_22[,6])
observed_cases_ce <- matrix(observed_cases_ce, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_mg <- as.matrix(bra_mg_22[,6])
observed_cases_mg <- matrix(observed_cases_mg, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pn <- as.matrix(bra_pn_22[,6])
observed_cases_pn <- matrix(observed_cases_pn, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pa <- as.matrix(bra_pa_22[,6])
observed_cases_pa <- matrix(observed_cases_pa, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_rg <- as.matrix(bra_rg_22[,6])
observed_cases_rg <- matrix(observed_cases_rg, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pi <- as.matrix(bra_pi_22[,6])
observed_cases_pi <- matrix(observed_cases_pi, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_ag <- as.matrix(bra_ag_22[,6])
observed_cases_ag <- matrix(observed_cases_ag, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_tc <- as.matrix(bra_tc_22[,6])
observed_cases_tc <- matrix(observed_cases_tc, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_se <- as.matrix(bra_se_22[,6])
observed_cases_se <- matrix(observed_cases_se, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_go <- as.matrix(bra_go_22[,6])
observed_cases_go <- matrix(observed_cases_go, nrow = 20, ncol = 52, byrow = FALSE)

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_bra_seir_m2.stan")

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
age_groups <- c(mean(0:1),
                mean(1:4),
                mean(5:9),
                mean(10:11),
                mean(12:17),
                mean(18:19),
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
prior_I0_bahia <- pmax(
  as.numeric(
    observed_cases_bahia[, 1]
  ),
  1
)

prior_sd_I0_bahia <- 100

bahia_foi <- bra_foi_state_summ$avg_foi[
  bra_foi_state_summ$NAME_1 == "Bahia"
]

# Define the first year of chikungunya circulation.
introduction_year_bahia <- 2014

# Define the year corresponding to the start of the fitted epidemic.
fitting_year_bahia <- 2022

# Calculate the maximum possible duration of prior exposure.
years_since_introduction_bahia <-
  fitting_year_bahia - introduction_year_bahia

years_since_introduction_bahia

exposure_years_bahia_finite <- pmin(
  age_groups,
  years_since_introduction_bahia
)

sero_bahia_finite <- 1 - exp(
  - bahia_foi * exposure_years_bahia_finite
)

stan_data_prevacc_bh <- list(
  
  # Number of observed weeks.
  T = 52,
  
  # Number of unobserved initialization weeks.
  B = 2,
  
  # Number of age groups.
  A = 20,
  
  # Reported cases arranged as age group by week.
  observed_cases_by_age =
    observed_cases_bahia,
  
  # Population size by age group.
  N =
    as.numeric(N_bahia$Bahia),
  
  # No aging over the one-year fitting period.
  r =
    rep(0, 20),
  
  # Prior mean for the initial infectious population.
  prior_I0 =
    prior_I0_bahia,
  
  # Prior standard deviation for the initial infectious population.
  prior_sd_I0 =
    prior_sd_I0_bahia,
  
  # Prespecified baseline seroprevalence.
  sero =
    1 - exp(
      - bra_foi_state_summ$avg_foi[
        bra_foi_state_summ$NAME_1 == "Bahia"
      ] * age_groups
    )
)

stan_data_prevacc_ce <- list(
  T = 52,
  A = 20,
  N = as.vector(N_ceara$Ceará), 
  observed_cases_by_age = round(observed_cases_ce),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ce[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)
)

stan_data_prevacc_mg <- list(
  T = 52,
  A = 20,
  N = as.vector(N_mg$`Minas Gerais`), 
  observed_cases_by_age = round(observed_cases_mg),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_mg[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Minas Gerais"] * age_groups)
)

stan_data_prevacc_pn <- list(
  T = 52,
  A = 20,
  N = as.vector(N_pemam$Pernambuco), 
  observed_cases_by_age = round(observed_cases_pn),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pn[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Pernambuco"] * age_groups)
)

stan_data_prevacc_pa <- list(
  T = 52,
  A = 20,
  N = as.vector(N_pa$Paraíba), 
  observed_cases_by_age = round(observed_cases_pa),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pa[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Paraíba"] * age_groups)
)

stan_data_prevacc_rg <- list(
  T = 52,
  A = 20,
  N = as.vector(N_rg$`Rio Grande do Norte`), 
  observed_cases_by_age = round(observed_cases_rg),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_rg[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Rio Grande do Norte"] * age_groups)
)

stan_data_prevacc_pi <- list(
  T = 52,
  A = 20,
  N = as.vector(N_pi$Piauí), 
  observed_cases_by_age = round(observed_cases_pi),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_pi[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Piauí"] * age_groups)
)

stan_data_prevacc_ag <- list(
  T = 52,
  A = 20,
  N = as.vector(N_ag$Alagoas), 
  observed_cases_by_age = round(observed_cases_ag),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_ag[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Alagoas"] * age_groups)
)

stan_data_prevacc_tc <- list(
  T = 52,
  A = 20,
  N = as.vector(N_tc$Tocantins), 
  observed_cases_by_age = round(observed_cases_tc),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_tc[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Tocantins"] * age_groups)
)

stan_data_prevacc_se <- list(
  T = 52,
  A = 20,
  N = as.vector(N_se$Sergipe), 
  observed_cases_by_age = round(observed_cases_se),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_se[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Sergipe"] * age_groups)
)

stan_data_prevacc_go <- list(
  T = 52,
  A = 20,
  N = as.vector(N_go$Goiás), 
  observed_cases_by_age = round(observed_cases_go),
  r = rep(0, 20),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,20),
  indexP = rep(0,52),
  eta = 0,
  prior_I0 = round(observed_cases_go[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Goiás"] * age_groups)
)

### fitting 
fit_prevacc_bh <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_bh,
  iter = 2000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 12    # Reduced tree depth
  )
)

fit_prevacc_ce <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ce,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 1,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 12    # Reduced tree depth
  )
)

fit_prevacc_mg <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mg,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 1,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pn <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pn,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 1,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pa <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pa,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.8,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_rg <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_rg,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_pi <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pi,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_ag <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ag,
  iter = 5000,            # Reduced from 40009
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_tc <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_tc,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_se <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_se,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

fit_prevacc_go <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_go,
  iter = 5000,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 1000,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 15    # Reduced tree depth
  )
)

save(fit_prevacc_bh, file = "00_Data/0_2_Processed/fit_prevacc_bh.RData")
save(fit_prevacc_ce, file = "00_Data/0_2_Processed/fit_prevacc_ce.RData")
save(fit_prevacc_mg, file = "00_Data/0_2_Processed/fit_prevacc_mg.RData")
save(fit_prevacc_pn, file = "00_Data/0_2_Processed/fit_prevacc_pn.RData")
save(fit_prevacc_pa, file = "00_Data/0_2_Processed/fit_prevacc_pa.RData")
save(fit_prevacc_rg, file = "00_Data/0_2_Processed/fit_prevacc_rg.RData")
save(fit_prevacc_pi, file = "00_Data/0_2_Processed/fit_prevacc_pi.RData")
save(fit_prevacc_ag, file = "00_Data/0_2_Processed/fit_prevacc_ag.RData")
save(fit_prevacc_tc, file = "00_Data/0_2_Processed/fit_prevacc_tc.RData")
save(fit_prevacc_se, file = "00_Data/0_2_Processed/fit_prevacc_se.RData")
save(fit_prevacc_go, file = "00_Data/0_2_Processed/fit_prevacc_go.RData")

## post processing 
list_bh <- create_summary_df(fit_prevacc_bh,
                             bra_sum_bh,
                             region = "Bahia")

list_ce <- create_summary_df(fit_prevacc_ce,
                             bra_sum_ce,
                             region = "Ceará")

list_mg <- create_summary_df(fit_prevacc_mg,
                             bra_sum_mg,
                             region = "Minas Gerais")

list_pn <- create_summary_df(fit_prevacc_pn,
                             bra_sum_pn,
                             region = "Pernambuco")

list_pa <- create_summary_df(fit_prevacc_pa,
                             bra_sum_pa,
                             region = "Paraíba")

list_rg <- create_summary_df(fit_prevacc_rg,
                             bra_sum_rg,
                             region = "Rio Grande do Norte")

list_pi <- create_summary_df(fit_prevacc_pi,
                             bra_sum_pi,
                             region = "Piauí")

list_ag <- create_summary_df(fit_prevacc_ag,
                             bra_sum_ag,
                             region = "Alagoas")

list_tc <- create_summary_df(fit_prevacc_tc,
                             bra_sum_tc,
                             region = "Tocantins")

list_se <- create_summary_df(fit_prevacc_se,
                                bra_sum_se,
                                region = "Sergipe")

list_go <- create_summary_df(fit_prevacc_go,
                                bra_sum_go,
                                region = "Goiás")

## bahia
df_bh <- list_bh$df_out

df_bh_summ <- list_bh$df_summ

observed_bh <- list_bh$observed

overall_fit_gg(observed_bh, df_bh_summ)


## ceara
df_ce <- list_ce$df_out

df_ce_summ <- list_ce$df_summ

observed_ce <- list_ce$observed


overall_fit_gg(observed_ce, df_ce_summ)


## minas gerais
df_mg <- list_mg$df_out

df_mg_summ <- list_mg$df_summ

observed_mg <- list_mg$observed

age_strat_gg(df_mg)

overall_fit_gg(observed_mg, df_mg_summ)


## pernambuco
df_pn <- list_pn$df_out

df_pn_summ <- list_pn$df_summ

observed_pn <- list_pn$observed

age_strat_gg(df_pn)

overall_fit_gg(observed_pn, df_pn_summ)


## paraiba
df_pa <- list_pa$df_out

df_pa_summ <- list_pa$df_summ

observed_pa <- list_pa$observed


overall_fit_gg(observed_pa, df_pa_summ)


## rio grande norte
df_rg <- list_rg$df_out

df_rg_summ <- list_rg$df_summ

observed_rg <- list_rg$observed

age_strat_gg(df_rg)

overall_fit_gg(observed_rg, df_rg_summ)

## Piuai
df_pi <- list_pi$df_out

df_pi_summ <- list_pi$df_summ

observed_pi <- list_pi$observed

age_strat_gg(df_pi)

overall_fit_gg(observed_pi, df_pi_summ)

## Alagoas
df_ag <- list_ag$df_out

df_ag_summ <- list_ag$df_summ

observed_ag <- list_ag$observed

age_strat_gg(df_ag)

overall_fit_gg(observed_ag, df_ag_summ)

## tocantins
df_tc <- list_tc$df_out

df_tc_summ <- list_tc$df_summ

observed_tc <- list_tc$observed

age_strat_gg(df_tc)

overall_fit_gg(observed_tc, df_tc_summ)

## Sergipe
df_se_22 <- list_se$df_out

df_se_summ_22 <- list_se$df_summ

observed_se <- list_se$observed


overall_fit_gg(observed_se, df_se_summ_22)

## goias
df_go_22 <- list_go$df_out

df_go_summ_22 <- list_go$df_summ

observed_go <- list_go$observed

age_strat_gg(df_go_22)

overall_fit_gg(observed_go, df_go_summ_22)

## overall graph 
observed_ce$region <- "Ceará" 
observed_bh$region <- "Bahia" 
observed_ag$region <- "Alagoas" 
observed_mg$region <- "Minas Gerais" 
observed_tc$region <- "Tocantins"
observed_pa$region <- "Paraíba" 
observed_pi$region <- "Piauí"
observed_pn$region <- "Pernambuco"
observed_rg$region <- "Rio Grande do Norte"
observed_se$region <- "Sergipe"
observed_go$region <- "Goiás"

observed_all <- bind_rows(observed_ce, observed_bh, observed_ag,
                          observed_mg, observed_tc, observed_pa,
                          observed_pi, observed_pn, observed_rg,
                          observed_se, observed_go)

pred_all <- bind_rows(
  df_ag_summ, df_bh_summ, df_ce_summ,
  df_mg_summ, df_pa_summ, df_pi_summ,
  df_rg_summ, df_pn_summ, df_tc_summ,
  df_go_summ_22, df_se_summ_22
)

model_fit <- 

ggplot()+
  geom_point(data = observed_all, aes(x = Week, y = Observed), size = 0.8)+
  facet_wrap(~region, scales = "free_y")+
  # Predicted line
  geom_line(data = pred_all, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2)+
  theme_pubclean()+
  theme(legend.position = "right")+
  ylab("Predicted and observed reported symptomatic cases")+
  scale_y_continuous(labels = comma)

model_fit_2 <- 
ggarrange(model_fit, 
          ncol = 1,
          labels = c("D"),
          common.legend = TRUE,
          legend = "bottom",
          align = "none")

ggsave(filename = "02_Outputs/2_1_Figures/figs5.jpg", model_fit, width = 12, height = 8, dpi = 1200)


save("observed_ce", "observed_bh", "observed_ag", 
     "observed_mg", "observed_tc", "observed_pa",
     "observed_se", "observed_go",
     "observed_pi", "observed_pn", "observed_rg", "observed_all", "bra_all_sum",
     file = "00_Data/0_2_Processed/observed_2022.RData")


