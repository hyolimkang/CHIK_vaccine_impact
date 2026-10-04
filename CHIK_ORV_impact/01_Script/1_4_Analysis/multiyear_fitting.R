# multi year data prep
bra_ce_multi <- bind_rows(
  bra_ce_22, bra_ce_23, bra_ce_24
)

bra_ce_multi$epi_week <- rep(1:156, each = 18)

bra_sum_ce_multi <- bind_rows(
  bra_sum_ce, bra_sum_ce_23, bra_sum_ce_24
)

bra_sum_ce_multi$epi_week <- c(1:156)

# multi year data prep
bra_bh_multi <- bind_rows(
  bra_bh_22, bra_bh_23, bra_bh_24
)

bra_bh_multi$epi_week <- rep(1:119, each = 18)

bra_sum_bh_multi <- bind_rows(
  bra_sum_bh, bra_sum_bh_23, bra_sum_bh_24
)

bra_sum_bh_multi$epi_week <- c(1:119)

# multi year data prep
bra_mg_multi <- bind_rows(
  bra_mg_22, bra_mg_23, bra_mg_24
)

bra_mg_multi$epi_week <- rep(1:96, each = 18)

bra_sum_mg_multi <- bind_rows(
  bra_sum_mg, bra_sum_mg_23, bra_sum_mg_24
)

bra_sum_mg_multi$epi_week <- c(1:96)

# multi year data prep
bra_sp_multi <- bind_rows(
  bra_sp_23, bra_sp_24
)

bra_sp_multi$epi_week <- rep(53:118, each = 18)

bra_sum_sp_multi <- bind_rows(
  bra_sum_sp_23, bra_sum_sp_24
)

bra_sum_sp_multi$epi_week <- c(53:118)

# multi year data prep
bra_tc_multi <- bind_rows(
  bra_tc_22, bra_tc_23
)

bra_tc_multi$epi_week <- rep(1:82, each = 18)

bra_sum_tc_multi <- bind_rows(
  bra_sum_tc, bra_sum_tc_23
)

bra_sum_tc_multi$epi_week <- c(1:82)

# multi year data prep
bra_ms_multi <- bind_rows(
  bra_ms_23, bra_ms_24
)

bra_ms_multi$epi_week <- rep(53:118, each = 18)

bra_sum_ms_multi <- bind_rows(
  bra_sum_ms_23, bra_sum_ms_24
)

bra_sum_ms_multi$epi_week <- c(53:118)

bra_sum_mt_24$epi_week <- c(105:144)

# multi year data prep
bra_pn_multi <- bind_rows(
  bra_pn_22, bra_pn_22, bra_pn_24
)

bra_pn_multi$epi_week <- rep(1:156, each = 18)

bra_sum_pn_multi <- bind_rows(
  bra_sum_pn, bra_sum_pn, bra_sum_pn_24
)

bra_sum_pn_multi$epi_week <- c(1:156)
bra_sum_pn_multi$tot_cases[53:108] <- 0

combined_multi_year <- bind_rows(bra_sum_bh_multi, bra_sum_mg_multi, bra_sum_ce_multi,
                                 bra_sum_sp_multi, bra_sum_tc_multi, bra_sum_ms_multi,
                                 bra_sum_pn_multi, bra_sum_mt_24, bra_sum_pa, bra_sum_pi,
                                 bra_sum_ag, bra_sum_rg)

p <- 
ggplot(combined_multi_year)+
  geom_line(aes(x= epi_week, y= tot_cases))+
  facet_wrap(~region, scales = "free_y")+
  theme_light()

ggsave(filename = "02_Outputs/2_1_Figures/fig_3yr_outbreaks.jpg", p, width = 10, height = 8, dpi = 1200)

# pre-process obs data 
observed_cases_ce_multi <- as.matrix(bra_ce_multi[,4])
observed_cases_ce_multi <- matrix(observed_cases_ce_multi, nrow = 18, ncol = 156, byrow = FALSE)

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

# indexP
indexP_ce <- indexP_matrix_long %>% filter(state == "Ceará")
indexP_ce <- rep(indexP_ce[[3]], times = 3)

## prevacc data 
stan_data_prevacc_ce <- list(
  T = 156,
  A = 18,
  N = as.vector(N_ceara), 
  observed_cases_by_age = round(observed_cases_ce_multi),
  r = rep(0.25/1000, 18),
  delay = 53,
  VE_block = 0,
  vaccination_rate = 0,
  vaccine_target_age = rep(0,18),
  indexP = indexP_ce,
  eta = 0,
  prior_I0 = round(observed_cases_ce[,1]),
  prior_sd_I0 = 100,
  sero = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)
)


### fitting 
stan_model_multi <- stan_model("01_Script/1_2_SIR_models/age_struc_vacc_bra_v4.1.stan")

fit_prevacc_ce_multi <- sampling(
  object = stan_model_multi,
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


# post-processing
list_ce_multi <- create_summary_df(fit_prevacc_ce_multi,
                                   bra_sum_ce_multi,
                                   region = "Ceará")

## ceara
df_ce_multi <- list_ce_multi$df_out

df_ce_summ_multi <- list_ce_multi$df_summ

observed_ce_multi <- list_ce_multi$observed

age_strat_gg(df_ce_multi)

a <- overall_fit_gg(observed_ce_multi, df_ce_summ_multi)

combined_ce <- age_struc_fit(fit_prevacc_ce_multi,
                             model_age_bins,
                             bra_ce_multi,
                             bra_sum_ce_multi)
agestrat_ui_gg(combined_ce)

