rho_estim(posterior_ce)
rho_estim(posterior_bh)
rho_estim(posterior_pa)
rho_estim(posterior_pn)
rho_estim(posterior_rg)
rho_estim(posterior_pi)
rho_estim(posterior_ag)
rho_estim(posterior_tc)
rho_estim(posterior_mg)
rho_estim(posterior_se)
rho_estim(posterior_go)

posterior_list <- list(
  ce = posterior_ce,
  bh = posterior_bh,
  pa = posterior_pa,
  pn = posterior_pn,
  rg = posterior_rg,
  pi = posterior_pi,
  ag = posterior_ag,
  tc = posterior_tc,
  mg = posterior_mg,
  se = posterior_se,
  go = posterior_go
)

rho_pool <- purrr::imap_dfr(posterior_list, function(post, region_name){
  tibble(
    Region = region_name,
    rho_id = seq_along(post$rho),
    rho    = as.numeric(post$rho)
  )
})


beta_df <- map_dfr(names(posterior_list), function(region) {
  est <- beta_estim(posterior_list[[region]])  # est is 3 x N_week matrix
  est_df <- as.data.frame(t(est))              # transpose: rows = weeks
  colnames(est_df) <- rownames(est)            # name columns as "50%", "2.5%", etc.
  est_df$week <- seq_len(nrow(est_df))         # add week number
  est_df$region <- region                      # add region name
  est_df
}) %>%
  dplyr::rename(
    beta_p50 = `50%`,
    beta_p2.5 = `2.5%`,
    beta_p97.5 = `97.5%`
  ) %>%
  select(region, week, beta_p50, beta_p2.5, beta_p97.5)

# Combine each rho_estim() result into a row with region label
rho_df <- map_dfr(names(posterior_list), function(region) {
  est <- rho_estim(posterior_list[[region]])
  # If est is a matrix or vector with rownames like "50%", "2.5%", "97.5%"
  tibble(
    region = region,
    rho_p50 = est["50%", 1],
    rho_p2.5 = est["2.5%", 1],
    rho_p97.5 = est["97.5%", 1]
  )
})

region_names <- c("Ceará", "Bahia","Paraíba","Pernambuco" , "Rio Grande do Norte",
                  "Piauí", "Alagoas", "Tocantins", "Minas Gerais", "Sergipe", "Goiás")
rho_df$region_full <- region_names

rho_df <- rho_df[,c(2:5)]
colnames(rho_df)[4] <- "region"

save(rho_df, file = "00_Data/0_2_Processed/rho_df.RData")

gamma_df <- map_dfr(names(posterior_list), function(region) {
  est <- gamma_estim(posterior_list[[region]])
  
  tibble(
    region = region,
    gamma_p50 = est["50%"],
    gamma_p2.5 = est["2.5%"],
    gamma_p97.5 = est["97.5%"]
  )
}) 

param_df <- beta_df %>%
  left_join(gamma_df, by = "region")

param_df <- param_df %>%
  left_join(rho_df, by = "region")

param_df <- param_df %>%
  mutate(region = dplyr::recode(region,
                         "ce" = "Ceará",
                         "ag" = "Alagoas",
                         "bh" = "Bahia",
                         "pa" = "Paraíba",
                         "pn" = "Pernambuco",
                         "tc" = "Tocantins",
                         "rg" = "Rio Grande do Norte",
                         "mg" = "Minas Gerais",
                         "pi" = "Piauí",
                         "se" = "Sergipe",
                         "go" = "Goiás"
                           
  ))

param_df <- param_df %>% mutate(
  r0_p50 = beta_p50 / gamma_p50,
  r0_lo = beta_p2.5 / gamma_p2.5,
  r0_hi = beta_p97.5 / gamma_p97.5
)

beta_state <- 
ggplot(param_df) +
  geom_line(aes(x = week, y = beta_p50))+
  geom_ribbon(aes(x = week, ymin = beta_p2.5, ymax = beta_p97.5, fill = region),
              alpha = 0.3)+
  facet_wrap(~region)+
  theme_pubclean()+
  theme(legend.position = "right")+
  ylab("Predicted transmission rates by state")

ggsave(filename = "02_Outputs/2_1_Figures/figs4.jpg", beta_state, width = 10, height = 7, dpi = 1200)

## 

library(bayesplot)

ppc_dens_overlay(
  y = bra_sum_ce$tot_cases,
  yrep = posterior_ce$pred_cases[1:100, ]
) +
  ggtitle("Posterior Predictive Check: Observed vs. Predicted Cases")

ppc_intervals(bra_sum_ce$tot_cases,posterior_ce$pred_cases[1:100, ]) +
  ggtitle("Posterior Predictive Check: Predictive Intervals of Cases")

print(fit_prevacc_ce, pars = c("base_beta[3]", "base_beta[5]"))  # Check Rhat and ESS
mcmc_trace(as.array(fit_prevacc_ce), pars = "rho")

y_pred_median <- apply(posterior_ce$pred_cases, 2, median)
yrep <- posterior_ce$pred_cases
# Observed data (should be same length)
y_obs <- as.numeric(bra_sum_ce$tot_cases)  # or by time if temporal

ggplot(data.frame(observed = y_obs, predicted = y_pred_median), aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Observed Cases",
    y = "Predicted Cases (Posterior Median)",
    title = "Observed vs. Predicted Cases"
  ) +
  theme_minimal()

bayes_R2 <- function(y, yrep_draws) {
  var_fit <- apply(yrep_draws, 1, var)                   # variance of predictions (across data points) per draw
  var_resid <- apply(yrep_draws, 1, function(yp) var(y - yp))  # variance of residuals
  var_fit / (var_fit + var_resid)
}

r2_samples <- bayes_R2(y_obs, yrep)

# 3. Summarize
mean(r2_samples)       # Posterior mean R²
quantile(r2_samples, c(0.025, 0.975))


obs_list <- list(
  "Ceará" = bra_sum_ce$cases,
  "Pernambuco" = bra_sum_pn$cases,
  "Minas Gerais" = bra_sum_mg$cases,
  "Bahia" = bra_sum_bh$cases,
  "Paraíba" = bra_sum_pa$cases,
  "Rio Grande do Norte" = bra_sum_rg$cases,
  "Piauí" = bra_sum_pi$cases,
  "Alagoas" = bra_sum_ag$cases,
  "Tocantins" = bra_sum_tc$cases,
  "Sergipe"  = bra_sum_se$cases,
  "Goiás"   = bra_sum_go$cases
)

yrep_list <- list(
  "Ceará" = posterior_ce$pred_cases,
  "Pernambuco" = posterior_pn$pred_cases,
  "Minas Gerais" = posterior_mg$pred_cases,
  "Bahia" = posterior_bh$pred_cases,
  "Paraíba" = posterior_pa$pred_cases,
  "Rio Grande do Norte" = posterior_rg$pred_cases,
  "Piauí" = posterior_pi$pred_cases,
  "Alagoas" = posterior_ag$pred_cases,
  "Tocantins" = posterior_tc$pred_cases,
  "Sergipe"  = posterior_se$pred_cases,
  "Goiás"   = posterior_go$cases
)

# Preprocess
ppc_df <- prep_ppc_interval_data(obs_list, yrep_list)
obs_pred_df <- prep_obs_pred_data(obs_list, yrep_list)
r2_df <- compute_all_r2(obs_list, yrep_list)

# Plot
bayes_pred_plot <- plot_obs_pred_facet(obs_pred_df)


ggsave(filename = "02_Outputs/2_1_Figures/bayes_pred_plot.jpg", bayes_pred_plot, width = 10, height = 7, dpi = 1200)



#-------------------------------------------------------------------------------
fit_ag_summ <- as.data.frame(rstan::summary(fit_prevacc_ag)$summary)
fit_ag_summ$parameter <- rownames(fit_ag_summ)
fit_bh_summ <- as.data.frame(rstan::summary(fit_prevacc_bh)$summary)
fit_bh_summ$parameter <- rownames(fit_bh_summ)
fit_ce_summ <- as.data.frame(rstan::summary(fit_prevacc_ce)$summary)
fit_ce_summ$parameter <- rownames(fit_ce_summ)
fit_mg_summ <- as.data.frame(rstan::summary(fit_prevacc_mg)$summary)
fit_mg_summ$parameter <- rownames(fit_mg_summ)
fit_pn_summ <- as.data.frame(rstan::summary(fit_prevacc_pn)$summary)
fit_pn_summ$parameter <- rownames(fit_pn_summ)
fit_pi_summ <- as.data.frame(rstan::summary(fit_prevacc_pi)$summary)
fit_pi_summ$parameter <- rownames(fit_pi_summ)
fit_go_summ <- as.data.frame(rstan::summary(fit_prevacc_go)$summary)
fit_go_summ$parameter <- rownames(fit_go_summ)
fit_se_summ <- as.data.frame(rstan::summary(fit_prevacc_se)$summary)
fit_se_summ$parameter <- rownames(fit_se_summ)
fit_rg_summ <- as.data.frame(rstan::summary(fit_prevacc_rg)$summary)
fit_rg_summ$parameter <- rownames(fit_rg_summ)
fit_tc_summ <- as.data.frame(rstan::summary(fit_prevacc_tc)$summary)
fit_tc_summ$parameter <- rownames(fit_tc_summ)
fit_pa_summ <- as.data.frame(rstan::summary(fit_prevacc_pa)$summary)
fit_pa_summ$parameter <- rownames(fit_pa_summ)

fit_df_func <- function(fit_summ, region){
  
  fit_df <- fit_summ %>%
    #filter(grepl("^base_beta\\[", parameter) | parameter %in% c("rho", "gamma", "sigma")) %>%
    filter(grepl("^base_beta\\[", parameter)) %>%
    select(parameter, `50%`, `2.5%`, `97.5%`, n_eff, Rhat) %>%
    dplyr::rename(
      median = `50%`,
      lower_95 = `2.5%`,
      upper_95 = `97.5%`,
      ESS = n_eff,
      Rhat = Rhat
    )%>%
    mutate(region = region)
  
  return(fit_df)
}

fit_ag_summ_df <- fit_df_func(fit_ag_summ, "Alagoas")
fit_bh_summ_df <- fit_df_func(fit_bh_summ, "Bahia")
fit_ce_summ_df <- fit_df_func(fit_ce_summ, "Ceará")
fit_mg_summ_df <- fit_df_func(fit_mg_summ, "Minas Gerais")
fit_pn_summ_df <- fit_df_func(fit_pn_summ, "Pernambuco")
fit_pi_summ_df <- fit_df_func(fit_pi_summ, "Piauí")
fit_go_summ_df <- fit_df_func(fit_go_summ, "Goiás")
fit_rg_summ_df <- fit_df_func(fit_rg_summ, "Rio Grande do Norte")
fit_tc_summ_df <- fit_df_func(fit_tc_summ, "Tocantins")
fit_pa_summ_df <- fit_df_func(fit_pa_summ, "Paraíba")
fit_se_summ_df <- fit_df_func(fit_se_summ, "Sergipe")

posterior_check_all <- rbind(fit_ag_summ_df,
                             fit_bh_summ_df,
                             fit_ce_summ_df,
                             fit_mg_summ_df,
                             fit_pn_summ_df,
                             fit_pi_summ_df,
                             fit_go_summ_df,
                             fit_rg_summ_df,
                             fit_tc_summ_df,
                             fit_pa_summ_df,
                             fit_se_summ_df
)

posterior_check_all <- posterior_check_all %>%
  mutate(parameter = str_replace(parameter, "^base_beta\\[(\\d+)\\]$", "beta_week\\1"))

write_xlsx(posterior_check_all, path = "02_Outputs/2_2_Tables/table_s3.xlsx")

calc_rho_beta_cor <- function(posterior_obj) {
  
  # 1) Posterior draws
  rho_draws  <- posterior_obj$rho            # vector of MCMC draws for rho (length ~4,000)
  beta_draws <- posterior_obj$base_beta      # matrix [draw x week]  e.g. 4000 x 52
  
  # 2) Beta summary per draw (ex: mean beta across weeks)
  beta_mean <- rowMeans(beta_draws)
  
  # 3) Compute correlation (Spearman or Pearson)
  cor_val <- cor(rho_draws, beta_mean, method = "spearman")
  
  return(cor_val)
}

# Apply function to each state posterior
posterior_list <- list(
  ce = posterior_ce,   # Ceará
  bh = posterior_bh,   # Bahia
  pa = posterior_pa,   # Paraíba
  pn = posterior_pn,   # Pernambuco
  rg = posterior_rg,   # Rio Grande do Norte
  pi = posterior_pi,   # Piauí
  ag = posterior_ag,   # Alagoas
  tc = posterior_tc,   # Tocantins
  mg = posterior_mg,   # Minas Gerais
  se = posterior_se,   # Sergipe
  go = posterior_go    # Goiás
)
rho_beta_cor_list <- map_dbl(posterior_list, calc_rho_beta_cor)

rho_beta_cor_list

rho_beta_cor_df <- as.data.frame(rho_beta_cor_list)

write.csv(rho_beta_cor_df, file = "rho_beta_cor_R1_comment3.csv")

rho_df_plot <- rho_df %>%
  mutate(region = factor(region, levels = region[order(rho_p50)]))

ggplot(rho_df_plot, aes(x = rho_p50, y = region)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = rho_p2.5, xmax = rho_p97.5), height = 0) +
  labs(x = "Predicted reporting rate (median and 95%UI)", y = NULL) +
  theme_bw()

symp_df <- data.frame(
  year         = c(2015:2024),
  reported     = c(16411, 558542, 195962, 87687, 178147,
                   98177, 132587, 265289, 266297, 425587),
  predicted_med = rep(883544, 10),
  predicted_lo  = rep(684889, 10),
  predicted_hi  = rep(1079146, 10)
)

rho_nat <- mean(symp_df$reported) / mean(symp_df$predicted_med)
rho_nat

symp_df <- symp_df %>%
  mutate(pred_reported_nat = predicted_med * rho_nat)

symp_df

write.csv(rho_df, file = "rho_df.csv")

