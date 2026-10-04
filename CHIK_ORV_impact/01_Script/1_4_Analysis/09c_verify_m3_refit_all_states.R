# =============================================================================
# 09c_verify_m3_refit_all_states.R
#
# Verification of the age_struc_bra_seir_m3.stan refit (rho = p_sym * rho_sym
# separation) across all 11 states, extending the single-state (Bahia) check
# to the full set. Two outputs:
#
#   1) bra_foi_state_summ + a new H_2022 column = median_over_draws(sum_t
#      phi_pred[t]) from the NEW (m3) fit -- the 2022 cumulative infection
#      hazard implied by the refit, alongside the existing serosurvey-derived
#      long-term avg_foi column (a genuinely different quantity: avg_foi is
#      lifetime-average annual FOI from the catalytic model; H_2022 is the
#      2022 outbreak's own cumulative hazard from the case-fitting model).
#
#   2) Old (m2) vs New (m3) comparison table, one row per state: rho_sym,
#      effective rho (=p_sym*rho_sym), total new_exposed, cumulative
#      phi_pred, AR (=1-exp(-sum phi)), and a PPC quality check (correlation
#      + RMSE of predicted vs observed weekly reported cases).
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
suppressMessages({library(rstan); library(dplyr)})

p_sym <- 0.524

state_abbr <- c(
  "Ceará" = "ce", "Bahia" = "bh", "Paraíba" = "pa", "Pernambuco" = "pn",
  "Rio Grande do Norte" = "rg", "Piauí" = "pi", "Alagoas" = "ag",
  "Tocantins" = "tc", "Minas Gerais" = "mg", "Sergipe" = "se", "Goiás" = "go"
)

load("00_Data/0_2_Processed/fits_prevacc_finite.RData")   # fit_prevacc_XX_finite x11
load("00_Data/0_2_Processed/posterior_finite_all.RData")  # posterior_XX (old, m2) x11
load("00_Data/0_2_Processed/observed_2022.RData")         # observed_XX (Week, Observed) x11
load("00_Data/0_2_Processed/bra_foi_state_summ.RData")    # bra_foi_state_summ

results <- list()

for (region in names(state_abbr)) {
  ab <- state_abbr[[region]]
  fit_obj  <- get(paste0("fit_prevacc_", ab, "_finite"))
  post_old <- get(paste0("posterior_", ab))
  obs_df   <- get(paste0("observed_", ab))

  post_new <- rstan::extract(fit_obj)

  rho_sym_med   <- median(post_new$rho_sym)
  rho_eff_new   <- p_sym * post_new$rho_sym
  rho_eff_med   <- median(rho_eff_new)
  rho_old_med   <- median(post_old$rho)

  tot_inf_old <- rowSums(post_old$incident_infections)
  tot_inf_new <- rowSums(post_new$incident_infections)

  sum_phi_old <- rowSums(post_old$phi_pred)
  sum_phi_new <- rowSums(post_new$phi_pred)

  AR_old <- 1 - exp(-sum_phi_old)
  AR_new <- 1 - exp(-sum_phi_new)

  pred_old_weekly_med <- apply(post_old$pred_cases, 2, median)
  pred_new_weekly_med <- apply(post_new$pred_cases, 2, median)
  obs_weekly <- obs_df$Observed

  results[[region]] <- data.frame(
    region              = region,
    rho_sym_med         = rho_sym_med,
    rho_old_med         = rho_old_med,
    rho_eff_new_med     = rho_eff_med,
    tot_new_exposed_old = median(tot_inf_old),
    tot_new_exposed_new = median(tot_inf_new),
    ratio_infections    = median(tot_inf_new) / median(tot_inf_old),
    sum_phi_old         = median(sum_phi_old),
    sum_phi_new         = median(sum_phi_new),
    H_2022_new          = median(sum_phi_new),
    AR_old              = median(AR_old),
    AR_new              = median(AR_new),
    AR_pct_change       = 100 * (median(AR_new) - median(AR_old)) / median(AR_old),
    ppc_cor_old         = cor(pred_old_weekly_med, obs_weekly),
    ppc_cor_new         = cor(pred_new_weekly_med, obs_weekly),
    ppc_rmse_old        = sqrt(mean((pred_old_weekly_med - obs_weekly)^2)),
    ppc_rmse_new        = sqrt(mean((pred_new_weekly_med - obs_weekly)^2))
  )
  cat("Done:", region, "\n")
}

comparison_all <- dplyr::bind_rows(results) %>% dplyr::arrange(dplyr::desc(AR_pct_change))

cat("\n=== Old (m2) vs New (m3) comparison, all 11 states ===\n")
print(as.data.frame(comparison_all %>%
  dplyr::select(region, rho_old_med, rho_sym_med, rho_eff_new_med,
                tot_new_exposed_old, tot_new_exposed_new, ratio_infections,
                AR_old, AR_new, AR_pct_change)), digits = 4)

cat("\n=== PPC quality check (should stay similar old vs new) ===\n")
print(as.data.frame(comparison_all %>%
  dplyr::select(region, ppc_cor_old, ppc_cor_new, ppc_rmse_old, ppc_rmse_new)), digits = 4)

# ---- Add H_2022 column to bra_foi_state_summ (sf object -- avoid dplyr verbs
# on the geometry column, use base R indexing instead) ----
H_2022_lookup <- setNames(comparison_all$H_2022_new, comparison_all$region)
bra_foi_state_summ$H_2022 <- H_2022_lookup[bra_foi_state_summ$NAME_1]

foi_tbl <- data.frame(
  NAME_1  = bra_foi_state_summ[["NAME_1"]],
  avg_foi = bra_foi_state_summ[["avg_foi"]],
  H_2022  = bra_foi_state_summ[["H_2022"]]
)
foi_tbl <- foi_tbl[!is.na(foi_tbl$H_2022), ]
foi_tbl <- foi_tbl[order(-foi_tbl$H_2022), ]

cat("\n=== bra_foi_state_summ with new H_2022 column (2022 cumulative hazard, m3 refit) ===\n")
print(foi_tbl, digits = 4, row.names = FALSE)

save(comparison_all, bra_foi_state_summ,
     file = "00_Data/0_2_Processed/m3_refit_verification_all_states.RData")
cat("\nSaved: 00_Data/0_2_Processed/m3_refit_verification_all_states.RData\n")
