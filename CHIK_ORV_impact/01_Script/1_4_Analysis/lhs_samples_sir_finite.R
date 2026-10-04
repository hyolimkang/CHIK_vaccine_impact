## LHS sampling for the finite-history pipeline -----------------------------
## Mirrors the lhs_combined / lhs_idx_list / posterior_idx_list section of
## lhs_samples_sir.R, but keeps that original (long-term FOI) pipeline
## untouched. This version uses:
##   - posterior_ce...posterior_go from the finite-history fits
##     (fit_prevacc_*_finite, extracted via extract_params() in
##     age_struc_fitting_region_func_updated.R)
##   - foi_draws_list (long-term average FOI, from make_region_foi_draws() in
##     age_struc_fitting_region_func_updated.R) -- matches the FOI source
##     make_finite_prevacc_data() now uses to build sero_finite when fitting
##     (bra_foi_state_summ$avg_foi, capped at 8 years). foi_draws_list_shortterm
##     (bra_state_short_term_foi$H_median-based) is no longer used here since
##     fitting reverted from short-term to long-term FOI.
##
## NOTE: simulate_pre_ui_age()/run_simulation_scenarios_ui_ixchiq() look up
## posterior_idx_list/lhs_idx_list as globals (they are not function
## parameters). To actually run the finite pipeline, point those globals at
## the _finite versions right before calling, e.g.:
##   posterior_idx_list <- posterior_idx_list_finite
##   lhs_idx_list        <- lhs_idx_list_finite
## and pass lhs_sample = lhs_combined_finite[[region]] via the regions list.
## Restore posterior_idx_list/lhs_idx_list from the originals afterwards if
## you need to re-run the long-term-FOI pipeline in the same session.

library(MASS)
library(truncnorm)

load("00_Data/0_2_Processed/posterior_finite_all.RData")

posterior_list_finite <- list(
  "Ceará" = posterior_ce,
  "Bahia"= posterior_bh,
  "Paraíba" = posterior_pa,
  "Pernambuco" = posterior_pn,
  "Rio Grande do Norte" = posterior_rg,
  "Piauí" = posterior_pi,
  "Alagoas" = posterior_ag,
  "Tocantins" = posterior_tc,
  "Minas Gerais" = posterior_mg,
  "Sergipe" = posterior_se,
  "Goiás" = posterior_go
)

# New (finite-history) fits store the weekly transmission trajectory as
# beta_observed rather than base_beta -- same rename handled by
# extract_params()/simulate_pre_ui_age() in
# age_struc_fitting_region_func_updated.R.
get_beta_draws <- function(posterior) {
  if ("beta_observed" %in% names(posterior)) {
    posterior$beta_observed
  } else if ("base_beta" %in% names(posterior)) {
    posterior$base_beta
  } else {
    stop("Neither beta_observed nor base_beta exists.")
  }
}

set.seed(123)
n_draws <- 1000
D <- randomLHS(n_draws, 12)

lhs_combined_finite <- lapply(names(posterior_list_finite), function(reg) {

  posterior  <- posterior_list_finite[[reg]]

  # rho/gamma/sigma are taken DIRECTLY from the Stan posterior (no
  # refit-to-lognormal-then-independently-resample) so they stay jointly
  # paired, per Stan draw, with beta_observed/I0 -- run_simulation_scenarios_ui_ixchiq()
  # now indexes all five with the SAME idx_post, preserving the fitted joint
  # posterior correlation instead of combining independently-resampled
  # marginals drawn via a separate index (idx_lhs). That mismatch previously
  # made the resimulated FOI diverge from the Stan fit's own phi_pred -- see
  # chat record, diagnostic scripts 09_diagnose_ar_incidence_vs_hazard.R /
  # 09b_diagnose_rho_consistency.R. I0/base_beta are no longer built here:
  # run_simulation_scenarios_ui_ixchiq() already reads
  # posterior$I0[idx_post,]/posterior$beta_observed[idx_post,] directly and
  # never consumed lhs_sample$I0/base_beta.

  # LHS sampling (foi/ve/vc/wd are genuinely independent of the Stan fit, so
  # remain LHS-sampled from their own priors as before)
  data.frame(
    rho       = posterior$rho,
    gamma     = posterior$gamma,
    sigma     = posterior$sigma,

    foi       = foi_draws_list[[reg]],

    # Ixchiq VE
    ve_ix     = qtruncnorm(D[,6], a = 0.967, b = 0.998,
                           mean = 0.989,
                           sd   = (0.998 - 0.967) / (2 * qnorm(0.975))),

    # Vimkunya VE
    ve_vimkun = qtruncnorm(D[,6], a = 0.972, b = 0.983,
                           mean = 0.978,
                           sd   = (0.983 - 0.972) / (2 * qnorm(0.975))),

    vc   = qtruncnorm(D[,7], a = 0.4, b = 0.6,
                      mean = 0.5,
                      sd   = (0.6 - 0.4) / (2 * qnorm(0.975))),

    vc10 = qtruncnorm(D[,7], a = 0.05, b = 0.15, mean = 0.10,
                      sd = (0.15 - 0.05) / (2 * qnorm(0.975))),

    vc50 = qtruncnorm(D[,7], a = 0.40, b = 0.60, mean = 0.50,
                      sd = (0.60 - 0.40) / (2 * qnorm(0.975))),

    vc90 = qtruncnorm(D[,7], a = 0.85, b = 0.95, mean = 0.90,
                      sd = (0.95 - 0.85) / (2 * qnorm(0.975))),

    wd   = qtruncnorm(D[,8], a = 0.1*0.9, b = 0.1*1.1,
                      mean = 0.1,
                      sd   = (0.1*1.1 - 0.1*0.9) / (2 * qnorm(0.975)))
  )
})

names(lhs_combined_finite) <- names(posterior_list_finite)

# for reproducibility, create LHS sample ids
set.seed(123)
n_draws <- 1000

lhs_idx_list_finite <- lapply(lhs_combined_finite, function(df) {
  sample(seq_len(nrow(df)), size = n_draws, replace = FALSE)
})
names(lhs_idx_list_finite) <- names(lhs_combined_finite)

# Identity mapping: rho/gamma/sigma now come directly from the Stan
# posterior (same data.frame row order as the posterior itself), so draw d
# must index the posterior's own draw d -- a random resample-with-replacement
# here would re-break the beta_observed/I0 <-> rho/gamma/sigma pairing this
# fix exists to restore.
posterior_idx_list_finite <- lapply(posterior_list_finite, function(posterior) {
  seq_len(nrow(posterior$I0))
})
names(posterior_idx_list_finite) <- names(posterior_list_finite)

save(
  posterior_list_finite, lhs_combined_finite,
  lhs_idx_list_finite, posterior_idx_list_finite,
  file = "00_Data/0_2_Processed/lhs_combined_finite.RData"
)
