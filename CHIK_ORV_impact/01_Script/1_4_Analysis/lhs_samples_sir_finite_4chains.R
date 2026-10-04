# =============================================================================
# lhs_samples_sir_finite_4chains.R
#
# 4-chain counterpart to lhs_samples_sir_finite.R (chat record 2026-09-05).
# Pairs the 4-chain posterior draws (16,000/state -- see
# age_struc_fitting_region_finite_2022_4chains.R) with independently-sampled
# LHS draws (foi/ve/vc/wd), same joint-posterior-preserving design as the
# 1-chain version (rho/gamma/sigma taken directly from the Stan posterior,
# same idx_post as beta_observed/I0 -- NOT an independently resampled
# marginal). Does NOT touch lhs_samples_sir_finite.R, posterior_finite_all.
# RData, or lhs_combined_finite.RData -- the 1-chain production pipeline is
# unaffected.
#
# Fixes the two things that would otherwise silently break with more draws
# than 1,000 (chat record 2026-09-05):
#   1) n_draws <- 1000 (both places in the 1-chain script) -> must equal the
#      actual number of posterior draws used downstream, else data.frame()
#      recycles the 1000-row LHS block instead of giving independent draws
#      (silently, no error/warning, whenever the draw count is an exact
#      multiple of 1000).
#   2) foi_draws_list <- make_region_foi_draws(regions, bra_foi_states) is a
#      top-level side effect of sourcing age_struc_fitting_region_func_
#      updated.R, using that function's n_draws=1000 DEFAULT -- regenerated
#      here at the correct size (function itself is untouched; only the call
#      is redone, after sourcing).
#
# SUBSAMPLE, not full 16,000 (chat record 2026-09-05): propagating all 16,000
# 4-chain draws through the downstream simulation (preui -> model_0803 ->
# all_draws) made sim_results_vc_ixchiq_model_finite_4chains.RData /
# postsim_vc_ixchiq_model_finite_4chains.RData balloon to an estimated
# ~23GB/~46GB (x2 -- both repos save a copy -- ~138GB combined), against the
# 1-chain pipeline's 1.42GB/2.89GB at 1,000 draws. The 4-chain refit's actual
# purpose was Stan convergence (Rhat 0.9998-1.0015, ESS 4,245-27,843 across
# all 11 states x 6 hyperparameters -- confirmed clean), not a 16x inflation
# of the PSA/simulation draw count the rest of the pipeline was built around.
# N_SUBSAMPLE draws are randomly selected (without replacement, seeded) from
# the 16,000 per state -- posterior_idx_list_finite stores which of the
# original 16,000 rows each of the N_SUBSAMPLE downstream "draws" maps to,
# so simulate_pre_ui_age()/run_simulation_scenarios_ui_ixchiq() (which index
# posterior$I0[idx_post,] etc. via posterior_idx_list) pull from the full,
# better-converged posterior without ever materializing more than
# N_SUBSAMPLE rows through the expensive simulation loop.
#
# bra_foi_states.RData (00_Data/0_2_Processed/, 345MB) is ALREADY the cached
# per-pixel FOI/state join -- make_region_foi_draws() resamples from its
# foi1..foi100 columns per state. No need to reprocess the 3GB allfoi_s1.
# RData + shapefile join again; that heavy step already ran once (whenever
# bra_foi_states.RData was first built) and its result is just loaded here.
#
# Output: 00_Data/0_2_Processed/lhs_combined_finite_4chains.RData
#   -> posterior_list_finite (full 16,000-draw posteriors, unchanged --
#      kept only because finite_weeksweep_brr.R/finite_weeksweep_fulldraws.R
#      etc. reference this object in the 1-chain file; a _4chains
#      counterpart for those hasn't been built yet), lhs_combined_finite,
#      lhs_idx_list_finite, posterior_idx_list_finite (all N_SUBSAMPLE rows
#      -- same shape as lhs_combined_finite.RData, just a subsample)
# =============================================================================

N_SUBSAMPLE <- 4000

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
source("01_Script/1_1_Functions/library.R")

message("Loading bra_foi_states.RData (cached, 345MB -- no need to rebuild) ...")
load("00_Data/0_2_Processed/bra_foi_states.RData")  # -> bra_foi_states

message("Sourcing age_struc_fitting_region_func_updated.R (functions + regions + foi_draws_list@1000, to be overridden below) ...")
source("01_Script/1_1_Functions/age_struc_fitting_region_func_updated.R")

message("Loading posterior_finite_all_4chains.RData ...")
load("00_Data/0_2_Processed/posterior_finite_all_4chains.RData")

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

draw_counts <- vapply(posterior_list_finite, function(p) nrow(p$I0), integer(1))
stopifnot("all 11 states must have the same number of posterior draws" =
            length(unique(draw_counts)) == 1)
n_draws_full <- unique(draw_counts)
cat(sprintf("[%s] Posterior draws/state: %d -- subsampling to %d for downstream simulation\n",
            format(Sys.time(), "%H:%M:%S"), n_draws_full, N_SUBSAMPLE))

# Which of the 16,000 posterior rows each of the N_SUBSAMPLE downstream
# "draws" maps to, per state -- NOT identity (unlike the full-16,000-draw
# version this replaces). simulate_pre_ui_age()/run_simulation_scenarios_ui_
# ixchiq() index posterior$I0[idx_post,] etc. via posterior_idx_list, so the
# full posterior arrays never need to be subset in memory -- only these
# index vectors are shorter.
set.seed(123)
posterior_idx_list_finite <- lapply(posterior_list_finite, function(posterior) {
  sort(sample(seq_len(nrow(posterior$I0)), size = N_SUBSAMPLE, replace = FALSE))
})
names(posterior_idx_list_finite) <- names(posterior_list_finite)

n_draws <- N_SUBSAMPLE

# Override the n_draws=1000-default foi_draws_list that sourcing the function
# file just created, at the correct (subsampled) size.
foi_draws_list <- make_region_foi_draws(regions, bra_foi_states, n_draws = n_draws)

library(MASS)
library(truncnorm)

set.seed(123)
D <- randomLHS(n_draws, 12)

lhs_combined_finite <- lapply(names(posterior_list_finite), function(reg) {

  posterior  <- posterior_list_finite[[reg]]
  idx        <- posterior_idx_list_finite[[reg]]

  data.frame(
    rho       = posterior$rho[idx],
    gamma     = posterior$gamma[idx],
    sigma     = posterior$sigma[idx],

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

# for reproducibility, create LHS sample ids -- n_draws == nrow(df) here, so
# this is a full permutation, not a subsample (matches the 1-chain script's
# semantics: only meaningful as a subsample when n_draws < nrow(df)).
set.seed(123)
lhs_idx_list_finite <- lapply(lhs_combined_finite, function(df) {
  sample(seq_len(nrow(df)), size = n_draws, replace = FALSE)
})
names(lhs_idx_list_finite) <- names(lhs_combined_finite)

# posterior_idx_list_finite was already built above (the N_SUBSAMPLE-of-
# n_draws_full subsample) -- rho/gamma/sigma in lhs_combined_finite were
# indexed by that same subsample, so draw d's rho/gamma/sigma (row d of
# lhs_combined_finite) and draw d's I0/beta_observed (posterior_idx_list_
# finite[[reg]][d]) refer to the SAME original posterior row.

save(
  posterior_list_finite, lhs_combined_finite,
  lhs_idx_list_finite, posterior_idx_list_finite,
  file = "00_Data/0_2_Processed/lhs_combined_finite_4chains.RData"
)
cat(sprintf("[%s] Saved: 00_Data/0_2_Processed/lhs_combined_finite_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))
