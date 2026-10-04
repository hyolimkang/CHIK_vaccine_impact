# =============================================================================
# age_struc_fitting_region_finite_2022_4chains.R
#
# Re-fit of the 11-state finite-history age-structured Stan model, matching
# the manuscript Methods text (chat record 2026-09-02):
#   "Four Markov chains were run for 5,000 iterations each, including 1,000
#    warm-up iterations, yielding 4,000 post-warm-up draws per chain."
#   -> chains = 4, iter = 5000, warmup = 1000 (4000 post-warmup draws/chain,
#      16000 total across chains -- vs the ORIGINAL 1-chain/2000-iter/1000-
#      warmup fit in age_struc_fitting_region_finite_2022.R, which yields
#      only 1000 post-warmup draws total).
#
# REDESIGN (chat record 2026-09-04) after the first attempt at this run was
# killed ~14h in, stuck on state 10/11 (Sergipe), losing all 9 already-
# completed states' work -- because the original version only called save()
# once at the very end. Two problems fixed here:
#
#   1) LOST WORK ON INTERRUPT. Each state is now checkpointed to its own
#      .rds file (00_Data/0_2_Processed/checkpoints_4chains/) IMMEDIATELY
#      after fitting, and the per-state loop SKIPS any state whose
#      checkpoint already exists. Re-running this script after a kill/crash
#      resumes from the next unfinished state instead of restarting from
#      scratch.
#
#   2) MEMORY / DISK BLOW-UP. 4 chains x 5000 iter = 16,000 post-warmup
#      draws/state vs. 1,000 before (16x). The original version's OOM
#      ("cannot allocate vector of size 67.7 Mb") happened DURING sampling,
#      not while saving -- rstan accumulates a full draws-by-parameter array
#      per iteration for every returned parameter, so an unrestricted
#      sampling() call keeps growing that accumulator for every large
#      per-age-per-week matrix (S/E/I/R, S_pred/E_pred/I_pred/R_pred,
#      log_lik, age_stratified_cases, expected_reported_cases, ...) even
#      though nothing downstream reads them. KEEP_PARS below restricts
#      sampling() to the ~14 scalar/short-vector quantities actually used
#      (by extract_params() here, by 45_figure2B_remaining_risk_brr.R's
#      posterior_ce$lambda, and by posterior_diagnostics_finite_4chains.R's
#      trace/ESS/PPC), which is what actually caps both the peak RAM during
#      sampling AND the final file size -- back-of-envelope, the excluded
#      A x T matrices alone accounted for an estimated ~20GB across 11
#      states at 16,000 draws/state. The raw stanfit object (which ALSO
#      carries warmup iterations and per-iteration sampler diagnostics on
#      top of every parameter) is no longer saved to disk at all -- the 3
#      things posterior_diagnostics_finite.R needs from it (trace array,
#      ESS/Rhat summary, pred_cases for the PPC) are computed immediately
#      after each state's fit, while the stanfit is still in memory, and
#      checkpointed instead; see posterior_diagnostics_finite_4chains.R.
#
# Chains are parallelised via cores = 4 (already the case in the original
# version of this script) -- states are still fit SEQUENTIALLY, not in
# parallel with each other: a single state's 4 parallel chains already used
# ~26GB+ RAM in the first attempt (this machine: 65GB total, 28 logical
# cores), so running 2+ states' chains concurrently risked another OOM for
# no benefit (total wall time is dominated by the slowest state either way).
#
# IMPORTANT: saves under NEW filenames (*_4chains suffix), does NOT
# overwrite fits_prevacc_finite.RData / posterior_finite_all.RData. Every
# downstream script in both this project and CHIK_benefit_risk currently
# reads the un-suffixed (1-chain) files, and there is a known latent bug
# (01_Script/1_4_Analysis/lhs_samples_sir_finite.R hardcodes n_draws <- 1000
# when pairing LHS draws to posterior draws in simulate_pre_ui_age()) that
# will silently propagate NAs if fed a >1000-draw posterior without a fix
# first. Swapping the production files for the 4-chain fit is a deliberate
# follow-up step, not part of this run.
#
# Output:
#   00_Data/0_2_Processed/checkpoints_4chains/state_<code>_4chains.rds
#     -- one per state, written as soon as that state finishes (resume point)
#   00_Data/0_2_Processed/posterior_finite_all_4chains.RData
#     -- posterior_ce, posterior_bh, ... (drop-in compatible with the
#        1-chain posterior_finite_all.RData's list shape/field names, minus
#        the large unused matrices listed above)
#   00_Data/0_2_Processed/diagnostics_finite_4chains.RData
#     -- trace_arr / diag_summary / pred_cases per state, consumed by
#        posterior_diagnostics_finite_4chains.R (run automatically at the
#        end of this script -- see bottom)
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
source("01_Script/1_1_Functions/library.R")

# NOTE: deliberately NOT sourcing 01_Script/1_1_Functions/age_struc_fitting_
# region_func_updated.R -- that file is 8716 lines and mixes function defs
# with top-level executed code throughout (verified: line 333 alone runs
# make_region_foi_draws(regions, bra_foi_states) unconditionally at source
# time, requiring bra_foi_states to already exist in the calling env, and
# there is no guarantee that's the ONLY such ordering dependency in a file
# that size). Inlining just the 2 functions actually needed here avoids
# executing thousands of lines of unrelated/unverified side-effect code
# during an unattended background run.

N_CHAINS <- 4
N_ITER   <- 5000
N_WARMUP <- 1000
N_CORES  <- 4

# Scalar/short-vector quantities actually consumed downstream -- everything
# else (S/E/I/R and their *_pred twins, new_exposed/new_infectious/
# new_recovered, expected_reported_cases, age_stratified_cases, log_lik,
# incident_infections/incident_infectious_onsets/aggregated_infections,
# infectious_prevalence, R_eff, phi_pred, N_pred, susceptible_fraction,
# p_infection, beta_rw_raw/beta_rw_centered, z_beta, pop_total) is an A x T
# or A x (T+1) matrix (A=20, T~52-54) that nothing downstream reads -- these
# are the ones responsible for the ~20GB estimated blow-up at 16,000
# draws/state and are excluded on purpose.
KEEP_PARS <- c(
  "I0", "alpha_log_beta", "sigma_beta_rw", "shape", "gamma", "sigma",
  "rho_sym", "rho", "beta", "beta_observed", "lambda",
  "pred_cases", "lp__"
)
DIAG_PARS <- c("rho", "gamma", "sigma", "alpha_log_beta", "sigma_beta_rw", "lp__")

CKPT_DIR <- "00_Data/0_2_Processed/checkpoints_4chains"
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[%s] Starting 4-chain finite-history refit (chains=%d, iter=%d, warmup=%d)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), N_CHAINS, N_ITER, N_WARMUP))

# ---- 0) Data (verbatim from age_struc_fitting_region_finite_2022.R) --------
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

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_bra_seir_m3.stan")

cat(sprintf("[%s] Loading allfoi_s1.RData (~3GB, this is slow)...\n", format(Sys.time(), "%H:%M:%S")))
load("00_Data/0_2_Processed/allfoi_s1.RData")  # -> allfoi
bra_foi <- allfoi %>% filter(country == "Brazil")
rm(allfoi); gc()
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))
bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)

bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = mean(foi_mid, na.rm = TRUE)) %>%
  filter(!is.na(NAME_1))

age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))

make_finite_prevacc_data <- function(
    old_data, state_name, fitting_year, introduction_year = 2014L,
    age_midpoints = age_groups, B = 2L
) {
  observed_cases <- round(as.matrix(old_data$observed_cases_by_age))
  storage.mode(observed_cases) <- "integer"
  A <- nrow(observed_cases)
  T <- ncol(observed_cases)
  stopifnot(length(age_midpoints) == A)

  state_foi <- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == state_name]
  stopifnot("state_foi not found for state_name in bra_foi_state_summ" = length(state_foi) == 1)

  years_since_introduction <- fitting_year - introduction_year
  exposure_years_finite <- pmin(age_midpoints, years_since_introduction)
  sero_finite <- 1 - exp(-state_foi * exposure_years_finite)
  prior_I0 <- pmax(as.numeric(old_data$prior_I0), 1)

  list(
    T = as.integer(T), B = as.integer(B), A = as.integer(A),
    observed_cases_by_age = observed_cases,
    N = as.numeric(old_data$N),
    r = rep(0, A),
    prior_I0 = prior_I0,
    prior_sd_I0 = as.numeric(old_data$prior_sd_I0),
    sero = as.numeric(sero_finite),
    p_sym = 0.524
  )
}

fitting_years <- c(
  "Bahia" = 2022L, "Ceará" = 2022L, "Minas Gerais" = 2022L, "Pernambuco" = 2022L,
  "Paraíba" = 2022L, "Rio Grande do Norte" = 2022L, "Piauí" = 2022L,
  "Alagoas" = 2022L, "Tocantins" = 2022L, "Sergipe" = 2022L, "Goiás" = 2022L
)

prevacc_longterm_data <- list(
  "Bahia" = list(observed_cases_by_age = observed_cases_bahia, N = as.numeric(N_bahia$Bahia),
                 prior_I0 = observed_cases_bahia[, 1], prior_sd_I0 = 100),
  "Ceará" = list(observed_cases_by_age = observed_cases_ce, N = as.numeric(N_ceara$Ceará),
                 prior_I0 = observed_cases_ce[, 1], prior_sd_I0 = 100),
  "Minas Gerais" = list(observed_cases_by_age = observed_cases_mg, N = as.numeric(N_mg$`Minas Gerais`),
                        prior_I0 = observed_cases_mg[, 1], prior_sd_I0 = 100),
  "Pernambuco" = list(observed_cases_by_age = observed_cases_pn, N = as.numeric(N_pemam$Pernambuco),
                      prior_I0 = observed_cases_pn[, 1], prior_sd_I0 = 100),
  "Paraíba" = list(observed_cases_by_age = observed_cases_pa, N = as.numeric(N_pa$Paraíba),
                   prior_I0 = observed_cases_pa[, 1], prior_sd_I0 = 100),
  "Rio Grande do Norte" = list(observed_cases_by_age = observed_cases_rg, N = as.numeric(N_rg$`Rio Grande do Norte`),
                               prior_I0 = observed_cases_rg[, 1], prior_sd_I0 = 100),
  "Piauí" = list(observed_cases_by_age = observed_cases_pi, N = as.numeric(N_pi$Piauí),
                 prior_I0 = observed_cases_pi[, 1], prior_sd_I0 = 100),
  "Alagoas" = list(observed_cases_by_age = observed_cases_ag, N = as.numeric(N_ag$Alagoas),
                   prior_I0 = observed_cases_ag[, 1], prior_sd_I0 = 100),
  "Tocantins" = list(observed_cases_by_age = observed_cases_tc, N = as.numeric(N_tc$Tocantins),
                     prior_I0 = observed_cases_tc[, 1], prior_sd_I0 = 100),
  "Sergipe" = list(observed_cases_by_age = observed_cases_se, N = as.numeric(N_se$Sergipe),
                   prior_I0 = observed_cases_se[, 1], prior_sd_I0 = 100),
  "Goiás" = list(observed_cases_by_age = observed_cases_go, N = as.numeric(N_go$Goiás),
                 prior_I0 = observed_cases_go[, 1], prior_sd_I0 = 100)
)

prevacc_finite_data <- lapply(names(prevacc_longterm_data), function(state_name) {
  make_finite_prevacc_data(
    old_data = prevacc_longterm_data[[state_name]], state_name = state_name,
    fitting_year = fitting_years[[state_name]], introduction_year = 2014L,
    age_midpoints = age_groups, B = 2L
  )
})
names(prevacc_finite_data) <- names(prevacc_longterm_data)

# ---- 1) Fitting: 4 chains x 5000 iter x 1000 warmup, per state, checkpointed
# adapt_delta/max_treedepth kept identical to the original 1-chain tuning per
# state (those were set to address that state's specific sampling geometry;
# more chains/iterations doesn't remove the need for them).
state_specs <- list(
  list(label = "Bahia",               code = "bh", adapt_delta = 0.95, max_treedepth = 12),
  list(label = "Ceará",                code = "ce", adapt_delta = 0.95, max_treedepth = 12),
  list(label = "Minas Gerais",         code = "mg", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Pernambuco",           code = "pn", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Paraíba",              code = "pa", adapt_delta = 0.8,  max_treedepth = 15),
  list(label = "Rio Grande do Norte",  code = "rg", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Piauí",                code = "pi", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Alagoas",              code = "ag", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Tocantins",            code = "tc", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Sergipe",              code = "se", adapt_delta = 0.95, max_treedepth = 15),
  list(label = "Goiás",                code = "go", adapt_delta = 0.95, max_treedepth = 15)
)

fit_and_checkpoint_state <- function(spec) {
  ckpt_path <- file.path(CKPT_DIR, sprintf("state_%s_4chains.rds", spec$code))
  if (file.exists(ckpt_path)) {
    cat(sprintf("[%s] %s already checkpointed -- skipping (resume).\n",
                format(Sys.time(), "%H:%M:%S"), spec$label))
    return(invisible(NULL))
  }

  data <- prevacc_finite_data[[spec$label]]
  cat(sprintf("[%s] Fitting %s (finite, 4 chains)...\n", format(Sys.time(), "%H:%M:%S"), spec$label))
  fit <- sampling(
    object = stan_model_age, data = data,
    iter = N_ITER, chains = N_CHAINS, warmup = N_WARMUP, thin = 1, seed = 123,
    cores = N_CORES, pars = KEEP_PARS,
    control = list(adapt_delta = spec$adapt_delta, max_treedepth = spec$max_treedepth)
  )
  cat(sprintf("[%s] Done sampling %s -- extracting + computing diagnostics before discarding raw fit...\n",
              format(Sys.time(), "%H:%M:%S"), spec$label))

  post <- rstan::extract(fit)
  beta_draws <- if ("beta_observed" %in% names(post)) post$beta_observed else post$base_beta
  result <- list(
    label = spec$label,
    posterior_prevacc = post,
    base_beta = apply(beta_draws, 2, median),
    I0 = apply(post$I0, 2, median),
    gamma = median(post$gamma),
    rho = median(post$rho),
    trace_arr = as.array(fit, pars = DIAG_PARS),
    diag_summary = rstan::summary(fit, pars = DIAG_PARS)$summary
  )
  saveRDS(result, ckpt_path, compress = "xz")
  cat(sprintf("[%s] Checkpointed %s -> %s (%.0f MB)\n", format(Sys.time(), "%H:%M:%S"),
              spec$label, ckpt_path, file.info(ckpt_path)$size / 1e6))

  rm(fit, post, result); gc()
  invisible(NULL)
}

for (spec in state_specs) {
  fit_and_checkpoint_state(spec)
}

# ---- 2) Combine checkpoints into the final drop-in-compatible outputs ------
cat(sprintf("[%s] All states checkpointed -- combining...\n", format(Sys.time(), "%H:%M:%S")))

all_ckpts <- lapply(state_specs, function(spec) {
  readRDS(file.path(CKPT_DIR, sprintf("state_%s_4chains.rds", spec$code)))
})
names(all_ckpts) <- vapply(state_specs, `[[`, character(1), "code")

code_to_var <- c(ce = "posterior_ce", bh = "posterior_bh", pa = "posterior_pa",
                  pn = "posterior_pn", rg = "posterior_rg", pi = "posterior_pi",
                  ag = "posterior_ag", tc = "posterior_tc", mg = "posterior_mg",
                  se = "posterior_se", go = "posterior_go")

posterior_env <- new.env()
for (code in names(code_to_var)) {
  assign(code_to_var[[code]], all_ckpts[[code]]$posterior_prevacc, envir = posterior_env)
}
save(list = unname(code_to_var), envir = posterior_env,
     file = "00_Data/0_2_Processed/posterior_finite_all_4chains.RData")
cat(sprintf("[%s] Saved: 00_Data/0_2_Processed/posterior_finite_all_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))

trace_arr_list    <- lapply(all_ckpts, `[[`, "trace_arr")
diag_summary_list <- lapply(all_ckpts, `[[`, "diag_summary")
pred_cases_list   <- lapply(all_ckpts, function(x) x$posterior_prevacc$pred_cases)
names(trace_arr_list) <- names(diag_summary_list) <- names(pred_cases_list) <-
  vapply(state_specs, `[[`, character(1), "label")

save(trace_arr_list, diag_summary_list, pred_cases_list, DIAG_PARS,
     file = "00_Data/0_2_Processed/diagnostics_finite_4chains.RData")
cat(sprintf("[%s] Saved: 00_Data/0_2_Processed/diagnostics_finite_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))

cat(sprintf("[%s] ALL DONE -- 4-chain finite-history refit complete. Case data (observed_2022.RData) is unchanged by refitting -- reuse it directly, no need to reconstruct.\n", format(Sys.time(), "%H:%M:%S")))

# ---- 3) Auto-run diagnostics (trace plots / ESS+Rhat table / PPC) ----------
cat(sprintf("[%s] Running posterior_diagnostics_finite_4chains.R ...\n", format(Sys.time(), "%H:%M:%S")))
source("01_Script/1_4_Analysis/posterior_diagnostics_finite_4chains.R")
