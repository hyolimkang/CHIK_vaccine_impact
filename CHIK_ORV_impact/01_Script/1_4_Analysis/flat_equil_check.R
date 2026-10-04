# =============================================================================
# flat_equil_check.R
# -----------------------------------------------------------------------------
# WHY prevaccination burden (and therefore BRR) differs by orders of magnitude
# between the three immunity assumptions.
#
# The ONLY structural difference between the three runs is the INITIAL immunity
# profile fed into the SIR as  R0 = 1 - exp(-FOI * x)  ->  R_[,1] = R0 * N,
#                                                          S[,1]  = N - I0 - R0*N
#
#   1) "equil" : long-term average FOI, age-dependent
#                immune(age) = 1 - exp(-FOI_avg * age)      (older = more immune)
#   2) "flat"  : flat immunity, 12-year exposure
#                immune      = 1 - exp(-FOI_avg * 12)        (same for ALL ages)
#   3) "bahia" : single epidemic, huge one-off FOI
#                immune      = 1 - exp(-0.85 * 1)            (same for ALL ages)
#
# beta (transmissibility) comes from the pre-vacc posterior fit and is, in
# principle, shared. The susceptible POOL (S) is what changes, and it scales
# the outbreak size -> prevacc burden -> averted burden -> BRR.
#
# This script needs NO large sim objects for the S / beta diagnostics (they are
# analytic / from the small posterior). The optional realized-burden section
# (SECTION 4) only runs if the (very large) sim/postsim objects are already in
# memory.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# small helper: return object if it exists in the session, else NULL
g <- function(name) if (exists(name, envir = .GlobalEnv)) get(name, envir = .GlobalEnv) else NULL

# -----------------------------------------------------------------------------
# SECTION 0. Inputs (age midpoints, FOI, population by age)
# -----------------------------------------------------------------------------
REGION    <- "Bahia"
BAHIA_FOI <- 0.85   # single-epidemic one-off FOI used in the Bahia run
FLAT_YEARS <- 12    # exposure years used in the flat-immunity run

# 20 age-bin midpoints (same as age_struc_fitting_region_func_updated.R)
age_mid <- g("age_groups")
if (is.null(age_mid)) {
  age_mid <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
               mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
               mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
               mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))
}

age_lab <- g("age_gr_levels")
if (is.null(age_lab) || length(age_lab) != length(age_mid)) {
  age_lab <- c("<1","1-4","5-9","10-11","12-17","18-19","20-24","25-29",
               "30-34","35-39","40-44","45-49","50-54","55-59","60-64",
               "65-69","70-74","75-79","80-84","85+")
}

# average FOI for the region (long-term equilibrium assumption)
foi_avg <- NA_real_
fss <- g("bra_foi_state_summ")
if (!is.null(fss)) {
  fss_df <- tryCatch(sf::st_drop_geometry(fss), error = function(e) as.data.frame(fss))
  hit <- fss_df$avg_foi[fss_df$NAME_1 == REGION]
  if (length(hit) == 1) foi_avg <- as.numeric(hit)
}
if (is.na(foi_avg)) foi_avg <- 0.009060327  # Bahia avg_foi fallback
message(sprintf("[check] %s long-term average FOI = %.6f", REGION, foi_avg))

# population by age (optional, for population-weighted susceptible pool)
N_age <- g("N_bahia")
if (is.list(N_age) && !is.null(N_age[[REGION]])) N_age <- N_age[[REGION]]
N_age <- tryCatch(as.numeric(N_age), error = function(e) NULL)
if (!is.null(N_age) && length(N_age) != length(age_mid)) N_age <- NULL

# -----------------------------------------------------------------------------
# SECTION 1. Initial immunity and susceptible proportion by age
# -----------------------------------------------------------------------------
imm_df <- tibble(
  age_bin = factor(age_lab, levels = age_lab),
  age_mid = age_mid,
  # immune (recovered) fraction = R0 input
  imm_equil = 1 - exp(-foi_avg   * age_mid),
  imm_flat  = 1 - exp(-foi_avg   * FLAT_YEARS),
  imm_bahia = 1 - exp(-BAHIA_FOI * 1)
) %>%
  mutate(
    # susceptible fraction = 1 - immune  (ignoring the tiny I0 seed)
    S_equil = 1 - imm_equil,
    S_flat  = 1 - imm_flat,
    S_bahia = 1 - imm_bahia
  )

cat("\n=== Susceptible proportion by age (S = exp(-FOI * x)) ===\n")
print(
  imm_df %>%
    transmute(age_bin, age_mid,
              S_equil = round(S_equil, 3),
              S_flat  = round(S_flat,  3),
              S_bahia = round(S_bahia, 3),
              `flat/equil` = round(S_flat  / S_equil, 2),
              `bahia/equil` = round(S_bahia / S_equil, 2)),
  n = nrow(imm_df)
)

# long format for plotting
S_long <- imm_df %>%
  select(age_bin, age_mid, S_equil, S_flat, S_bahia) %>%
  pivot_longer(starts_with("S_"), names_to = "assumption",
               values_to = "S_prop", names_prefix = "S_")

p_S <- ggplot(S_long, aes(age_mid, S_prop, colour = assumption)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = paste0(REGION, ": initial susceptible proportion by age"),
       subtitle = "S(age) = exp(-FOI * x); the pool that drives the outbreak",
       x = "Age (midpoint)", y = "Susceptible proportion", colour = "Assumption") +
  theme_bw()
print(p_S)

# -----------------------------------------------------------------------------
# SECTION 2. Population-weighted susceptible pool (the headline driver)
# -----------------------------------------------------------------------------
if (!is.null(N_age)) {
  pool <- imm_df %>%
    summarise(
      S_equil = sum(S_equil * N_age),
      S_flat  = sum(S_flat  * N_age),
      S_bahia = sum(S_bahia * N_age)
    )
  totN <- sum(N_age)
  cat("\n=== Population-weighted susceptible pool ===\n")
  cat(sprintf("  total population        : %s\n", format(round(totN), big.mark = ",")))
  cat(sprintf("  susceptibles  equil     : %s  (%.1f%%)\n",
              format(round(pool$S_equil), big.mark = ","), 100 * pool$S_equil / totN))
  cat(sprintf("  susceptibles  flat      : %s  (%.1f%%)\n",
              format(round(pool$S_flat), big.mark = ","), 100 * pool$S_flat / totN))
  cat(sprintf("  susceptibles  bahia     : %s  (%.1f%%)\n",
              format(round(pool$S_bahia), big.mark = ","), 100 * pool$S_bahia / totN))
  cat(sprintf("\n  susceptible-pool ratio  flat / equil  = %.2fx\n", pool$S_flat  / pool$S_equil))
  cat(sprintf("  susceptible-pool ratio  bahia / equil = %.2fx\n", pool$S_bahia / pool$S_equil))
  cat("  (outbreak size / prevacc burden scales strongly with this pool)\n")
} else {
  cat("\n[note] N_bahia (population by age) not in session -> skipping pool weighting.\n")
  cat("       unweighted mean S: equil = ", round(mean(imm_df$S_equil), 3),
      " | flat = ", round(mean(imm_df$S_flat), 3),
      " | bahia = ", round(mean(imm_df$S_bahia), 3), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# SECTION 3. beta / transmissibility comparison (from the small posterior)
# -----------------------------------------------------------------------------
# R0_transmission = base_beta / gamma. If the same pre-vacc posterior is used,
# beta is identical across assumptions and the burden gap is purely from S.
summ_beta <- function(post, tag) {
  if (is.null(post)) return(NULL)
  bb <- post$base_beta            # draws x T  (or vector)
  gm <- post$gamma
  beta_med <- if (is.matrix(bb)) median(apply(bb, 1, mean), na.rm = TRUE) else median(bb, na.rm = TRUE)
  gam_med  <- median(gm, na.rm = TRUE)
  tibble(source = tag,
         beta_mean_med = beta_med,
         gamma_med     = gam_med,
         R0_transmission = beta_med / gam_med)
}

beta_tbl <- bind_rows(
  summ_beta(g("posterior_bh"), "posterior_bh (in session)"),
  summ_beta(g("posterior_equil"), "posterior_equil"),
  summ_beta(g("posterior_flat"),  "posterior_flat")
)
if (nrow(beta_tbl)) {
  cat("\n=== beta / transmissibility (median across draws) ===\n")
  print(beta_tbl)
} else {
  cat("\n[note] no posterior_* objects in session -> skipping beta comparison.\n")
}

# -----------------------------------------------------------------------------
# SECTION 4. Realized burden comparison (OPTIONAL - needs large objects loaded)
# -----------------------------------------------------------------------------
# Total symptomatic infections summed over age x week, median across draws,
# for a chosen VE / coverage / scenario slot. Compares the two models directly.
# Only runs if the objects are ALREADY in memory (these files are multi-GB).
total_symp <- function(model, region = REGION, ve = 1, cov = 1, scen = 1) {
  if (is.null(model)) return(NA_real_)
  node <- tryCatch(model[[region]][[ve]][[cov]]$scenario_result[[scen]]$sim_result$age_array_raw_symp,
                   error = function(e) NULL)
  if (is.null(node)) return(NA_real_)
  d <- length(dim(node))
  per_draw <- apply(node, d, sum, na.rm = TRUE)   # sum over age & week, per draw
  median(per_draw, na.rm = TRUE)
}

ps_equil <- g("postsim_equil"); ps_flat <- g("postsim_flat")
if (!is.null(ps_equil) || !is.null(ps_flat)) {
  cat("\n=== Realized total symptomatic burden (median draw; VE slot 1, cov 1, scenario 1) ===\n")
  b_eq <- total_symp(ps_equil)
  b_fl <- total_symp(ps_flat)
  cat(sprintf("  equil : %s\n", ifelse(is.na(b_eq), "n/a", format(round(b_eq), big.mark = ","))))
  cat(sprintf("  flat  : %s\n", ifelse(is.na(b_fl), "n/a", format(round(b_fl), big.mark = ","))))
  if (!is.na(b_eq) && !is.na(b_fl) && b_eq > 0)
    cat(sprintf("  flat / equil burden ratio = %.1fx\n", b_fl / b_eq))
  cat("  NOTE: adjust ve/cov/scenario indices to match the no-vaccine baseline slot.\n")
} else {
  cat("\n[note] postsim_equil / postsim_flat not in session -> skipping burden check.\n")
  cat("       (these objects are multi-GB; load one at a time if needed.)\n")
}

# -----------------------------------------------------------------------------
# SECTION 5. pre/post CONSISTENCY + transmissibility check (the real bug hunt)
# -----------------------------------------------------------------------------
# For a NO-real-effect baseline, post-vacc burden should be the SAME order of
# magnitude as pre-vacc burden. A huge pre vs tiny post (or vice versa) means
# preui and postsim were generated with DIFFERENT base_beta / FOI -> the
# pre - post "averted" is an artifact and BRR explodes.
#
# Run this with the objects already in your session:
#   preui_bh_equil, preui_bh_flat (= Bahia), postsim_equil, postsim_flat (= Bahia)

ve_slot <- "VE0"; cov_slot <- "cov50"; scen_slot <- 1

get_node <- function(model) {
  tryCatch(model[["Bahia"]][[ve_slot]][[cov_slot]]$scenario_result[[scen_slot]],
           error = function(e) NULL)
}

# pull transmissibility diagnostics from a sim_out (if present)
beta_diag <- function(node, tag) {
  if (is.null(node) || is.null(node$sim_out)) {
    cat(sprintf("  [%s] no sim_out\n", tag)); return(invisible())
  }
  so <- node$sim_out
  bb <- so$base_beta
  beta1 <- if (!is.null(bb)) bb[1] else NA_real_
  R0v  <- if (!is.null(so$R0_vec))    so$R0_vec[1]    else NA_real_
  Reff <- if (!is.null(so$R_eff_vec)) so$R_eff_vec[1] else NA_real_
  S1   <- if (!is.null(so$S)) sum(so$S[, 1]) else NA_real_
  Ntot <- if (!is.null(so$S)) sum(so$S[, ncol(so$S)]) else NA_real_  # rough
  cat(sprintf("  [%s] base_beta[1]=%.4g  R0_vec[1]=%.3f  R_eff[1]=%.3f  S0=%s\n",
              tag, beta1, R0v, Reff,
              ifelse(is.na(S1), "n/a", format(round(S1), big.mark = ","))))
}

cat("\n=== pre vs post total symptomatic (Bahia) ===\n")
pe <- g("preui_bh_equil"); pf <- g("preui_bh_flat")
if (!is.null(pe)) cat(sprintf("  pre  equil : %s\n", format(round(sum(pe$weekly_cases_median_rawsymp)), big.mark = ",")))
if (!is.null(pf)) cat(sprintf("  pre  bahia : %s\n", format(round(sum(pf$weekly_cases_median_rawsymp)), big.mark = ",")))
ne <- get_node(g("postsim_equil")); nf <- get_node(g("postsim_flat"))
if (!is.null(ne)) cat(sprintf("  post equil : %s\n", format(round(sum(ne$sim_result$weekly_cases_median_rawsymp)), big.mark = ",")))
if (!is.null(nf)) cat(sprintf("  post bahia : %s\n", format(round(sum(nf$sim_result$weekly_cases_median_rawsymp)), big.mark = ",")))

cat("\n=== transmissibility at t=1 (postsim sim_out) ===\n")
beta_diag(ne, "post equil")
beta_diag(nf, "post bahia")
cat("\nINTERPRETATION:\n")
cat("  * If post-bahia R_eff[1] < 1 while pre-bahia is a huge epidemic,\n")
cat("    preui and postsim used DIFFERENT base_beta -> averted is an artifact.\n")
cat("  * Fix: regenerate preui AND postsim for Bahia with the SAME posterior/\n")
cat("    base_beta and the SAME immunity init (R0 = 1 - exp(-0.85*1)).\n")

cat("\n[done] flat_equil_check.R\n")
