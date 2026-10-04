# =============================================================================
# preui_flat_run_4chains.R
#
# 4-chain counterpart to preui_flat_run.R (chat record 2026-09-05). Builds
# the pre-vaccination-uptake simulation ("preui") per state from the 4-chain
# posterior + LHS pairing, using the SAME simulate_pre_ui_age() function as
# the 1-chain pipeline (that function already derives n_draws dynamically
# from posterior_idx_list/lhs_idx_list -- length(posterior_idx) -- no
# hardcoded 1000 inside it; see age_struc_fitting_region_func_updated.R).
#
# observed_<code> (per-state weekly case totals) does NOT depend on the fit
# at all -- create_summary_df()'s $observed is built purely from bra_sum_<code>
# $cases (raw case data), never from the stanfit/posterior. bra_sum_ce etc.
# are already in bra_cases_2022_cleaned.RData (same file the 4-chain refit
# itself loads), so they're reconstructed directly here instead of calling
# create_summary_df() (which needs the raw stanfit's age_stratified_cases --
# deliberately not saved for the 4-chain run; see age_struc_fitting_region_
# finite_2022_4chains.R's header for why).
#
# Does NOT touch preui_flat_run.R or either of its outputs (00_Data/0_2_
# Processed/preui_flat_all.RData, CHIK_benefit_risk/01_Data/preui_flat_all.
# RData) -- the 1-chain production pipeline is unaffected.
#
# Output (*_4chains suffix, both copies matching preui_flat_run.R's
# convention of saving into both repos):
#   00_Data/0_2_Processed/preui_flat_all_4chains.RData
#   C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/preui_flat_all_4chains.RData
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
source("01_Script/1_1_Functions/library.R")

message("Loading bra_foi_states.RData ...")
load("00_Data/0_2_Processed/bra_foi_states.RData")  # -> bra_foi_states

message("Sourcing age_struc_fitting_region_func_updated.R (functions + regions + age_groups + bra_foi_state_summ) ...")
source("01_Script/1_1_Functions/age_struc_fitting_region_func_updated.R")

message("Sourcing sim_functions_final.R (sirv_sim_coverageSwitch(), used by simulate_pre_ui_age()) ...")
source("01_Script/1_1_Functions/sim_functions_final.R")

message("Loading posterior_finite_all_4chains.RData, lhs_combined_finite_4chains.RData ...")
load("00_Data/0_2_Processed/posterior_finite_all_4chains.RData")
load("00_Data/0_2_Processed/lhs_combined_finite_4chains.RData")
# -> lhs_combined_finite, lhs_idx_list_finite, posterior_idx_list_finite

# simulate_pre_ui_age() looks these two up as GLOBALS (not function args) --
# same requirement as model_0803_finite_history.R's setup.
posterior_idx_list <- posterior_idx_list_finite
lhs_idx_list        <- lhs_idx_list_finite

message("Loading bra_cases_2022_cleaned.RData, bra_pop_2022_cleaned.RData (case totals + population, unchanged by refit) ...")
load("00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")   # -> bra_sum_ce, bra_sum_bh, ...
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")     # -> N_ceara, N_bahia, ...

make_observed <- function(bra_sum) {
  data.frame(Week = seq_len(nrow(bra_sum)), Observed = bra_sum$cases, Type = "Observed")
}
observed_ce <- make_observed(bra_sum_ce)
observed_bh <- make_observed(bra_sum_bh)
observed_pa <- make_observed(bra_sum_pa)
observed_pn <- make_observed(bra_sum_pn)
observed_rg <- make_observed(bra_sum_rg)
observed_pi <- make_observed(bra_sum_pi)
observed_ag <- make_observed(bra_sum_ag)
observed_tc <- make_observed(bra_sum_tc)
observed_mg <- make_observed(bra_sum_mg)
observed_se <- make_observed(bra_sum_se)
observed_go <- make_observed(bra_sum_go)

cat(sprintf("[%s] Running simulate_pre_ui_age() per state (%d draws/state) ...\n", format(Sys.time(), "%H:%M:%S"), length(posterior_idx_list[[1]])))

preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups,
                                N = N_ceara$Ceará, region = "Ceará",
                                observed = observed_ce, lhs_sample = lhs_combined_finite$Ceará)
cat(sprintf("[%s] Done: Ceará\n", format(Sys.time(), "%H:%M:%S")))

preui_bh <- simulate_pre_ui_age(posterior = posterior_bh, bra_foi_state_summ, age_groups,
                                N = N_bahia$Bahia, region = "Bahia",
                                observed = observed_bh, lhs_sample = lhs_combined_finite$Bahia)
cat(sprintf("[%s] Done: Bahia\n", format(Sys.time(), "%H:%M:%S")))

preui_pa <- simulate_pre_ui_age(posterior = posterior_pa, bra_foi_state_summ, age_groups,
                                N = N_pa$Paraíba, region = "Paraíba",
                                observed = observed_pa, lhs_sample = lhs_combined_finite$Paraíba)
cat(sprintf("[%s] Done: Paraíba\n", format(Sys.time(), "%H:%M:%S")))

preui_pn <- simulate_pre_ui_age(posterior = posterior_pn, bra_foi_state_summ, age_groups,
                                N = N_pemam$Pernambuco, region = "Pernambuco",
                                observed = observed_pn, lhs_sample = lhs_combined_finite$Pernambuco)
cat(sprintf("[%s] Done: Pernambuco\n", format(Sys.time(), "%H:%M:%S")))

preui_rg <- simulate_pre_ui_age(posterior = posterior_rg, bra_foi_state_summ, age_groups,
                                N = N_rg$`Rio Grande do Norte`, region = "Rio Grande do Norte",
                                observed = observed_rg, lhs_sample = lhs_combined_finite$`Rio Grande do Norte`)
cat(sprintf("[%s] Done: Rio Grande do Norte\n", format(Sys.time(), "%H:%M:%S")))

preui_pi <- simulate_pre_ui_age(posterior = posterior_pi, bra_foi_state_summ, age_groups,
                                N = N_pi$Piauí, region = "Piauí",
                                observed = observed_pi, lhs_sample = lhs_combined_finite$Piauí)
cat(sprintf("[%s] Done: Piauí\n", format(Sys.time(), "%H:%M:%S")))

preui_ag <- simulate_pre_ui_age(posterior = posterior_ag, bra_foi_state_summ, age_groups,
                                N = N_ag$Alagoas, region = "Alagoas",
                                observed = observed_ag, lhs_sample = lhs_combined_finite$Alagoas)
cat(sprintf("[%s] Done: Alagoas\n", format(Sys.time(), "%H:%M:%S")))

preui_tc <- simulate_pre_ui_age(posterior = posterior_tc, bra_foi_state_summ, age_groups,
                                N = N_tc$Tocantins, region = "Tocantins",
                                observed = observed_tc, lhs_sample = lhs_combined_finite$Tocantins)
cat(sprintf("[%s] Done: Tocantins\n", format(Sys.time(), "%H:%M:%S")))

preui_mg <- simulate_pre_ui_age(posterior = posterior_mg, bra_foi_state_summ, age_groups,
                                N = N_mg$`Minas Gerais`, region = "Minas Gerais",
                                observed = observed_mg, lhs_sample = lhs_combined_finite$`Minas Gerais`)
cat(sprintf("[%s] Done: Minas Gerais\n", format(Sys.time(), "%H:%M:%S")))

preui_se <- simulate_pre_ui_age(posterior = posterior_se, bra_foi_state_summ, age_groups,
                                N = N_se$Sergipe, region = "Sergipe",
                                observed = observed_se, lhs_sample = lhs_combined_finite$Sergipe)
cat(sprintf("[%s] Done: Sergipe\n", format(Sys.time(), "%H:%M:%S")))

preui_go <- simulate_pre_ui_age(posterior = posterior_go, bra_foi_state_summ, age_groups,
                                N = N_go$Goiás, region = "Goiás",
                                observed = observed_go, lhs_sample = lhs_combined_finite$Goiás)
cat(sprintf("[%s] Done: Goiás\n", format(Sys.time(), "%H:%M:%S")))

save(
  preui_ce, preui_bh, preui_pa, preui_pn, preui_rg, preui_pi,
  preui_ag, preui_tc, preui_mg, preui_se, preui_go,
  file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/preui_flat_all_4chains.RData"
)
cat(sprintf("[%s] Saved: C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/preui_flat_all_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))

save(
  preui_ce, preui_bh, preui_pa, preui_pn, preui_rg, preui_pi,
  preui_ag, preui_tc, preui_mg, preui_se, preui_go,
  file = "00_Data/0_2_Processed/preui_flat_all_4chains.RData"
)
cat(sprintf("[%s] Saved: 00_Data/0_2_Processed/preui_flat_all_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))
