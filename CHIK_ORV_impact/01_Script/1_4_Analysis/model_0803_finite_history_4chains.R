# =============================================================================
# model_0803_finite_history_4chains.R
#
# 4-chain counterpart to model_0803_finite_history.R (chat record 2026-09-05).
# Same pipeline (prevacc simulation -> post-vaccination scenario simulation
# -> NNV), pointed at the 4-chain posterior/LHS/preui instead of the 1-chain
# ones. Differences from the 1-chain script, all deliberate:
#
#   1) Does NOT reload fits_prevacc_finite.RData + re-run extract_params()
#      to rebuild posterior_ce/etc. -- that step re-derives the exact same
#      posterior_finite_all.RData that's already sitting on disk from the
#      fitting script; for the 4-chain run there IS no saved raw stanfit to
#      re-extract from (deliberately -- see age_struc_fitting_region_finite_
#      2022_4chains.R's header), so posterior_finite_all_4chains.RData is
#      loaded directly instead.
#
#   2) finite_sero_list (1-chain script lines 129-141) and rho_df_finite.RData
#      (line 12) are both dead code in the original -- finite_sero_list is
#      built but never passed to make_prevacc() (sero_vec stays NULL at
#      every call site), and rho_df is loaded but never referenced again.
#      Neither is reproduced here.
#
#   3) The "debug test" block (1-chain script lines 241-259) that overwrites
#      regions/regions_base down to Ceará-only before the post-vaccination
#      step is NOT reproduced (chat record 2026-09-05: confirmed with user
#      to run all 11 states, matching what postsim_vc_ixchiq_model_finite.
#      RData / finite_weeksweep_brr.R actually consume).
#
# VE/coverage grid kept IDENTICAL to the 1-chain script (ve_inf_set = c(0,
# 0.989), coverage_set = 0.5) -- only the posterior/LHS/preui inputs change.
#
# Output (*_4chains suffix -- does not overwrite the 1-chain production
# files in either repo):
#   00_Data/0_2_Processed/prevacc_sets_finite_4chains.RData
#   00_Data/0_2_Processed/sim_results_vc_ixchiq_model_finite_4chains.RData
#   00_Data/0_2_Processed/postsim_vc_ixchiq_model_finite_4chains.RData
#   C:/.../CHIK_benefit_risk/01_Data/sim_results_vc_ixchiq_model_finite_4chains.RData
#   C:/.../CHIK_benefit_risk/01_Data/postsim_vc_ixchiq_model_finite_4chains.RData
#   00_Data/0_2_Processed/combined_nnv_national_age_ixchiq_4chains.RData
#   C:/.../CHIK_benefit_risk/01_Data/combined_nnv_df_region_coverage_model_finite_4chains.RData
#   02_Outputs/2_2_Tables/combined_nnv_national_age_ixchiq_4chains.xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
source("01_Script/1_1_Functions/library.R")

message("Loading bra_foi_states.RData ...")
load("00_Data/0_2_Processed/bra_foi_states.RData")  # -> bra_foi_states

message("Sourcing age_struc_fitting_region_func_updated.R + sim_functions_final.R ...")
source("01_Script/1_1_Functions/age_struc_fitting_region_func_updated.R")
source("01_Script/1_1_Functions/sim_functions_final.R")

message("Loading 4-chain posterior + LHS pairing ...")
load("00_Data/0_2_Processed/posterior_finite_all_4chains.RData")
load("00_Data/0_2_Processed/lhs_combined_finite_4chains.RData")
# -> lhs_combined_finite, lhs_idx_list_finite, posterior_idx_list_finite

message("Loading disability/mortality/life-expectancy inputs (unchanged by refit) ...")
lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")
le_sample <- readRDS("00_Data/0_2_Processed/le_sample.RDS")
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")  # -> hosp, fatal, nh_fatal

message("Loading case totals + population (unchanged by refit) ...")
load("00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")   # -> bra_sum_ce, ...
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")     # -> N_ceara, ...

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

# scenarios
target_age_list <- list(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 1-11 years
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 12-17 years
                        c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0), # 18-64
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1) # 65+
)

regions <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 posterior = posterior_ce, lhs_sample = lhs_combined_finite$Ceará),
  "Bahia" = list(observed = observed_bh, N = N_bahia$Bahia,
                 posterior = posterior_bh, lhs_sample = lhs_combined_finite$Bahia),
  "Paraíba" = list(observed = observed_pa, N = N_pa$Paraíba,
                   posterior = posterior_pa, lhs_sample = lhs_combined_finite$Paraíba),
  "Pernambuco" = list(observed = observed_pn, N = N_pemam$Pernambuco,
                      posterior = posterior_pn, lhs_sample = lhs_combined_finite$Pernambuco),
  "Rio Grande do Norte" = list(observed = observed_rg, N = N_rg$`Rio Grande do Norte`,
                               posterior = posterior_rg, lhs_sample = lhs_combined_finite$`Rio Grande do Norte`),
  "Piauí" = list(observed = observed_pi, N = N_pi$Piauí,
                 posterior = posterior_pi, lhs_sample = lhs_combined_finite$Piauí),
  "Tocantins" = list(observed = observed_tc, N = N_tc$Tocantins,
                     posterior = posterior_tc, lhs_sample = lhs_combined_finite$Tocantins),
  "Alagoas" = list(observed = observed_ag, N = N_ag$Alagoas,
                   posterior = posterior_ag, lhs_sample = lhs_combined_finite$Alagoas),
  "Minas Gerais" = list(observed = observed_mg, N = N_mg$`Minas Gerais`,
                        posterior = posterior_mg, lhs_sample = lhs_combined_finite$`Minas Gerais`),
  "Sergipe" = list(observed = observed_se, N = N_se$Sergipe,
                   posterior = posterior_se, lhs_sample = lhs_combined_finite$Sergipe),
  "Goiás" = list(observed = observed_go, N = N_go$Goiás,
                 posterior = posterior_go, lhs_sample = lhs_combined_finite$Goiás)
)

# simulate_pre_ui_age()/run_simulation_scenarios_ui_ixchiq() look up these
# two as globals (not function arguments).
posterior_idx_list <- posterior_idx_list_finite
lhs_idx_list        <- lhs_idx_list_finite

make_prevacc <- function(
    posterior, observed, N, region, lhs_sample, sero_vec = NULL,
    finite_history = TRUE
) {
  pre_ui <- simulate_pre_ui_age(
    posterior = posterior,
    bra_foi_state_summ = bra_foi_state_summ,
    age_groups = age_groups,
    N = N,
    region = region,
    observed = observed,
    lhs_sample = lhs_sample,
    sero_vec = sero_vec,
    finite_history = finite_history
  )

  presum <- summarise_presim_ui(
    sim_result = pre_ui,
    observed = observed,
    age_gr_levels = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old = lhs_old,
    le_sample = le_sample,
    hosp = hosp,
    fatal = fatal,
    nh_fatal = nh_fatal,
    region = region
  )

  list(
    prevacc_ui = presum$summary_cases_pre,
    prevacc_ui_all = presum$summary_cases_pre_all,
    pre_summary_age = presum$summary_cases_pre_age
  )
}

cat(sprintf("[%s] Running make_prevacc() (simulate_pre_ui_age + summarise_presim_ui) per state, all 11 states ...\n", format(Sys.time(), "%H:%M:%S")))
prevacc_sets_finite <- purrr::imap(
  regions,
  function(reg_args, region_name) {
    result <- make_prevacc(
      posterior = reg_args$posterior,
      observed = reg_args$observed,
      N = reg_args$N,
      region = region_name,
      lhs_sample = reg_args$lhs_sample,
      finite_history = TRUE
    )
    cat(sprintf("[%s] Done prevacc: %s\n", format(Sys.time(), "%H:%M:%S"), region_name))
    result
  }
)

save("prevacc_sets_finite", file = "00_Data/0_2_Processed/prevacc_sets_finite_4chains.RData")
cat(sprintf("[%s] Saved: 00_Data/0_2_Processed/prevacc_sets_finite_4chains.RData\n", format(Sys.time(), "%H:%M:%S")))

## 1) Post vacc simulation --------------------------------------------------------
regions <- purrr::imap(
  regions,
  ~ c(.x, list(prevacc_ui = prevacc_sets_finite[[.y]]$prevacc_ui))
)

ve_inf_set     <- c(0, 0.989)
coverage_set   <- 0.5

cat(sprintf("[%s] Running run_simulation_scenarios_ui_ixchiq() per state x VE x coverage, all 11 states ...\n", format(Sys.time(), "%H:%M:%S")))
sim_results_vc_ixchiq_model  <- imap(regions, function(reg_args, reg_name) {

  n_draws <- length(reg_args$posterior$gamma)

  ## VE loop
  result <- imap(set_names(ve_inf_set, paste0("VE", ve_inf_set*100)), function(ve_val, ve_tag) {

    ## coverage loop
    imap(set_names(coverage_set, paste0("cov", coverage_set*100)),
         function(cov_val, cov_tag) {

           ve_input <- if (abs(ve_val - 0.989) < 1e-6) NA else ve_val

           run_simulation_scenarios_ui_ixchiq(
             target_age_list    = target_age_list,
             observed           = reg_args$observed,
             N                  = reg_args$N,
             bra_foi_state_summ = bra_foi_state_summ,
             age_groups         = age_groups,
             region_name        = reg_name,
             hosp               = hosp,
             fatal              = fatal,
             nh_fatal           = nh_fatal,
             lhs_sample_young   = lhs_sample_young,
             lhs_old            = lhs_old,
             le_sample          = le_sample,
             age_gr_levels      = age_gr_levels,
             prevacc_ui         = reg_args$prevacc_ui,
             posterior          = reg_args$posterior,
             ve_inf             = ve_input,
             total_coverage     = cov_val,
             lhs_sample         = reg_args$lhs_sample,
             finite_history     = TRUE
           )
         })
  })
  cat(sprintf("[%s] Done post-vacc sim: %s\n", format(Sys.time(), "%H:%M:%S"), reg_name))
  result
})

postsim_vc_ixchiq_model  <- imap(regions, function(reg_args, region_name) {

  ve_cov_list <- sim_results_vc_ixchiq_model[[region_name]]
  prevacc      <- prevacc_sets_finite[[region_name]]

  imap(ve_cov_list, function(cov_list, ve_tag) {
    imap(cov_list, function(scenario_result, cov_tag) {
      postsim_all_ui(
        scenario_result       = scenario_result,
        observed              = reg_args$observed,
        age_gr_levels         = age_gr_levels,
        pre_summary_cases_age = prevacc$pre_summary_age,
        pre_summary_cases     = prevacc$prevacc_ui,
        pre_summary_cases_all = prevacc$prevacc_ui_all,
        region                = region_name
      )
    }) %>% set_names(names(cov_list))
  }) %>% set_names(names(ve_cov_list))
})

save("sim_results_vc_ixchiq_model", file = "00_Data/0_2_Processed/sim_results_vc_ixchiq_model_finite_4chains.RData")
save("postsim_vc_ixchiq_model", file = "00_Data/0_2_Processed/postsim_vc_ixchiq_model_finite_4chains.RData")

save(sim_results_vc_ixchiq_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/sim_results_vc_ixchiq_model_finite_4chains.RData")
save(postsim_vc_ixchiq_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/postsim_vc_ixchiq_model_finite_4chains.RData")
cat(sprintf("[%s] Saved sim_results/postsim (both repos)\n", format(Sys.time(), "%H:%M:%S")))

# NNV----------------------------------------------------------------------------
n_scenarios <- 4
nnv_results_coverage_model <- imap(postsim_vc_ixchiq_model, function(ve_cov_list, region_name) {

  reg_args <- regions[[region_name]]

  imap(ve_cov_list, function(cov_list, ve_tag) {
    imap(cov_list, function(postsim_ui, cov_tag) {
      nnv_list(
        vacc_alloc     = vacc_allocation(postsim_ui, reg_args$observed, region_name),
        postsim_all_ui = postsim_ui,
        N              = reg_args$N,
        region         = region_name,
        observed       = reg_args$observed
      )
    }) |> set_names(names(cov_list))
  }) |> set_names(names(ve_cov_list))
})

setting_key <- c(
  "Ceará"             = "High",
  "Bahia"             = "Low",
  "Paraíba"           = "High",
  "Pernambuco"        = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"             = "High",
  "Tocantins"         = "Moderate",
  "Alagoas"           = "High",
  "Minas Gerais"      = "Low",
  "Sergipe"           = "Low",
  "Goiás"             = "Low"
)

age_seq <- rep(age_gr[1:length(age_gr_levels)], n_scenarios)

combined_nnv_df_region_coverage_model <- imap_dfr(nnv_results_coverage_model,
  function(ve_cov_list, region_name) {
    imap_dfr(ve_cov_list,
      function(cov_list, ve_tag) {
        imap_dfr(cov_list,
          function(nnv, cov_tag) {
            nnv$final_summ_df %>%
              mutate(
                region   = region_name,
                VE       = ve_tag,
                VC       = cov_tag,
                setting  = setting_key[region_name],
                age_gr   = factor(age_seq, levels = age_gr_levels)
              )
          }
        )
      }
    )
  }
)

combined_nnv_national_age_ixchiq <- combined_nnv_df_region_coverage_model %>%
  group_by(VE, VC, scenario, AgeGroup) %>%
  summarise(
    across(c(tot_vacc, pre_inf:diff_hosp_hi), sum, na.rm = TRUE),
    .groups = "drop"
  )

write_xlsx(combined_nnv_national_age_ixchiq, path = "02_Outputs/2_2_Tables/combined_nnv_national_age_ixchiq_4chains.xlsx")
save(combined_nnv_national_age_ixchiq, file = "00_Data/0_2_Processed/combined_nnv_national_age_ixchiq_4chains.RData")
save(combined_nnv_df_region_coverage_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/combined_nnv_df_region_coverage_model_finite_4chains.RData")

cat(sprintf("[%s] ALL DONE -- model_0803_finite_history_4chains.R complete.\n", format(Sys.time(), "%H:%M:%S")))
