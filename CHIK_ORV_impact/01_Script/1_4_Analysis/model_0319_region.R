load("00_Data/0_2_Processed/allfoi_s1.RData")
load("00_Data/0_2_Processed/bra_foi_states.RData")
load("00_Data/0_2_Processed/fit_prevacc_bh.RData")
load("00_Data/0_2_Processed/fit_prevacc_ce.RData")
load("00_Data/0_2_Processed/fit_prevacc_mg.RData")
load("00_Data/0_2_Processed/fit_prevacc_pn.RData")
load("00_Data/0_2_Processed/fit_prevacc_pa.RData")
load("00_Data/0_2_Processed/fit_prevacc_rg.RData")
load("00_Data/0_2_Processed/fit_prevacc_pi.RData")
load("00_Data/0_2_Processed/fit_prevacc_ag.RData")
load("00_Data/0_2_Processed/fit_prevacc_tc.RData")
load("00_Data/0_2_Processed/fit_prevacc_se.RData")
load("00_Data/0_2_Processed/fit_prevacc_go.RData")
load("00_Data/0_2_Processed/lhs_orv.RData")
load("00_Data/0_2_Processed/region_coverage.RData")
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")
load("00_Data/0_2_Processed/observed_2022.RData")
load("00_Data/0_2_Processed/bra_foi_state_summ.RData")
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")
load("00_Data/0_2_Processed/combined_contour.RData")
load("00_Data/0_2_Processed/combined_heatmap.RData")
load("00_Data/0_2_Processed/tornado_results_all.RData")
load("00_Data/0_2_Processed/rho_df.RData")
load("00_Data/0_2_Processed/sim_results_ve.RData")
load("00_Data/0_2_Processed/presim_all_regions.RData")
load("00_Data/0_2_Processed/prevacc_sets_ve.RData")

lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")
le_sample <- readRDS("00_Data/0_2_Processed/le_sample.RDS")
source("01_Script/1_1_Functions/sim_functions_final.R")
source("01_Script/1_1_Functions/age_struc_fitting_region_func.R")
source("01_Script/1_1_Functions/library.R")

# extract params-----------------------------------------------------------------
param_ce <- extract_params(fit_prevacc_ce)
param_pa <- extract_params(fit_prevacc_pa)
param_pn <- extract_params(fit_prevacc_pn)
param_bh <- extract_params(fit_prevacc_bh)
param_rg <- extract_params(fit_prevacc_rg)
param_pi <- extract_params(fit_prevacc_pi)
param_ag <- extract_params(fit_prevacc_ag)
param_tc <- extract_params(fit_prevacc_tc)
param_mg <- extract_params(fit_prevacc_mg)
param_se <- extract_params(fit_prevacc_se)
param_go <- extract_params(fit_prevacc_go)

posterior_ce <- param_ce$posterior_prevacc
posterior_bh <- param_bh$posterior_prevacc
posterior_pa <- param_pa$posterior_prevacc
posterior_pn <- param_pn$posterior_prevacc
posterior_rg <- param_rg$posterior_prevacc
posterior_pi <- param_pi$posterior_prevacc
posterior_ag <- param_ag$posterior_prevacc
posterior_tc <- param_tc$posterior_prevacc
posterior_mg <- param_mg$posterior_prevacc
posterior_se <- param_se$posterior_prevacc
posterior_go <- param_go$posterior_prevacc

# scenarios
target_age_list <- list(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 1-11 years
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 12-17 years
                        c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0), # 18-64
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1) # 65+
                       
                        ) 


# prevacc 95%UI simulation ----------------------------------------------------
preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups, 
                                N = N_ceara$Ceará, 
                                region = "Ceará",
                                observed = observed_ce,
                                lhs_sample = lhs_combined$Ceará
)

preui_bh <- simulate_pre_ui_age(posterior = posterior_bh, bra_foi_state_summ, age_groups, 
                                N = N_bahia$Bahia, 
                                region = "Bahia",
                                observed = observed_bh,
                                lhs_sample = lhs_combined$Bahia
)

preui_pa <- simulate_pre_ui_age(posterior = posterior_pa, bra_foi_state_summ, age_groups, 
                                N = N_pa$Paraíba, 
                                region = "Paraíba",
                                observed = observed_pa,
                                lhs_sample = lhs_combined$Paraíba
)

preui_pn <- simulate_pre_ui_age(posterior = posterior_pn, bra_foi_state_summ, age_groups, 
                                N = N_pemam$Pernambuco, 
                                region = "Pernambuco",
                                observed = observed_pn,
                                lhs_sample = lhs_combined$Pernambuco
)

preui_rg <- simulate_pre_ui_age(posterior = posterior_rg, bra_foi_state_summ, age_groups, 
                                N = N_rg$`Rio Grande do Norte`, 
                                region = "Rio Grande do Norte",
                                observed = observed_rg,
                                lhs_sample = lhs_combined$`Rio Grande do Norte`
)

preui_pi <- simulate_pre_ui_age(posterior = posterior_pi, bra_foi_state_summ, age_groups, 
                                N = N_pi$Piauí, 
                                region = "Piauí",
                                observed = observed_pi,
                                lhs_sample = lhs_combined$Piauí
)

preui_tc <- simulate_pre_ui_age(posterior = posterior_tc, bra_foi_state_summ, age_groups, 
                                N = N_tc$Tocantins, 
                                region = "Tocantins",
                                observed = observed_tc,
                                lhs_sample = lhs_combined$Tocantins
)

preui_ag <- simulate_pre_ui_age(posterior = posterior_ag, bra_foi_state_summ, age_groups, 
                                N = N_ag$Alagoas, 
                                region = "Alagoas",
                                observed = observed_ag,
                                lhs_sample = lhs_combined$Alagoas
)

preui_mg <- simulate_pre_ui_age(posterior = posterior_mg, bra_foi_state_summ, age_groups, 
                                N = N_mg$`Minas Gerais`, 
                                region = "Minas Gerais",
                                observed = observed_mg,
                                lhs_sample = lhs_combined$`Minas Gerais`
)

preui_se <- simulate_pre_ui_age(posterior = posterior_se, bra_foi_state_summ, age_groups, 
                                N = N_se$Sergipe, 
                                region = "Sergipe",
                                observed = observed_se,
                                lhs_sample = lhs_combined$Sergipe
)

preui_go <- simulate_pre_ui_age(posterior = posterior_go, bra_foi_state_summ, age_groups, 
                                N = N_go$Goiás, 
                                region = "Goiás",
                                observed = observed_go,
                                lhs_sample = lhs_combined$Goiás
)

## summarise ----------------------------------------------------
pre_results_ce_ui <- summarise_presim_ui(sim_result     = preui_ce, 
                                          observed         = observed_ce,
                                          age_gr_levels    = age_gr_levels,
                                          lhs_sample_young = lhs_sample_young,
                                          lhs_old          = lhs_old,
                                          le_sample        = le_sample,
                                          hosp             = hosp,
                                          fatal            = fatal,
                                          nh_fatal         = nh_fatal,
                                          region           = "Ceará"
)

pre_results_bh_ui <- summarise_presim_ui(sim_result       = preui_bh, 
                                         observed         = observed_bh,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Bahia"
)

pre_results_pa_ui <- summarise_presim_ui(sim_result       = preui_pa, 
                                         observed         = observed_pa,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Paraíba"
)

pre_results_pn_ui <- summarise_presim_ui(sim_result       = preui_pn, 
                                         observed         = observed_pn,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Pernambuco"
)

pre_results_rg_ui <- summarise_presim_ui(sim_result       = preui_rg, 
                                         observed         = observed_rg,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Rio Grande do Norte"
)

pre_results_pi_ui <- summarise_presim_ui(sim_result       = preui_pi, 
                                         observed         = observed_pi,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Piauí"
)

pre_results_tc_ui <- summarise_presim_ui(sim_result       = preui_tc, 
                                         observed         = observed_tc,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Tocantins"
)

pre_results_ag_ui <- summarise_presim_ui(sim_result       = preui_ag, 
                                         observed         = observed_ag,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Alagoas"
)

pre_results_mg_ui <- summarise_presim_ui(sim_result       = preui_mg, 
                                         observed         = observed_mg,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Minas Gerais"
)

pre_results_se_ui <- summarise_presim_ui(sim_result       = preui_se, 
                                         observed         = observed_se,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Sergipe"
)

pre_results_go_ui <- summarise_presim_ui(sim_result       = preui_go, 
                                         observed         = observed_go,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Goiás"
)

## prevacc outputs ----------------------------------------------------
prevacc_ui_ce <- pre_results_ce_ui$summary_cases_pre
prevacc_ui_ce_all <- pre_results_ce_ui$summary_cases_pre_all
pre_summary_age_ce <- pre_results_ce_ui$summary_cases_pre_age

prevacc_ui_bh <- pre_results_bh_ui$summary_cases_pre
prevacc_ui_bh_all <- pre_results_bh_ui$summary_cases_pre_all
pre_summary_age_bh <- pre_results_bh_ui$summary_cases_pre_age

prevacc_ui_pa <- pre_results_pa_ui$summary_cases_pre
prevacc_ui_pa_all <- pre_results_pa_ui$summary_cases_pre_all
pre_summary_age_pa <- pre_results_pa_ui$summary_cases_pre_age

prevacc_ui_pn <- pre_results_pn_ui$summary_cases_pre
prevacc_ui_pn_all <- pre_results_pn_ui$summary_cases_pre_all
pre_summary_age_pn <- pre_results_pn_ui$summary_cases_pre_age

prevacc_ui_rg <- pre_results_rg_ui$summary_cases_pre
prevacc_ui_rg_all <- pre_results_rg_ui$summary_cases_pre_all
pre_summary_age_rg <- pre_results_rg_ui$summary_cases_pre_age

prevacc_ui_pi <- pre_results_pi_ui$summary_cases_pre
prevacc_ui_pi_all <- pre_results_pi_ui$summary_cases_pre_all
pre_summary_age_pi <- pre_results_pi_ui$summary_cases_pre_age

prevacc_ui_tc <- pre_results_tc_ui$summary_cases_pre
prevacc_ui_tc_all <- pre_results_tc_ui$summary_cases_pre_all
pre_summary_age_tc <- pre_results_tc_ui$summary_cases_pre_age

prevacc_ui_ag <- pre_results_ag_ui$summary_cases_pre
prevacc_ui_ag_all <- pre_results_ag_ui$summary_cases_pre_all
pre_summary_age_ag <- pre_results_ag_ui$summary_cases_pre_age

prevacc_ui_mg <- pre_results_mg_ui$summary_cases_pre
prevacc_ui_mg_all <- pre_results_mg_ui$summary_cases_pre_all
pre_summary_age_mg <- pre_results_mg_ui$summary_cases_pre_age

prevacc_ui_se <- pre_results_se_ui$summary_cases_pre
prevacc_ui_se_all <- pre_results_se_ui$summary_cases_pre_all
pre_summary_age_se <- pre_results_se_ui$summary_cases_pre_age

prevacc_ui_go <- pre_results_go_ui$summary_cases_pre
prevacc_ui_go_all <- pre_results_go_ui$summary_cases_pre_all
pre_summary_age_go <- pre_results_go_ui$summary_cases_pre_age

# postsim 95% UI run ----------------------------------------------------------
postsim_ce_ui <- run_simulation_scenarios_ui_ixchiq(
  target_age_list   = target_age_list,
  #supply            = supply_ce,
  observed          = observed_ce,
  N                 = N_ceara$Ceará,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Ceará",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_ce,
  posterior         = posterior_ce,
  lhs_sample        = lhs_combined$Ceará
)

postsim_bh_ui <- run_simulation_scenarios_ui_vimkun(
  target_age_list   = target_age_list,
  #supply            = supply_bh,
  observed           = observed_bh,
  N                 = N_bahia$Bahia,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Bahia",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_bh,
  posterior         = posterior_bh,
  lhs_sample        = lhs_combined$Bahia
)

postsim_pa_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_pa,
  observed           = observed_pa,
  N                 = N_pa$Paraíba,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Paraíba",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_pa,
  posterior         = posterior_pa,
  lhs_sample        = lhs_combined$Paraíba
)

postsim_pn_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_pn,
  observed          = observed_pn,
  N                 = N_pemam$Pernambuco,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Pernambuco",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_pn,
  posterior         = posterior_pn,
  lhs_sample        = lhs_combined$Pernambuco
)

postsim_rg_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_rg,
  observed           = observed_rg,
  N                 = N_rg$`Rio Grande do Norte`,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Rio Grande do Norte",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_rg,
  posterior         = posterior_rg,
  lhs_sample        = lhs_combined$`Rio Grande do Norte`
)

postsim_pi_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_pi,
  observed          = observed_pi,
  N                 = N_pi$Piauí,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Piauí",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_pi,
  posterior         = posterior_pi,
  lhs_sample        = lhs_combined$Piauí
)

postsim_tc_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_tc,
  observed          = observed_tc,
  N                 = N_tc$Tocantins,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Tocantins",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_tc,
  posterior         = posterior_tc,
  lhs_sample        = lhs_combined$Tocantins
)

postsim_ag_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_ag,
  observed          = observed_ag,
  N                 = N_ag$Alagoas,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Alagoas",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_ag,
  posterior         = posterior_ag,
  lhs_sample        = lhs_combined$Alagoas
)

postsim_mg_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_mg,
  N                 = N_mg$`Minas Gerais`,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Minas Gerais",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_mg,
  posterior         = posterior_mg,
  lhs_sample        = lhs_combined$`Minas Gerais`
)

postsim_se_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_se,
  N                 = N_se$Sergipe,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Sergipe",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_se,
  posterior         = posterior_se,
  lhs_sample        = lhs_combined$Sergipe
)

postsim_go_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_go,
  N                 = N_go$Goiás,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Goiás",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_go,
  posterior         = posterior_go,
  lhs_sample        = lhs_combined$Goiás
)

# save pre-postsim data-------------------
save(
  postsim_ce_ui, postsim_bh_ui, postsim_pa_ui, postsim_pn_ui,
  postsim_rg_ui, postsim_pi_ui, postsim_tc_ui, postsim_ag_ui,
  postsim_mg_ui, postsim_se_ui, postsim_go_ui,
  file = "00_Data/0_2_Processed/postsim_all_regions_ixchiqcov.RData"
)

save(
  preui_ce, preui_bh, preui_pa, preui_pn,
  preui_rg, preui_pi, preui_tc, preui_ag,
  preui_mg, preui_se, preui_go,
  file = "00_Data/0_2_Processed/presim_all_regions.RData"
)

# summarise for post sim ------------------------------------------------------
postsim_all_ce_ui <- postsim_all_ui(
  scenario_result      = postsim_ce_ui,
  observed             = observed_ce,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_ce,
  pre_summary_cases     = prevacc_ui_ce,
  pre_summary_cases_all     = prevacc_ui_ce_all,
  region                = "Ceará"
)

postsim_all_bh_ui <- postsim_all_ui(
  scenario_result      = postsim_bh_ui,
  observed             = observed_bh,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_bh,
  pre_summary_cases     = prevacc_ui_bh,
  pre_summary_cases_all     = prevacc_ui_bh_all,
  region                = "Bahia"
)

postsim_all_pa_ui <- postsim_all_ui(
  scenario_result      = postsim_pa_ui,
  observed             = observed_pa,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_pa,
  pre_summary_cases     = prevacc_ui_pa,
  pre_summary_cases_all     = prevacc_ui_pa_all,
  region                = "Paraíba"
)

postsim_all_pn_ui <- postsim_all_ui(
  scenario_result      = postsim_pn_ui,
  observed             = observed_pn,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_pn,
  pre_summary_cases     = prevacc_ui_pn,
  pre_summary_cases_all     = prevacc_ui_pn_all,
  region                = "Pernambuco"
)

postsim_all_rg_ui <- postsim_all_ui(
  scenario_result      = postsim_rg_ui,
  observed             = observed_rg,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_rg,
  pre_summary_cases     = prevacc_ui_rg,
  pre_summary_cases_all     = prevacc_ui_rg_all,
  region                = "Rio Grande do Norte"
)

postsim_all_pi_ui <- postsim_all_ui(
  scenario_result      = postsim_pi_ui,
  observed             = observed_pi,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_pi,
  pre_summary_cases     = prevacc_ui_pi,
  pre_summary_cases_all     = prevacc_ui_pi_all,
  region                = "Piauí"
)

postsim_all_tc_ui <- postsim_all_ui(
  scenario_result      = postsim_tc_ui,
  observed             = observed_tc,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_tc,
  pre_summary_cases     = prevacc_ui_tc,
  pre_summary_cases_all     = prevacc_ui_tc_all,
  region                = "Tocantins"
)

postsim_all_ag_ui <- postsim_all_ui(
  scenario_result      = postsim_ag_ui,
  observed             = observed_ag,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_ag,
  pre_summary_cases     = prevacc_ui_ag,
  pre_summary_cases_all     = prevacc_ui_ag_all,
  region                = "Alagoas"
)

postsim_all_mg_ui <- postsim_all_ui(
  scenario_result      = postsim_mg_ui,
  observed             = observed_mg,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_mg,
  pre_summary_cases     = prevacc_ui_mg,
  pre_summary_cases_all     = prevacc_ui_mg_all,
  region                = "Minas Gerais"
)

postsim_all_se_ui <- postsim_all_ui(
  scenario_result      = postsim_se_ui,
  observed             = observed_se,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_se,
  pre_summary_cases     = prevacc_ui_se,
  pre_summary_cases_all     = prevacc_ui_se_all,
  region                = "Sergipe"
)

postsim_all_go_ui <- postsim_all_ui(
  scenario_result      = postsim_go_ui,
  observed             = observed_go,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_go,
  pre_summary_cases     = prevacc_ui_go,
  pre_summary_cases_all     = prevacc_ui_go_all,
  region                = "Goiás"
)

# pre-post graphs -------------------------------------------------------------
p1 <- epi_graph(postsim_all_ce_ui,
                observed_ce)

p2 <- epi_graph(postsim_all_bh_ui,
                observed_bh)

p3 <- epi_graph(postsim_all_pa_ui,
                observed_pa)

p4 <- epi_graph(postsim_all_pn_ui,
                observed_pn)

p5 <- epi_graph(postsim_all_rg_ui,
                observed_rg)

p6 <- epi_graph(postsim_all_pi_ui,
                observed_pi)

p7 <- epi_graph(postsim_all_tc_ui,
                observed_tc)

p8 <- epi_graph(postsim_all_ag_ui,
                observed_ag)

p9 <- epi_graph(postsim_all_mg_ui,
                observed_mg)

p10 <- epi_graph(postsim_all_se_ui,
                observed_se)

p11 <- epi_graph(postsim_all_go_ui,
                observed_go)

combined_prepost_case <- bind_rows(
  postsim_all_ce_ui$summary_week_df,
  postsim_all_bh_ui$summary_week_df,
  postsim_all_pa_ui$summary_week_df,
  postsim_all_pn_ui$summary_week_df,
  postsim_all_rg_ui$summary_week_df,
  postsim_all_pi_ui$summary_week_df,
  postsim_all_tc_ui$summary_week_df,
  postsim_all_ag_ui$summary_week_df,
  postsim_all_mg_ui$summary_week_df,
  postsim_all_se_ui$summary_week_df,
  postsim_all_go_ui$summary_week_df
)

global_impact_all <- combined_prepost_case %>%
  group_by(Scenario, region) %>%
  summarise(
    total_post_cases = sum(post_cases),
    total_pre_case   = sum(pre_cases),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_case - total_post_cases,
    impact = diff / total_pre_case * 100
  )

annotation_df <- global_impact_all %>%
  mutate(Strategy = dplyr::recode(Scenario,
                           "Scenario_1" = "Strategy 1",
                           "Scenario_2" = "Strategy 2",
                           "Scenario_3" = "Strategy 3"))  %>%
  group_by(region) %>%
  summarise(ann = paste0(Strategy, ": ", round(impact, 1), "%", collapse = "\n"),
            .groups = "drop")

overall_max <- max(combined_prepost_case$hi95, na.rm = TRUE)

# Create an annotation data frame that has one row per region
vaccine_ann_df <- combined_prepost_case %>%
  distinct(region) %>%
  mutate(x = 4, 
         y = overall_max * 0.95, 
         label = "<----- Vaccine impact start (+ 2 wks after initiation)")

max_cases <- max(combined_prepost_case$hi95, na.rm = TRUE)
vacc_start_week_s1 <- postsim_all_ce_ui$vacc_weeks$scenario1$start
vacc_end_week_s1   <- postsim_all_ce_ui$vacc_weeks$scenario1$end
vacc_start_week_s2 <- postsim_all_ce_ui$vacc_weeks$scenario2$start
vacc_end_week_s2   <- postsim_all_ce_ui$vacc_weeks$scenario2$end
vacc_start_week_s3 <- postsim_all_ce_ui$vacc_weeks$scenario3$start
vacc_end_week_s3   <- postsim_all_ce_ui$vacc_weeks$scenario3$end

pre_post_graph <- 
ggplot(combined_prepost_case) +
  # Scenario shading ribbons (existing)
  geom_ribbon(data = combined_prepost_case %>% filter(Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 0, ymax = Inf), fill = "grey70",
              alpha = 0.3) +
  #geom_ribbon(data = combined_prepost_case %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
   #           aes(x = Week, ymin = 1/3 * max_cases, ymax = 2/3 * max_cases, fill = Scenario),
  #            alpha = 0.4) +
  #geom_ribbon(data = combined_prepost_case %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
   #           aes(x = Week, ymin = 0, ymax = 1/3 * max_cases, fill = Scenario),
   #           alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
              alpha = 0.3) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = 0.6) +
  # Post-vaccination case line (colored by scenario)
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "With vaccination"), size = 0.2) +
  # Pre-vaccination case line (dashed, in black)
  geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Without vaccination"), color = "black", size = 0.3) +
  scale_linetype_manual(values = c("With vaccination" = "solid", "Without vaccination" = "dashed")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Scenario_1" = "1-11 years only", 
                               "Scenario_2" = "12-17 years only", 
                               "Scenario_3" = "18-59 years only",
                               "Scenario_4" = "60+ years only")) +
  scale_color_brewer(palette = "Set1",
                     labels = c("Scenario_1" = "1-11 years only", 
                                "Scenario_2" = "12-17 years only", 
                                "Scenario_3" = "18-59 years only",
                                "Scenario_4" = "60+ years only")) +
  scale_y_continuous(labels = scales::comma) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       #title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
       x = "Week", y = "Predicted symptomatic cases") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(5, 30, 5, 5)
  ) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = Inf, label = ann), 
            inherit.aes = FALSE,
            hjust = 1.1, 
            vjust = 1.1, 
            size = 3) +
  #geom_text(data = vaccine_ann_df,
  #          aes(x = x, y = y, label = label),
  #          inherit.aes = FALSE,
  #          hjust = 0, vjust = 1, size = 2)+
  #coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  geom_point(data = observed_all, aes(x = Week, y = Observed),size = 0.05) + 
  facet_wrap(~region, scales = "free_y")

ggsave(filename = "02_Outputs/2_1_Figures/pre_post_region.jpg", p, width = 15, height = 7, dpi = 1200)

#### 
n_scenarios = 4
# vacc allocations -------------------------------------------------------------
vacc_alloc_ce <- vacc_allocation(postsim_all_ce_ui, observed_ce, "Ceará")

vacc_alloc_bh <- vacc_allocation(postsim_all_bh_ui, observed_bh, "Bahia")

vacc_alloc_pa <- vacc_allocation(postsim_all_pa_ui, observed_pa, "Paraíba")

vacc_alloc_pn <- vacc_allocation(postsim_all_pn_ui, observed_pn, "Pernambuco")

vacc_alloc_rg <- vacc_allocation(postsim_all_rg_ui, observed_rg, "Rio Grande do Norte")

vacc_alloc_pi <- vacc_allocation(postsim_all_pi_ui, observed_pi, "Piauí")

vacc_alloc_tc <- vacc_allocation(postsim_all_tc_ui, observed_tc, "Tocantins" )

vacc_alloc_ag <- vacc_allocation(postsim_all_ag_ui, observed_ag, "Alagoas")

vacc_alloc_mg <- vacc_allocation(postsim_all_mg_ui, observed_mg, "Minas Gerais" )

vacc_alloc_se <- vacc_allocation(postsim_all_se_ui, observed_se, "Sergipe" )

vacc_alloc_go <- vacc_allocation(postsim_all_go_ui, observed_go, "Goiás" )


# vacc allocations graph -------------------------------------------------------------
vacc_alloc_graph(vacc_alloc_ce)
vacc_alloc_graph(vacc_alloc_bh)
vacc_alloc_graph(vacc_alloc_pa)
vacc_alloc_graph(vacc_alloc_pn)
vacc_alloc_graph(vacc_alloc_rg)
vacc_alloc_graph(vacc_alloc_pi)
vacc_alloc_graph(vacc_alloc_tc)
vacc_alloc_graph(vacc_alloc_ag)
vacc_alloc_graph(vacc_alloc_mg)
vacc_alloc_graph(vacc_alloc_se)
vacc_alloc_graph(vacc_alloc_go)

combined_vacc_alloc <- bind_rows(
  vacc_alloc_ce$weekly_allocation_long,
  vacc_alloc_bh$weekly_allocation_long,
  vacc_alloc_pa$weekly_allocation_long,
  vacc_alloc_pn$weekly_allocation_long,
  vacc_alloc_rg$weekly_allocation_long,
  vacc_alloc_pi$weekly_allocation_long,
  vacc_alloc_tc$weekly_allocation_long,
  vacc_alloc_ag$weekly_allocation_long,
  vacc_alloc_mg$weekly_allocation_long,
  vacc_alloc_se$weekly_allocation_long,
  vacc_alloc_go$weekly_allocation_long
)

combined_vacc_alloc_region <- combined_vacc_alloc %>% 
  group_by(Scenario, region) %>%
  summarise(
    Vaccinated = sum(Vaccinated, na.rm = TRUE)
  ) %>%
  ungroup()

colnames(combined_vacc_alloc_region) <- c("Scenario", "Region", "tot_vacc")
combined_vacc_alloc_region <- combined_vacc_alloc_region %>%
  mutate(
    Scenario = case_when(
      str_detect(Scenario, "^Scenario_") ~ Scenario,
      str_detect(Scenario, "^Scenario\\d+$") ~ str_replace(Scenario, "Scenario(\\d+)", "Scenario_\\1"),
      TRUE ~ Scenario
    )
  )

combined_vacc_alloc_summ <- combined_vacc_alloc %>%
  group_by(Scenario, AgeGroup, age_gr, Week) %>%
  summarise(
    Vaccinated = sum(Vaccinated, na.rm = TRUE),
    tot_sum = sum(tot_sum, na.rm = TRUE)
  ) %>%
  ungroup()

paired <- brewer.pal(12, "Paired")
extra <- brewer.pal(8, "Dark2")
full_palette <- c(paired, extra)

p <- 
ggplot(combined_vacc_alloc_summ, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group (millions)",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=scales::label_number(scale_cut = scales::cut_short_scale()),
                     breaks = seq(0, 3000000, by = 500000))+
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c("Scenario1" = "<20 years only", 
                                              "Scenario2" = "20-59 years only", 
                                              "Scenario3" = ">60 years only")))  +         # Facet by Scenario if you want separate plots for each
  theme_minimal()

p1 <- 
ggplot(combined_vacc_alloc, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=comma)+
  scale_y_continuous(label=scales::label_number(scale_cut = scales::cut_short_scale()),
                     breaks = seq(0, 3000000, by = 500000))+
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c("Scenario1" = "<20 years only", 
                                              "Scenario2" = "20-59 years only", 
                                              "Scenario3" = ">60 years only")))  +         # Facet by Scenario if you want separate plots for each
  theme_pubclean()+
  theme(
    legend.position = "right",         
    legend.direction = "vertical",    
    legend.box = "vertical",          
    legend.title = element_text(face = "bold")
  )+
  facet_wrap(~region)


ggsave(filename = "02_Outputs/2_1_Figures/figs7.jpg", p, width = 12, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/fig_vacc_alloc_region.jpg", p, width = 11, height = 8, dpi = 1200)

# nnv  ------------------------------------------------------------------------
nnv_ce <- nnv_list(vacc_alloc_ce, postsim_all_ce_ui, N_ceara$Ceará, "Ceará", observed_ce)

nnv_bh <- nnv_list(vacc_alloc_bh, postsim_all_bh_ui, N_bahia$Bahia, "Bahia", observed_bh)

nnv_pa <- nnv_list(vacc_alloc_pa, postsim_all_pa_ui, N_pa$Paraíba, "Paraíba", observed_pa)

nnv_pn <- nnv_list(vacc_alloc_pn, postsim_all_pn_ui, N_pemam$Pernambuco, "Pernambuco", observed_pn)

nnv_rg <- nnv_list(vacc_alloc_rg, postsim_all_rg_ui, N_rg$`Rio Grande do Norte`, "Rio Grande do Norte", observed_rg)

nnv_pi <- nnv_list(vacc_alloc_pi, postsim_all_pi_ui, N_pi$Piauí, "Piauí", observed_pi)

nnv_tc <- nnv_list(vacc_alloc_tc, postsim_all_tc_ui, N_tc$Tocantins, "Tocantins", observed_tc)

nnv_ag <- nnv_list(vacc_alloc_ag, postsim_all_ag_ui, N_ag$Alagoas, "Alagoas", observed_ag)

nnv_mg <- nnv_list(vacc_alloc_mg, postsim_all_mg_ui, N_mg$`Minas Gerais`, "Minas Gerais", observed_mg)

nnv_se <- nnv_list(vacc_alloc_se, postsim_all_se_ui, N_se$Sergipe, "Sergipe", observed_se)

nnv_go <- nnv_list(vacc_alloc_go, postsim_all_go_ui, N_go$Goiás, "Goiás", observed_go)

# global total -----------------------------------------------------------------
## total pre-post infection graph (national level) ------------------------------
df_ce <- postsim_all_ce_ui$summary_week_df
df_bh <- postsim_all_bh_ui$summary_week_df
df_pa <- postsim_all_pa_ui$summary_week_df
df_pn <- postsim_all_pn_ui$summary_week_df
df_rg <- postsim_all_rg_ui$summary_week_df
df_pi <- postsim_all_pi_ui$summary_week_df
df_tc <- postsim_all_tc_ui$summary_week_df
df_ag <- postsim_all_ag_ui$summary_week_df
df_mg <- postsim_all_mg_ui$summary_week_df
df_se <- postsim_all_se_ui$summary_week_df
df_go <- postsim_all_go_ui$summary_week_df

# Combine the two data frames
combined_df <- bind_rows(df_ce, df_bh, df_pa,
                         df_pn, df_rg, df_pi,
                         df_tc, df_ag, df_mg,
                         df_se, df_go) %>%
  group_by(Scenario, Week) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

# (2) Calculate global impact from your combined_df
global_impact_week <- combined_df %>%
  group_by(Scenario) %>%
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_case   = sum(pre_cases, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_case - total_post_cases,
    impact = diff / total_pre_case * 100
  )

# (3) Create annotation text from the global impact
annotation_text <- paste0(
  dplyr::recode(
    global_impact_week$Scenario,
    "Scenario_1" = "Strategy 1",
    "Scenario_2" = "Strategy 2",
    "Scenario_3" = "Strategy 3",
    "Scenario_4" = "Strategy 4"
  ),
  ": ",
  round(global_impact_week$impact, 1),
  "%",
  collapse = "\n"
)

postsim_all_global_ui <- list(
  postsim_all_ce_ui = postsim_all_ce_ui,
  postsim_all_bh_ui = postsim_all_bh_ui,
  summary_week_df   = combined_df,
  global_impact_week = global_impact_week,
  annotation_text = annotation_text,
  vacc_start_week_s1 = postsim_all_ce_ui$vacc_weeks$scenario1,
  vacc_start_week_s2 = postsim_all_ce_ui$vacc_weeks$scenario2,
  vacc_start_week_s3 = postsim_all_ce_ui$vacc_weeks$scenario3
  
)

epi_graph_all <- 
epi_graph_nat(postsim_all_global_ui, bra_all_sum)

ggsave(filename = "02_Outputs/2_1_Figures/epi_graph_all.jpg", p, width = 8, height = 6, dpi = 1200)


## combined graph
pre_post_graph_comb <- pre_post_graph +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) 

# Remove legend from epi_graph_all (Scenario A)
epi_graph_all <- epi_graph_all + theme(legend.position = "none")

# Arrange plots side-by-side with a common legend
combined_graph <- ggarrange(epi_graph_all, pre_post_graph_comb, 
                           ncol = 2, nrow = 1,
                           labels = c("A", "B"),
                           common.legend = TRUE,
                           legend = "bottom",
                           align = "none")
combined_graph <- combined_graph +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil_0617.jpg", combined_graph, width = 16, height = 6, dpi = 1500)
ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil.pdf", combined_graph, width = 23, height = 10)

### national level scenario graph -------------------------------
target_week_ce <- target_week(postsim_all_ce_ui, prevacc_ui_ce, "Ceará")
target_week_bh <- target_week(postsim_all_bh_ui, prevacc_ui_bh, "Bahia")
target_week_pa <- target_week(postsim_all_pa_ui, prevacc_ui_pa,  "Paraíba")
target_week_pn <- target_week(postsim_all_pn_ui, prevacc_ui_pn, "Pernambuco")
target_week_rg <- target_week(postsim_all_rg_ui, prevacc_ui_rg, "Rio Grande do Norte")
target_week_pi <- target_week(postsim_all_pi_ui, prevacc_ui_pi, "Piauí")
target_week_tc <- target_week(postsim_all_tc_ui, prevacc_ui_tc, "Tocantins")
target_week_ag <- target_week(postsim_all_ag_ui, prevacc_ui_ag, "Alagoas")
target_week_mg <- target_week(postsim_all_mg_ui, prevacc_ui_mg, "Minas Gerais")
target_week_sg <- target_week(postsim_all_se_ui, prevacc_ui_se, "Sergipe")
target_week_go <- target_week(postsim_all_go_ui, prevacc_ui_go, "Goiás")

combined_target_week <- bind_rows(
  target_week_ce,
  target_week_bh,
  target_week_pa,
  target_week_pn,
  target_week_rg,
  target_week_pi,
  target_week_tc,
  target_week_ag,
  target_week_mg,
  target_week_sg,
  target_week_go
)

combined_target_week <- combined_target_week %>%
  left_join(scenario_max_df, by = "Scenario") %>%
  mutate(
    ribbon_ymin = case_when(
      Scenario == "Scenario_1" ~ 2/3 * max_cases,
      Scenario == "Scenario_2" ~ 1/3 * max_cases,
      Scenario == "Scenario_3" ~ 0
    ),
    ribbon_ymax = case_when(
      Scenario == "Scenario_1" ~ max_cases,
      Scenario == "Scenario_2" ~ 2/3 * max_cases,
      Scenario == "Scenario_3" ~ 1/3 * max_cases
    ),
    annotation_y = 0.9 * max_cases
  )

summ_target_df <- combined_target_week %>%
  group_by(Scenario, Week) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

global_impact_week <- combined_target_week %>%
  group_by(Scenario) %>%
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_case   = sum(pre_cases, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_case - total_post_cases,
    impact = diff / total_pre_case * 100
  )

# (3) Create annotation text from the global impact

annotation_df <- global_impact_week %>%
  mutate(
    Strategy = dplyr::recode(Scenario,
                             "Scenario_1" = "Vaccine impact",
                             "Scenario_2" = "Vaccine impact",
                             "Scenario_3" = "Vaccine impact"),
    ann = paste0(Strategy, ": ", round(impact, 1), "%")
  ) %>%
  select(Scenario, ann)

facet_labels <- c(
  "Scenario_1" = "Direct effect in 1-19 years",
  "Scenario_2" = "Direct effect in 20-59 years",
  "Scenario_3" = "Direct effect in ≥60 years"
)

legend_labels <- c(
  "Scenario_1" = "1-19 years",
  "Scenario_2" = "20-59 years",
  "Scenario_3" = "≥60 years"
)

y_pos <- max(summ_target_df$post_cases_hi, na.rm = TRUE) 

target_p <- 
ggplot(summ_target_df) +
  geom_ribbon(data = summ_target_df %>% filter(Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 0, ymax = Inf), fill = "grey70",
              alpha = 0.3)  +
  geom_ribbon(aes(x = Week, ymin = post_cases_lo, ymax = post_cases_hi, fill = Scenario),
              alpha = 0.3) +
  geom_line(aes(x = Week, y = pre_cases, linetype = "Without vaccination"), color = "black", linewidth = 0.8) +
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "With vaccination"), linewidth = 0.8) + 
  geom_point(data = bra_all_sum_target, aes(x = Week, y = Observed, shape = "Observed"),size = 0.6) + 
  geom_ribbon(aes(x = Week, ymin = pre_cases_lo, ymax = pre_cases_hi),
              fill = "lightgrey", alpha = 0.6) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  #geom_text(
  #  data = summ_target_df %>%
  #    filter(Week == 4) %>%
  #    distinct(Scenario, annotation_y),
  #  mapping = aes(x = 4, y = y_pos, label = "<----- Vaccine impact start"),
  #  inherit.aes = FALSE,
  #  hjust = 0, vjust = 1, size = 3
  #) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = Inf, label = ann), 
            inherit.aes = FALSE,
            hjust = 1.1, 
            vjust = 1.1, 
            size = 3) +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  #geom_point(data = observed_all, aes(x = Week, y = Observed),size = 0.05) + 
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = facet_labels))+
  theme_pubclean()+
  theme(legend.position = "right")+
  scale_y_continuous(labels = comma)+
  scale_fill_discrete(
    name   = "Vaccination strategy",
    labels = legend_labels
  )+
  scale_color_discrete(
    name   = "Vaccination strategy",
    labels = legend_labels
  )+
  labs(
    y = "Predicted symptomatic cases",
    linetype = "Type"
  ) 

target_p <- target_p + 
  theme(
    axis.title.y = element_blank()
  )

combined_nat <- ggarrange(epi_graph_all, target_p, 
                            ncol = 2, nrow = 1,
                            labels = c("A", "B"),
                            common.legend = TRUE,
                            legend = "bottom",
                            align = "none")


subnat <- ggarrange(pre_post_graph, 
                    ncol = 1, nrow = 1,
                    labels = c("C"),
                    common.legend = TRUE,
                    legend = "bottom",
                    align = "none")


subnat <- annotate_figure(
  subnat,
  bottom = text_grob(
    "Note: Grey dashed line shows predicted symptomatic cases adjusted for underreporting; black dots show predicted reported symptomatic cases", 
    hjust = 0, x = 0,           
    face  = "italic",           
    size  = 10
  )
)

ggsave(filename = "02_Outputs/2_1_Figures/combined_nat_0617.jpg", combined_nat, width = 16, height = 6, dpi = 1500)
ggsave(filename = "02_Outputs/2_1_Figures/combined_subnat_0617.jpg", subnat, width = 12, height = 6)

## total nnv graph (national level) ------------------------------
n_scenarios = 4

df_nnv_ce <- nnv_ce$final_summ_df
df_nnv_bh <- nnv_bh$final_summ_df
df_nnv_pa <- nnv_pa$final_summ_df
df_nnv_pn <- nnv_pn$final_summ_df
df_nnv_rg <- nnv_rg$final_summ_df
df_nnv_pi <- nnv_pi$final_summ_df
df_nnv_tc <- nnv_tc$final_summ_df
df_nnv_ag <- nnv_ag$final_summ_df
df_nnv_mg <- nnv_mg$final_summ_df
df_nnv_se <- nnv_se$final_summ_df
df_nnv_go <- nnv_go$final_summ_df

df_nnv_ce$setting <- "High"
df_nnv_bh$setting <- "Low"
df_nnv_pa$setting <- "High"
df_nnv_pn$setting <- "Moderate"
df_nnv_rg$setting <- "Low"
df_nnv_pi$setting <- "High"
df_nnv_tc$setting <- "Moderate"
df_nnv_ag$setting <- "High"
df_nnv_mg$setting <- "Low"
df_nnv_se$setting <- "Low"
df_nnv_go$setting <- "Low"

df_nnv_ce$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_ce$age_gr <- factor(df_nnv_ce$age_gr, levels = age_gr_levels)

df_nnv_bh$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_bh$age_gr <- factor(df_nnv_bh$age_gr, levels = age_gr_levels)

df_nnv_pa$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_pa$age_gr <- factor(df_nnv_pa$age_gr, levels = age_gr_levels)

df_nnv_pn$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_pn$age_gr <- factor(df_nnv_pn$age_gr, levels = age_gr_levels)

df_nnv_rg$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_rg$age_gr <- factor(df_nnv_rg$age_gr, levels = age_gr_levels)

df_nnv_pi$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_pi$age_gr <- factor(df_nnv_pi$age_gr, levels = age_gr_levels)

df_nnv_tc$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_tc$age_gr <- factor(df_nnv_tc$age_gr, levels = age_gr_levels)

df_nnv_ag$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_ag$age_gr <- factor(df_nnv_ag$age_gr, levels = age_gr_levels)

df_nnv_mg$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_mg$age_gr <- factor(df_nnv_mg$age_gr, levels = age_gr_levels)

df_nnv_se$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_se$age_gr <- factor(df_nnv_se$age_gr, levels = age_gr_levels)

df_nnv_go$age_gr <- rep(age_gr[1:length(age_gr_levels)],n_scenarios)
df_nnv_go$age_gr <- factor(df_nnv_go$age_gr, levels = age_gr_levels)


combined_nnv_df_region <- bind_rows(df_nnv_ce, df_nnv_bh, df_nnv_pa,
                                    df_nnv_pn, df_nnv_rg, df_nnv_pi,
                                    df_nnv_tc, df_nnv_ag, df_nnv_mg,
                                    df_nnv_go, df_nnv_se)


combined_nnv_df_region_all <- combined_nnv_df_region %>%
  group_by(region) %>% 
  summarise(
    tot_vacc       = sum(tot_vacc,      na.rm = TRUE),
    tot_pop        = sum(tot_pop),
    # raw nums 
    pre_vacc      = sum(pre_vacc),
    pre_vacc_lo   = sum(pre_vacc_lo),
    pre_vacc_hi   = sum(pre_vacc_hi),
    post_vacc     = sum(post_vacc),
    post_vacc_lo  = sum(post_vacc_lo),
    post_vacc_hi  = sum(post_vacc_hi),
    pre_fatal     = sum(pre_fatal),
    pre_fatal_lo  = sum(pre_fatal_lo),
    pre_fatal_hi  = sum(pre_fatal_hi),
    post_fatal    = sum(post_fatal),
    post_fatal_lo = sum(post_fatal_lo),
    post_fatal_hi = sum(post_fatal_hi),
    pre_daly      = sum(pre_daly),
    pre_daly_lo   = sum(pre_daly_lo),
    pre_daly_hi   = sum(pre_daly_hi),
    post_daly     = sum(post_daly),
    post_daly_lo  = sum(post_daly_lo),
    post_daly_hi  = sum(post_daly_hi),
    pre_vacc_per100k = pre_vacc / tot_pop * 1e5,
    pre_vacc_per100k_lo = pre_vacc_lo / tot_pop * 1e5, 
    pre_vacc_per100k_hi = pre_vacc_hi / tot_pop * 1e5,
    # diffs
    diff           = sum(diff,          na.rm = TRUE),
    diff_low       = sum(diff_low,      na.rm = TRUE),
    diff_hi        = sum(diff_hi,       na.rm = TRUE),
    diff_fatal     = sum(diff_fatal,    na.rm = TRUE),
    diff_fatal_low = sum(diff_fatal_low,na.rm = TRUE),
    diff_fatal_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    diff_daly      = sum(diff_daly,     na.rm = TRUE),
    diff_daly_low  = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi   = sum(diff_daly_hi,  na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(
    nnv         = tot_vacc / diff,
    nnv_lo      = tot_vacc / diff_hi,
    nnv_hi      = tot_vacc / diff_low,
    nnv_fatal       = tot_vacc / diff_fatal,
    nnv_fatal_lo    = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi    = tot_vacc / diff_fatal_low,
    nnv_daly        = tot_vacc / diff_daly,
    nnv_daly_lo     = tot_vacc / diff_daly_hi,
    nnv_daly_hi     = tot_vacc / diff_daly_low
  )


# 전체 시나리오별 총합 (세팅 모두 aggregated)
combined_nnv_national <- combined_nnv_df_region %>%
  group_by(scenario) %>% 
  summarise(
    tot_vacc       = sum(tot_vacc,      na.rm = TRUE),
    # raw nums 
    pre_vacc      = sum(pre_vacc),
    pre_vacc_lo   = sum(pre_vacc_lo),
    pre_vacc_hi   = sum(pre_vacc_hi),
    post_vacc     = sum(post_vacc),
    post_vacc_lo  = sum(post_vacc_lo),
    post_vacc_hi  = sum(post_vacc_hi),
    pre_fatal     = sum(pre_fatal),
    pre_fatal_lo  = sum(pre_fatal_lo),
    pre_fatal_hi  = sum(pre_fatal_hi),
    post_fatal    = sum(post_fatal),
    post_fatal_lo = sum(post_fatal_lo),
    post_fatal_hi = sum(post_fatal_hi),
    pre_daly      = sum(pre_daly),
    pre_daly_lo   = sum(pre_daly_lo),
    pre_daly_hi   = sum(pre_daly_hi),
    post_daly     = sum(post_daly),
    post_daly_lo  = sum(post_daly_lo),
    post_daly_hi  = sum(post_daly_hi),
    # diffs
    diff           = sum(diff,          na.rm = TRUE),
    diff_low       = sum(diff_low,      na.rm = TRUE),
    diff_hi        = sum(diff_hi,       na.rm = TRUE),
    diff_fatal     = sum(diff_fatal,    na.rm = TRUE),
    diff_fatal_low = sum(diff_fatal_low,na.rm = TRUE),
    diff_fatal_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    diff_daly      = sum(diff_daly,     na.rm = TRUE),
    diff_daly_low  = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi   = sum(diff_daly_hi,  na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(
    nnv         = tot_vacc / diff,
    nnv_lo      = tot_vacc / diff_hi,
    nnv_hi      = tot_vacc / diff_low,
    nnv_fatal       = tot_vacc / diff_fatal,
    nnv_fatal_lo    = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi    = tot_vacc / diff_fatal_low,
    nnv_daly        = tot_vacc / diff_daly,
    nnv_daly_lo     = tot_vacc / diff_daly_hi,
    nnv_daly_hi     = tot_vacc / diff_daly_low
  )

# setting 별, 시나리오별, 연령별 
combined_nnv_setting <- combined_nnv_df_region %>% 
  group_by(setting, scenario, target) %>% 
  summarise(
    tot_vacc       = sum(tot_vacc,      na.rm = TRUE),
    # raw nums 
    pre_vacc      = sum(pre_vacc),
    pre_vacc_lo   = sum(pre_vacc_lo),
    pre_vacc_hi   = sum(pre_vacc_hi),
    post_vacc     = sum(post_vacc),
    post_vacc_lo  = sum(post_vacc_lo),
    post_vacc_hi  = sum(post_vacc_hi),
    pre_fatal     = sum(pre_fatal),
    pre_fatal_lo  = sum(pre_fatal_lo),
    pre_fatal_hi  = sum(pre_fatal_hi),
    post_fatal    = sum(post_fatal),
    post_fatal_lo = sum(post_fatal_lo),
    post_fatal_hi = sum(post_fatal_hi),
    pre_daly      = sum(pre_daly),
    pre_daly_lo   = sum(pre_daly_lo),
    pre_daly_hi   = sum(pre_daly_hi),
    post_daly     = sum(post_daly),
    post_daly_lo  = sum(post_daly_lo),
    post_daly_hi  = sum(post_daly_hi),
    #diffs
    diff           = sum(diff,          na.rm = TRUE),
    diff_low       = sum(diff_low,      na.rm = TRUE),
    diff_hi        = sum(diff_hi,       na.rm = TRUE),
    diff_fatal     = sum(diff_fatal,    na.rm = TRUE),
    diff_fatal_low = sum(diff_fatal_low,na.rm = TRUE),
    diff_fatal_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    diff_daly      = sum(diff_daly,     na.rm = TRUE),
    diff_daly_low  = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi   = sum(diff_daly_hi,  na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  # 1) Scenario 전체 투여량 추가
  group_by(setting, scenario) %>%
  mutate(
    scenario_vacc = sum(tot_vacc, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    nnv_fatal       = scenario_vacc / diff_fatal,
    nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
    nnv_daly        = scenario_vacc / diff_daly,
    nnv_daly_lo     = scenario_vacc / diff_daly_hi,
    nnv_daly_hi     = scenario_vacc / diff_daly_low
  )

# scenario, target (national)
combined_nnv_setting_nat_target <- combined_nnv_df_region %>% 
  group_by(scenario, target) %>% 
  summarise(
    tot_vacc       = sum(tot_vacc,      na.rm = TRUE),
    # raw nums 
    pre_vacc      = sum(pre_vacc),
    pre_vacc_lo   = sum(pre_vacc_lo),
    pre_vacc_hi   = sum(pre_vacc_hi),
    post_vacc     = sum(post_vacc),
    post_vacc_lo  = sum(post_vacc_lo),
    post_vacc_hi  = sum(post_vacc_hi),
    pre_fatal     = sum(pre_fatal),
    pre_fatal_lo  = sum(pre_fatal_lo),
    pre_fatal_hi  = sum(pre_fatal_hi),
    post_fatal    = sum(post_fatal),
    post_fatal_lo = sum(post_fatal_lo),
    post_fatal_hi = sum(post_fatal_hi),
    pre_daly      = sum(pre_daly),
    pre_daly_lo   = sum(pre_daly_lo),
    pre_daly_hi   = sum(pre_daly_hi),
    post_daly     = sum(post_daly),
    post_daly_lo  = sum(post_daly_lo),
    post_daly_hi  = sum(post_daly_hi),
    
    diff           = sum(diff,          na.rm = TRUE),
    diff_low       = sum(diff_low,      na.rm = TRUE),
    diff_hi        = sum(diff_hi,       na.rm = TRUE),
    diff_fatal     = sum(diff_fatal,    na.rm = TRUE),
    diff_fatal_low = sum(diff_fatal_low,na.rm = TRUE),
    diff_fatal_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    diff_daly      = sum(diff_daly,     na.rm = TRUE),
    diff_daly_low  = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi   = sum(diff_daly_hi,  na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  ungroup() %>%
  # 시나리오별 전체 백신량(=scenario_vacc) 계산
  group_by(scenario) %>%
  mutate(
    scenario_vacc = sum(tot_vacc, na.rm = TRUE)
  )%>%
  ungroup() %>%
  mutate(
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    nnv_fatal       = scenario_vacc / diff_fatal,
    nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
    nnv_daly        = scenario_vacc / diff_daly,
    nnv_daly_lo     = scenario_vacc / diff_daly_hi,
    nnv_daly_hi     = scenario_vacc / diff_daly_low
  )

# 지역별, 시나리오별 총합
combined_nnv_setting_summ <- combined_nnv_setting %>% 
  group_by(setting, scenario) %>% 
  summarise(
    tot_vacc       = sum(tot_vacc,      na.rm = TRUE),
    # raw nums 
    pre_vacc      = sum(pre_vacc),
    pre_vacc_lo   = sum(pre_vacc_lo),
    pre_vacc_hi   = sum(pre_vacc_hi),
    post_vacc     = sum(post_vacc),
    post_vacc_lo  = sum(post_vacc_lo),
    post_vacc_hi  = sum(post_vacc_hi),
    pre_fatal     = sum(pre_fatal),
    pre_fatal_lo  = sum(pre_fatal_lo),
    pre_fatal_hi  = sum(pre_fatal_hi),
    post_fatal    = sum(post_fatal),
    post_fatal_lo = sum(post_fatal_lo),
    post_fatal_hi = sum(post_fatal_hi),
    pre_daly      = sum(pre_daly),
    pre_daly_lo   = sum(pre_daly_lo),
    pre_daly_hi   = sum(pre_daly_hi),
    post_daly     = sum(post_daly),
    post_daly_lo  = sum(post_daly_lo),
    post_daly_hi  = sum(post_daly_hi),
    
    diff           = sum(diff,          na.rm = TRUE),
    diff_low       = sum(diff_low,      na.rm = TRUE),
    diff_hi        = sum(diff_hi,       na.rm = TRUE),
    diff_fatal     = sum(diff_fatal,    na.rm = TRUE),
    diff_fatal_low = sum(diff_fatal_low,na.rm = TRUE),
    diff_fatal_hi  = sum(diff_fatal_hi, na.rm = TRUE),
    diff_daly      = sum(diff_daly,     na.rm = TRUE),
    diff_daly_low  = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi   = sum(diff_daly_hi,  na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(
    nnv         = tot_vacc / diff,
    nnv_lo      = tot_vacc / diff_hi,
    nnv_hi      = tot_vacc / diff_low,
    nnv_fatal       = tot_vacc / diff_fatal,
    nnv_fatal_lo    = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi    = tot_vacc / diff_fatal_low,
    nnv_daly        = tot_vacc / diff_daly,
    nnv_daly_lo     = tot_vacc / diff_daly_hi,
    nnv_daly_hi     = tot_vacc / diff_daly_low
  )




write_xlsx(combined_nnv_region_all, path = "02_Outputs/2_2_Tables/combined_nnv_region_all.xlsx")
write_xlsx(combined_nnv_df_region, path = "02_Outputs/2_2_Tables/combined_nnv_df_region.xlsx")


# combining dfs
combined_nnv_national$setting <- "National"
combined_nnv_setting_summ <- rbind(combined_nnv_national, combined_nnv_setting_summ)

combined_nnv_setting_nat_target$setting <- "National"
combined_nnv_setting <- rbind(combined_nnv_setting, combined_nnv_setting_nat_target)


## for faceting all subnat
nnv_nat_gg(combined_nnv_df_region, 
           y_var = "nnv",
           y_lab = "NNV to avert a single symptomatic case",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")
nnv_nat_gg(combined_nnv_df_region, 
           y_var = "nnv_fatal",
           y_lab = "NNV to avert a single fatal case",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")
nnv_nat_gg(combined_nnv_df_region, 
           y_var = "nnv_daly",
           y_lab = "NNV to avert a single DALY",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")



p1 <- 
  nnv_nat_setting_summ(combined_nnv_setting_summ, 
                y_var = "nnv",
                y_lab = "NNV to avert a symptomatic case",
                x_lab = "Age group",
                title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

p2 <- 
  nnv_nat_setting_summ(combined_nnv_setting_summ, 
                     y_var = "nnv_fatal",
                     y_lab = "NNV to avert a death",
                     x_lab = "Age group",
                     title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")
p3 <- 
  nnv_nat_setting_summ(combined_nnv_setting_summ, 
                     y_var = "nnv_daly",
                     y_lab = "NNV to avert a DALY",
                     x_lab = "Age group",
                     title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

p1 <- p1 +
  #scale_x_discrete(labels = function(x) gsub(" years", "", x)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a symptomatic case")

p2 <- p2 +
  #scale_x_discrete(labels = function(x) gsub(" years", "", x)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a death")

p3 <- p3 +
  #scale_x_discrete(labels = function(x) gsub(" years", "", x)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a DALY")

# Arrange plots side-by-side with a common legend
combined_nnv_graph <- ggarrange(p1, p2, p3,
                            ncol = 1, nrow = 3,
                            labels = c("A", "B", "C"),
                            common.legend = TRUE,
                            legend = "bottom",
                            align = "v")

ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_v2.jpg", combined_nnv_graph, width = 12, height = 10, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_all_v3.jpg", combined_nnv_graph, width = 12, height = 10, dpi = 1200)

write_xlsx(combined_nnv_setting, path = "02_Outputs/2_2_Tables/combined_nnv_setting.xlsx")
write_xlsx(combined_nnv_setting_summ, path = "02_Outputs/2_2_Tables/combined_nnv_setting_summ.xlsx")

# table for publication (nnv table) -------------------------------------------------------
library(glue)
library(writexl)
# remove all zero rows

combined_nnv_national <- combined_nnv_national %>% filter(if_all(everything(), ~. !=0))

metrics <- c("pre_vacc", "post_vacc", "pre_fatal",
             "pre_daly", "post_fatal", "post_daly",
             "nnv", "nnv_fatal", "nnv_daly")

table_nnv_ci <- combined_nnv_national %>%
  rowwise() %>%
  mutate(
    
    prevacc_ci = glue(
      "{formatC(pre_vacc,        format='f', digits=0, big.mark=',')} ",
      "({formatC(pre_vacc_lo,  format='f', digits=0, big.mark=',')}–",
      "{formatC(pre_vacc_hi,   format='f', digits=0, big.mark=',')})"
    ),
    pre_fatal_ci = glue(
      # fatal 계열만 소수점 2자리
      "{formatC(pre_fatal,       format='f', digits=0, big.mark=',')} ",
      "({formatC(pre_fatal_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(pre_fatal_hi,  format='f', digits=0, big.mark=',')})"
    ),
    pre_daly_ci = glue(
      "{formatC(pre_daly,       format='f', digits=0, big.mark=',')} ",
      "({formatC(pre_daly_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(pre_daly_hi,  format='f', digits=0, big.mark=',')})"
    ),
    postvacc_ci = glue(
      "{formatC(post_vacc,       format='f', digits=0, big.mark=',')} ",
      "({formatC(post_vacc_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(post_vacc_hi,  format='f', digits=0, big.mark=',')})"
    ),
    post_fatal_ci = glue(
      "{formatC(post_fatal,       format='f', digits=0, big.mark=',')} ",
      "({formatC(post_fatal_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(post_fatal_hi,  format='f', digits=0, big.mark=',')})"
    ),
    post_daly_ci = glue(
      "{formatC(post_daly,       format='f', digits=0, big.mark=',')} ",
      "({formatC(post_daly_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(post_daly_hi,  format='f', digits=0, big.mark=',')})"
    )
  ) %>%
  ungroup() %>%
  select(
    scenario, 
    ends_with("_ci")
  )

write_xlsx(table_nnv_ci, path = "02_Outputs/2_2_Tables/table_nnv_ci_veinf0.xlsx")
write_xlsx(table_nnv_ci, path = "02_Outputs/2_2_Tables/table_nnv_ci_veinf99.xlsx")


# make per1M dataframe -----------------------------------------------------------------
n_scenarios = 4

per1m_ce <- nnv_ce$per1M_summary
per1m_bh <- nnv_bh$per1M_summary
per1m_pa <- nnv_pa$per1M_summary
per1m_pn <- nnv_pn$per1M_summary
per1m_rg <- nnv_rg$per1M_summary
per1m_pi <- nnv_pi$per1M_summary
per1m_tc <- nnv_tc$per1M_summary
per1m_ag <- nnv_ag$per1M_summary
per1m_mg <- nnv_mg$per1M_summary
per1m_se <- nnv_se$per1M_summary
per1m_go <- nnv_go$per1M_summary

per1m_ce$setting <- "High"
per1m_bh$setting <- "Low"
per1m_pa$setting <- "High"
per1m_pn$setting <- "Moderate"
per1m_rg$setting <- "Low"
per1m_pi$setting <- "High"
per1m_tc$setting <- "Moderate"
per1m_ag$setting <- "High"
per1m_mg$setting <- "Low"
per1m_se$setting <- "Low"
per1m_go$setting <- "Low"

per1m_all <- bind_rows(per1m_ag, per1m_bh, per1m_ce,
                       per1m_go, per1m_mg, per1m_pi, 
                       per1m_pn, per1m_rg, per1m_se, per1m_tc)

# 1) setting × scenario별 합산 후 per1M 재계산
summ_by_setting <- per1m_all  %>%
  group_by(setting, scenario, target) %>%
  summarise(
    scenario_vacc    = sum(scenario_vacc, na.rm = TRUE),
    tot_pop          = sum(tot_pop),
    vacc_prop        = scenario_vacc / tot_pop,
    diff             = sum(diff,           na.rm = TRUE),
    diff_low         = sum(diff_low,       na.rm = TRUE),
    diff_hi          = sum(diff_hi,        na.rm = TRUE),
    diff_fatal       = sum(diff_fatal,     na.rm = TRUE),
    diff_fatal_low   = sum(diff_fatal_low, na.rm = TRUE),
    diff_fatal_hi    = sum(diff_fatal_hi,  na.rm = TRUE),
    diff_daly        = sum(diff_daly,      na.rm = TRUE),
    diff_daly_low    = sum(diff_daly_low,  na.rm = TRUE),
    diff_daly_hi     = sum(diff_daly_hi,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    per1M_diff      = diff           / scenario_vacc * 1e5,
    per1M_diff_lo   = diff_low       / scenario_vacc * 1e5,
    per1M_diff_hi   = diff_hi        / scenario_vacc * 1e5,
    per1M_fatal     = diff_fatal     / scenario_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low / scenario_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi  / scenario_vacc * 1e5,
    per1M_daly      = diff_daly      / scenario_vacc * 1e5,
    per1M_daly_lo   = diff_daly_low  / scenario_vacc * 1e5,
    per1M_daly_hi   = diff_daly_hi   / scenario_vacc * 1e5,
    
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    
    nnv_fatal       = scenario_vacc / diff_fatal,
    nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
    
    nnv_daly        = scenario_vacc / diff_daly,
    nnv_daly_lo     = scenario_vacc / diff_daly_hi,
    nnv_daly_hi     = scenario_vacc / diff_daly_low
  )

# 2) scenario별(모든 setting 합친) 합산 후 per1M 재계산, setting = "National"
summ_national <- per1m_all %>%
  group_by(scenario, target) %>%
  summarise(
    # case
    diff       = sum(diff,           na.rm = TRUE),
    diff_low    = sum(diff_low,       na.rm = TRUE),
    diff_hi    = sum(diff_hi,        na.rm = TRUE),
    scenario_vacc   = sum(scenario_vacc,  na.rm = TRUE),
    tot_pop          = sum(tot_pop),
    vacc_prop        = scenario_vacc / tot_pop,
    
    # fatal
    diff_fatal       = sum(diff_fatal,           na.rm = TRUE),
    diff_fatal_low    = sum(diff_fatal_low, na.rm = TRUE),
    diff_fatal_hi    = sum(diff_fatal_hi,  na.rm = TRUE),
    
    # daly
    diff_daly       = sum(diff_daly,           na.rm = TRUE),
    diff_daly_low    = sum(diff_daly_low, na.rm = TRUE),
    diff_daly_hi    = sum(diff_daly_hi,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    setting        = "National",
    per1M_diff     = diff       / scenario_vacc * 1e5,
    per1M_diff_lo  = diff_low    / scenario_vacc * 1e5,
    per1M_diff_hi  = diff_hi    / scenario_vacc * 1e5,
    
    per1M_fatal     = diff_fatal       / scenario_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low    / scenario_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi    / scenario_vacc * 1e5,
    
    per1M_daly     = diff_daly       / scenario_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low    / scenario_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi    / scenario_vacc * 1e5,
    
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    
    nnv_fatal       = scenario_vacc / diff_fatal,
    nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi    = scenario_vacc / diff_fatal_low,
    
    nnv_daly        = scenario_vacc / diff_daly,
    nnv_daly_lo     = scenario_vacc / diff_daly_hi,
    nnv_daly_hi     = scenario_vacc / diff_daly_low
  ) 

# 3) 둘을 합치기
per1m_summary_full <- bind_rows(summ_by_setting, summ_national)

plot_df <- per1m_summary_full %>%
  # 1) 시나리오별 실제 target age 지정
  mutate(
    target_of_scenario = case_when(
      scenario == "Scenario_1" ~ "1-11 years",
      scenario == "Scenario_2" ~ "12-17 years",
      scenario == "Scenario_3" ~ "18-59 years",
      scenario == "Scenario_4" ~ "60+ years"
    ),
    # 2) Direct vs Indirect + target 명명
    effect_group = if_else(
      target == target_of_scenario,
      "Total effect in target group",
      paste0("Indirect — ", target)
    ),
    # 3) 순서 고정
    scenario = factor(scenario, levels = c("Scenario_1","Scenario_2","Scenario_3", "Scenario_4")),
    setting  = factor(setting,   levels = c("National", "High", "Moderate","Low"))
  )

plot_df <- plot_df %>%
  mutate(
    effect_group = factor(effect_group,
                          levels = c(
                            "Total effect in target group",
                            "Indirect — <1 (novacc)",
                            "Indirect — 1-11 years",
                            "Indirect — 12-17 years",
                            "Indirect — 18-59 years",
                            "Indirect — 60+ years"  # ✅ 오타 수정
                          )
    )
  ) %>%
  mutate(
    scenario = factor(scenario,
                      levels = c("Scenario_1","Scenario_2","Scenario_3", "Scenario_4"),
                      labels = c(
                        "Strategy 1 (1–11 years)",
                        "Strategy 2 (12-17 years)",
                        "Strategy 3 (18-59 years)",
                        "Strategy 4 (60+ years)"
                      )
    )
  )

# combined nnv + per1m graphs -----------------------------------------------
# symp
per1m_symp <- make_per1m_plot(plot_df, 
                y_var = "per1M_diff",
                y_lo  = "per1M_diff_lo",     
                y_hi  = "per1M_diff_hi",
                y_label = "Cases averted per 100,000 doses") +
  theme(legend.position = "none") +
  guides(fill = "none", colour = "none")

a <- combine_nnv_per1m(p1, per1m_symp, tag = "A") +theme(axis.title.x = element_blank(),
                                                         axis.text.x  = element_blank(),
                                                         axis.ticks.x = element_blank()
                                                         ) +  guides(fill = "none", colour = "none")
# fatal
per1m_fatal <- make_per1m_plot(plot_df, 
                              y_var = "per1M_fatal",
                              y_lo  = "per1M_fatal_lo",     
                              y_hi  = "per1M_fatal_hi",
                              y_label = "Deaths averted per 100,000 doses")

b <- combine_nnv_per1m(p2, per1m_fatal, tag = "B")+theme(axis.title.x = element_blank(),
                                                         axis.text.x  = element_blank(),
                                                         axis.ticks.x = element_blank()
)

# daly 
per1m_daly <- make_per1m_plot(plot_df, 
                               y_var = "per1M_daly",
                               y_lo  = "per1M_daly_lo",     
                               y_hi  = "per1M_daly_hi",
                               y_label = "DALYs averted per 100,000 doses")+
  theme(legend.position = "none") +
  guides(fill = "none", colour = "none")

c <- combine_nnv_per1m(p3, per1m_daly, tag = "C")

combined_nnv_graph <-  wrap_plots(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "right"),
  c + theme(legend.position = "none"),
  ncol   = 1,
  guides = "collect"        # 남아 있는 a의 범례만 수집
) 

ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_per1m.jpg", combined_nnv_graph, width = 14, height = 16, dpi = 1200)

# table for publication -------------------------------------------------------
library(glue)
library(writexl)
# remove all zero rows

table_per1m <- plot_df %>% filter(if_all(everything(), ~. !=0))

metrics <- c("diff", "diff_fatal", "diff_daly",
             "per1M_diff", "per1M_fatal", "per1M_daly",
             "nnv", "nnv_fatal", "nnv_daly")

table_per1m_ci <- table_per1m %>%
  rowwise() %>%
  mutate(
    
    vacc_prop = glue("{formatC(vacc_prop * 100, format = 'f', digits = 1)}%"),
    scenario_vacc = formatC(scenario_vacc,
                            format   = "f",
                            digits   = 0,
                            big.mark = ","),
    diff_ci = glue(
      "{formatC(diff,        format='f', digits=0, big.mark=',')} ",
      "({formatC(diff_low,  format='f', digits=0, big.mark=',')}–",
      "{formatC(diff_hi,   format='f', digits=0, big.mark=',')})"
    ),
    diff_fatal_ci = glue(
      # fatal 계열만 소수점 2자리
      "{formatC(diff_fatal,       format='f', digits=0, big.mark=',')} ",
      "({formatC(diff_fatal_low, format='f', digits=0, big.mark=',')}–",
      "{formatC(diff_fatal_hi,  format='f', digits=0, big.mark=',')})"
    ),
    diff_daly_ci = glue(
      "{formatC(diff_daly,       format='f', digits=0, big.mark=',')} ",
      "({formatC(diff_daly_low, format='f', digits=0, big.mark=',')}–",
      "{formatC(diff_daly_hi,  format='f', digits=0, big.mark=',')})"
    ),
    per1M_diff_ci = glue(
      "{formatC(per1M_diff,       format='f', digits=0, big.mark=',')} ",
      "({formatC(per1M_diff_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(per1M_diff_hi,  format='f', digits=0, big.mark=',')})"
    ),
    per1M_fatal_ci = glue(
      "{formatC(per1M_fatal,       format='f', digits=2, big.mark=',')} ",
      "({formatC(per1M_fatal_lo, format='f', digits=2, big.mark=',')}–",
      "{formatC(per1M_fatal_hi,  format='f', digits=2, big.mark=',')})"
    ),
    per1M_daly_ci = glue(
      "{formatC(per1M_daly,       format='f', digits=0, big.mark=',')} ",
      "({formatC(per1M_daly_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(per1M_daly_hi,  format='f', digits=0, big.mark=',')})"
    ),
    nnv_ci = glue(
      "{formatC(nnv,       format='f', digits=0, big.mark=',')} ",
      "({formatC(nnv_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(nnv_hi,  format='f', digits=0, big.mark=',')})"
    ),
    nnv_fatal_ci = glue(
      "{formatC(nnv_fatal,       format='f', digits=0, big.mark=',')} ",
      "({formatC(nnv_fatal_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(nnv_fatal_hi,  format='f', digits=0, big.mark=',')})"
    ),
    nnv_daly_ci = glue(
      "{formatC(nnv_daly,       format='f', digits=0, big.mark=',')} ",
      "({formatC(nnv_daly_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(nnv_daly_hi,  format='f', digits=0, big.mark=',')})"
    )
  ) %>%
  ungroup() %>%
  select(
    setting, scenario, scenario_vacc, vacc_prop,
    ends_with("_ci")
  )

order_vec <- c("National", "High", "Moderate", "Low")
table_per1m_ci <- table_per1m_ci %>%
  arrange(match(setting, order_vec))

colnames(table_per1m_ci) <- c("Setting", "Vaccination strategy", "Total number of vaccines",
                              "Vaccine coverage", "Differences in cases (95%UI)",
                              "Differences in deaths (95%UI)", 
                              "Differences in DALYs (95%UI)",
                              "Cases averted per 100K doses (95%UI)",
                              "Deaths averted per 100K doses (95%UI)",
                              "DALYs averted per 100K doses (95%UI)",
                              "NNV to avert a symptomatic case (95%UI)",
                              "NNV to avert a death (95%UI)",
                              "NNV to avert a DALY (95%UI)"
                              )

write_xlsx(table_per1m_ci, path = "02_Outputs/2_2_Tables/table_per1m_ci_veinf0.xlsx")
write_xlsx(table_per1m_ci, path = "02_Outputs/2_2_Tables/table_per1m_ci_veinf99.xlsx")

# case reduction by region -----------------------------------------------------
global_impact_ce <- postsim_all_ce_ui$global_impact_week
global_impact_bh <- postsim_all_bh_ui$global_impact_week
global_impact_pn <- postsim_all_pn_ui$global_impact_week
global_impact_pi <- postsim_all_pi_ui$global_impact_week
global_impact_pa <- postsim_all_pa_ui$global_impact_week
global_impact_rg <- postsim_all_rg_ui$global_impact_week
global_impact_ag <- postsim_all_ag_ui$global_impact_week
global_impact_tc <- postsim_all_tc_ui$global_impact_week
global_impact_mg <- postsim_all_mg_ui$global_impact_week
global_impact_se <- postsim_all_se_ui$global_impact_week
global_impact_go <- postsim_all_go_ui$global_impact_week

combined_impact_df <- bind_rows(global_impact_ce, global_impact_bh, global_impact_pn,
                                global_impact_pi, global_impact_pa, global_impact_rg,
                                global_impact_ag, global_impact_tc, global_impact_mg,
                                global_impact_se, global_impact_go) 

combined_impact_s1 <- combined_impact_df %>% filter(Scenario == "Scenario_1")
combined_impact_s2 <- combined_impact_df %>% filter(Scenario == "Scenario_2")
combined_impact_s3 <- combined_impact_df %>% filter(Scenario == "Scenario_3")

colnames(combined_impact_s1)[2] <-"post_cases_s1"
colnames(combined_impact_s2)[2] <-"post_cases_s2"
colnames(combined_impact_s3)[2] <-"post_cases_s3"

br_states <- br_states %>% 
  left_join(combined_impact_s1[, c("region", colnames(combined_impact_s1)[2])], by = "region") %>%
  left_join(combined_impact_s2[, c("region", colnames(combined_impact_s2)[2])], by = "region") %>%
  left_join(combined_impact_s3[, c("region", colnames(combined_impact_s3)[2])], by = "region") 

br_states <- br_states %>% 
  mutate(post_cases_s1 = coalesce(post_cases_s1, case_22))%>%
  mutate(post_cases_s2 = coalesce(post_cases_s2, case_22))%>%
  mutate(post_cases_s3 = coalesce(post_cases_s3, case_22))

p1 <- 
region_prepost_map("post_cases_s1",
                 title = "Predicted post-vaccination symptomatic cases Scenario 1 (2022)")+
  theme(plot.margin = margin(1, 1, 1, 1, "pt"))

p2 <- 
region_prepost_map("post_cases_s2",
                   title = "Predicted post-vaccination symptomatic cases Scenario 2 (2022)")+
  theme(plot.margin = margin(1, 1, 1, 1, "pt"))

p3 <- 
region_prepost_map("post_cases_s3",
                   title = "Predicted post-vaccination symptomatic cases Scenario 3 (2022)")

ggsave(filename = "02_Outputs/2_1_Figures/map_s1.jpg", p1, width = 15, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/map_s2.jpg", p2, width = 15, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/map_s3.jpg", p3, width = 15, height = 8, dpi = 1200)


# R0 estimation ----------------------------------------------------------------
r0_list_ce <- r0_estim(posterior_ce)
r0_list_bh <- r0_estim(posterior_bh)
r0_list_ag <- r0_estim(posterior_ag)
r0_list_mg <- r0_estim(posterior_mg)
r0_list_rg <- r0_estim(posterior_rg)
r0_list_tc <- r0_estim(posterior_tc)
r0_list_pi <- r0_estim(posterior_pi)
r0_list_pa <- r0_estim(posterior_pa)
r0_list_pn <- r0_estim(posterior_pn)

ggplot(r0_list_ce$r0_df)+
  geom_line(aes(x = time, y = median))


# delay heatmap -----------------------------------------------------------
heatmap_ce <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_ce,
  posterior            = posterior_ce,
  N                    = N_ceara$Ceará,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Ceará",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ce
)

heatmap_bh <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_bh,
  posterior            = posterior_bh,
  N                    = N_bahia$Bahia,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Bahia",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_bh
)

heatmap_pn <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pn,
  posterior            = posterior_pn,
  N                    = N_pemam$Pernambuco,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Pernambuco",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pn
)

heatmap_rg <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_rg,
  posterior            = posterior_rg,
  N                    = N_rg$`Rio Grande do Norte`,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Rio Grande do Norte",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_rg
)

heatmap_pa <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pa,
  posterior            = posterior_pa,
  N                    = N_pa$Paraíba,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Paraíba",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pa
)

heatmap_pi <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pi,
  posterior            = posterior_pi,
  N                    = N_pi$Piauí,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Piauí",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pi
)

heatmap_tc <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_tc,
  posterior            = posterior_tc,
  N                    = N_tc$Tocantins,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Tocantins",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_tc
)

heatmap_ag <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_ag,
  posterior            = posterior_ag,
  N                    = N_ag$Alagoas,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Alagoas",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ag
)

heatmap_mg <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_mg,
  posterior            = posterior_mg,
  N                    = N_mg$`Minas Gerais`,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Minas Gerais",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_mg
)

heatmap_se <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_se,
  posterior            = posterior_se,
  N                    = N_se$Sergipe,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Sergipe",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_se
)

heatmap_go <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_go,
  posterior            = posterior_go,
  N                    = N_go$Goiás,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Goiás",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_go
)

combined_heatmap <- bind_rows(
  heatmap_ag,
  heatmap_bh,
  heatmap_ce,
  heatmap_mg,
  heatmap_pa,
  heatmap_pi,
  heatmap_pn,
  heatmap_rg,
  heatmap_tc,
  heatmap_se,
  heatmap_go
)

save("combined_heatmap", file = "00_Data/0_2_Processed/combined_heatmap_coverage_ixchiq.RData")
save("combined_heatmap", file = "00_Data/0_2_Processed/combined_heatmap_coverage_vimkun.RData")

heatmap_summ <- combined_heatmap %>%
  mutate(
    scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
  )

vacc_df <- combined_vacc_alloc_region %>%
  dplyr::rename(region = Region) %>%
  select(region, Scenario, tot_vacc)

vacc_wide <- combined_vacc_alloc_region %>%
  dplyr::rename(region = Region) %>%
  select(region, Scenario, tot_vacc) %>%
  pivot_wider(
    names_from   = Scenario,
    values_from  = tot_vacc,
    names_prefix = "totvacc_"
  )

heatmap_summ_2 <- heatmap_summ %>%
  left_join(vacc_wide, by = "region") %>%
  select(-starts_with("impact_")) 

# updated heatmap summ
heatmap_summ_2 <- heatmap_summ_2 %>%
  mutate(
    # vaccine parameters
    impact_cases_agegr1 = (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) / totvacc_Scenario_1 * 1e5,
  
    impact_fatal_agegr1 = (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) / totvacc_Scenario_1 * 1e5,
    
    impact_daly_agegr1  = (diff_daly_agegr1  + diff_daly_agegr2  + diff_daly_agegr3)  / totvacc_Scenario_1 * 1e5,
    
    # — Age group 2 direct/indirect/total —
    impact_cases_agegr2 = (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) / totvacc_Scenario_2 * 1e5,
    
    impact_fatal_agegr2 = (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) / totvacc_Scenario_2 * 1e5,
    
    impact_daly_agegr2  = (diff_daly_agegr1  + diff_daly_agegr2  + diff_daly_agegr3)  / totvacc_Scenario_2 * 1e5,
    
    # — Age group 3 direct/indirect/total —
    impact_cases_agegr3 = (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) / totvacc_Scenario_3 * 1e5,
    
    impact_fatal_agegr3 = (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) / totvacc_Scenario_3 * 1e5,
    
    impact_daly_agegr3  = (diff_daly_agegr1  + diff_daly_agegr2  + diff_daly_agegr3)   / totvacc_Scenario_3 * 1e5
  )

collapsed_heatmap <- heatmap_summ_2 %>%
  group_by(Scenario, Supply, Delay) %>%
  summarise(
    across(
      .cols = where(is.numeric),  
      .fns = mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_xlsx(collapsed_heatmap, path = "02_Outputs/2_2_Tables/collapsed_heatmap_v2.xlsx")

library(grid)  
common_margin <- theme(
  plot.margin = unit(c(5, 5, 5, 5), "pt")
)

p_cases <- 
  plot_heatmap(
    collapsed_heatmap,
  "impact_cases_agegr",
  "Cases averted per 100,000 doses"
)

p_fatal <-
  plot_heatmap(
    collapsed_heatmap,
  "impact_fatal_agegr",
  "Deaths averted per 100,000 doses"
)

p_daly <-
  plot_heatmap(
    collapsed_heatmap,
    "impact_daly_agegr",
    "DALYs averted per 100,000 doses"
  )

p_cases <- p_cases + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

p_fatal <- p_fatal + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    strip.text      = element_blank(),      
    strip.background = element_blank()
  )

p_daly <- 
  plot_heatmap(
    collapsed_heatmap,
    "impact_daly_agegr",
    "DALYs averted per 100,000 doses"
  )+
  theme(
    strip.text      = element_blank(),
    strip.background = element_blank()
  )
p_cases  <- p_cases  + common_margin
p_fatal  <- p_fatal  + common_margin
p_daly   <- p_daly   + common_margin

combined_sens_g_heatmap <- 
  ggarrange(p_cases, 
            p_fatal,
            p_daly,
            ncol = 1, nrow = 3,
            labels = c("A", "B", "C"),
            common.legend = F,
            legend = "right",
            align = "v",
            heights      = c(1, 1, 1.2) 
            )

ggsave(filename = "02_Outputs/2_1_Figures/fig_tot_heatmap_0617.jpg", combined_sens_g_heatmap, width = 10, height = 8, dpi = 1200)


## fatal summary----------------------------------------------------------------
cum_fatal_ce <- fatal_summ_func(postsim_all_ce_ui,
                                N_ceara$Ceará,
                                observed_ce,
                                "Ceará")
cum_fatal_bh <- fatal_summ_func(postsim_all_bh_ui,
                                N_bahia$Bahia,
                                observed_bh,
                                "Bahia")
cum_fatal_ag <- fatal_summ_func(postsim_all_ag_ui,
                                N_ag$Alagoas,
                                observed_ag,
                                "Alagoas")
cum_fatal_mg <- fatal_summ_func(postsim_all_mg_ui,
                                N_mg$`Minas Gerais`,
                                observed_mg,
                                "Minas Gerais")
cum_fatal_pa <- fatal_summ_func(postsim_all_pa_ui,
                                N_pa$Paraíba,
                                observed_pa,
                                "Paraíba")
cum_fatal_pi <- fatal_summ_func(postsim_all_pi_ui,
                                N_pi,
                                observed_pi,
                                "Piauí")
cum_fatal_pn <- fatal_summ_func(postsim_all_pn_ui,
                                N_pemam$Pernambuco,
                                observed_pn,
                                "Pernambuco")
cum_fatal_rg <- fatal_summ_func(postsim_all_rg_ui,
                                N_rg$`Rio Grande do Norte`,
                                observed_rg,
                                "Rio Grande do Norte")
cum_fatal_tc <- fatal_summ_func(postsim_all_tc_ui,
                                N_tc$Tocantins,
                                observed_tc,
                                "Tocantins")

cum_fatal_se <- fatal_summ_func(postsim_all_se_ui,
                                N_se$Sergipe,
                                observed_se,
                                "Sergipe")

cum_fatal_go <- fatal_summ_func(postsim_all_go_ui,
                                N_go$Goiás,
                                observed_go,
                                "Goiás")

cum_fatal_ce_df <- cum_fatal_ce$fatal_agegrp3_sum
cum_fatal_bh_df <- cum_fatal_bh$fatal_agegrp3_sum
cum_fatal_pa_df <- cum_fatal_pa$fatal_agegrp3_sum
cum_fatal_pn_df <- cum_fatal_pn$fatal_agegrp3_sum
cum_fatal_rg_df <- cum_fatal_rg$fatal_agegrp3_sum
cum_fatal_pi_df <- cum_fatal_pi$fatal_agegrp3_sum
cum_fatal_tc_df <- cum_fatal_tc$fatal_agegrp3_sum
cum_fatal_ag_df <- cum_fatal_ag$fatal_agegrp3_sum
cum_fatal_mg_df <- cum_fatal_mg$fatal_agegrp3_sum
cum_fatal_se_df <- cum_fatal_se$fatal_agegrp3_sum
cum_fatal_go_df <- cum_fatal_go$fatal_agegrp3_sum

cum_fatal_ce_df$setting <- "High"
cum_fatal_bh_df$setting <- "Low"
cum_fatal_pa_df$setting <- "High"
cum_fatal_pn_df$setting <- "Moderate"
cum_fatal_rg_df$setting <- "Low"
cum_fatal_pi_df$setting <- "High"
cum_fatal_tc_df$setting <- "Moderate"
cum_fatal_ag_df$setting <- "High"
cum_fatal_mg_df$setting <- "Low"
cum_fatal_se_df$setting <- "Low"
cum_fatal_go_df$setting <- "Low"


combined_prepost_fatal <- bind_rows(
  cum_fatal_ce_df,
  cum_fatal_bh_df,
  cum_fatal_ag_df,
  cum_fatal_mg_df,
  cum_fatal_pa_df,
  cum_fatal_pi_df,
  cum_fatal_pn_df,
  cum_fatal_rg_df,
  cum_fatal_tc_df,
  cum_fatal_se_df,
  cum_fatal_go_df
)

combined_prepost_fatal_nat <- combined_prepost_fatal %>% group_by(Scenario) %>% summarise(
  diff_per_million = mean(diff_per_million),
  diff_per_million_lo = mean(diff_per_million_lo),
  diff_per_million_hi = mean(diff_per_million_hi),
)

combined_prepost_fatal_nat$setting <- "National"

combined_prepost_fatal_summ <- combined_prepost_fatal %>% group_by(Scenario, setting) %>% summarise(
  diff_per_million = mean(diff_per_million),
  diff_per_million_lo = mean(diff_per_million_lo),
  diff_per_million_hi = mean(diff_per_million_hi),
)

combined_prepost_fatal_summ <- rbind(combined_prepost_fatal_nat, combined_prepost_fatal_summ)

combined_prepost_fatal_summ$setting <- factor(combined_prepost_fatal_summ$setting, levels = c("National", "Low", "Moderate", "High"))

# 2) grouped bar chart: x=setting, y=diff_per_million, fill=Scenario
ggplot(combined_prepost_fatal_summ, 
       aes(x = setting, y = diff_per_million, fill = Scenario)) +
  # 1) 반투명 막대 (alpha = 0.5)
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.7,
           alpha    = 0.5) +
  # 2) 95% UI 추가
  geom_errorbar(aes(ymin = diff_per_million_lo, ymax = diff_per_million_hi),
                position = position_dodge(width = 0.8),
                width    = 0.2) +
  # 3) ColorBrewer Set1 팔레트 적용
  scale_fill_brewer(palette = "Set1") +
  # 레이블
  labs(
    x    = "Setting",
    y    = "Deaths averted per million vaccinated",
    fill = "Scenario"
  ) +
  # y축 여유, 테마 조정
  scale_y_continuous(expand = expansion(c(0, .05))) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank()
  )




global_impact_fatal <- combined_prepost_fatal %>%
  group_by(Scenario, region) %>%
  summarise(
    total_post_fatal = sum(fatal),
    total_pre_fatal  = sum(pre_fatal),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_fatal - total_post_fatal,
    impact = diff / total_pre_fatal * 100
  )

annotation_df <- global_impact_fatal %>%
  group_by(region) %>%
  summarise(ann = paste0(Scenario, ": ", round(impact, 1), "%", collapse = "\n"),
            .groups = "drop")

overall_max <- max(combined_prepost_fatal$pre_fatal_prop_hi, na.rm = TRUE)

# Create an annotation data frame that has one row per region
vaccine_ann_df <- combined_prepost_fatal %>%
  distinct(region) %>%
  mutate(x = 4, 
         y = overall_max * 0.95, 
         label = "<----- Vaccine impact start (+ 2 wks after initiation)")

pre_post_fatal_graph <- 
  
  ggplot(combined_prepost_fatal) +
  # Scenario shading ribbons (existing)
  geom_ribbon(data = combined_prepost_fatal %>% filter(Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 0, ymax = Inf), fill = "grey90",
              alpha = 0.4) +
  #geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
  #            aes(x = Week, ymin = 1/3 * overall_max, ymax = 2/3 * overall_max, fill = Scenario),
  #            alpha = 0.4) +
  #geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
  #            aes(x = Week, ymin = 0, ymax = 1/3 * overall_max, fill = Scenario),
  #            alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_fatal_prop_lo, ymax = post_fatal_prop_hi, fill = Scenario),
             alpha = 0.2) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = pre_fatal_prop_lo, ymax = pre_fatal_prop_hi),
              fill = "#D8BFD8", alpha = 0.2) +
  # Post-vaccination case line (colored by scenario)
  geom_line(aes(x = Week, y = post_fatal_prop, color = Scenario, linetype = "With vaccination"), size = 0.5) +
  # Pre-vaccination case line (dashed, in black)
  geom_line(aes(x = Week, y = pre_fatal_prop, group = Scenario, linetype = "Without vaccination"), color = "black", size = 0.5) +
  scale_linetype_manual(values = c("With vaccination" = "solid", "Without vaccination" = "dashed")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Scenario_1" = "<20 years only", 
                               "Scenario_2" = "20-59 years only", 
                               "Scenario_3" = ">60 years only")) +
  scale_color_brewer(palette = "Set1",
                     labels = c("Scenario_1" = "<20 years only", 
                                "Scenario_2" = "20-59 years only", 
                                "Scenario_3" = ">60 years only")) +
  #scale_y_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) + 
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       title = "Vaccine impact at sub-national level (reduction in fatal cases)",
       x = "Week", y = "Predicted cumulative reported fatal cases (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(5, 30, 5, 5)
  ) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = Inf, label = ann), 
            inherit.aes = FALSE,
            hjust = 1.1, 
            vjust = 1.1, 
            size = 3) +
  #geom_text(data = vaccine_ann_df,
  #          aes(x = x, y = y, label = label),
  #          inherit.aes = FALSE,
  #          hjust = 0, vjust = 1, size = 2)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  #geom_point(data = observed_all, aes(x = Week, y = Observed), size = 0.5) + 
  facet_wrap(~region)


combined_prepost_fatal_nat <- combined_prepost_fatal %>%
  group_by(Scenario, Week) %>%
  summarise(
    fatal    = sum(fatal),
    fatal_lo = sum(fatal_lo),
    fatal_hi = sum(fatal_hi),
    pre_fatal = sum(pre_fatal),
    pre_fatal_lo = sum(pre_fatal_lo),
    pre_fatal_hi = sum(pre_fatal_hi),
    tot_pop          = sum(tot_pop),
    pre_cases = sum(pre_cases),
    pre_cases_lo = sum(pre_cases_lo),
    pre_cases_hi = sum(pre_cases_hi),
    post_cases = sum(post_cases),
    post_cases_lo = sum(post_cases_lo),
    post_cases_hi = sum(post_cases_hi),
    .groups = "drop"
  ) %>% arrange(Scenario, Week) %>%
  group_by(Scenario) %>%
  mutate(
    # Now do cumsum within each Scenario in ascending Week order
    cum_fatal         = cumsum(fatal),
    cum_fatal_lo      = cumsum(fatal_lo),
    cum_fatal_hi      = cumsum(fatal_hi),
    cum_fatal_pre     = cumsum(pre_fatal),
    cum_fatal_pre_lo  = cumsum(pre_fatal_lo),
    cum_fatal_pre_hi  = cumsum(pre_fatal_hi)
  ) %>%
  ungroup() %>%
  mutate(
    pre_fatal_prop      = (cum_fatal_pre / pre_cases) * 100,
    pre_fatal_prop_lo   = (cum_fatal_pre_lo / post_cases_lo) * 100,
    pre_fatal_prop_hi   = (cum_fatal_pre_hi / post_cases_hi) * 100,
    post_fatal_prop     = (cum_fatal / post_cases) * 100,
    post_fatal_prop_lo  = (cum_fatal_lo / post_cases_lo) * 100,
    post_fatal_prop_hi  = (cum_fatal_hi / post_cases_hi) * 100
  )

global_impact_fatal_nat <- combined_prepost_fatal_nat %>%
  group_by(Scenario) %>%
  summarise(
    total_post_fatal = sum(fatal),
    total_pre_fatal  = sum(pre_fatal),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_fatal - total_post_fatal,
    impact = diff / total_pre_fatal * 100
  )

global_impact_fatal_nat <- combined_prepost_fatal %>%
  group_by(Scenario) %>%
  summarise(
    total_post_fatal = sum(cum_fatal),  
    total_post_fatal_lo = sum(cum_fatal_lo),
    total_post_fatal_hi = sum(cum_fatal_hi),
    total_pre_fatal   = sum(cum_fatal_pre),
    total_pre_fatal_lo = sum(cum_fatal_pre_lo),
    total_pre_fatal_hi = sum(cum_fatal_pre_hi),
    .groups = "drop"
  ) %>%
  mutate(
    diff_fatal   = `total_pre_fatal` - `total_post_fatal`,
    impact_fatal = diff_fatal / `total_pre_fatal`  * 100,
    diff_fatal_lo = `total_pre_fatal_lo` - `total_post_fatal_lo`,
    impact_fatal_lo = diff_fatal_lo / `total_pre_fatal_lo`  * 100,
    diff_fatal_hi = `total_pre_fatal_hi` - `total_post_fatal_hi`,
    impact_fatal_hi = diff_fatal_hi / `total_pre_fatal_hi`  * 100
  )

annotation_text_fatal <- paste0(
  global_impact_fatal_nat$Scenario, ": ", 
  round(global_impact_fatal_nat$impact_fatal, 1), "%",
  collapse = "\n"
)

cum_fatal_nat <- list(
  fatal_week_df   = combined_prepost_fatal_nat,
  annotation_text_fatal = annotation_text_fatal
)


fatal_graph_nat <- 
cum_fatal_plot(cum_fatal_nat, postsim_all_ag_ui)


# combining
pre_post_fatal_comb <- pre_post_fatal_graph +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Vaccine impact at sub-national level (reduction in fatal cases)")

# Remove legend from epi_graph_all (Scenario A)
fatal_graph_nat <- fatal_graph_nat + theme(legend.position = "none") + 
  ggtitle("Vaccine impact at national level (reduction in fatal cases)")

# Arrange plots side-by-side with a common legend
combined_fatal_graph <- ggarrange(fatal_graph_nat, pre_post_fatal_comb,
                            ncol = 2, nrow = 1,
                            labels = c("C", "D"),
                            common.legend = TRUE,
                            legend = "bottom",
                            align = "none")

combined_fatal_graph <- combined_fatal_graph +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil_fatal.jpg", combined_fatal_graph, width = 12, height = 6, dpi = 1200)


combined_national <- ggarrange(epi_graph_all, fatal_graph_nat, 
                            ncol = 2, nrow = 1,
                            labels = c("A", "B"),
                            common.legend = TRUE,
                            legend = "bottom",
                            align = "none")

combined_subnat <- ggarrange(pre_post_graph_comb, pre_post_fatal_graph, 
                               ncol = 2, nrow = 1,
                               labels = c("C", "D"),
                               common.legend = TRUE,
                               legend = "bottom",
                               align = "none")
ggsave(filename = "02_Outputs/2_1_Figures/combined_national.jpg", combined_national, width = 15, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/combined_subnat.jpg", combined_subnat, width = 15, height = 6, dpi = 1200)


# bar graphs 
subnat_fatal_bar   <- combined_prepost_fatal %>% select("Scenario", "Week", "setting",
                                                        "post_fatal_prop", "post_fatal_prop_lo", "post_fatal_prop_hi",
                                                        "pre_fatal_prop", "pre_fatal_prop_lo", "pre_fatal_prop_hi")
nat_fatal_bar      <- cum_fatal_nat$fatal_week_df
nat_fatal_bar$setting <- "National"

nat_fatal_bar <-  nat_fatal_bar %>% select("Scenario", "Week", "setting",
                                           "post_fatal_prop", "post_fatal_prop_lo", "post_fatal_prop_hi",
                                           "pre_fatal_prop", "pre_fatal_prop_lo", "pre_fatal_prop_hi")

subnat_fatal_bar <- subnat_fatal_bar %>% filter(Week == 52)
subnat_fatal_bar <- subnat_fatal_bar %>% group_by(Scenario, Week, setting) %>%  summarise(
  post_fatal_prop    = mean(post_fatal_prop,    na.rm = TRUE),
  post_fatal_prop_lo = mean(post_fatal_prop_lo, na.rm = TRUE),
  post_fatal_prop_hi = mean(post_fatal_prop_hi, na.rm = TRUE),
  pre_fatal_prop    = mean(pre_fatal_prop,    na.rm = TRUE),
  pre_fatal_prop_lo = mean(pre_fatal_prop_lo, na.rm = TRUE),
  pre_fatal_prop_hi = mean(pre_fatal_prop_hi, na.rm = TRUE),
  .groups = "drop"
)
  
nat_fatal_bar <- nat_fatal_bar %>% filter(Week == 52)
fatal_bar <- rbind(nat_fatal_bar, subnat_fatal_bar)
fatal_bar <- fatal_bar %>%
  mutate(setting = factor(setting, levels = c("National", "High", "Moderate", "Low")))

fatal_bar <- fatal_bar %>%
  mutate(
    reduction_fatal_percent = 100 * (pre_fatal_prop - post_fatal_prop) / pre_fatal_prop,
    reduction_fatal_percent_lo = 100 * (pre_fatal_prop_lo - post_fatal_prop_hi) / pre_fatal_prop_lo,
    reduction_fatal_percent_hi = 100 * (pre_fatal_prop_hi - post_fatal_prop_lo) / pre_fatal_prop_hi
  )

write_xlsx(fatal_bar, path = "02_Outputs/2_2_Tables/fatal_bar.xlsx")


fatal_bar_pre <- fatal_bar %>%
  transmute(
    Scenario,
    Week,
    setting,
    Type = "Without vaccination",
    fatal_prop    = pre_fatal_prop,
    fatal_prop_lo = pre_fatal_prop_lo,
    fatal_prop_hi = pre_fatal_prop_hi
  )

fatal_bar_post <- fatal_bar %>%
  transmute(
    Scenario,
    Week,
    setting,
    Type = "With vaccination",
    fatal_prop    = post_fatal_prop,
    fatal_prop_lo = post_fatal_prop_lo,
    fatal_prop_hi = post_fatal_prop_hi
  )

fatal_bar_long <- bind_rows(fatal_bar_pre, fatal_bar_post)
fatal_bar_long$Type <- factor(
  fatal_bar_long$Type,
  levels = c("Without vaccination", "With vaccination")
)

p <- 
ggplot(fatal_bar_long, 
       aes(x = Type, 
           y = fatal_prop, 
           fill = Type)) +
  geom_col(
    position = position_dodge2(width = 0.8, preserve = "single"), 
    width    = 0.7, 
    alpha    = 0.6
  ) +
  geom_errorbar(
    aes(ymin = fatal_prop_lo, ymax = fatal_prop_hi),
    position = position_dodge2(width = 0.8, preserve = "single"),
    width    = 0.2
  ) +
  facet_grid(
    Scenario ~ setting,
    labeller = labeller(
      Scenario = c("Scenario_1" = "<20 years",
                   "Scenario_2" = "20–59 years",
                   "Scenario_3" = "≥60 years")
    ),
    scales = "free_y",  
    space  = "free_y"
  ) +
  scale_fill_brewer(
    palette = "Set1",
    name    = "Vaccination",
    labels  = c("Without vaccination", "With vaccination")
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x     = NULL,
    y     = "Cumulative fatal cases (%)"
  ) +
  theme_pubclean(base_size = 11) +
  theme(
    axis.text.x  = element_text(hjust = 1),
    legend.position = "bottom"
  )

p <- 
ggarrange(p, 
          ncol = 1, nrow = 1,
          labels = c("C"),
          common.legend = TRUE,
          legend = "bottom",
          align = "v")
ggsave(filename = "02_Outputs/2_1_Figures/fig_cumfatal.jpg", p, width = 10, height = 5, dpi = 1200)

# contour map for VE ----------------------------------------------------------
contour_ce <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_ce,
  posterior            = posterior_ce,
  N                    = N_ceara$Ceará,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Ceará",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ce
)

contour_bh <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_bh,
  posterior            = posterior_bh,
  N                    = N_bahia$Bahia,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Bahia",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_bh
)

contour_pn <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_pn,
  posterior            = posterior_pn,
  N                    = N_pemam$Pernambuco,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Pernambuco",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pn
)

contour_rg <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_rg,
  posterior            = posterior_rg,
  N                    = N_rg$`Rio Grande do Norte`,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Rio Grande do Norte",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_rg
)

contour_pa <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_pa,
  posterior            = posterior_pa,
  N                    = N_pa$Paraíba,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Paraíba",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pa
)

contour_pi <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_pi,
  posterior            = posterior_pi,
  N                    = N_pi$Piauí,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Piauí",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pi
)

contour_tc <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_tc,
  posterior            = posterior_tc,
  N                    = N_tc$Tocantins,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Tocantins",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_tc
)

contour_ag <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_ag,
  posterior            = posterior_ag,
  N                    = N_ag$Alagoas,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Alagoas",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ag
)

contour_mg <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_mg,
  posterior            = posterior_mg,
  N                    = N_mg$`Minas Gerais`,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Minas Gerais",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_mg
)

contour_se <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_se,
  posterior            = posterior_se,
  N                    = N_se$Sergipe,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Sergipe",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_se
)

contour_go <- ve_contour(
  target_age_list      = target_age_list,
  observed             = observed_go,
  posterior            = posterior_go,
  N                    = N_go$Goiás,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Goiás",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_go
)

combined_contour <- bind_rows(
  contour_ag,
  contour_bh,
  contour_ce,
  contour_mg,
  contour_pa,
  contour_pi,
  contour_pn,
  contour_rg,
  contour_tc,
  contour_se,
  contour_go,
) %>%
  mutate(
    # 1) 직접효과 = 시나리오에 대응하는 agegr의 diff_cases
    direct = case_when(
      Scenario == "Scenario_1" ~ diff_cases_agegr1,
      Scenario == "Scenario_2" ~ diff_cases_agegr2,
      Scenario == "Scenario_3" ~ diff_cases_agegr3,
      TRUE                     ~ NA_real_
    ),
    # 2) 간접효과 = 나머지 두 개 agegr diff_cases 합
    indirect = case_when(
      Scenario == "Scenario_1" ~ diff_cases_agegr2 + diff_cases_agegr3,
      Scenario == "Scenario_2" ~ diff_cases_agegr1 + diff_cases_agegr3,
      Scenario == "Scenario_3" ~ diff_cases_agegr1 + diff_cases_agegr2,
      TRUE                     ~ NA_real_
    ),
    # 3) 간접/직접 비율
    ind_d_ratio = indirect / direct,
    indirect_frac = indirect / (direct + indirect)
  )

contour_summ <- combined_contour %>%
  mutate(
    scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
  )

contour_summ <- contour_summ %>%
  left_join(vacc_wide, by = "region") %>%
  select(-starts_with("impact_")) 

collapsed_contour <- contour_summ %>%
  group_by(Scenario, VE_inf, VE_block) %>%
  summarise(
    across(
      .cols = where(is.numeric),  
      .fns = mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

collapsed_contour <- collapsed_contour %>%
  mutate(
    total_vacc = case_when(
      Scenario == "Scenario_1" ~ `totvacc_Scenario_1`,
      Scenario == "Scenario_2" ~ `totvacc_Scenario_2`,
      Scenario == "Scenario_3" ~ `totvacc_Scenario_3`,
      TRUE                     ~ NA_real_
    ),
    tot_diff_cases     = diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3,
    impact_cases_total = tot_diff_cases / total_vacc * 1e6
  )%>%
  mutate(
    # 1) 연령대별 pre/post 합쳐서 전체 사례수
    total_pre_cases  = tot_pre_cases_agegr1 +
      tot_pre_cases_agegr2 +
      tot_pre_cases_agegr3,
    total_post_cases = tot_post_cases_agegr1 +
      tot_post_cases_agegr2 +
      tot_post_cases_agegr3,
    # 2) (pre − post) / pre 로 “전체 케이스 감소율”
    reduction_rate   = (total_pre_cases - total_post_cases) / total_pre_cases
  )


collapsed_contour <- collapsed_contour %>% 
  ## 1. 시나리오별 총 접종량 지정 ----
mutate(
  total_vacc = case_when(
    Scenario == "Scenario_1" ~ totvacc_Scenario_1,
    Scenario == "Scenario_2" ~ totvacc_Scenario_2,
    Scenario == "Scenario_3" ~ totvacc_Scenario_3,
    TRUE                     ~ NA_real_
  )
) %>% 
  
  ## 2. 연령대 합계(직접+간접) ----
mutate(
  tot_diff_cases = diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3,
  tot_diff_fatal = diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3,
  tot_diff_daly  = diff_daly_agegr1  + diff_daly_agegr2  + diff_daly_agegr3,
  
  total_pre_cases  = tot_pre_cases_agegr1 + tot_pre_cases_agegr2 + tot_pre_cases_agegr3,
  total_post_cases = tot_post_cases_agegr1 + tot_post_cases_agegr2 + tot_post_cases_agegr3,
  
  total_pre_fatal  = tot_pre_fatal_agegr1 + tot_pre_fatal_agegr2 + tot_pre_fatal_agegr3,
  total_post_fatal = tot_post_fatal_agegr1 + tot_post_fatal_agegr2 + tot_post_fatal_agegr3,
  
  total_pre_daly   = tot_pre_daly_agegr1  + tot_pre_daly_agegr2  + tot_pre_daly_agegr3,
  total_post_daly  = tot_post_daly_agegr1 + tot_post_daly_agegr2 + tot_post_daly_agegr3
) %>% 
  
  ## 3. 100만 도즈당 효과 및 감소율 ----
mutate(
  impact_cases_total = tot_diff_cases / total_vacc * 1e5,
  impact_fatal_total = tot_diff_fatal / total_vacc * 1e5,
  impact_daly_total  = tot_diff_daly  / total_vacc * 1e5,
  
  reduction_cases = (total_pre_cases - total_post_cases) / total_pre_cases,
  reduction_fatal = (total_pre_fatal - total_post_fatal) / total_pre_fatal,
  reduction_daly  = (total_pre_daly  - total_post_daly)  / total_pre_daly
)

# contour plots  for indirect frac
ggplot(collapsed_contour, aes(x = VE_inf, y = VE_block, z = indirect_frac)) +
  geom_contour_filled(breaks = seq(0, 1, by = 0.1), show.legend = TRUE) +
  # → discrete viridis scale
  scale_fill_viridis_d(
    option = "magma",
    name   = "Indirect\nfraction",
    labels = function(breaks) {
      mid <- (as.numeric(sub("\\((.+),.*", "\\1", breaks)) +
                as.numeric(sub(".*,(.+)\\]", "\\1", breaks))) / 2
      percent(mid, accuracy = 1)
    }
  ) +
  labs(
    x     = "Infection-blocking VE",
    y     = "Disease-blocking VE",
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  )+
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    aspect.ratio    = 1
  ) +
  facet_wrap(~ Scenario, labeller = labeller(
    Scenario = c(
      Scenario_1 = "<20 years target",
      Scenario_2 = "20–59 years target",
      Scenario_3 = "≥60 years target"
    )
  ))

# contour plot for reduction rate-----------------------------------------------
d <- 
ggplot(collapsed_contour, aes(x = VE_inf, y = VE_block, z = impact_cases_total)) +
  geom_contour_filled() +
  facet_wrap(~ Scenario)  +
  theme_minimal()+
  # → discrete viridis scale
  scale_fill_viridis_d(
    option = "magma",
    name   = "Cases averted per 100,000 doses",
    labels = function(breaks) {
      mid <- (as.numeric(sub("\\((.+),.*", "\\1", breaks)) +
                as.numeric(sub(".*,(.+)\\]", "\\1", breaks))) / 2
    }
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    aspect.ratio    = 1
  ) +
  facet_wrap(~ Scenario, labeller = labeller(
    Scenario = c(
      Scenario_1 = "<20 years target",
      Scenario_2 = "20–59 years target",
      Scenario_3 = "≥60 years target"
    ),
    drop = FALSE
  ))+
  labs(
    x     = "Infection-blocking VE",
    y     = "Disease-blocking VE",
  ) + scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    #labels = percent_format(accuracy = 1),
    labels = function(x) x * 100,
    expand = c(0, 0)
  )

e <- 
ggplot(collapsed_contour, aes(x = VE_inf, y = VE_block, z = impact_fatal_total)) +
  geom_contour_filled() +
  facet_wrap(~ Scenario)  +
  theme_minimal()+
  # → discrete viridis scale
  scale_fill_viridis_d(
    option = "magma",
    name   = "Deaths averted per 100,000 doses",
    labels = function(breaks) {
      mid <- (as.numeric(sub("\\((.+),.*", "\\1", breaks)) +
                as.numeric(sub(".*,(.+)\\]", "\\1", breaks))) / 2
    }
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    aspect.ratio    = 1
  ) +
  facet_wrap(~ Scenario, labeller = labeller(
    Scenario = c(
      Scenario_1 = "<20 years target",
      Scenario_2 = "20–59 years target",
      Scenario_3 = "≥60 years target"
    ),
    drop = FALSE
  ))+
  labs(
    x     = "Infection-blocking VE",
    y     = "Disease-blocking VE",
  ) + scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    #labels = percent_format(accuracy = 1),
    labels = function(x) x * 100,
    expand = c(0, 0)
  )

f <- 
ggplot(collapsed_contour, aes(x = VE_inf, y = VE_block, z = impact_daly_total)) +
  geom_contour_filled() +
  facet_wrap(~ Scenario)  +
  theme_minimal()+
  # → discrete viridis scale
  scale_fill_viridis_d(
    option = "magma",
    name   = "DALYs averted per 100,000 doses",
    labels = function(breaks) {
      mid <- (as.numeric(sub("\\((.+),.*", "\\1", breaks)) +
                as.numeric(sub(".*,(.+)\\]", "\\1", breaks))) / 2
    }
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    aspect.ratio    = 1
  ) +
  facet_wrap(~ Scenario, labeller = labeller(
    Scenario = c(
      Scenario_1 = "<20 years target",
      Scenario_2 = "20–59 years target",
      Scenario_3 = "≥60 years target"
    ),
    drop = FALSE
  ))+
  labs(
    x     = "Infection-blocking VE",
    y     = "Disease-blocking VE",
  ) + scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    #labels = percent_format(accuracy = 0.01),
    #labels = scales::number_format(accuracy = 1),
    labels = function(x) x * 100,
    expand = c(0, 0)
  )

d <- d + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

e <- e + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    strip.text      = element_blank(),      
    strip.background = element_blank()
  )

f <- f +
  theme(
    strip.text      = element_blank(),
    strip.background = element_blank()
  )

d <- d + 
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank()
  )

e <- e +
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c(
               Scenario_1 = "",
               Scenario_2 = "",
               Scenario_3 = ""
             ))
  )
f <- f +
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c(
               Scenario_1 = "",
               Scenario_2 = "",
               Scenario_3 = ""
             ))
  )
d <- d + labs(tag = "D")
e <- e + labs(tag = "E")
f <- f + labs(tag = "F")

combined_sens_g_contour <- 
  (d / e / f) + plot_layout(ncol = 1, heights = c(1, 1, 1))

ggsave(filename = "02_Outputs/2_1_Figures/combined_sens_g_contour.jpg", combined_sens_g_contour, width = 10, height = 8, dpi = 1200)


###
pcrit <- postsim_bh_ui[[1]][["sim_out"]][["pcrit_vec"]]
c_eff   <- postsim_bh_ui[[1]]$sim_out$coverage_target * 0.989
plot(pcrit, type="l", ylim=c(0,1), ylab="Coverage", xlab="Week")
lines(c_eff, col=2, lwd=2)
legend("topright", c("p_crit(t)", "actual coverage×VE"), lty=1, col=c(1,2))
