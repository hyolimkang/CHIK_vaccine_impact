load("00_Data/0_2_Processed/fit_prevacc_bh.RData")
load("00_Data/0_2_Processed/fit_prevacc_ce.RData")
load("00_Data/0_2_Processed/fit_prevacc_mg.RData")
load("00_Data/0_2_Processed/fit_prevacc_pn.RData")
load("00_Data/0_2_Processed/fit_prevacc_pa.RData")
load("00_Data/0_2_Processed/fit_prevacc_rg.RData")
load("00_Data/0_2_Processed/fit_prevacc_pi.RData")
load("00_Data/0_2_Processed/fit_prevacc_ag.RData")
load("00_Data/0_2_Processed/fit_prevacc_tc.RData")
load("00_Data/0_2_Processed/brazil_pop_transposed.RData")
load("00_Data/0_2_Processed/observed_2022.RData")
load("00_Data/0_2_Processed/bra_foi_state_summ.RData")
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")
load("00_Data/0_2_Processed/combined_contour.RData")
load("00_Data/0_2_Processed/combined_heatmap.RData")
load("00_Data/0_2_Processed/tornado_results_all.RData")

lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")
le_sample <- readRDS("00_Data/0_2_Processed/le_sample.RDS")
source("01_Script/1_1_Functions/sim_function_final.R")
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

posterior_ce <- param_ce$posterior_prevacc
posterior_bh <- param_bh$posterior_prevacc
posterior_pa <- param_pa$posterior_prevacc
posterior_pn <- param_pn$posterior_prevacc
posterior_rg <- param_rg$posterior_prevacc
posterior_pi <- param_pi$posterior_prevacc
posterior_ag <- param_ag$posterior_prevacc
posterior_tc <- param_tc$posterior_prevacc
posterior_mg <- param_mg$posterior_prevacc

# calculate national total supply 
N_bahia <- brazil_pop_transposed[,5]
N_ceara <- brazil_pop_transposed[,6]
N_mg    <- brazil_pop_transposed[,13]
N_pemam <- brazil_pop_transposed[,17]
N_pa    <- brazil_pop_transposed[,15]
N_rg    <- brazil_pop_transposed[,19]
N_pi    <- brazil_pop_transposed[,18]
N_ag    <- brazil_pop_transposed[,2]
N_tc    <- brazil_pop_transposed[,27]

# scenarios
target_age_list <- list(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)) # >60 yrs old) # <20 yrs old

# prevacc 95%UI simulation ----------------------------------------------------

## pre-processing ----------------------------------------------------
preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups, 
                            N = N_ceara$Ceará, 
                            region = "Ceará",
                            observed = observed_ce
                            )

preui_bh <- simulate_pre_ui_age(posterior = posterior_bh, bra_foi_state_summ, age_groups, 
                                N = N_bahia$Bahia, 
                                region = "Bahia",
                                observed = observed_bh
)

preui_pa <- simulate_pre_ui_age(posterior = posterior_pa, bra_foi_state_summ, age_groups, 
                                N = N_pa$Paraíba, 
                                region = "Paraíba",
                                observed = observed_pa
)

preui_pn <- simulate_pre_ui_age(posterior = posterior_pn, bra_foi_state_summ, age_groups, 
                                N = N_pemam$Pernambuco, 
                                region = "Pernambuco",
                                observed = observed_pn
)

preui_rg <- simulate_pre_ui_age(posterior = posterior_rg, bra_foi_state_summ, age_groups, 
                                N = N_rg$`Rio Grande do Norte`, 
                                region = "Rio Grande do Norte",
                                observed = observed_rg
)

preui_pi <- simulate_pre_ui_age(posterior = posterior_pi, bra_foi_state_summ, age_groups, 
                                N = N_pi$Piauí, 
                                region = "Piauí",
                                observed = observed_pi
)

preui_tc <- simulate_pre_ui_age(posterior = posterior_tc, bra_foi_state_summ, age_groups, 
                                N = N_tc$Tocantins, 
                                region = "Tocantins",
                                observed = observed_tc
)

preui_ag <- simulate_pre_ui_age(posterior = posterior_ag, bra_foi_state_summ, age_groups, 
                                N = N_ag$Alagoas, 
                                region = "Alagoas",
                                observed = observed_ag
)

preui_mg <- simulate_pre_ui_age(posterior = posterior_mg, bra_foi_state_summ, age_groups, 
                                N = N_mg$`Minas Gerais`, 
                                region = "Minas Gerais",
                                observed = observed_mg
)

## summarise ----------------------------------------------------
pre_results_ce_ui <- summarise_presim_ui(sim_result        = preui_ce, 
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

# postsim 95% UI run ----------------------------------------------------------
postsim_ce_ui <- run_simulation_scenarios_ui(
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
  posterior         = posterior_ce
)

postsim_bh_ui <- run_simulation_scenarios_ui(
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
  posterior         = posterior_bh
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
  posterior         = posterior_pa
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
  posterior         = posterior_pn
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
  posterior         = posterior_rg
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
  posterior         = posterior_pi
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
  posterior         = posterior_tc
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
  posterior         = posterior_ag
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
  posterior         = posterior_mg
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

combined_prepost_case <- bind_rows(
  postsim_all_ce_ui$summary_week_df,
  postsim_all_bh_ui$summary_week_df,
  postsim_all_pa_ui$summary_week_df,
  postsim_all_pn_ui$summary_week_df,
  postsim_all_rg_ui$summary_week_df,
  postsim_all_pi_ui$summary_week_df,
  postsim_all_tc_ui$summary_week_df,
  postsim_all_ag_ui$summary_week_df,
  postsim_all_mg_ui$summary_week_df
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
  group_by(region) %>%
  summarise(ann = paste0(Scenario, ": ", round(impact, 1), "%", collapse = "\n"),
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
  geom_ribbon(data = combined_prepost_case %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 2/3 * max_cases, ymax = max_cases, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_case %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
              aes(x = Week, ymin = 1/3 * max_cases, ymax = 2/3 * max_cases, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_case %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
              aes(x = Week, ymin = 0, ymax = 1/3 * max_cases, fill = Scenario),
              alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
              alpha = 0.3) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
              fill = "gray", alpha = 0.6) +
  # Post-vaccination case line (colored by scenario)
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "With vaccination"), size = 0.2) +
  # Pre-vaccination case line (dashed, in black)
  geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Without vaccination"), color = "black", size = 0.3) +
  scale_linetype_manual(values = c("With vaccination" = "solid", "Without vaccination" = "dashed")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Scenario_1" = "<20 years only", 
                               "Scenario_2" = "20-59 years only", 
                               "Scenario_3" = ">60 years only")) +
  scale_color_brewer(palette = "Set1",
                     labels = c("Scenario_1" = "<20 years only", 
                                "Scenario_2" = "20-59 years only", 
                                "Scenario_3" = ">60 years only")) +
  scale_y_continuous(labels = scales::comma) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
       x = "Week", y = "Predicted symptomatic reported cases") +
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
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  geom_point(data = observed_all, aes(x = Week, y = Observed),size = 0.00000005) + 
  facet_wrap(~region)

ggsave(filename = "02_Outputs/2_1_Figures/pre_post_region.jpg", p, width = 15, height = 7, dpi = 1200)


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

combined_vacc_alloc <- bind_rows(
  vacc_alloc_ce$weekly_allocation_long,
  vacc_alloc_bh$weekly_allocation_long,
  vacc_alloc_pa$weekly_allocation_long,
  vacc_alloc_pn$weekly_allocation_long,
  vacc_alloc_rg$weekly_allocation_long,
  vacc_alloc_pi$weekly_allocation_long,
  vacc_alloc_tc$weekly_allocation_long,
  vacc_alloc_ag$weekly_allocation_long,
  vacc_alloc_mg$weekly_allocation_long
)

combined_vacc_alloc_summ <- combined_vacc_alloc %>%
  group_by(Scenario, AgeGroup, age_gr, Week) %>%
  summarise(
    Vaccinated = sum(Vaccinated, na.rm = TRUE),
    tot_sum = sum(tot_sum, na.rm = TRUE)
  ) %>%
  ungroup()

paired <- brewer.pal(12, "Paired")
extra <- brewer.pal(8, "Dark2")[1:6]
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

p <- 
ggplot(combined_vacc_alloc, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=comma)+
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c("Scenario1" = "<20 years only", 
                                              "Scenario2" = "20-59 years only", 
                                              "Scenario3" = ">60 years only")))  +         # Facet by Scenario if you want separate plots for each
  theme_light()+
  facet_wrap(~region)


ggsave(filename = "02_Outputs/2_1_Figures/fig_vacc_alloc_tot.jpg", p, width = 12, height = 7, dpi = 1200)
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

# Combine the two data frames
combined_df <- bind_rows(df_ce, df_bh, df_pa,
                         df_pn, df_rg, df_pi,
                         df_tc, df_ag, df_mg) %>%
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
  global_impact_week$Scenario, ": ",
  round(global_impact_week$impact, 1), "%",
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
  ) + 
  ggtitle("Vaccine impact at sub-national level (reduction in symptomatic cases)")

# Remove legend from epi_graph_all (Scenario A)
epi_graph_all <- epi_graph_all + theme(legend.position = "none") + 
                 ggtitle("Vaccine impact at national level (reduction in symptomatic cases)")

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

ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil.jpg", combined_graph, width = 12, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil.pdf", combined_graph, width = 23, height = 10)
## total nnv graph (national level) ------------------------------
n_scenarios = 3

df_nnv_ce <- nnv_ce$final_summ_df
df_nnv_bh <- nnv_bh$final_summ_df
df_nnv_pa <- nnv_pa$final_summ_df
df_nnv_pn <- nnv_pn$final_summ_df
df_nnv_rg <- nnv_rg$final_summ_df
df_nnv_pi <- nnv_pi$final_summ_df
df_nnv_tc <- nnv_tc$final_summ_df
df_nnv_ag <- nnv_ag$final_summ_df
df_nnv_mg <- nnv_mg$final_summ_df

df_nnv_ce$setting <- "high"
df_nnv_bh$setting <- "moderate"
df_nnv_pa$setting <- "moderate"
df_nnv_pn$setting <- "high"
df_nnv_rg$setting <- "moderate"
df_nnv_pi$setting <- "moderate"
df_nnv_tc$setting <- "low"
df_nnv_ag$setting <- "moderate"
df_nnv_mg$setting <- "moderate"


df_nnv_ce$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_ce$age_gr <- factor(df_nnv_ce$age_gr, levels = age_gr_levels)

df_nnv_bh$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_bh$age_gr <- factor(df_nnv_bh$age_gr, levels = age_gr_levels)

df_nnv_pa$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_pa$age_gr <- factor(df_nnv_pa$age_gr, levels = age_gr_levels)

df_nnv_pn$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_pn$age_gr <- factor(df_nnv_pn$age_gr, levels = age_gr_levels)

df_nnv_rg$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_rg$age_gr <- factor(df_nnv_rg$age_gr, levels = age_gr_levels)

df_nnv_pi$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_pi$age_gr <- factor(df_nnv_pi$age_gr, levels = age_gr_levels)

df_nnv_tc$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_tc$age_gr <- factor(df_nnv_tc$age_gr, levels = age_gr_levels)

df_nnv_ag$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_ag$age_gr <- factor(df_nnv_ag$age_gr, levels = age_gr_levels)

df_nnv_mg$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_mg$age_gr <- factor(df_nnv_mg$age_gr, levels = age_gr_levels)

combined_nnv_df_region <- bind_rows(df_nnv_ce, df_nnv_bh, df_nnv_pa,
                                    df_nnv_pn, df_nnv_rg, df_nnv_pi,
                                    df_nnv_tc, df_nnv_ag, df_nnv_mg) 


combined_nnv_region_all <- combined_nnv_df_region %>%
  group_by(scenario, region) %>%
  summarise(
    Median            = sum(Median, na.rm = TRUE),
    low95             = sum(low95, na.rm = TRUE),
    hi95              = sum(hi95, na.rm = TRUE),
    hospitalised      = sum(hospitalised, na.rm = TRUE),
    hospitalised_lo   = sum(hospitalised_lo, na.rm = TRUE),
    hospitalised_hi   = sum(hospitalised_hi, na.rm = TRUE),
    yld_acute         = sum(yld_acute, na.rm = TRUE),
    yld_acute_lo      = sum(yld_acute_lo, na.rm = TRUE),
    yld_acute_hi      = sum(yld_acute_hi, na.rm = TRUE),
    yld_subacute      = sum(yld_subacute, na.rm = TRUE),
    yld_subacute_lo   = sum(yld_subacute_lo, na.rm = TRUE),
    yld_subacute_hi   = sum(yld_subacute_hi, na.rm = TRUE),
    yld_chronic       = sum(yld_chronic, na.rm = TRUE),
    yld_chronic_lo    = sum(yld_chronic_lo, na.rm = TRUE),
    yld_chronic_hi    = sum(yld_chronic_hi, na.rm = TRUE),
    yld_total         = sum(yld_total, na.rm = TRUE),
    yld_total_lo      = sum(yld_total_lo, na.rm = TRUE),
    yld_total_hi      = sum(yld_total_hi, na.rm = TRUE),
    yll               = sum(yll, na.rm = TRUE),
    yll_lo            = sum(yll_lo, na.rm = TRUE),
    yll_hi            = sum(yll_hi, na.rm = TRUE),
    daly_tot          = sum(daly_tot, na.rm = TRUE),
    daly_tot_lo       = sum(daly_tot_lo, na.rm = TRUE),
    daly_tot_hi       = sum(daly_tot_hi, na.rm = TRUE),
    fatal             = sum(fatal, na.rm = TRUE),
    fatal_lo          = sum(fatal_lo, na.rm = TRUE),
    fatal_hi          = sum(fatal_hi, na.rm = TRUE),
    pre_vacc          = sum(pre_vacc, na.rm = TRUE),
    pre_vacc_low95    = sum(pre_vacc_low95, na.rm = TRUE),
    pre_vacc_hi95     = sum(pre_vacc_hi95, na.rm = TRUE),
    pre_yld_acute     = sum(pre_yld_acute, na.rm = TRUE),
    pre_yld_acute_low95 = sum(pre_yld_acute_low95, na.rm = TRUE),
    pre_yld_acute_hi  = sum(pre_yld_acute_hi, na.rm = TRUE),
    pre_yld_subacute  = sum(pre_yld_subacute, na.rm = TRUE),
    pre_yld_subacute_low95 = sum(pre_yld_subacute_low95, na.rm = TRUE),
    pre_yld_subacute_hi = sum(pre_yld_subacute_hi, na.rm = TRUE),
    pre_yld_chronic   = sum(pre_yld_chronic, na.rm = TRUE),
    pre_yld_chronic_low95 = sum(pre_yld_chronic_low95, na.rm = TRUE),
    pre_yld_chronic_hi = sum(pre_yld_chronic_hi, na.rm = TRUE),
    pre_yld_total     = sum(pre_yld_total, na.rm = TRUE),
    pre_yld_total_low95 = sum(pre_yld_total_low95, na.rm = TRUE),
    pre_yld_total_hi  = sum(pre_yld_total_hi, na.rm = TRUE),
    pre_yll           = sum(pre_yll, na.rm = TRUE),
    pre_yll_low95     = sum(pre_yll_low95, na.rm = TRUE),
    pre_yll_hi        = sum(pre_yll_hi, na.rm = TRUE),
    pre_daly          = sum(pre_daly, na.rm = TRUE),
    pre_daly_low95    = sum(pre_daly_low95, na.rm = TRUE),
    pre_daly_hi       = sum(pre_daly_hi, na.rm = TRUE),
    pre_fatal         = sum(pre_fatal, na.rm = TRUE),
    pre_fatal_low95   = sum(pre_fatal_low95, na.rm = TRUE),
    pre_fatal_hi      = sum(pre_fatal_hi, na.rm = TRUE),
    tot_vacc          = sum(tot_vacc, na.rm = TRUE)
  ) %>% 
  ungroup()


write_xlsx(combined_nnv_region_all, path = "02_Outputs/2_2_Tables/combined_nnv_region_all.xlsx")
write_xlsx(combined_nnv_df_region, path = "02_Outputs/2_2_Tables/combined_nnv_df_region.xlsx")

combined_nnv_df <- combined_nnv_df_region %>%
  group_by(scenario, AgeGroup) %>%
  summarise(
    Median            = sum(Median, na.rm = TRUE),
    low95             = sum(low95, na.rm = TRUE),
    hi95              = sum(hi95, na.rm = TRUE),
    hospitalised      = sum(hospitalised, na.rm = TRUE),
    hospitalised_lo   = sum(hospitalised_lo, na.rm = TRUE),
    hospitalised_hi   = sum(hospitalised_hi, na.rm = TRUE),
    yld_acute         = sum(yld_acute, na.rm = TRUE),
    yld_acute_lo      = sum(yld_acute_lo, na.rm = TRUE),
    yld_acute_hi      = sum(yld_acute_hi, na.rm = TRUE),
    yld_subacute      = sum(yld_subacute, na.rm = TRUE),
    yld_subacute_lo   = sum(yld_subacute_lo, na.rm = TRUE),
    yld_subacute_hi   = sum(yld_subacute_hi, na.rm = TRUE),
    yld_chronic       = sum(yld_chronic, na.rm = TRUE),
    yld_chronic_lo    = sum(yld_chronic_lo, na.rm = TRUE),
    yld_chronic_hi    = sum(yld_chronic_hi, na.rm = TRUE),
    yld_total         = sum(yld_total, na.rm = TRUE),
    yld_total_lo      = sum(yld_total_lo, na.rm = TRUE),
    yld_total_hi      = sum(yld_total_hi, na.rm = TRUE),
    yll               = sum(yll, na.rm = TRUE),
    yll_lo            = sum(yll_lo, na.rm = TRUE),
    yll_hi            = sum(yll_hi, na.rm = TRUE),
    daly_tot          = sum(daly_tot, na.rm = TRUE),
    daly_tot_lo       = sum(daly_tot_lo, na.rm = TRUE),
    daly_tot_hi       = sum(daly_tot_hi, na.rm = TRUE),
    fatal             = sum(fatal, na.rm = TRUE),
    fatal_lo          = sum(fatal_lo, na.rm = TRUE),
    fatal_hi          = sum(fatal_hi, na.rm = TRUE),
    pre_vacc          = sum(pre_vacc, na.rm = TRUE),
    pre_vacc_low95    = sum(pre_vacc_low95, na.rm = TRUE),
    pre_vacc_hi95     = sum(pre_vacc_hi95, na.rm = TRUE),
    pre_yld_acute     = sum(pre_yld_acute, na.rm = TRUE),
    pre_yld_acute_low95 = sum(pre_yld_acute_low95, na.rm = TRUE),
    pre_yld_acute_hi  = sum(pre_yld_acute_hi, na.rm = TRUE),
    pre_yld_subacute  = sum(pre_yld_subacute, na.rm = TRUE),
    pre_yld_subacute_low95 = sum(pre_yld_subacute_low95, na.rm = TRUE),
    pre_yld_subacute_hi = sum(pre_yld_subacute_hi, na.rm = TRUE),
    pre_yld_chronic   = sum(pre_yld_chronic, na.rm = TRUE),
    pre_yld_chronic_low95 = sum(pre_yld_chronic_low95, na.rm = TRUE),
    pre_yld_chronic_hi = sum(pre_yld_chronic_hi, na.rm = TRUE),
    pre_yld_total     = sum(pre_yld_total, na.rm = TRUE),
    pre_yld_total_low95 = sum(pre_yld_total_low95, na.rm = TRUE),
    pre_yld_total_hi  = sum(pre_yld_total_hi, na.rm = TRUE),
    pre_yll           = sum(pre_yll, na.rm = TRUE),
    pre_yll_low95     = sum(pre_yll_low95, na.rm = TRUE),
    pre_yll_hi        = sum(pre_yll_hi, na.rm = TRUE),
    pre_daly          = sum(pre_daly, na.rm = TRUE),
    pre_daly_low95    = sum(pre_daly_low95, na.rm = TRUE),
    pre_daly_hi       = sum(pre_daly_hi, na.rm = TRUE),
    pre_fatal         = sum(pre_fatal, na.rm = TRUE),
    pre_fatal_low95   = sum(pre_fatal_low95, na.rm = TRUE),
    pre_fatal_hi      = sum(pre_fatal_hi, na.rm = TRUE),
    tot_vacc          = sum(tot_vacc, na.rm = TRUE)
  ) %>% 
  ungroup()

# Next, compute differences, impacts, and NNV values (with UI) for each metric.
final_nnv_df <- combined_nnv_df %>%
  mutate(
    # For Cases:
    diff           = pre_vacc - Median,
    diff_low       = pre_vacc_low95 - hi95,
    diff_hi        = pre_vacc_hi95 - low95,
    impact         = diff / pre_vacc * 100,
    impact_low     = diff_low / pre_vacc_hi95 * 100,
    impact_hi      = diff_hi / pre_vacc_low95 * 100,
    nnv            = tot_vacc / diff,
    nnv_lo         = tot_vacc / diff_hi,
    nnv_hi         = tot_vacc / diff_low,
    
    # For Fatal outcomes:
    diff_fatal     = pre_fatal - fatal,
    diff_fatal_low = pre_fatal_low95 - fatal_hi,
    diff_fatal_hi  = pre_fatal_hi - fatal_lo,
    impact_fatal   = diff_fatal / pre_fatal * 100,
    impact_fatal_low = diff_fatal_low / pre_fatal_hi * 100,
    impact_fatal_hi  = diff_fatal_hi / pre_fatal_low95 * 100,
    nnv_fatal      = tot_vacc / diff_fatal,
    nnv_fatal_lo   = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi   = tot_vacc / diff_fatal_low,
    
    # For YLD Acute:
    diff_yld_acute       = pre_yld_acute - yld_acute,
    diff_yld_acute_low   = pre_yld_acute_low95 - yld_acute_hi,
    diff_yld_acute_hi    = pre_yld_acute_hi - yld_acute_lo,
    impact_yld_acute     = diff_yld_acute / pre_yld_acute * 100,
    impact_yld_acute_low = diff_yld_acute_low / pre_yld_acute_hi * 100,
    impact_yld_acute_hi  = diff_yld_acute_hi / pre_yld_acute_low95 * 100,
    nnv_yld_acute        = tot_vacc / diff_yld_acute,
    nnv_yld_acute_lo     = tot_vacc / diff_yld_acute_hi,
    nnv_yld_acute_hi     = tot_vacc / diff_yld_acute_low,
    
    # For YLD Subacute:
    diff_yld_subacute       = pre_yld_subacute - yld_subacute,
    diff_yld_subacute_low   = pre_yld_subacute_low95 - yld_subacute_hi,
    diff_yld_subacute_hi    = pre_yld_subacute_hi - yld_subacute_lo,
    impact_yld_subacute     = diff_yld_subacute / pre_yld_subacute * 100,
    impact_yld_subacute_low = diff_yld_subacute_low / pre_yld_subacute_hi * 100,
    impact_yld_subacute_hi  = diff_yld_subacute_hi / pre_yld_subacute_low95 * 100,
    nnv_yld_subacute        = tot_vacc / diff_yld_subacute,
    nnv_yld_subacute_lo     = tot_vacc / diff_yld_subacute_hi,
    nnv_yld_subacute_hi     = tot_vacc / diff_yld_subacute_low,
    
    # For YLD Chronic:
    diff_yld_chronic       = pre_yld_chronic - yld_chronic,
    diff_yld_chronic_low   = pre_yld_chronic_low95 - yld_chronic_hi,
    diff_yld_chronic_hi    = pre_yld_chronic_hi - yld_chronic_lo,
    impact_yld_chronic     = diff_yld_chronic / pre_yld_chronic * 100,
    impact_yld_chronic_low = diff_yld_chronic_low / pre_yld_chronic_hi * 100,
    impact_yld_chronic_hi  = diff_yld_chronic_hi / pre_yld_chronic_low95 * 100,
    nnv_yld_chronic        = tot_vacc / diff_yld_chronic,
    nnv_yld_chronic_lo     = tot_vacc / diff_yld_chronic_hi,
    nnv_yld_chronic_hi     = tot_vacc / diff_yld_chronic_low,
    
    # For YLD Total:
    diff_yld_total       = pre_yld_total - yld_total,
    diff_yld_total_low   = pre_yld_total_low95 - yld_total_hi,
    diff_yld_total_hi    = pre_yld_total_hi - yld_total_lo,
    impact_yld_total     = diff_yld_total / pre_yld_total * 100,
    impact_yld_total_low = diff_yld_total_low / pre_yld_total_hi * 100,
    impact_yld_total_hi  = diff_yld_total_hi / pre_yld_total_low95 * 100,
    nnv_yld_total        = tot_vacc / diff_yld_total,
    nnv_yld_total_lo     = tot_vacc / diff_yld_total_hi,
    nnv_yld_total_hi     = tot_vacc / diff_yld_total_low,
    
    # For YLL:
    diff_yll       = pre_yll - yll,
    diff_yll_low   = pre_yll_low95 - yll_hi,
    diff_yll_hi    = pre_yll_hi - yll_lo,
    impact_yll     = diff_yll / pre_yll * 100,
    impact_yll_low = diff_yll_low / pre_yll_hi * 100,
    impact_yll_hi  = diff_yll_hi / pre_yll_low95 * 100,
    nnv_yll        = tot_vacc / diff_yll,
    nnv_yll_lo     = tot_vacc / diff_yll_hi,
    nnv_yll_hi     = tot_vacc / diff_yll_low,
    
    # For DALY:
    diff_daly       = pre_daly - daly_tot,
    diff_daly_low   = pre_daly_low95 - daly_tot_hi,
    diff_daly_hi    = pre_daly_hi - daly_tot_lo,
    impact_daly     = diff_daly / pre_daly * 100,
    impact_daly_low = diff_daly_low / pre_daly_hi * 100,
    impact_daly_hi  = diff_daly_hi / pre_daly_low95 * 100,
    nnv_daly        = tot_vacc / diff_daly,
    nnv_daly_lo     = tot_vacc / diff_daly_hi,
    nnv_daly_hi     = tot_vacc / diff_daly_low
  )


final_nnv_df$age_gr <- rep(age_gr[1:18],n_scenarios)
final_nnv_df$age_gr <- factor(final_nnv_df$age_gr, levels = age_gr_levels)


save("final_nnv_df", file = "02_Outputs/2_2_Tables/final_nnv_df.RData")
write_xlsx(final_nnv_df, path = "02_Outputs/2_2_Tables/final_nnv_df.xlsx")


## nnv graphs
p1 <- 
nnv_gg(
  final_nnv_df,
  y_var = "nnv",  # name of y variable as a string
  y_lab = "NNV to avert a single symptomatic case",
  x_lab = "Age group",
  title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
)

p2 <- 
nnv_gg(
  final_nnv_df,
  y_var = "nnv_fatal",  # name of y variable as a string
  y_lab = "NNV to avert a single fatal case",
  x_lab = "Age group",
  title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
)

p3 <- 
nnv_gg(
  final_nnv_df,
  y_var = "nnv_daly",  # name of y variable as a string
  y_lab = "NNV to avert a single DALY",
  x_lab = "Age group",
  title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
)

p1 <- p1 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single symptomatic case")

p2 <- p2 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single fatal case")

p3 <- p3 + theme(
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single DALY")


combined_nnv <- ggarrange(p1, p2, p3,
                          ncol = 1, nrow = 3,
                          labels = c("A", "B", "C"),
                          common.legend = TRUE,
                          legend = "bottom",
                          align = "v")
combined_nnv <- annotate_figure(combined_nnv,
                left = text_grob("Number needed to vaccinate to avert each outcome", rot = 90, size = 12),
                bottom = text_grob("Age group", size = 12))


ggsave(filename = "02_Outputs/2_1_Figures/nnv_symp_nat.jpg", p1, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_fatal_nat.jpg", p2, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_daly_nat.jpg", p3, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv.jpg", combined_nnv, width = 9, height = 14, dpi = 1200)

## for faceting
g1 <- 
nnv_nat_gg(combined_nnv_df_region, 
           y_var = "nnv",
           y_lab = "NNV to avert a single symptomatic case",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")
g2 <-
  nnv_nat_gg(combined_nnv_df_region, 
             y_var = "nnv_fatal",
             y_lab = "NNV to avert a single fatal case",
             x_lab = "Age group",
             title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

g3 <- 
nnv_nat_gg(combined_nnv_df_region, 
           y_var = "nnv_daly",
           y_lab = "NNV to avert a single DALY",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

g1 <- g1 +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single symptomatic case")

g2 <- g2 +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single fatal case")

g3 <- g3 +
  scale_x_discrete(labels = function(x) gsub(" years", "", x)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single DALY")

# Arrange plots side-by-side with a common legend
combined_nnv_graph <- ggarrange(g1, g2, g3,
                            ncol = 1, nrow = 3,
                            labels = c("D", "E", "F"),
                            common.legend = TRUE,
                            legend = "bottom",
                            align = "v")

combined_nnv_graph <- annotate_figure(combined_nnv_graph,
                left = text_grob("Number needed to vaccinate to avert each outcome", rot = 90, size = 12),
                bottom = text_grob("Age group", size = 12))

ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_region.jpg", combined_nnv_graph, width = 13, height = 9, dpi = 1200)

# table: pre-post cases -----------------------------------------------------
scenario_result_ce <- attach_target_column(df_nnv_ce,
                                           observed = observed_ce,
                                           region   = "Ceará")
scenario_result_bh <- attach_target_column(df_nnv_bh,
                                           observed = observed_bh,
                                           region   = "Bahia")
scenario_result_pa <- attach_target_column(df_nnv_pa,
                                           observed = observed_pa,
                                           region   = "Paraíba")
scenario_result_pn <- attach_target_column(df_nnv_pn,
                                           observed = observed_pn,
                                           region   = "Pernambuco")
scenario_result_pi <- attach_target_column(df_nnv_pi,
                                           observed = observed_pi,
                                           region   = "Piauí")
scenario_result_rg <- attach_target_column(df_nnv_rg,
                                           observed = observed_rg,
                                           region   = "Rio Grande do Norte")
scenario_result_tc <- attach_target_column(df_nnv_tc,
                                           observed = observed_tc,
                                           region   = "Tocantins")
scenario_result_ag <- attach_target_column(df_nnv_ag,
                                           observed = observed_ag,
                                           region   = "Alagoas")
scenario_result_mg <- attach_target_column(df_nnv_mg,
                                           observed = observed_mg,
                                           region   = "Minas Gerais")

combined_scenario <- bind_rows(scenario_result_ce, scenario_result_bh, scenario_result_pa,
                               scenario_result_pn, scenario_result_pi, scenario_result_rg,
                               scenario_result_tc, scenario_result_ag, scenario_result_mg)

combined_scenario_summ <- combined_scenario %>%
  group_by(target, scenario) %>%
  summarise(
    pre_cases    = sum(pre_cases, na.rm = TRUE),
    post_cases   = sum(post_cases, na.rm = TRUE),
    post_low95   = sum(post_low95, na.rm = TRUE),
    post_hi95    = sum(post_hi95, na.rm = TRUE),
    pre_low95    = sum(pre_low95, na.rm = TRUE),
    pre_hi95     = sum(pre_hi95, na.rm = TRUE),
    post_fatal   = sum(post_fatal, na.rm = TRUE),
    pre_fatal    = sum(pre_fatal, na.rm = TRUE),
    post_fatal_lo = sum(post_fatal_lo, na.rm = TRUE),
    post_fatal_hi = sum(post_fatal_hi, na.rm = TRUE),
    pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
    pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
    
    tot_vacc     = sum(tot_vacc),
    
    post_daly      = sum(post_daly, na.rm = TRUE),
    pre_daly       = sum(pre_daly, na.rm = TRUE),
    post_daly_lo   = sum(post_daly_lo, na.rm = TRUE),
    post_daly_hi   = sum(post_daly_hi, na.rm = TRUE),
    pre_daly_lo    = sum(pre_daly_lo, na.rm = TRUE),
    pre_daly_hi    = sum(pre_daly_hi, na.rm = TRUE),

    .groups      = "drop"
  ) %>%
  mutate(
    ## For Cases:
    diff_case      = pre_cases - post_cases,
    diff_case_low  = pre_low95 - post_low95,
    diff_case_hi   = pre_hi95 - post_hi95,
    impact         = diff_case / pre_cases * 100,
    impact_low     = diff_case_low / pre_low95 * 100,
    impact_hi      = diff_case_hi / pre_hi95 * 100,
    
    ## For Fatal outcomes:
    diff_fatal      = pre_fatal - post_fatal,
    diff_fatal_low  = pre_fatal_lo - post_fatal_lo,
    diff_fatal_hi   = pre_fatal_hi - post_fatal_hi,
    impact_fatal    = diff_fatal / pre_fatal * 100,
    impact_fatal_low = diff_fatal_low / pre_fatal_lo * 100,
    impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,
    
    ## For DALY:
    diff_daly       = pre_daly - post_daly,
    diff_daly_low   = pre_daly_lo - post_daly_lo,
    diff_daly_hi    = pre_daly_hi - post_daly_hi,
    impact_daly     = diff_daly / pre_daly * 100,
    impact_daly_low = diff_daly_low / pre_daly_lo * 100,
    impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100
  )

write_xlsx(combined_scenario_summ, path = "00_Data/0_2_Processed/combined_scenario_summ_v2.xlsx")


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

combined_impact_df <- bind_rows(global_impact_ce, global_impact_bh, global_impact_pn,
                                global_impact_pi, global_impact_pa, global_impact_rg,
                                global_impact_ag, global_impact_tc, global_impact_mg) 

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

combined_heatmap <- bind_rows(
  heatmap_ag,
  heatmap_bh,
  heatmap_ce,
  heatmap_mg,
  heatmap_pa,
  heatmap_pi,
  heatmap_pn,
  heatmap_rg,
  heatmap_tc
)

save("combined_heatmap", file = "00_Data/0_2_Processed/combined_heatmap.RData")

heatmap_summ <- combined_heatmap %>% group_by(scenario_id, Scenario, Supply, Delay) %>% summarise(
  tot_pre_cases  = sum(tot_pre_cases),
  tot_post_cases = sum(tot_post_cases),
  tot_pre_fatal  = sum(tot_pre_fatal),
  tot_post_fatal = sum(tot_post_fatal)
) %>% mutate(
  tot_diff = tot_pre_cases - tot_post_cases,
  impact = tot_diff / tot_pre_cases * 100,
  tot_pre_fatal = sum(tot_pre_fatal, na.rm = TRUE),
  tot_post_fatal = sum(tot_post_fatal, na.rm = TRUE),
  diff_fatal = tot_pre_fatal - tot_post_fatal,
  impact_fatal = diff_fatal / tot_pre_fatal * 100
)

write_xlsx(heatmap_summ, path = "02_Outputs/2_2_Tables/heatmap_summ.xlsx")

combined_heatmap_g <- 
ggplot(heatmap_summ, aes(x = Supply, y = Delay, fill = impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% Cases averted") +  # Add color legend breaks
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~ Scenario, 
             labeller = as_labeller(c("Scenario_1" = "<20 years only", 
                                      "Scenario_2" = "20-59 years only", 
                                      "Scenario_3" = ">60 years only")))+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig_tot_heatmap.jpg", combined_heatmap_g, width = 10, height = 4, dpi = 1200)


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

combined_prepost_fatal <- bind_rows(
  cum_fatal_ce$fatal_week_df,
  cum_fatal_bh$fatal_week_df,
  cum_fatal_ag$fatal_week_df,
  cum_fatal_mg$fatal_week_df,
  cum_fatal_pa$fatal_week_df,
  cum_fatal_pi$fatal_week_df,
  cum_fatal_pn$fatal_week_df,
  cum_fatal_rg$fatal_week_df,
  cum_fatal_tc$fatal_week_df
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
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 2/3 * overall_max, ymax = overall_max, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
              aes(x = Week, ymin = 1/3 * overall_max, ymax = 2/3 * overall_max, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
              aes(x = Week, ymin = 0, ymax = 1/3 * overall_max, fill = Scenario),
              alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_fatal_prop_lo, ymax = post_fatal_prop_hi, fill = Scenario),
              alpha = 0.2) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = pre_fatal_prop_lo, ymax = pre_fatal_prop_hi),
              fill = "gray", alpha = 0.2) +
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
       title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
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
    pre_fatal_prop      = (cum_fatal_pre / tot_pop) * 100,
    pre_fatal_prop_lo   = (cum_fatal_pre_lo / tot_pop) * 100,
    pre_fatal_prop_hi   = (cum_fatal_pre_hi / tot_pop) * 100,
    post_fatal_prop     = (cum_fatal / tot_pop) * 100,
    post_fatal_prop_lo  = (cum_fatal_lo / tot_pop) * 100,
    post_fatal_prop_hi  = (cum_fatal_hi / tot_pop) * 100
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

# contour map for VE ----------------------------------------------------------
contour_ce <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_ce,
  posterior            = posterior_ce,
  N                    = N_ceara,
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

contour_bh <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_bh,
  posterior            = posterior_bh,
  N                    = N_bahia,
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

contour_pn <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_pn,
  posterior            = posterior_pn,
  N                    = N_pemam,
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

contour_rg <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_rg,
  posterior            = posterior_rg,
  N                    = N_rg,
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

contour_pa <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_pa,
  posterior            = posterior_pa,
  N                    = N_pa,
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

contour_pi <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_pi,
  posterior            = posterior_pi,
  N                    = N_pi,
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

contour_tc <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_tc,
  posterior            = posterior_tc,
  N                    = N_tc,
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

contour_ag <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_ag,
  posterior            = posterior_ag,
  N                    = N_ag,
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

contour_mg <- contour_ve(
  target_age_list      = target_age_list,
  observed             = observed_mg,
  posterior            = posterior_mg,
  N                    = N_mg,
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

combined_contour <- bind_rows(
  contour_ag,
  contour_bh,
  contour_ce,
  contour_mg,
  contour_pa,
  contour_pi,
  contour_pn,
  contour_rg,
  contour_tc
)

save("combined_contour", file = "00_Data/0_2_Processed/combined_contour.RData")


contour_summ <- combined_contour %>% group_by(scenario_id, Scenario, Supply, VE) %>% summarise(
  tot_pre_cases  = sum(tot_pre_cases),
  tot_post_cases = sum(tot_post_cases),
  tot_pre_fatal  = sum(tot_pre_fatal),
  tot_post_fatal = sum(tot_post_fatal)
) %>% mutate(
  tot_diff = tot_pre_cases - tot_post_cases,
  impact = tot_diff / tot_pre_cases * 100,
  tot_pre_fatal = sum(tot_pre_fatal, na.rm = TRUE),
  tot_post_fatal = sum(tot_post_fatal, na.rm = TRUE),
  diff_fatal = tot_pre_fatal - tot_post_fatal,
  impact_fatal = diff_fatal / tot_pre_fatal * 100
)

conmbined_contour <- 
ggplot(data = contour_summ, aes(x = Supply, y = VE, z = impact)) +
  geom_contour_filled(alpha = 0.9) +  # Create filled contour plot
  scale_fill_viridis_d(name = "Impact (%)") +  # Use a nice color scale
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_blank()
  )+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black", size = 0.2) +
  facet_wrap(~Scenario)+
  facet_wrap(~ Scenario, 
             labeller = as_labeller(c("Scenario_1" = "<20 years only", 
                                      "Scenario_2" = "20-59 years only", 
                                      "Scenario_3" = ">60 years only")))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))+
  xlab("Vaccine coverage")+
  ylab("Vaccine efficacy")

ggsave(filename = "02_Outputs/2_1_Figures/fig_tot_contour.jpg", conmbined_contour, width = 10, height = 4, dpi = 1200)

combined_heatmap_g <- combined_heatmap_g + 
theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank()
)

combined_sens_g_heatmap <- 
ggarrange(combined_heatmap_g, 
          ncol = 1, nrow = 1,
          labels = c("A"),
          common.legend = TRUE,
          legend = "right",
          align = "v")
ggsave(filename = "02_Outputs/2_1_Figures/fig_combined_sens_g.jpg", combined_sens_g_heatmap, width = 10, height = 4, dpi = 1200)

owsa_g <- 
ggarrange(owsa_all,
          ncol = 1, nrow = 1,
          labels = c("B"),
          common.legend = TRUE,
          legend = "right",
          align = "v")
ggsave(filename = "02_Outputs/2_1_Figures/owsa_g.jpg", owsa_g, width = 10, height = 4, dpi = 1200)
