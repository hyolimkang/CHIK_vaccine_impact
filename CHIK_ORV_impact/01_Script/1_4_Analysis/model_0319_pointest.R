## point estimate results
### prevacc sims 
sim_result_novacc_ce <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_ceara,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = param_ce$base_beta,
  I0_draw = param_ce$I0,
  R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Ceará"] * age_groups)),
  #R0 = rep(0, 18),
  rho = param_ce$rho,
  gamma = param_ce$gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  total_supply = 0,
  weekly_delivery_speed = 0
)

sim_result_novacc_pn <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_pemam,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = param_pn$base_beta,
  I0_draw = param_pn$I0,
  R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Pernambuco"] * age_groups)),
  #R0 = rep(0, 18),
  rho = param_pn$rho,
  gamma = param_pn$gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  total_supply = 0,
  weekly_delivery_speed = 0
)

sim_result_novacc_bh <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_bahia,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = param_bh$base_beta,
  I0_draw = param_bh$I0,
  R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Bahia"] * age_groups)),
  #R0 = rep(0, 18),
  rho = param_bh$rho,
  gamma = param_bh$gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  total_supply = 0,
  weekly_delivery_speed = 0
)

sim_result_novacc_pa <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_pa,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = param_pa$base_beta,
  I0_draw = param_pa$I0,
  R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Paraíba"] * age_groups)),
  #R0 = rep(0, 18),
  rho = param_pa$rho,
  gamma = param_pa$gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  total_supply = 0,
  weekly_delivery_speed = 0
)

sim_result_novacc_rg <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_rg,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = param_rg$base_beta,
  I0_draw = param_rg$I0,
  R0  = (1 - exp(-bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == "Rio Grande do Norte"] * age_groups)),
  #R0 = rep(0, 18),
  rho = param_rg$rho,
  gamma = param_rg$gamma,
  
  delay = 53,                   # Vaccination starts from Week x
  VE_block = 0,                 # Vaccine efficacy
  target_age = c(rep(0,18)),    # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  weekly_delivery_speed = 0
)

pre_results_ce <- summarize_simulation(
  sim_result       = sim_result_novacc_ce,
  age_gr_levels    = age_gr_levels,
  age_gr           = age_gr,
  lhs_sample_young = lhs_sample_young,
  lhs_old          = lhs_old,
  le_sample        = le_sample,
  hosp             = hosp,
  fatal            = fatal,
  nh_fatal         = nh_fatal
)

pre_results_pn <- summarize_simulation(
  sim_result       = sim_result_novacc_pn,
  age_gr_levels    = age_gr_levels,
  age_gr           = age_gr,
  lhs_sample_young = lhs_sample_young,
  lhs_old          = lhs_old,
  le_sample        = le_sample,
  hosp             = hosp,
  fatal            = fatal,
  nh_fatal         = nh_fatal
)

pre_results_bh <- summarize_simulation(
  sim_result       = sim_result_novacc_bh,
  age_gr_levels    = age_gr_levels,
  age_gr           = age_gr,
  lhs_sample_young = lhs_sample_young,
  lhs_old          = lhs_old,
  le_sample        = le_sample,
  hosp             = hosp,
  fatal            = fatal,
  nh_fatal         = nh_fatal
)

pre_results_pa <- summarize_simulation(
  sim_result       = sim_result_novacc_pa,
  age_gr_levels    = age_gr_levels,
  age_gr           = age_gr,
  lhs_sample_young = lhs_sample_young,
  lhs_old          = lhs_old,
  le_sample        = le_sample,
  hosp             = hosp,
  fatal            = fatal,
  nh_fatal         = nh_fatal
)

pre_results_rg <- summarize_simulation(
  sim_result       = sim_result_novacc_rg,
  age_gr_levels    = age_gr_levels,
  age_gr           = age_gr,
  lhs_sample_young = lhs_sample_young,
  lhs_old          = lhs_old,
  le_sample        = le_sample,
  hosp             = hosp,
  fatal            = fatal,
  nh_fatal         = nh_fatal
)

# dfs for pre-age
pre_summary_age_ce <- pre_results_ce$summary_cases_pre_age
pre_summary_ce     <- pre_results_ce$summary_cases_pre

pre_summary_age_bh <- pre_results_bh$summary_cases_pre_age
pre_summary_bh     <- pre_results_bh$summary_cases_pre
pre_summary_bh_all     <- pre_results_bh$summary_cases_pre_all

pre_summary_age_rg <- pre_results_bh$summary_cases_pre_age
pre_summary_rg    <- pre_results_bh$summary_cases_pre
pre_summary_rg_all     <- pre_results_rg$summary_cases_pre_all

ggplot(pre_summary_rg_all)+
  geom_line(aes(x = Week, y = Median))+
  geom_point(data = observed_rg, aes(x = Week, y = Observed))
  
## postvacc sims
postsim_ce <- run_simulation_scenarios(
  target_age_list   = target_age_list,
  supply            = supply_ce,
  N                 = N_pa,
  param             = param_pa,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Ceará",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr            = age_gr,
  age_gr_levels     = age_gr_levels
)

postsim_pn <- run_simulation_scenarios(
  target_age_list   = target_age_list,
  supply            = supply_pn,
  N                 = N_pemam,
  param             = param_pn,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Pernambuco",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr            = age_gr,
  age_gr_levels     = age_gr_levels
)

postsim_bh <- run_simulation_scenarios(
  target_age_list   = target_age_list,
  supply            = supply_bh,
  N                 = N_bahia,
  param             = param_bh,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Bahia",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr            = age_gr,
  age_gr_levels     = age_gr_levels
)

postsim_pa <- run_simulation_scenarios(
  target_age_list   = target_age_list,
  supply            = supply_pa,
  N                 = N_pa,
  param             = param_pa,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Paraíba",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr            = age_gr,
  age_gr_levels     = age_gr_levels
)

postsim_rg <- run_simulation_scenarios(
  target_age_list   = target_age_list,
  supply            = supply_rg,
  N                 = N_rg,
  param             = param_rg,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Rio Grande do Norte",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr            = age_gr,
  age_gr_levels     = age_gr_levels
)

postsim_all_ce <- postsim_all(
  target_age_list = target_age_list,
  supply          = supply_ce,
  N               = N_ceara, 
  param           = param_ce, 
  age_groups      = age_groups,
  lhs_sample_young = lhs_sample_young, 
  lhs_old          = lhs_old, 
  le_sample        = le_sample,
  hosp             = hosp, 
  fatal            = fatal, 
  nh_fatal         = nh_fatal,
  age_gr           = age_gr, 
  age_gr_levels    = age_gr_levels,
  bra_foi_state_summ = bra_foi_state_summ, 
  region          = "Ceará",
  pre_summary_cases_age = pre_summary_age_ce,
  pre_summary_cases     = pre_summary_ce
)

postsim_all_bh <- postsim_all(
  target_age_list = target_age_list,
  supply          = supply_bh,
  N               = N_bahia, 
  param           = param_bh, 
  age_groups      = age_groups,
  lhs_sample_young = lhs_sample_young, 
  lhs_old          = lhs_old, 
  le_sample        = le_sample,
  hosp             = hosp, 
  fatal            = fatal, 
  nh_fatal         = nh_fatal,
  age_gr           = age_gr, 
  age_gr_levels    = age_gr_levels,
  bra_foi_state_summ = bra_foi_state_summ, 
  region          = "Bahia",
  pre_summary_cases_age = pre_summary_age_bh,
  pre_summary_cases     = pre_summary_bh
)
