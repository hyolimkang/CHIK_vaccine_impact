### Load fitted posteriors
load("00_Data/0_2_Processed/fits_prevacc_finite.RData")

regions <- list(
  "Ceará" = list(
    observed = observed_ce,
    N = N_ceara$Ceará,
    posterior = posterior_ce,
    lhs_sample = lhs_combined$Ceará
  ),
  "Bahia" = list(
    observed = observed_bh,
    N = N_bahia$Bahia,
    posterior = posterior_bh,
    lhs_sample = lhs_combined$Bahia
  ),
  "Paraíba" = list(
    observed = observed_pa,
    N = N_pa$Paraíba,
    posterior = posterior_pa,
    lhs_sample = lhs_combined$Paraíba
  ),
  "Pernambuco" = list(
    observed = observed_pn,
    N = N_pemam$Pernambuco,
    posterior = posterior_pn,
    lhs_sample = lhs_combined$Pernambuco
  ),
  "Rio Grande do Norte" = list(
    observed = observed_rg,
    N = N_rg$`Rio Grande do Norte`,
    posterior = posterior_rg,
    lhs_sample = lhs_combined$`Rio Grande do Norte`
  ),
  "Piauí" = list(
    observed = observed_pi,
    N = N_pi$Piauí,
    posterior = posterior_pi,
    lhs_sample = lhs_combined$Piauí
  ),
  "Tocantins" = list(
    observed = observed_tc,
    N = N_tc$Tocantins,
    posterior = posterior_tc,
    lhs_sample = lhs_combined$Tocantins
  ),
  "Alagoas" = list(
    observed = observed_ag,
    N = N_ag$Alagoas,
    posterior = posterior_ag,
    lhs_sample = lhs_combined$Alagoas
  ),
  "Minas Gerais" = list(
    observed = observed_mg,
    N = N_mg$`Minas Gerais`,
    posterior = posterior_mg,
    lhs_sample = lhs_combined$`Minas Gerais`
  ),
  "Sergipe" = list(
    observed = observed_se,
    N = N_se$Sergipe,
    posterior = posterior_se,
    lhs_sample = lhs_combined$Sergipe
  ),
  "Goiás" = list(
    observed = observed_go,
    N = N_go$Goiás,
    posterior = posterior_go,
    lhs_sample = lhs_combined$Goiás
  )
)

finite_sero_list <- list(
  "Ceará" = stan_data_prevacc_ce_finite$sero,
  "Bahia" = stan_data_prevacc_bh_finite$sero,
  "Paraíba" = stan_data_prevacc_pa_finite$sero,
  "Pernambuco" = stan_data_prevacc_pn_finite$sero,
  "Rio Grande do Norte" = stan_data_prevacc_rg_finite$sero,
  "Piauí" = stan_data_prevacc_pi_finite$sero,
  "Tocantins" = stan_data_prevacc_tc_finite$sero,
  "Alagoas" = stan_data_prevacc_ag_finite$sero,
  "Minas Gerais" = stan_data_prevacc_mg_finite$sero,
  "Sergipe" = stan_data_prevacc_se_finite$sero,
  "Goiás" = stan_data_prevacc_go_finite$sero
)

# extract params-----------------------------------------------------------------
param_ce <- extract_params(fit_prevacc_ce_finite)
param_pa <- extract_params(fit_prevacc_pa_finite)
param_pn <- extract_params(fit_prevacc_pn_finite)
param_bh <- extract_params(fit_prevacc_bh_finite)
param_rg <- extract_params(fit_prevacc_rg_finite)
param_pi <- extract_params(fit_prevacc_pi_finite)
param_ag <- extract_params(fit_prevacc_ag_finite)
param_tc <- extract_params(fit_prevacc_tc_finite)
param_mg <- extract_params(fit_prevacc_mg_finite)
param_se <- extract_params(fit_prevacc_se_finite)
param_go <- extract_params(fit_prevacc_go_finite)

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

prevacc_sets_finite <- purrr::imap(
  regions,
  function(reg_args, region_name) {
    make_prevacc(
      posterior = reg_args$posterior,
      observed = reg_args$observed,
      N = reg_args$N,
      region = region_name,
      lhs_sample = reg_args$lhs_sample
    )
  }
)

save("prevacc_sets_finite", file = "00_Data/0_2_Processed/prevacc_sets_finite.RData")

## 1) Post vacc simulation --------------------------------------------------------

regions_base <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 posterior = posterior_ce, lhs_sample = lhs_combined$Ceará),
  "Bahia" = list(observed = observed_bh, N = N_bahia$Bahia,
                 posterior = posterior_bh, lhs_sample = lhs_combined$Bahia),
  "Paraíba" = list(observed = observed_pa, N = N_pa$Paraíba,
                   posterior = posterior_pa, lhs_sample = lhs_combined$Paraíba),
  "Pernambuco" = list(observed = observed_pn, N = N_pemam$Pernambuco,
                      posterior = posterior_pn, lhs_sample = lhs_combined$Pernambuco),
  "Rio Grande do Norte" = list(observed = observed_rg, N = N_rg$`Rio Grande do Norte`,
                               posterior = posterior_rg,
                               lhs_sample = lhs_combined$`Rio Grande do Norte`),
  "Piauí" = list(observed = observed_pi, N = N_pi$Piauí,
                 posterior = posterior_pi, lhs_sample = lhs_combined$Piauí),
  "Tocantins" = list(observed = observed_tc, N = N_tc$Tocantins,
                     posterior = posterior_tc, lhs_sample = lhs_combined$Tocantins),
  "Alagoas" = list(observed = observed_ag, N = N_ag$Alagoas,
                   posterior = posterior_ag, lhs_sample = lhs_combined$Alagoas),
  "Minas Gerais" = list(observed = observed_mg, N = N_mg$`Minas Gerais`,
                        posterior = posterior_mg,
                        lhs_sample = lhs_combined$`Minas Gerais`),
  "Sergipe" = list(observed = observed_se, N = N_se$Sergipe,
                   posterior = posterior_se, lhs_sample = lhs_combined$Sergipe),
  "Goiás" = list(observed = observed_go, N = N_go$Goiás,
                 posterior = posterior_go, lhs_sample = lhs_combined$Goiás)
)

regions <- purrr::imap(
  regions_base,
  ~ c(.x, list(prevacc_ui = prevacc_sets_finite[[.y]]$prevacc_ui))
)

# debug test
regions_base <- list(
  "Ceará" = list(
    observed = observed_ce,
    N = N_ceara$Ceará,
    posterior = posterior_ce,
    lhs_sample = lhs_combined$Ceará
  )
)

regions <- purrr::imap(
  regions_base,
  function(reg_args, region_name) {
    reg_args$prevacc_ui <-
      prevacc_sets_finite[[region_name]]$prevacc_ui
    
    reg_args
  }
)
coverage_set  <- 0.5  

## 0) set ------------------------------------------------------------------
ve_inf_set     <- c(0, 0.989)         
coverage_set   <- c(0, 0.10, 0.50, 0.9) 


sim_results_vc_ixchiq_model  <- imap(regions, function(reg_args, reg_name) {
  
  n_draws <- length(reg_args$posterior$gamma)
  
  ## ─ VE loop  ──────────────────────────────────────────────────────────
  imap(set_names(ve_inf_set, paste0("VE", ve_inf_set*100)), function(ve_val, ve_tag) {
    
    ## ─ coverage loop ───────────────────────────────────────────────────────
    imap(set_names(coverage_set, paste0("cov", coverage_set*100)),
         function(cov_val, cov_tag) {
           
           # if !ve=0, then use ve_val 
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
             lhs_sample         = reg_args$lhs_sample        
           )
         })
  })
})

postsim_vc_ixchiq_model  <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- sim_results_vc_ixchiq_model [[region_name]]    
  prevacc      <- prevacc_sets_finite[[region_name]]
  
  # VE 수준 루프
  imap(ve_cov_list, function(cov_list, ve_tag) {
    
    # coverage 수준 루프
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
      
    }) %>% set_names(names(cov_list))  # cov 태그로 네임 설정
    
  }) %>% set_names(names(ve_cov_list))  # VE 태그로 네임 설정
  
})


