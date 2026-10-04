### Load fitted posteriors (original / long-term-FOI fits, from
### age_struc_fitting_region_2022.R -- fit_prevacc_bh_longterm,
### fit_prevacc_ce_longterm, etc., bundled the same way as fits_prevacc_finite.RData)
load("00_Data/0_2_Processed/fits_prevacc_longterm.RData")
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")

# LHS samples built from the original fits + long-term average FOI --
# see lhs_samples_sir.R. Loads posterior_idx_list, lhs_idx_list, lhs_combined.
load("00_Data/0_2_Processed/lhs_orv.RData")

load("00_Data/0_2_Processed/region_coverage.RData")

# extract params-----------------------------------------------------------------
param_ce <- extract_params(fit_prevacc_ce_longterm)
param_pa <- extract_params(fit_prevacc_pa_longterm)
param_pn <- extract_params(fit_prevacc_pn_longterm)
param_bh <- extract_params(fit_prevacc_bh_longterm)
param_rg <- extract_params(fit_prevacc_rg_longterm)
param_pi <- extract_params(fit_prevacc_pi_longterm)
param_ag <- extract_params(fit_prevacc_ag_longterm)
param_tc <- extract_params(fit_prevacc_tc_longterm)
param_mg <- extract_params(fit_prevacc_mg_longterm)
param_se <- extract_params(fit_prevacc_se_longterm)
param_go <- extract_params(fit_prevacc_go_longterm)

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

save(
  posterior_ce, posterior_bh, posterior_pa, posterior_pn, posterior_rg,
  posterior_pi, posterior_ag, posterior_tc, posterior_mg, posterior_se,
  posterior_go,
  file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/posterior_longterm_all.RData"
)

# scenarios
target_age_list <- list(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 1-11 years
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 12-17 years
                        c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0), # 18-64
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1) # 65+

)

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

longterm_sero_list <- list(
  "Ceará" = stan_data_prevacc_ce_longterm$sero,
  "Bahia" = stan_data_prevacc_bh_longterm$sero,
  "Paraíba" = stan_data_prevacc_pa_longterm$sero,
  "Pernambuco" = stan_data_prevacc_pn_longterm$sero,
  "Rio Grande do Norte" = stan_data_prevacc_rg_longterm$sero,
  "Piauí" = stan_data_prevacc_pi_longterm$sero,
  "Tocantins" = stan_data_prevacc_tc_longterm$sero,
  "Alagoas" = stan_data_prevacc_ag_longterm$sero,
  "Minas Gerais" = stan_data_prevacc_mg_longterm$sero,
  "Sergipe" = stan_data_prevacc_se_longterm$sero,
  "Goiás" = stan_data_prevacc_go_longterm$sero
)


# simulate_pre_ui_age()/run_simulation_scenarios_ui_ixchiq() look up these
# two as globals (not function arguments). lhs_orv.RData already loads them
# under these exact names, so no renaming/reassignment is needed here (unlike
# the finite-history script, where the saved objects carry a _finite suffix).

make_prevacc <- function(
    posterior, observed, N, region, lhs_sample, sero_vec = NULL,
    # TRUE: lhs_sample must be from lhs_combined_finite (short-term FOI).
    # FALSE: lhs_sample must be from lhs_combined (long-term FOI).
    finite_history = FALSE
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

prevacc_sets_longterm <- purrr::imap(
  regions,
  function(reg_args, region_name) {
    make_prevacc(
      posterior = reg_args$posterior,
      observed = reg_args$observed,
      N = reg_args$N,
      region = region_name,
      lhs_sample = reg_args$lhs_sample,
      finite_history = FALSE  # reg_args$lhs_sample is from lhs_combined (long-term)
    )
  }
)

save("prevacc_sets_longterm", file = "00_Data/0_2_Processed/prevacc_sets_longterm.RData")


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
  ~ c(.x, list(prevacc_ui = prevacc_sets_longterm[[.y]]$prevacc_ui))
)

# debug test----------------------------
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
      prevacc_sets_longterm[[region_name]]$prevacc_ui

    reg_args
  }
)
coverage_set  <- 0.5

## 0) set ------------------------------------------------------------------
ve_inf_set     <- c(0, 0.989)
coverage_set   <- c(0.10, 0.50, 0.9)
coverage_set = 0.5

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
             lhs_sample         = reg_args$lhs_sample,
             finite_history     = FALSE  # reg_args$lhs_sample is from lhs_combined (long-term)
           )
         })
  })
})

postsim_vc_ixchiq_model  <- imap(regions, function(reg_args, region_name) {

  ve_cov_list <- sim_results_vc_ixchiq_model [[region_name]]
  prevacc      <- prevacc_sets_longterm[[region_name]]

  # VE
  imap(ve_cov_list, function(cov_list, ve_tag) {

    # coverage
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

save("sim_results_vc_ixchiq_model", file = "00_Data/0_2_Processed/sim_results_vc_ixchiq_model_longterm.RData")
save("postsim_vc_ixchiq_model", file = "00_Data/0_2_Processed/postsim_vc_ixchiq_model_longterm.RData")

save(sim_results_vc_ixchiq_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/sim_results_vc_ixchiq_model_longterm.RData")
save(postsim_vc_ixchiq_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/postsim_vc_ixchiq_model_longterm.RData")

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

    }) |>
      set_names(names(cov_list))         # "cov1" "cov5" …
  }) |>
    set_names(names(ve_cov_list))        # "VE0" "VE98.9" …
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

## 1) nnv_results --------------------------------------------------------------
combined_nnv_df_region_coverage_model <- imap_dfr(nnv_results_coverage_model,
                                                  function(ve_cov_list, region_name) {
                                                    imap_dfr(ve_cov_list,                     # VE loop
                                                             function(cov_list, ve_tag) {
                                                               imap_dfr(cov_list,                    # coverage loop
                                                                        function(nnv, cov_tag) {
                                                                          nnv$final_summ_df %>%            # (age_seq × n_scenarios 행)
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


write_xlsx(combined_nnv_national_age_ixchiq, path = "02_Outputs/2_2_Tables/combined_nnv_national_age_ixchiq_longterm.xlsx")
save(combined_nnv_national_age_ixchiq, file = "00_Data/0_2_Processed/combined_nnv_national_age_ixchiq_longterm.RData")
save(combined_nnv_df_region_coverage_model, file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/combined_nnv_df_region_coverage_model_longterm.RData")
