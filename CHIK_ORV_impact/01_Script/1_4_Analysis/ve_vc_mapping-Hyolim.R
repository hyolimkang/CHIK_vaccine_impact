
regions <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 prevacc_ui = prevacc_ui_ce, posterior = posterior_ce,
                 lhs_sample = lhs_combined$Ceará),
  "Bahia" = list(observed = observed_bh, N = N_bahia$Bahia,
                 prevacc_ui = prevacc_ui_bh, posterior = posterior_bh,
                 lhs_sample = lhs_combined$Bahia),
  "Paraíba" = list(observed = observed_pa, N = N_pa$Paraíba,
                   prevacc_ui = prevacc_ui_pa, posterior = posterior_pa,
                   lhs_sample = lhs_combined$Paraíba),
  "Pernambuco" = list(observed = observed_pn, N = N_pemam$Pernambuco,
                      prevacc_ui = prevacc_ui_pn, posterior = posterior_pn,
                      lhs_sample = lhs_combined$Pernambuco),
  "Rio Grande do Norte" = list(observed = observed_rg, N = N_rg$`Rio Grande do Norte`,
                               prevacc_ui = prevacc_ui_rg, posterior = posterior_rg,
                               lhs_sample = lhs_combined$`Rio Grande do Norte`),
  "Piauí" = list(observed = observed_pi, N = N_pi$Piauí,
                 prevacc_ui = prevacc_ui_pi, posterior = posterior_pi,
                 lhs_sample = lhs_combined$Piauí),
  "Tocantins" = list(observed = observed_tc, N = N_tc$Tocantins,
                     prevacc_ui = prevacc_ui_tc, posterior = posterior_tc,
                     lhs_sample = lhs_combined$Tocantins),
  "Alagoas" = list(observed = observed_ag, N = N_ag$Alagoas,
                   prevacc_ui = prevacc_ui_ag, posterior = posterior_ag,
                   lhs_sample = lhs_combined$Alagoas),
  "Minas Gerais" = list(observed = observed_mg, N = N_mg$`Minas Gerais`,
                        prevacc_ui = prevacc_ui_mg, posterior = posterior_mg,
                        lhs_sample = lhs_combined$`Minas Gerais`),
  "Sergipe" = list(observed = observed_se, N = N_se$Sergipe,
                   prevacc_ui = prevacc_ui_se, posterior = posterior_se,
                   lhs_sample = lhs_combined$Sergipe),
  "Goiás" = list(observed = observed_go, N = N_go$Goiás,
                 prevacc_ui = prevacc_ui_go, posterior = posterior_go,
                 lhs_sample = lhs_combined$Goiás)
)

# debug test
regions <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 prevacc_ui = prevacc_ui_ce, posterior = posterior_ce,
                 lhs_sample = lhs_combined$Ceará)
)
coverage_set  <- 0.5 # test 
## 0) set ------------------------------------------------------------------
ve_inf_set     <- c(0, 0.989)         
coverage_set   <- c(0.10, 0.50, 0.9) 
doses <- 2e6                       

## 1) simulation --------------------------------------------------------
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
sim_results_dose_ixchiq_model <- imap(regions, function(reg_args, reg_name) {
  
  ## ─ 지역별 고정 coverage 계산 ────────────────────────────
  doses_this  <- doses_by_region[[reg_name]]          # 해당 지역 도즈
  total_cov_i <- doses_this / sum(reg_args$N)         # 인구 대비 비율 (0–1)
  
  ## ─ VE 수준 루프만 돌리기 ────────────────────────────────
  imap(set_names(ve_inf_set, paste0("VE", ve_inf_set * 100)),
       function(ve_val, ve_tag) {
         
         run_simulation_scenarios_ui(
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
           ve_inf             = ve_val,          # ← VE 값
           total_coverage     = total_cov_i      # ← 고정된 coverage
         )
       })
})

postsim_vc_ixchiq_model  <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- sim_results_vc_ixchiq_model [[region_name]]    # VE×cov 결과 리스트
  prevacc      <- prevacc_sets[[region_name]]
  
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
postsim_dose_ixchiq_model  <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- sim_results_dose_ixchiq_model [[region_name]]    # VE×cov 결과 리스트
  prevacc      <- prevacc_sets[[region_name]]
  
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

save("sim_results_vc_ixchiq_model", file = "00_Data/0_2_Processed/sim_results_vc_ixchiq_model.RData")
save("sim_results_dose_ixchiq_model", file = "00_Data/0_2_Processed/sim_results_dose_ixchiq_model.RData")

save("postsim_vc_ixchiq_model", file = "00_Data/0_2_Processed/postsim_vc_ixchiq_model.RData")
save("postsim_dose_ixchiq_model", file = "00_Data/0_2_Processed/postsim_dose_ixchiq_model.RData")

regions <- list(
  "Ceará" = list(
    vacc_alloc = vacc_alloc_ce,
    observed   = observed_ce,
    N          = N_ceara$Ceará,
    prevacc_ui = prevacc_ui_ce,
    posterior  = posterior_ce
  ),
  "Bahia" = list(
    vacc_alloc = vacc_alloc_bh,
    observed   = observed_bh,
    N          = N_bahia$Bahia,
    prevacc_ui = prevacc_ui_bh,
    posterior  = posterior_bh
  ),
  "Paraíba" = list(
    vacc_alloc = vacc_alloc_pa,
    observed   = observed_pa,
    N          = N_pa$Paraíba,
    prevacc_ui = prevacc_ui_pa,
    posterior  = posterior_pa
  ),
  "Pernambuco" = list(
    vacc_alloc = vacc_alloc_pn,
    observed   = observed_pn,
    N          = N_pemam$Pernambuco,
    prevacc_ui = prevacc_ui_pn,
    posterior  = posterior_pn
  ),
  "Rio Grande do Norte" = list(
    vacc_alloc = vacc_alloc_rg,
    observed   = observed_rg,
    N          = N_rg$`Rio Grande do Norte`,
    prevacc_ui = prevacc_ui_rg,
    posterior  = posterior_rg
  ),
  "Piauí" = list(
    vacc_alloc = vacc_alloc_pi,
    observed   = observed_pi,
    N          = N_pi$Piauí,
    prevacc_ui = prevacc_ui_pi,
    posterior  = posterior_pi
  ),
  "Tocantins" = list(
    vacc_alloc = vacc_alloc_tc,
    observed   = observed_tc,
    N          = N_tc$Tocantins,
    prevacc_ui = prevacc_ui_tc,
    posterior  = posterior_tc
  ),
  "Alagoas" = list(
    vacc_alloc = vacc_alloc_ag,
    observed   = observed_ag,
    N          = N_ag$Alagoas,
    prevacc_ui = prevacc_ui_ag,
    posterior  = posterior_ag
  ),
  "Minas Gerais" = list(
    vacc_alloc = vacc_alloc_mg,
    observed   = observed_mg,
    N          = N_mg$`Minas Gerais`,
    prevacc_ui = prevacc_ui_mg,
    posterior  = posterior_mg
  ),
  "Sergipe" = list(
    vacc_alloc = vacc_alloc_se,
    observed   = observed_se,
    N          = N_se$Sergipe,
    prevacc_ui = prevacc_ui_se,
    posterior  = posterior_se
  ),
  "Goiás" = list(
    vacc_alloc = vacc_alloc_go,
    observed   = observed_go,
    N          = N_go$Goiás,
    prevacc_ui = prevacc_ui_go,
    posterior  = posterior_go
  )
)


nnv_results_coverage_model <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- postsim_vc_coverage_model[[region_name]]
  
  imap(ve_cov_list, function(cov_list, ve_tag) {
    
    imap(cov_list, function(postsim_ui, cov_tag) {
      
      nnv_list(
        vacc_alloc     = vacc_allocation(postsim_ui, reg_args$observed, region_name),
        postsim_all_ui = postsim_ui,
        N              = reg_args$N,
        region         = region_name,
        observed       = reg_args$observed
      )
      
    }) %>% set_names(names(cov_list))
    
  }) %>% set_names(names(ve_cov_list))
  
})


# upated 
postsim_vc_ixchiq_model <- purrr::map(
  postsim_vc_ixchiq_model,
  ~ purrr::map(.x, ~ .x[names(.x) != "cov0"])
)

n_scenarios <- 4
nnv_results_coverage_model <- imap(postsim_vc_ixchiq_model, function(ve_cov_list, region_name) {
  
  # 원본 관측·인구 데이터
  reg_args <- regions[[region_name]]
  
  ## ― VE 수준 루프 ―────────────────────────────
  imap(ve_cov_list, function(cov_list, ve_tag) {
    
    ## ― coverage 수준 루프 ―────────────────────
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

## 1) nnv_results (지역 → VE → VC → nnv_list 결과) 결합 --------------------------
combined_nnv_df_region_coverage_model <- imap_dfr(nnv_results_coverage_model,      # 지역 loop
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
  group_by(VE, VC, scenario, AgeGroup) %>%            # VE·VC·scenario 기준으로 그룹화
  summarise(
    across(c(tot_vacc, pre_inf:diff_hosp_hi), sum, na.rm = TRUE),
    .groups = "drop"
  )

write_xlsx(combined_nnv_national_age_ixchiq, path = "02_Outputs/2_2_Tables/combined_nnv_national_age_ixchiq.xlsx")
save(combined_nnv_national_age_ixchiq, file = "00_Data/0_2_Processed/combined_nnv_national_age_ixchiq.RData")
save(combined_nnv_df_region_coverage_model, file = "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_benefit_risk/01_Data/combined_nnv_df_region_coverage_model.RData")

combined_nnv_national_ve_ixchiq <- combined_nnv_df_region_coverage_model %>%
  group_by(VE, VC, scenario) %>%            # VE·VC·scenario 기준으로 그룹화
  summarise(
    across(c(tot_vacc, pre_inf:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ## NNV
    nnv_inf      = tot_vacc / diff_inf,
    nnv_inf_lo   = tot_vacc / diff_inf_hi,
    nnv_inf_hi   = tot_vacc / diff_inf_lo,
    nnv          = tot_vacc / diff,
    nnv_lo       = tot_vacc / diff_hi,
    nnv_hi       = tot_vacc / diff_low,
    nnv_fatal    = tot_vacc / diff_fatal,
    nnv_fatal_lo = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi = tot_vacc / diff_fatal_low,
    nnv_daly     = tot_vacc / diff_daly,
    nnv_daly_lo  = tot_vacc / diff_daly_hi,
    nnv_daly_hi  = tot_vacc / diff_daly_low,
    
    ## per 100k 
    per1M_inf     = diff_inf       / tot_vacc * 1e5,
    per1M_inf_lo  = diff_inf_lo    / tot_vacc * 1e5,
    per1M_inf_hi  = diff_inf_hi    / tot_vacc * 1e5,
    
    per1M_diff     = diff       / tot_vacc * 1e5,
    per1M_diff_lo  = diff_low    / tot_vacc * 1e5,
    per1M_diff_hi  = diff_hi    / tot_vacc * 1e5,
    
    per1M_fatal     = diff_fatal       / tot_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low    / tot_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi    / tot_vacc * 1e5,
    
    per1M_daly     = diff_daly       / tot_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low    / tot_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi    / tot_vacc * 1e5
  )

combined_nnv_setting_ve_ixchiq <- combined_nnv_df_region_coverage_model %>% 
  group_by(VE, VC, setting, scenario) %>%           # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_inf:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(VE, setting, scenario) %>%                   # 전체 투여량(시나리오) 계산
  mutate(scenario_vacc = sum(tot_vacc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    nnv_inf      = scenario_vacc / diff_inf,
    nnv_inf_lo   = scenario_vacc / diff_inf_hi,
    nnv_inf_hi   = scenario_vacc / diff_inf_lo,
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    nnv_fatal   = scenario_vacc / diff_fatal,
    nnv_fatal_lo= scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi= scenario_vacc / diff_fatal_low,
    nnv_daly    = scenario_vacc / diff_daly,
    nnv_daly_lo = scenario_vacc / diff_daly_hi,
    nnv_daly_hi = scenario_vacc / diff_daly_low,
    
    per1M_inf     = diff_inf       / scenario_vacc * 1e5,
    per1M_inf_lo  = diff_inf_lo    / scenario_vacc * 1e5,
    per1M_inf_hi  = diff_inf_hi    / scenario_vacc * 1e5,
    
    per1M_diff     = diff       / scenario_vacc * 1e5,
    per1M_diff_lo  = diff_low    / scenario_vacc * 1e5,
    per1M_diff_hi  = diff_hi    / scenario_vacc * 1e5,
    
    per1M_fatal     = diff_fatal       / scenario_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low    / scenario_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi    / scenario_vacc * 1e5,
    
    per1M_daly     = diff_daly       / scenario_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low    / scenario_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi    / scenario_vacc * 1e5,
    
  )

combined_nnv_setting_summ_ve_ixchiq <- combined_nnv_setting_ve_ixchiq %>% 
  group_by(VE, VC, setting, scenario) %>%                   # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_inf:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    nnv_inf      = tot_vacc / diff_inf,
    nnv_inf_lo   = tot_vacc / diff_inf_hi,
    nnv_inf_hi   = tot_vacc / diff_inf_lo,
    nnv         = tot_vacc / diff,
    nnv_lo      = tot_vacc / diff_hi,
    nnv_hi      = tot_vacc / diff_low,
    nnv_fatal   = tot_vacc / diff_fatal,
    nnv_fatal_lo= tot_vacc / diff_fatal_hi,
    nnv_fatal_hi= tot_vacc / diff_fatal_low,
    nnv_daly    = tot_vacc / diff_daly,
    nnv_daly_lo = tot_vacc / diff_daly_hi,
    nnv_daly_hi = tot_vacc / diff_daly_low,
    
    per1M_inf     = diff_inf       / tot_vacc * 1e5,
    per1M_inf_lo  = diff_inf_lo    / tot_vacc * 1e5,
    per1M_inf_hi  = diff_inf_hi    / tot_vacc * 1e5,
    
    per1M_diff     = diff       / tot_vacc * 1e5,
    per1M_diff_lo  = diff_low    / tot_vacc * 1e5,
    per1M_diff_hi  = diff_hi    / tot_vacc * 1e5,
    
    per1M_fatal     = diff_fatal       / tot_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low    / tot_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi    / tot_vacc * 1e5,
    
    per1M_daly     = diff_daly       / tot_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low    / tot_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi    / tot_vacc * 1e5,
    
  )

combined_nnv_national_ve_ixchiq$setting <- "National"
combined_nnv_setting_summ_ve_ixchiq <- rbind(combined_nnv_national_ve_ixchiq, combined_nnv_setting_summ_ve_ixchiq)

fix_factors <- function(df){
  df %>% 
    mutate(
      ## ── 공통 범주형 ───────────────────────────────
      setting  = factor(setting,
                        levels = c("National","High","Moderate","Low")),
      scenario = factor(scenario,
                        levels = c("Scenario_1","Scenario_2","Scenario_3","Scenario_4"),
                        labels = c("1–11 years","12–17 years","18–59 years","60+ years")),
      ## ── VE / Coverage ────────────────────────────
      VE_lab   = factor(VE, levels = c("VE0","VE98.9"),
                        labels = c("0 %","98·9 %")),
      VC_lab   = factor(case_when(
        VC == "cov10"  ~ "10 %",
        VC == "cov50"  ~ "50 %",
        VC == "cov90" ~ "90 %"),
        levels = c("10 %","50 %","90 %")),
      ## ── VE × Coverage 조합 (fill‧범례용) ────────
      VEVC     = factor(paste0(VE_lab, "_", VC_lab),
                        levels = c("0 %_10 %","0 %_50 %","0 %_90 %",
                                   "98·9 %_10 %","98·9 %_50 %","98·9 %_90 %"))
    )
}

combined_nnv_setting_summ_ve_ixchiq <- combined_nnv_setting_summ_ve_ixchiq %>% 
  fix_factors()

nnv_program_scale <- combined_nnv_national_ve_ixchiq %>% filter(
  VC == "cov50"
)

save(combined_nnv_setting_summ_ve_ixchiq, file = "00_Data/0_2_Processed/combined_nnv_setting_summ_ve_ixchiq.RData")
write.csv(combined_nnv_df_region_coverage_model,
          file = "00_Data/0_2_Processed/combined_nnv_df_region_coverage_model.csv",
          row.names = FALSE)

make_nnv_overlay_vc <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0","VE98.9"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_10 %"    = "#A6CEE3",
    "0 %_50 %"    = "#B2DF8A",
    "0 %_90 %"   = "#FB9A99",
    "98·9 %_10 %" = "#1F78B4",
    "98·9 %_50 %" = "#33A02C",
    "98·9 %_90 %"= "#E31A1C"
  )
  
  #cb_okabe <- c(
  #  "0 %_1 %"    = "#56B4E9",  # sky blue
  ##  "0 %_5 %"    = "#009E73",  # bluish green
  ##  "0 %_10 %"   = "#F0E442",  # yellow
  #  "98·9 %_1 %" = "#E69F00",  # orange
  #  "98·9 %_5 %" = "#0072B2",  # blue
  #  "98·9 %_10 %"= "#CC79A7"   # reddish purple
  #)
  dodge_cov <- position_dodge(width = 0.9)
  
  ggplot(df2, aes(x = scenario)) +
    
    ## ── VE 0 % (넓은 막대) ─────────────────────────
    geom_col(data = filter(df2, VE_lab == "0 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.8) +
    geom_errorbar(data = filter(df2, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.25) +
    
    ## ── VE 98·9 % (좁은 막대) ──────────────────────
    geom_col(data = filter(df2, VE_lab == "98·9 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "98·9 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.15) +
    
    facet_grid(~ setting) +
    
    ## ── 색·범례 설정 ───────────────────────────────
    scale_fill_manual(
      name   = "VE × Coverage",
      values = cols6,
      breaks = names(cols6)          # ← 범례 순서 고정 (2→5→10 %)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +  # ← 행 우선 두 줄
    
    ## ── 축·테마 ───────────────────────────────────
    scale_y_continuous(labels = short_si,
                       expand = expansion(mult = c(0, .05))) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position     = "bottom",
      axis.text.x         = element_text(angle = 45, hjust = 1),
      strip.background    = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.x  = element_blank()
    )
}

make_per1M_overlay_vc <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0", "VE98.9"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_2 %"    = "#A6CEE3",
    "0 %_5 %"    = "#B2DF8A",
    "0 %_10 %"   = "#FB9A99",
    "98·9 %_2 %" = "#1F78B4",
    "98·9 %_5 %" = "#33A02C",
    "98·9 %_10 %"= "#E31A1C"
  )
  
  dodge_cov <- position_dodge(width = 0.9)
  
  ggplot(df2, aes(x = scenario)) +
    
    ## ── VE 98.9 % (넓은 막대) ──────────────────────
    geom_col(data = filter(df2, VE_lab == "98·9 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.8) +
    geom_errorbar(data = filter(df2, VE_lab == "98·9 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.25) +
    
    ## ── VE 0 % (좁은 막대) ─────────────────────────
    geom_col(data = filter(df2, VE_lab == "0 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.15) +
    
    facet_grid(~ setting) +
    
    ## ── 색·범례 설정 ───────────────────────────────
    scale_fill_manual(
      name   = "VE × Supply",
      values = cols6,
      breaks = names(cols6)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    
    ## ── 축·테마 ───────────────────────────────────
    scale_y_continuous(labels = short_si,
                       expand = expansion(mult = c(0, .05))) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position     = "bottom",
      axis.text.x         = element_text(angle = 45, hjust = 1),
      strip.background    = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.x  = element_blank()
    )
}

library(colorspace)
short_si <- function(x) {
  dplyr::case_when(
    abs(x) >= 1e6 ~ paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
    abs(x) >= 1e3 ~ paste0(formatC(x / 1e3, format = "f", digits = 1), "K"),
    TRUE          ~ as.character(x)
  )
}

p1 <- 
make_nnv_overlay_vc(combined_nnv_setting_summ_ve_ixchiq,
                 "nnv",       "NNV to avert a symptomatic case") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a symptomatic case",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) + theme(axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.y  = element_text(size = 10),
           legend.position = "none") 

p2 <-
make_nnv_overlay_vc(combined_nnv_setting_summ_ve_vc,
                    "nnv_fatal",       "NNV to avert fatal case") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a fatal case",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  )+ theme(axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.y  = element_text(size = 10),
           legend.position = "none") 

p3 <- 
make_nnv_overlay_vc(combined_nnv_setting_summ_ve_vc,
                    "nnv_daly",       "NNV to avert a DALY") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a DALY",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) 

nnv_coverage <- 
(p1 / p2 / p3) +                      
  plot_layout(guides = "collect") &                 
  theme(legend.position = "bottom") 

ggsave(filename = "02_Outputs/2_1_Figures/nnv_coverage.jpg", nnv_coverage, width = 9, height = 7, dpi = 1200)

nnv_coverage <- 
  (p1 / p2 / p3) +                      
  plot_layout(guides = "collect") &                 
  theme(legend.position = "bottom") 

### per 1m graph ---------------------------------

region_setting_lut <- tribble(
  ~region , ~setting ,
  "Ceará"     , "High",
  "Bahia"     , "Low",
  "Paraíba"     , "High",
  "Pernambuco"     , "Moderate",
  "Rio Grande do Norte"     , "Low",
  "Piauí"     , "High",
  "Tocantins"     , "Moderate",
  "Alagoas"     , "High",
  "Minas Gerais"     , "Low",
  "Sergipe"     , "Low",
  "Goiás"     , "Low"
)

## 1. nnv_results_coverage_model  ➜  per-1M summary 한데 모으기 -----
per1m_all <- imap_dfr(
  nnv_results_coverage_model,                # ─ region 루프 ─
  function(ve_list, region_nm) {
    imap_dfr(ve_list, function(cov_list, ve_tag) {          # ─ VE 루프 ─
      imap_dfr(cov_list, function(nnv_obj, cov_tag) {       # ─ coverage 루프 ─
        nnv_obj[["per1M_summary"]] %>%                      # ✨ 여기!
          mutate(
            VE = ve_tag,
            VC = cov_tag
          )
      })
    }) %>%
      mutate(region = region_nm)
  }
) %>%
  left_join(region_setting_lut, by = "region")

summ_by_setting <- per1m_all %>% 
  group_by(setting, scenario, target, VE, VC) %>% 
  summarise(
    scenario_vacc   = sum(scenario_vacc, na.rm = TRUE),
    tot_pop         = sum(tot_pop),
    vacc_prop       = scenario_vacc / tot_pop,
    across(starts_with("diff"),  sum, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    per1M_diff     = diff        / scenario_vacc * 1e5,
    per1M_diff_lo  = diff_low    / scenario_vacc * 1e5,
    per1M_diff_hi  = diff_hi     / scenario_vacc * 1e5,
    per1M_fatal    = diff_fatal  / scenario_vacc * 1e5,
    per1M_fatal_lo = diff_fatal_low / scenario_vacc * 1e5,
    per1M_fatal_hi = diff_fatal_hi  / scenario_vacc * 1e5,
    per1M_daly     = diff_daly   / scenario_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low  / scenario_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi   / scenario_vacc * 1e5,
    
    nnv          = scenario_vacc / diff,
    nnv_lo       = scenario_vacc / diff_hi,
    nnv_hi       = scenario_vacc / diff_low,
    
    nnv_fatal    = scenario_vacc / diff_fatal,
    nnv_fatal_lo = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi = scenario_vacc / diff_fatal_low,
    
    nnv_daly     = scenario_vacc / diff_daly,
    nnv_daly_lo  = scenario_vacc / diff_daly_hi,
    nnv_daly_hi  = scenario_vacc / diff_daly_low
  )

## 3. National 합산(모든 setting)------------------------------------
summ_national <- per1m_all %>% 
  group_by(scenario, target, VE, VC) %>% 
  summarise(
    scenario_vacc   = sum(scenario_vacc, na.rm = TRUE),
    tot_pop         = sum(tot_pop),
    vacc_prop       = scenario_vacc / tot_pop,
    across(starts_with("diff"),  sum, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    setting        = "National",
    per1M_diff     = diff        / scenario_vacc * 1e5,
    per1M_diff_lo  = diff_low    / scenario_vacc * 1e5,
    per1M_diff_hi  = diff_hi     / scenario_vacc * 1e5,
    
    per1M_fatal    = diff_fatal  / scenario_vacc * 1e5,
    per1M_fatal_lo = diff_fatal_low / scenario_vacc * 1e5,
    per1M_fatal_hi = diff_fatal_hi  / scenario_vacc * 1e5,
    
    per1M_daly     = diff_daly   / scenario_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low  / scenario_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi   / scenario_vacc * 1e5,
    
    nnv          = scenario_vacc / diff,
    nnv_lo       = scenario_vacc / diff_hi,
    nnv_hi       = scenario_vacc / diff_low,
    
    nnv_fatal    = scenario_vacc / diff_fatal,
    nnv_fatal_lo = scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi = scenario_vacc / diff_fatal_low,
    
    nnv_daly     = scenario_vacc / diff_daly,
    nnv_daly_lo  = scenario_vacc / diff_daly_hi,
    nnv_daly_hi  = scenario_vacc / diff_daly_low
  )

## 4. setting 요약 + National 합치기 -------------------------------
per1m_summary_full <- bind_rows(summ_by_setting, summ_national)

## 5. plot_df 변환 --------------------------------------------------
plot_df <- per1m_summary_full %>% 
  mutate(
    target_of_scenario = case_when(
      scenario == "Scenario_1" ~ "1-11 years",
      scenario == "Scenario_2" ~ "12-17 years",
      scenario == "Scenario_3" ~ "18-59 years",
      scenario == "Scenario_4" ~ "60+ years"
    ),
    effect_group = if_else(
      target == target_of_scenario,
      "Total effect in target group",
      paste0("Indirect — ", target)
    ),
    # factor 레벨 고정
    effect_group = factor(effect_group,
                          levels = c(
                            "Total effect in target group",
                            "Indirect — <1 (novacc)",
                            "Indirect — 1-11 years",
                            "Indirect — 12-17 years",
                            "Indirect — 18-59 years",
                            "Indirect — 60+ years"
                          )
    ),
    scenario = factor(scenario,
                      levels = c("Scenario_1","Scenario_2","Scenario_3","Scenario_4"),
                      labels = c(
                        "Strategy 1 (1–11 years)",
                        "Strategy 2 (12–17 years)",
                        "Strategy 3 (18–59 years)",
                        "Strategy 4 (60+ years)"
                      )
    ),
    setting  = factor(setting, levels = c("National","High","Moderate","Low"))
  )

plot_df_ve989_vc50 <- plot_df %>% filter(VC == "cov50" & VE == "VE98.9")

make_per1m_plot <- function(df, y_var, y_lo, y_hi,
                            y_label,
                            x_label      = "Vaccination strategy",
                            legend_title = "Effect type") {
  
  ## ── 준비 ────────────────────────────────────────────────
  y_sym   <- rlang::sym(y_var)
  ylo_sym <- rlang::sym(y_lo)
  yhi_sym <- rlang::sym(y_hi)
  
  level_order <- c(
    "Total effect in target group",
    "Indirect — <1 (novacc)",
    "Indirect — 1-11 years",
    "Indirect — 12-17 years",
    "Indirect — 18-59 years",
    "Indirect — 60+ years"
  )
  
  df <- df %>% 
    mutate(
      ## 스택 순서 고정
      effect_group = factor(effect_group, levels = level_order),
      
      ## 보기 좋은 레이블 (원하면 여기서 수정)
      VE_lab = factor(VE, levels = c("VE0","VE98.9"),
                      labels = c("VE 0 %", "VE 98·9 %")),
      VC_lab = factor(VC, levels = c("cov10","cov50","cov90"),
                      labels = c("Coverage 10 %",
                                 "Coverage 50 %",
                                 "Coverage 90 %"))
    )
  
  ## ── 합계·UI 계산 ─────────────────────────────────────────
  plot_summary <- df %>%
    group_by(setting, VE_lab, VC_lab, scenario, effect_group) %>%
    summarise(
      mid_val = sum(!!y_sym,   na.rm = TRUE),
      lo_val  = sum(!!ylo_sym, na.rm = TRUE),
      hi_val  = sum(!!yhi_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  error_df <- plot_summary %>%
    group_by(setting, VE_lab, VC_lab, scenario) %>%
    summarise(
      y    = sum(mid_val),
      ymin = sum(lo_val),
      ymax = sum(hi_val),
      .groups = "drop"
    )
  
  ## ── 그래프 ──────────────────────────────────────────────
  ggplot(plot_summary,
         aes(x = scenario, y = mid_val, fill = effect_group)) +
    geom_col(position = "stack", width = 0.8) +
    #geom_errorbar(
    #  data        = error_df,
    #  inherit.aes = FALSE,
    #  aes(x = scenario, y = y, ymin = ymin, ymax = ymax),
    #  width  = 0.25, colour = "black") +
    facet_grid(VC_lab ~ setting + VE_lab, switch = "y") +   
    scale_fill_brewer(
      palette = "Set2",
      name    = legend_title,
      breaks  = level_order,
      labels  = level_order
    ) +
    scale_y_continuous(
      labels = function(x) ifelse(abs(x) >= 1000,
                                  paste0(round(x/1000), "k"), x)) +
    labs(x = x_label, y = y_label) +
    theme_pubclean() +
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text.y.left = element_text(angle = 0)           
    )
}

base_strip <- theme(
  strip.text.x       = element_blank(),
  strip.background.x = element_blank(),
  strip.text.y.left  = element_text(size = 10),
  strip.background.y = element_rect(fill = "grey95", colour = NA) 
)


g1 <- 
make_per1m_plot(
  plot_df_ve989_vc50,
  y_var  = "per1M_diff",
  y_lo   = "per1M_diff_lo",
  y_hi   = "per1M_diff_hi",
  y_label = "Averted symptomatic cases per 100 k doses"
) + theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 10),
          legend.position = "none")
  
g2 <- 
make_per1m_plot(
  plot_df_ve989_vc50,
  y_var  = "per1M_fatal",
  y_lo   = "per1M_fatal_lo",
  y_hi   = "per1M_fatal_hi",
  y_label = "Averted deaths per 100 k doses"
) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 10),
    legend.position = "none"
  ) +
  base_strip   

g3 <- 
make_per1m_plot(
  plot_df_ve989_vc50,
  y_var  = "per1M_daly",
  y_lo   = "per1M_daly_lo",
  y_hi   = "per1M_daly_hi",
  y_label = "Averted DALYs per 100 k doses"
)+ base_strip

per1m_coverage <- 
  #(g1 / g2 / g3) + 
  g3 + 
  plot_layout(guides = "collect") &                 
  theme(legend.position = "bottom") 

p <- 
make_per1m_plot(
  plot_df_ve989_vc50,
  y_var  = "per1M_daly",
  y_lo   = "per1M_daly_lo",
  y_hi   = "per1M_daly_hi",
  y_label = "Averted DALYs per 100 k doses"
)

ggsave(filename = "02_Outputs/2_1_Figures/per1m_coverage.jpg", per1m_coverage, width = 9, height = 13, dpi = 1200)



prevacc_df <- combined_nnv_national_age_ixchiq[1:20,]

prevacc_df$tot_pop <- pop_age_nat$tot_pop
prevacc_df$age_gr_levels <- age_gr_levels

prevacc_df <- prevacc_df %>% mutate(
  preinf_perM = pre_inf / tot_pop * 1e6,
  preinf_perM_lo = pre_inf_lo / tot_pop * 1e6,
  preinf_perM_hi = pre_inf_hi / tot_pop * 1e6,
  prevacc_perM = pre_vacc / tot_pop * 1e6,
  prevacc_perM_lo = pre_vacc_lo / tot_pop * 1e6,
  prevacc_perM_hi = pre_vacc_hi / tot_pop * 1e6,
  
  age_gr_levels = factor(age_gr_levels, levels = age_gr_levels)
)


ggplot(prevacc_df, 
       aes(x = age_gr_levels, y = preinf_perM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(
    aes(ymin = preinf_perM_lo, ymax = preinf_perM_hi),
    width = 0.25, colour = "black"
  ) +
  scale_y_log10()


ggplot(prevacc_df, 
       aes(x = age_gr_levels, y = prevacc_perM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(
    aes(ymin = prevacc_perM_lo, ymax = prevacc_perM_hi),
    width = 0.25, colour = "black"
  ) + labs(
    x = "Agr group",
    y = "Symptomatic cases per 1M"
    ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)   
  )

## dalys by component 
extract_final_summ_df <- function(model) {
  model %>%
    imap(function(ve_cov_list, region_name) {
      ve_cov_list %>%
        imap(function(cov_list, ve_tag) {
          cov_list %>%
            imap(function(postsim_ui, cov_tag) {
              # scenario list 펼치기
              postsim_ui$final_summ %>%
                imap_dfr(function(df, scen_idx) {
                  df %>%
                    mutate(
                      Region    = region_name,
                      VE_tag    = ve_tag,
                      Cov_tag   = cov_tag,
                      Scenario  = as.integer(scen_idx)
                    )
                })
            }) %>% bind_rows()
        }) %>% bind_rows()
    }) %>% bind_rows() %>%
    relocate(Region, VE_tag, Cov_tag, Scenario, AgeGroup)
}

final_summ_df <- extract_final_summ_df(postsim_vc_ixchiq_model)

national_totals <- final_summ_df %>%
  group_by(VE_tag, Cov_tag, Scenario) %>%
  summarise(
    pre_yld_acute     = sum(pre_yld_acute,     na.rm = TRUE),
    yld_acute         = sum(yld_acute,         na.rm = TRUE),
    pre_yld_subacute  = sum(pre_yld_subacute,  na.rm = TRUE),
    yld_subacute      = sum(yld_subacute,      na.rm = TRUE),
    pre_yld_chronic   = sum(pre_yld_chronic,   na.rm = TRUE),
    yld_chronic       = sum(yld_chronic,       na.rm = TRUE),
    pre_yll           = sum(pre_yll,           na.rm = TRUE),
    yll               = sum(yll,               na.rm = TRUE),
    .groups = "drop"
  )

national_components <- national_totals %>%
  transmute(
    VE_tag, Cov_tag, Scenario,
    `YLD (acute)`    = safe_pct(pre_yld_acute,    yld_acute),
    `YLD (subacute)` = safe_pct(pre_yld_subacute, yld_subacute),
    `YLD (chronic)`  = safe_pct(pre_yld_chronic,  yld_chronic),
    YLL              = safe_pct(pre_yll,          yll)
  ) %>%
  pivot_longer(cols = c(`YLD (acute)`,`YLD (subacute)`,`YLD (chronic)`, YLL),
               names_to = "Component", values_to = "PctReduction")


scenario_labels <- c(
  "1" = "Strategy1: 1-11 years old",
  "2" = "Strategy2: 12-17 years old",
  "3" = "Strategy3: 18-59 years old",
  "4" = "Strategy4: 60+ years old"
)

national_components <- national_components %>%
  mutate(
    Strategy = factor(scenario_labels[as.character(Scenario)],
                      levels = scenario_labels),   # facet 순서도 고정
    Component = factor(Component,
                       levels = c("YLD (acute)",
                                  "YLD (subacute)",
                                  "YLD (chronic)",
                                  "YLL"))
  )

plot_national_component_reduction <- function(df,
                                              ve_sel = "VE98.9",
                                              cov_sel = "cov50") {
  df %>%
    filter(VE_tag == ve_sel, Cov_tag == cov_sel) %>%
    ggplot(aes(x = Component, y = PctReduction)) +
    geom_col(fill = "grey60", width = 0.8) +
    geom_text(aes(label = ifelse(is.na(PctReduction), "NA",
                                 sprintf("%.1f%%", PctReduction))),
              vjust = -0.3, size = 3) +
    facet_wrap(~ Strategy, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = NULL,
      y = "% reduction vs pre (national total)",
      title = sprintf("National-level DALY component reductions (VE=%s, Coverage=%s)", ve_sel, cov_sel)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)  # x축 45도 회전
    )
}

plot_national_component_reduction(national_components,
                                  ve_sel = "VE98.9",
                                  cov_sel = "cov50")

totvacc <- combined_nnv_national_age_ixchiq %>% group_by(scenario) %>% summarise(
  tot_vacc = sum(tot_vacc)
)
  
national_totals <- final_summ_df %>%
  group_by(VE_tag, Cov_tag, Scenario) %>%
  summarise(
    pre_yld_acute     = sum(pre_yld_acute,     na.rm = TRUE),
    yld_acute         = sum(yld_acute,         na.rm = TRUE),
    pre_yld_subacute  = sum(pre_yld_subacute,  na.rm = TRUE),
    yld_subacute      = sum(yld_subacute,      na.rm = TRUE),
    pre_yld_chronic   = sum(pre_yld_chronic,   na.rm = TRUE),
    yld_chronic       = sum(yld_chronic,       na.rm = TRUE),
    pre_yll           = sum(pre_yll,           na.rm = TRUE),
    yll               = sum(yll,               na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = paste0("Scenario_", Scenario))

national_per_vaccM <- national_totals %>%
  left_join(totvacc, by = "scenario") %>%
  mutate(denom_perM = ifelse(tot_vacc > 0, tot_vacc/1e6, NA_real_)) %>%
  transmute(
    VE_tag, Cov_tag, Scenario = as.character(Scenario),
    denom_perM,
    `YLD (acute)`    = (pre_yld_acute    - yld_acute)    / denom_perM,
    `YLD (subacute)` = (pre_yld_subacute - yld_subacute) / denom_perM,
    `YLD (chronic)`  = (pre_yld_chronic  - yld_chronic)  / denom_perM,
    YLL              = (pre_yll          - yll)          / denom_perM
  ) %>%
  pivot_longer(cols = c(`YLD (acute)`,`YLD (subacute)`,`YLD (chronic)`, YLL),
               names_to = "Component", values_to = "Value")
scenario_labels <- c(
  "1" = "Strategy1: 1-11 years old",
  "2" = "Strategy2: 12-17 years old",
  "3" = "Strategy3: 18-59 years old",
  "4" = "Strategy4: 60+ years old"
)

national_per_vaccM <- national_per_vaccM %>%
  mutate(
    Strategy  = factor(scenario_labels[Scenario], levels = scenario_labels),
    Component = factor(Component,
                       levels = c("YLD (chronic)","YLD (subacute)","YLD (acute)","YLL"))
  )

plot_national_components_perVaccM <- function(df, ve_sel = "VE98.9", cov_sel = "cov50") {
  df %>%
    filter(VE_tag == ve_sel, Cov_tag == cov_sel) %>%
    ggplot(aes(y = Component, x = Value)) +  # 가로 막대
    geom_col(fill = "grey60", width = 0.72) +
    geom_text(aes(label = scales::comma(round(Value, 1))),
              hjust = -0.1, size = 3) +
    facet_wrap(~ Strategy, nrow = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "Averted DALYs per 1,000,000 vaccinated (national aggregate)",
      y = NULL,
      title = sprintf("DALY component reductions per vaccinated million (VE=%s, Cov=%s)",
                      ve_sel, cov_sel)
    ) +
    theme_classic() +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))
}

# 사용 예
plot_national_components_perVaccM(national_per_vaccM,
                                  ve_sel = "VE98.9",
                                  cov_sel = "cov50")

national_components_pct <- national_totals %>%
  transmute(
    VE_tag, Cov_tag, Scenario = as.character(Scenario),
    `YLD (acute)`    = safe_pct(pre_yld_acute,    yld_acute),
    `YLD (subacute)` = safe_pct(pre_yld_subacute, yld_subacute),
    `YLD (chronic)`  = safe_pct(pre_yld_chronic,  yld_chronic),
    YLL              = safe_pct(pre_yll,          yll)
  ) %>%
  pivot_longer(
    cols = c(`YLD (acute)`,`YLD (subacute)`,`YLD (chronic)`, YLL),
    names_to = "Component", values_to = "PctReduction"
  )

# 3) 전략 라벨 고정
scenario_labels <- c(
  "1" = "Strategy1: 1-11 years old",
  "2" = "Strategy2: 12-17 years old",
  "3" = "Strategy3: 18-59 years old",
  "4" = "Strategy4: 60+ years old"
)

national_components_pct <- national_components_pct %>%
  mutate(
    Strategy  = factor(scenario_labels[Scenario], levels = scenario_labels),
    Component = factor(Component,
                       levels = c("YLD (chronic)", "YLD (subacute)", "YLD (acute)", "YLL"))
  )

# 4) 가로 막대 그래프 (% 감소)
plot_national_components_pct <- function(df, ve_sel = "VE98.9", cov_sel = "cov50") {
  df %>%
    filter(VE_tag == ve_sel, Cov_tag == cov_sel) %>%
    ggplot(aes(y = Component, x = PctReduction)) +
    geom_col(fill = "grey60", width = 0.72) +
    geom_text(aes(label = ifelse(is.na(PctReduction), "NA",
                                 sprintf("%.1f%%", PctReduction))),
              hjust = -0.1, size = 3) +
    facet_wrap(~ Strategy, nrow = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "% reduction vs pre (national total basis)",
      y = NULL,
      title = sprintf("National %% reduction by DALY component (VE=%s, Cov=%s)",
                      ve_sel, cov_sel)
    ) +
    theme_classic() +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))
}

# 실행 예
plot_national_components_pct(national_components_pct, ve_sel = "VE98.9", cov_sel = "cov50")