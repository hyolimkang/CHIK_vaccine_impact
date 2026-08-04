library(writexl)
#vimkunya

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

# for debug 
regions <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 prevacc_ui = prevacc_ui_ce, posterior = posterior_ce,
                 lhs_sample = lhs_vimkunya$Ceará)
)
## 0) 준비 ------------------------------------------------------------------
ve_inf_set     <- c(0, 0.978)          
coverage_set   <- c(0.10, 0.50, 0.9)              
coverage_set   <- 0.5 # test
coverage_set   <- c(0.01, 0.05, 0.07)            
doses <- 2e6                       
region_prop <- as.data.frame(t(region_totals[1,])) %>%
  tibble::rownames_to_column("region") %>%              
  dplyr::rename(pop = "1") %>%                                  
  mutate(
    prop  = pop / sum(pop),                            
    doses = round(prop * doses)                         
  )

doses_by_region <- setNames(region_prop$doses, region_prop$region)


## 1) 시뮬레이션 실행 --------------------------------------------------------
sim_results_vc_coverage_vimkun  <- imap(regions, function(reg_args, reg_name) {
  
  ## ─ VE 수준 루프 ──────────────────────────────────────────────────────────
  imap(set_names(ve_inf_set, paste0("VE", ve_inf_set*100)), function(ve_val, ve_tag) {
    
    ## ─ coverage 루프 ───────────────────────────────────────────────────────
    imap(set_names(coverage_set, paste0("cov", coverage_set*100)),
         function(cov_val, cov_tag) {
           
           # if !ve=0, then use ve_val 
           ve_input <- if (abs(ve_val - 0.978) < 1e-6) NA else ve_val 
           
           run_simulation_scenarios_ui_vimkun(
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
sim_results_vc_dose_vimkun <- imap(regions, function(reg_args, reg_name) {
  
  ## ─ 지역별 고정 coverage 계산 ────────────────────────────
  doses_this  <- doses_by_region[[reg_name]]         
  total_cov_i <- doses_this / sum(reg_args$N)         
  
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

save("sim_results_vc_coverage_vimkun", file = "00_Data/0_2_Processed/sim_results_vc_coverage_vimkun.RData")
save("sim_results_vc_dose_vimkun", file = "00_Data/0_2_Processed/sim_results_vc_dose_vimkun.RData")

postsim_vc_vimkun  <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- sim_results_vc_coverage_vimkun [[region_name]]    # VE×cov 결과 리스트
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
postsim_dose_vimkun  <- imap(regions, function(reg_args, region_name) {
  
  ve_cov_list <- sim_results_vc_dose_vimkun [[region_name]]    # VE×cov 결과 리스트
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

save("postsim_vc_vimkun", file = "00_Data/0_2_Processed/postsim_vc_vimkun.RData")
save("postsim_dose_vimkun", file = "00_Data/0_2_Processed/postsim_dose_vimkun.RData")

n_scenarios = 4
nnv_results_vimkun_coverage <- imap(postsim_vc_vimkun, function(ve_cov_list, region_name) {
  
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
nnv_results_vimkun_dose <- imap(postsim_dose_vimkun, function(ve_cov_list, region_name) {
  
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

n_scenarios <- 4
age_seq <- rep(age_gr[1:length(age_gr_levels)], n_scenarios) 

## 1) nnv_results (지역 → VE → VC → nnv_list 결과) 결합 --------------------------
combined_nnv_df_vimkun_coverage    <- imap_dfr(nnv_results_vimkun_coverage,      # 지역 loop
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

write_xlsx(combined_nnv_df_vimkun_coverage, path = "02_Outputs/2_2_Tables/combined_nnv_df_vimkun_coverage.xlsx")


combined_nnv_national_vimkun_cov <- combined_nnv_df_vimkun_coverage %>%
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
    
    ## per 100k 예방 효과
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

combined_nnv_setting_vimkun_vc <- combined_nnv_df_vimkun_coverage %>% 
  group_by(VE, VC, setting, scenario) %>%           # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_inf:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(VE, setting, scenario) %>%                   # 전체 투여량(시나리오) 계산
  mutate(scenario_vacc = sum(tot_vacc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    nnv_inf      = tot_vacc / diff_inf,
    nnv_inf_lo   = tot_vacc / diff_inf_hi,
    nnv_inf_hi   = tot_vacc / diff_inf_lo,
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

combined_nnv_setting_sum_vimkun_vc <- combined_nnv_setting_vimkun_vc %>% 
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

combined_nnv_national_vimkun_cov$setting <- "National"
combined_nnv_setting_sum_vimkun_vc <- rbind(combined_nnv_national_vimkun_cov, combined_nnv_setting_sum_vimkun_vc)

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
      VE_lab   = factor(VE, levels = c("VE0","VE97.8"),
                        labels = c("0 %","97·8 %")),
      VC_lab   = factor(case_when(
        VC == "cov10"  ~ "10 %",
        VC == "cov50"  ~ "50 %",
        VC == "cov90" ~ "90 %"),
        levels = c("10 %","50 %","90 %")),
      ## ── VE × Coverage 조합 (fill‧범례용) ────────
      VEVC     = factor(paste0(VE_lab, "_", VC_lab),
                        levels = c("0 %_10 %","0 %_50 %","0 %_90 %",
                                   "97·8 %_10 %","97·8 %_50 %","97·8 %_90 %"))
    )
}

## (2) rbind 후 factor 고정
combined_nnv_setting_sum_vimkun_vc <- combined_nnv_setting_sum_vimkun_vc %>%  fix_factors()

save(combined_nnv_setting_sum_vimkun_vc, file = "00_Data/0_2_Processed/combined_nnv_setting_sum_vimkun_vc.RData")

short_si <- function(x) {
  dplyr::case_when(
    abs(x) >= 1e6 ~ paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
    abs(x) >= 1e3 ~ paste0(formatC(x / 1e3, format = "f", digits = 1), "K"),
    TRUE          ~ as.character(x)
  )
}

make_nnv_overlay_vc <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0","VE97.8"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_10 %"    = "#A6CEE3",
    "0 %_50 %"    = "#B2DF8A",
    "0 %_90 %"   = "#FB9A99",
    "97·8 %_10 %" = "#1F78B4",
    "97·8 %_50 %" = "#33A02C",
    "97·8 %_90 %"= "#E31A1C"
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
    geom_col(data = filter(df2, VE_lab == "97·8 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "97·8 %"),
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


make_nnv_overlay_vc(combined_nnv_setting_sum_vimkun_vc,
                    "nnv",       "NNV to avert a symptomatic case") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) + 
  #scale_y_log10()  +
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

#### all states
make_nnv_overlay_all <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0","VE98.9")) %>%
    mutate(
      setting = factor(setting, levels = c("High", "Moderate", "Low")),
      region_label = paste0(region, " (", setting, ")")
    ) %>%
    arrange(setting, region) %>%  # ← 정렬
    mutate(region_label = factor(region_label, levels = unique(region_label)))  # ← 순서 고정
  
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_10 %"    = "#A6CEE3",
    "0 %_50 %"    = "#B2DF8A",
    "0 %_90 %"    = "#FB9A99",
    "98·9 %_10 %" = "#1F78B4",
    "98·9 %_50 %" = "#33A02C",
    "98·9 %_90 %" = "#E31A1C"
  )
  
  dodge_cov <- position_dodge(width = 0.9)
  
  ggplot(df2, aes(x = scenario)) +
    geom_col(data = filter(df2, VE_lab == "0 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.8) +
    geom_errorbar(data = filter(df2, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.25) +
    
    geom_col(data = filter(df2, VE_lab == "98·9 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "98·9 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.15) +
    
    facet_wrap(~ region_label, nrow = 2) +  # ← strip 위에 setting 붙이기
    
    scale_fill_manual(
      name   = "VE × Coverage",
      values = cols6,
      breaks = names(cols6)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    
    scale_y_continuous(labels = short_si,
                       expand = expansion(mult = c(0, .05))) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position     = "bottom",
      axis.text.x         = element_text(angle = 45, hjust = 1),
      strip.background    = element_rect(fill = "grey95", colour = NA),
      strip.text          = element_text(face = "bold"),
      panel.grid.major.x  = element_blank()
    )
}

combined_nnv_df_vimkun_all <- combined_nnv_df_vimkun_coverage %>%
  group_by(VE, VC, scenario, region, setting) %>%            # VE·VC·scenario 기준으로 그룹화
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
    
    ## per 100k 예방 효과
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

combined_nnv_df_vimkun_all <- combined_nnv_df_vimkun_all %>% fix_factors()

combined_nnv_df_ixchiq_all <- combined_nnv_df_region_coverage_model %>%
  group_by(VE, VC, scenario, region, setting) %>%            # VE·VC·scenario 기준으로 그룹화
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
    
    ## per 100k 예방 효과
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

combined_nnv_df_ixchiq_all <- combined_nnv_df_ixchiq_all %>% fix_factors()


p1 <- 
make_nnv_overlay_all(combined_nnv_df_ixchiq_all,
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
  make_nnv_overlay_all(combined_nnv_df_ixchiq_all,
                       "nnv_fatal",       "NNV to avert a death") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a death",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) + theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y  = element_text(size = 10),
            legend.position = "none") 
p3 <- 
  make_nnv_overlay_all(combined_nnv_df_ixchiq_all,
                       "nnv_daly",       "NNV to avert a DALY") 

nnv_coverage_ixchiq_all <- 
  (p1 / p2 / p3) +                      # plot arrangement
  plot_layout(guides = "collect") +    # collect legends
  plot_annotation(tag_levels = list(c("A","B","C"))) & 
  theme(legend.position = "bottom")    # apply theme to all

ggsave(filename = "02_Outputs/2_1_Figures/figs8_4.jpg", p1, width = 10, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs8_5.jpg", p2, width = 10, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs8_6.jpg", p3, width = 10, height = 6, dpi = 1200)

# -------------------------------------------------------------------------------
# combining two vaccines
# -------------------------------------------------------------------------------

# 1. Ixchiq simulation results
load("00_Data/0_2_Processed/sim_results_vc_ixchiq_model.RData")
load("00_Data/0_2_Processed/sim_results_vc_coverage_vimkun.RData")
load("00_Data/0_2_Processed/combined_nnv_setting_summ_ve_ixchiq.RData")
load("00_Data/0_2_Processed/combined_nnv_setting_sum_vimkun_vc.RData")

combined_nnv_setting_summ_ve_ixchiq$vaccine <- "Ixchiq"
combined_nnv_setting_sum_vimkun_vc$vaccine <- "Vimkunya"

combined_nnv_two_vacc <- rbind(combined_nnv_setting_sum_vimkun_vc, combined_nnv_setting_summ_ve_ixchiq)

save(combined_nnv_two_vacc, file = "00_Data/0_2_Processed/combined_nnv_two_vacc.RData")

make_nnv_overlay_two_vacc <- function(df, y_var, y_lab) {
  
  library(dplyr); library(ggplot2); library(rlang); library(scales)
  
  ## ── 1. 준비: 사용할 행 ──────────────────────────────
  df2 <- df %>%                 # VE0 + “각 백신의 최고 VE“만 남김
    group_by(vaccine) %>%       
    filter(VE_lab == "0 %" | VE == max(VE)) %>% 
    ungroup()
  
  ## ── 2. 심볼 ─────────────────────────────────────────
  y_sym  <- sym(y_var)
  lo_sym <- sym(paste0(y_var, "_lo"))
  hi_sym <- sym(paste0(y_var, "_hi"))
  
  ## ── 3. 팔레트(동적으로) ────────────────────────────
  #   ▸ VE × coverage 조합별(fill) 팔레트
  fill_lvls <- sort(unique(df2$VEVC))
  fill_cols <- setNames(
    hue_pal()(length(fill_lvls)),  # 필요 개수만큼 자동 색상
    fill_lvls
  )
  #   ▸ 백신별(colour) 테두리 팔레트
  vac_lvls  <- sort(unique(df2$vaccine))
  vac_cols  <- setNames(c("#000000", "#565656", "#939393")[seq_along(vac_lvls)],
                        vac_lvls)
  
  ## ── 4. 가로 위치 dodge 정의 ───────────────────────
  dodge_cov_vac <- position_dodge(width = 0.9)
  
  ## ── 5. 플롯 ───────────────────────────────────────
  ggplot(df2, aes(x = scenario,
                  group  = interaction(VC_lab, vaccine))) +
    
    # (A) VE 0 %  ─ 넓은 막대
    geom_col(
      data      = filter(df2, VE_lab == "0 %"),
      aes(y     = !!y_sym,
          fill  = VEVC,
          colour= vaccine),
      width     = 0.8,
      position  = dodge_cov_vac
    ) +
    geom_errorbar(
      data      = filter(df2, VE_lab == "0 %"),
      aes(ymin  = !!lo_sym,
          ymax  = !!hi_sym,
          colour= vaccine),
      width     = 0.25,
      position  = dodge_cov_vac
    ) +
    
    # (B) 각 백신의 “최고 VE” ─ 좁은 막대
    geom_col(
      data      = filter(df2, VE_lab != "0 %"),
      aes(y     = !!y_sym,
          fill  = VEVC,
          colour= vaccine),
      width     = 0.5,
      position  = dodge_cov_vac
    ) +
    geom_errorbar(
      data      = filter(df2, VE_lab != "0 %"),
      aes(ymin  = !!lo_sym,
          ymax  = !!hi_sym,
          colour= vaccine),
      width     = 0.15,
      position  = dodge_cov_vac
    ) +
    
    facet_grid(~ setting) +
    
    scale_fill_manual(
      name   = "VE × Coverage",
      values = fill_cols,
      breaks = names(fill_cols)
    ) +
    scale_colour_manual(
      name   = "Vaccine",
      values = vac_cols
    ) +
    guides(
      fill   = guide_legend(nrow = 2, byrow = TRUE),
      colour = guide_legend(order = 1)
    ) +
    
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, .05))
    ) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position     = "bottom",
      axis.text.x         = element_text(angle = 45, hjust = 1),
      strip.background    = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.x  = element_blank()
    )
}
make_nnv_overlay_two_vacc(
  df     = combined_nnv_two_vacc,
  y_var  = "nnv",                 # 또는 "nnv_fatal", "nnv_daly"
  y_lab  = "NNV to avert a symptomatic case"
)


make_nnv_overlay_maxVE_nested_shadeVE <- function(
    df, y_var, y_lab,
    summary_fun  = mean,
    light_factor = 0.45,  # 높은 VE → 밝게
    dark_factor  = 0.25   # 낮은 VE → 짙게
) {
  library(dplyr); library(ggplot2); library(rlang); library(scales); library(stringr)
  
  ## 1 ─ ‘최고 VE’ 행만 -------------------------------------------------
  df_max <- df %>%
    group_by(vaccine) %>%
    filter(VE == max(VE, na.rm = TRUE)) %>%
    ungroup()
  
  ## 2 ─ 넓은↔좁은 막대용 백신 정렬 ------------------------------------
  vac_order <- df_max %>%
    group_by(vaccine) %>%
    summarise(score = summary_fun(.data[[y_var]], na.rm = TRUE)) %>%
    arrange(desc(score)) %>%
    pull(vaccine)
  df_max$vaccine <- factor(df_max$vaccine, levels = vac_order)
  
  ## 3 ─ coverage별 기본색 ---------------------------------------------
  base_cov_cols <- c("10 %" = "#1F78B4",   # 파랑
                     "50 %" = "#33A02C",   # 초록
                     "90 %" = "#E31A1C")   # 빨강
  
  lighten <- function(col, f = 0.4) {
    x <- col2rgb(col)/255
    rgb(x[1] + (1-x[1])*f, x[2] + (1-x[2])*f, x[3] + (1-x[3])*f)
  }
  darken  <- function(col, f = 0.25) {
    x <- col2rgb(col)/255
    rgb(x[1]*(1-f), x[2]*(1-f), x[3]*(1-f))
  }
  
  ## 4 ─ VE 숫자 추출 & 진하기 결정 ------------------------------------
  df_max <- df_max %>%
    mutate(
      # "97·8 %" -> "97.8" -> 97.8
      VE_num = readr::parse_number(str_replace_all(VE, "·", ".")),
    ) %>%
    group_by(VC_lab) %>%
    mutate(
      rank_desc = rank(-VE_num, ties.method = "first"), # 1 = 낮은 VE
      fill_col  = ifelse(
        rank_desc == 1,
        darken(base_cov_cols[VC_lab], dark_factor),     # 낮은 VE → 짙게
        lighten(base_cov_cols[VC_lab], light_factor)    # 높은 VE → 밝게
      )
    ) %>%
    ungroup()
  
  fill_cols <- setNames(df_max$fill_col, df_max$VEVC)
  
  ## 5 ─ 그래프 요소 ----------------------------------------------------
  vac_cols <- setNames(c("#555555", "#000000")[seq_along(vac_order)],
                       vac_order)
  
  y_sym  <- sym(y_var)
  lo_sym <- sym(paste0(y_var, "_lo"))
  hi_sym <- sym(paste0(y_var, "_hi"))
  dodge_cov <- position_dodge(width = 0.9)
  
  ## 6 ─ 플롯 -----------------------------------------------------------
  ggplot(df_max, aes(x = scenario,
                     y = !!y_sym,
                     fill = VEVC,
                     group = interaction(VC_lab, vaccine))) +
    
    # 넓은 막대 (값 큰 백신)
    geom_col(data = filter(df_max, vaccine == vac_order[1]),
             aes(colour = vaccine),
             width = 0.75, alpha = 0.9, position = dodge_cov) +
    geom_errorbar(data = filter(df_max, vaccine == vac_order[1]),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, colour = vaccine),
                  width = 0.18, position = dodge_cov) +
    
    # 좁은 막대 (값 작은 백신)
    geom_col(data = filter(df_max, vaccine == vac_order[2]),
             aes(colour = vaccine),
             width = 0.45, position = dodge_cov) +
    geom_errorbar(data = filter(df_max, vaccine == vac_order[2]),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, colour = vaccine),
                  width = 0.10, position = dodge_cov) +
    
    facet_grid(~ setting) +
    scale_fill_manual(name = "Coverage × VE", values = fill_cols,
                      breaks = names(fill_cols)) +
    scale_colour_manual(name = "Vaccine", values = vac_cols) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE),
           colour = guide_legend(order = 1)) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
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

make_nnv_overlay_maxVE_nested_shadeVE(
  df      = combined_nnv_two_vacc,  # 원본 데이터 프레임
  y_var   = "nnv",                            # y 축 변수
  y_lab   = "NNV to avert a symptomatic case"
) + scale_y_log10()


### both vimkunya + ixchiq
short_si <- function(x) {
  dplyr::case_when(
    abs(x) >= 1e6 ~ paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
    abs(x) >= 1e3 ~ paste0(formatC(x / 1e3, format = "f", digits = 1), "K"),
    TRUE          ~ as.character(x)
  )
}

make_nnv_overlay_two_vacc <- function(df, y_var, y_lab) {
  library(dplyr); library(ggplot2); library(rlang); library(scales)
  
  ## 1 ─ VE0 + 최고 VE만 ---------------------------------------------
  df2 <- df %>% 
    group_by(vaccine) %>% 
    filter(VE_lab == "0 %" | VE == max(VE)) %>% 
    ungroup() %>% 
    mutate(
      cov_key = factor(gsub("[^0-9]", "", VC_lab), levels = c("10","50","90")),
      VE_type = factor(
        ifelse(VE_lab == "0 %",
               "Disease blocking only",         
               "Disease and infection blocking"),      
        levels = c("Disease blocking only", "Disease and infection blocking")
      )
    )
  
  ## 2 ─ 심볼 ----------------------------------------------------------
  y_sym  <- sym(y_var)
  lo_sym <- sym(paste0(y_var, "_lo"))
  hi_sym <- sym(paste0(y_var, "_hi"))
  
  ## 3 ─ 팔레트 & 위치 -------------------------------------------------
  base_cov_cols <- c("10" = "#1F78B4", "50" = "#33A02C", "90" = "#E31A1C")
  dodge_cov <- position_dodge(width = 0.9)
  
  ## 4 ─ 플롯 ----------------------------------------------------------
  ggplot(df2, aes(x = scenario,
                  y = !!y_sym,
                  fill  = cov_key,      # 색 = Coverage
                  alpha = VE_type,      # 투명도 = VE
                  group = cov_key)) +
    
    # 넓은 막대 (No VE)
    geom_col(data = filter(df2, VE_type == "Disease blocking only"),
             width = 0.8, position = dodge_cov, colour = "black") +
    geom_errorbar(data = filter(df2, VE_type == "Disease blocking only"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym),
                  width = 0.22, position = dodge_cov) +
    
    # 좁은 막대 (Max VE)
    geom_col(data = filter(df2, VE_type == "Disease and infection blocking"),
             width = 0.5, position = dodge_cov, colour = "black") +
    geom_errorbar(data = filter(df2, VE_type == "Disease and infection blocking"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym),
                  width = 0.12, position = dodge_cov) +
    
    facet_grid(vaccine ~ setting, switch = "y") +
    
    ## 범례 ------------------------------------------------------------
  scale_fill_manual(
    name   = "Vaccine coverage",
    values = base_cov_cols,
    breaks = c("10","50","90"),
    labels = c("10 %","50 %","90 %")
  ) +
    scale_alpha_manual(
      name   = "Vaccine protection",
      values = c("Disease blocking only"  = 0.35,
                 "Disease and infection blocking" = 0.90),
      breaks = c("Disease blocking only", "Disease and infection blocking"),
      guide = "legend"
    ) +
    
    ## 축·테마 ---------------------------------------------------------
  scale_y_continuous(labels = short_si,
                     expand = expansion(mult = c(0, .05))) +
    labs(
      x       = NULL,
      y       = y_lab,
      caption = "Bar width: wide = Disease blocking only, narrow = Disease and infection blocking"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position     = "bottom",
      legend.box          = "vertical",
      #axis.text.x         = element_text(angle = 45, hjust = 1),
      strip.placement     = "inside",
      strip.background    = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.x  = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 13),
      legend.title        = element_text(size = 12), 
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.caption = element_text(size = 11, hjust = 1)
    )
}


combined_nat_twovacc <- combined_nnv_two_vacc %>% filter(setting == "National")

p1 <- 
make_nnv_overlay_two_vacc(
  df     = combined_nat_twovacc,
  y_var  = "nnv",
  y_lab  = "NNV to avert a symptomatic case"
)+ theme(axis.title.x = element_blank(),
         axis.text.x  = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.y  = element_text(size = 12),
         legend.position = "none")+
  labs(caption = NULL) + 
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  )

p2<-
make_nnv_overlay_two_vacc(
  df     = combined_nat_twovacc,
  y_var  = "nnv_fatal",
  y_lab  = "NNV to avert a death"
)+
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a death",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  )+ theme(axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text.y  = element_text(size = 12),
           legend.position = "none")+
  labs(caption = NULL) + 
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  )

p3<-
make_nnv_overlay_two_vacc(
  df     = combined_nat_twovacc,
  y_var  = "nnv_daly",
  y_lab  = "NNV to avert a DALY"
)+
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    x        = "Age specific vaccination strategy",
    y        = "NNV to avert a DALY",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) + theme(axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 14)) + 
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  )

legend_guides <- guides(
  fill  = guide_legend(title = "Vaccine coverage",
                       order = 1, nrow = 1, byrow = TRUE),
  alpha = guide_legend(title = "Vaccine protection",
                       order = 2, nrow = 1, byrow = TRUE)
)

nnv_coverage_twovacc <-
  (p1 / p2 / p3) +
  plot_annotation(
    tag_levels = list(c("A","B","C"))
    #,caption = "Bar width: wide = Disease blocking VE, "
    #~ "narrow = Disease and infection blocking VE"
  ) + 
  plot_layout(guides = "collect") &      
  legend_guides &                       
  theme(
    legend.position      = "bottom",
    legend.box           = "vertical",      
    legend.justification = "left",     
    legend.box.just      = "left",
    legend.text.align    = 0,
    legend.spacing.y     = unit(0, "pt"),
    legend.text          = element_text(size = 13) ,   
    legend.title         = element_text(size = 13)
  )  &
  theme(                              
    plot.caption = element_text(
      hjust  = 0,                  
      margin = margin(l = 6),        
      size   = 13                   
    ),
    plot.tag = element_text(      # ★ ABC 크기·굵기
      size = 13,                  # ← 원하는 크기로 (예: 20)
      face = "bold"
    )
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig2.jpg", nnv_coverage_twovacc, width = 9.5, height = 11, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/fig2.pdf", nnv_coverage_twovacc, width = 9.5, height = 11, dpi = 300)

# per 1m df for vimkunya -------------------------------------------------------

per1m_all_vimk <- imap_dfr(
  nnv_results_vimkun_coverage,                # ─ region 루프 ─
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

summ_by_setting_vimk <- per1m_all_vimk %>% 
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
summ_national_vimk <- per1m_all_vimk %>% 
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
per1m_summary_full_vimk <- bind_rows(summ_by_setting_vimk, summ_national_vimk)


plot_df_vimk <- per1m_summary_full_vimk %>% 
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

## combining vimk + ixchiq
plot_df$vaccine <- "Ixchiq"
plot_df_vimk$vaccine <- "Vimkunya"

plotdf_twovacc <- rbind(plot_df_vimk, plot_df)
plot_df_ve989 <- plotdf_twovacc %>% filter(VE!= "VE0" & VC == "cov50")

make_per1m_plot <- function(df, y_var, y_lo, y_hi,
                            y_label,
                            x_label      = "Vaccination strategy",
                            legend_title = "Effect type") {
  library(dplyr); library(ggplot2); library(rlang); library(scales)
  
  ## ── 심볼 준비 ───────────────────────────────────────
  y_sym   <- sym(y_var)
  ylo_sym <- sym(y_lo)
  yhi_sym <- sym(y_hi)
  
  ## ── 원하는 표시·순서 정의 ───────────────────────────
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
      ## 백신 이름‧순서(원하는 대로 수정)
      vaccine = factor(vaccine,
                       levels = c("Ixchiq", "Vimkunya")),
      
      ## 스택 순서 고정
      effect_group = factor(effect_group, levels = level_order),
      
      ## 보기 좋은 레이블
      VE_lab = factor(VE,
                      levels = c("VE0", "VE98.9"),
                      labels = c("VE 0 %", "VE 98·9 %")),
      VC_lab = factor(VC,
                      levels = c("cov10", "cov50", "cov90"),
                      labels = c("Coverage 10 %",
                                 "Coverage 50 %",
                                 "Coverage 90 %"))
    )
  
  ## ── 합계·UI 계산 ────────────────────────────────────
  plot_summary <- df %>%
    group_by(vaccine, setting, VE_lab, VC_lab, scenario, effect_group) %>%
    summarise(
      mid_val = sum(!!y_sym,   na.rm = TRUE),
      lo_val  = sum(!!ylo_sym, na.rm = TRUE),
      hi_val  = sum(!!yhi_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  error_df <- plot_summary %>%
    group_by(vaccine, setting, VE_lab, VC_lab, scenario) %>%
    summarise(
      y    = sum(mid_val),
      ymin = sum(lo_val),
      ymax = sum(hi_val),
      .groups = "drop"
    )
  
  ## ── 그래프 ─────────────────────────────────────────
  ggplot(plot_summary,
         aes(x = scenario, y = mid_val, fill = effect_group)) +
    geom_col(position = "stack", width = 0.8) +
    geom_errorbar(
      data        = error_df,
      inherit.aes = FALSE,
      aes(x = scenario, y = y, ymin = ymin, ymax = ymax),
      width  = 0.25, colour = "black") +
    
    ## ★ 행 = Coverage, 열 = 백신 → Setting → VE
    #facet_grid(VC_lab ~ vaccine + setting,
    #         switch = "y") +
    facet_grid(. ~ vaccine + setting)+
    #ggforce::facet_nested(VC_lab ~ vaccine + setting, switch = "y") +
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
      legend.position   = "bottom",
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text.y.left = element_text(angle = 0)   # 행 strip 가독성
    )
}

g1 <- 
make_per1m_plot(
  plot_df_ve989,
  y_var  = "per1M_diff",
  y_lo   = "per1M_diff_lo",
  y_hi   = "per1M_diff_hi",
  y_label = "Averted symptomatic cases per 100 k doses"
) + theme_pubclean(base_size = 8) + theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 8),
          legend.position = "none")
g2 <- 
make_per1m_plot(
  plot_df_ve989,
  y_var  = "per1M_fatal",
  y_lo   = "per1M_fatal_lo",
  y_hi   = "per1M_fatal_hi",
  y_label = "Averted deaths per 100 k doses"
) + theme_pubclean(base_size = 8) + theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 8),
          legend.position = "none",
          strip.text.x    = element_blank(),       # 상단 텍스트 숨기기
          strip.background = element_blank() )+
  labs(caption = NULL) 
g3 <- 
make_per1m_plot(
  plot_df_ve989,
  y_var  = "per1M_daly",
  y_lo   = "per1M_daly_lo",
  y_hi   = "per1M_daly_hi",
  y_label = "Averted DALYs per 100 k doses"
) +
  ggtitle(NULL) +
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "Averted DALYs per 100 k doses",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) +
  theme_pubclean(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.text.x     = element_blank(),
    strip.background = element_blank(),
    legend.position  = "none"
  )

per100k_coverage_twovacc <- 
  (g1 / g2 / g3) +                      
  plot_layout(guides = "collect") &                 
  theme(legend.position = "bottom") 

ggsave(filename = "02_Outputs/2_1_Figures/per100k_coverage_twovacc.jpg", per100k_coverage_twovacc, width = 11, height = 10, dpi = 1200)

# pub table ---------------------------------------------------------
library(glue)
metrics <- c("pre_inf", "post_inf",
             "pre_vacc", "post_vacc", "pre_fatal",
             "pre_daly", "post_fatal", "post_daly",
             "diff", "diff_fatal", "diff_daly",
             "nnv", "nnv_fatal", "nnv_daly")

table_nnv_ci_ve <- combined_nnv_two_vacc %>%
  rowwise() %>%
  mutate(
    diff_ci = glue(
      "{formatC(diff,        format='f', digits=0, big.mark=',')} ",
      "({formatC(diff_low,  format='f', digits=0, big.mark=',')}–",
      "{formatC(diff_hi,   format='f', digits=0, big.mark=',')})"
    ),
    preinf_ci = glue(
      "{formatC(pre_inf,        format='f', digits=0, big.mark=',')} ",
      "({formatC(pre_inf_lo,  format='f', digits=0, big.mark=',')}–",
      "{formatC(pre_inf_hi,   format='f', digits=0, big.mark=',')})"
    ),
    postinf_ci = glue(
      "{formatC(post_inf,        format='f', digits=0, big.mark=',')} ",
      "({formatC(post_inf_lo,  format='f', digits=0, big.mark=',')}–",
      "{formatC(post_inf_hi,   format='f', digits=0, big.mark=',')})"
    ),
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
    ),
    nnv_inf_ci = glue(
      "{formatC(nnv_inf,      format='f', digits=0, big.mark=',')} ",
      "({formatC(nnv_inf_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(nnv_inf_hi,  format='f', digits=0, big.mark=',')})"
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
    ),
    per1M_inf_ci = glue(
      "{formatC(per1M_inf,       format='f', digits=0, big.mark=',')} ",
      "({formatC(per1M_inf_lo, format='f', digits=0, big.mark=',')}–",
      "{formatC(per1M_inf_hi,  format='f', digits=0, big.mark=',')})"
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
    )
  ) %>%
  ungroup() %>%
  select(
    scenario, setting, VE, vaccine, VC,
    ends_with("_ci")
  )

table_nnv_ci_twovacc_v50 <- table_nnv_ci_ve %>% filter(VC == "cov50")
colnames(table_nnv_ci_twovacc_v50) <- c("Vaccination strategy",       
                                        "Setting",        
                                        "Vaccine efficacy",             
                                        "Vaccine type",        
                                        "Vaccine coverage",             
                                        "diff_ci",
                                        "pre_infection",
                                        "post_infection",
                                        "prevacc_ci",    
                                        "pre_fatal_ci",   
                                        "pre_daly_ci",    
                                        "postvacc_ci",    
                                        "post_fatal_ci",  
                                        "post_daly_ci", 
                                        "NNV (Infection)",
                                        "NNV (symptomatic)",         
                                        "NNV (fatal)",  
                                        "NNV (DALY)",    
                                        "Infections averted per 100,000 doses",
                                        "Cases averted per 100,000 doses",  
                                        "Deaths averted per 100,000 doses", 
                                        "DALYs averted per 100,000 doses" )

write_xlsx(table_nnv_ci_twovacc_v50, path = "02_Outputs/2_2_Tables/table_nnv_ci_twovacc_v50.xlsx")




