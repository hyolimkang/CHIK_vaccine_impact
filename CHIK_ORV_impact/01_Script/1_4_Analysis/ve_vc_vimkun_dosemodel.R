# vimkun dose model
## 1) nnv_results (지역 → VE → VC → nnv_list 결과) 결합 --------------------------
combined_nnv_df_vimkun_dose    <- imap_dfr(nnv_results_vimkun_dose,      # 지역 loop
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

combined_nnv_national_vimkun_dose <- combined_nnv_df_vimkun_dose %>%
  group_by(VE, VC, scenario) %>%            # VE·VC·scenario 기준으로 그룹화
  summarise(
    across(c(tot_vacc, pre_vacc:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ## NNV
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

combined_nnv_setting_vimkun_dose <- combined_nnv_df_vimkun_dose %>% 
  group_by(VE, VC, setting, scenario) %>%           # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_vacc:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(VE, setting, scenario) %>%                   # 전체 투여량(시나리오) 계산
  mutate(scenario_vacc = sum(tot_vacc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    nnv         = scenario_vacc / diff,
    nnv_lo      = scenario_vacc / diff_hi,
    nnv_hi      = scenario_vacc / diff_low,
    nnv_fatal   = scenario_vacc / diff_fatal,
    nnv_fatal_lo= scenario_vacc / diff_fatal_hi,
    nnv_fatal_hi= scenario_vacc / diff_fatal_low,
    nnv_daly    = scenario_vacc / diff_daly,
    nnv_daly_lo = scenario_vacc / diff_daly_hi,
    nnv_daly_hi = scenario_vacc / diff_daly_low,
    
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

combined_nnv_setting_sum_vimkun_dose <- combined_nnv_setting_vimkun_dose %>% 
  group_by(VE, VC, setting, scenario) %>%                   # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_vacc:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    nnv         = tot_vacc / diff,
    nnv_lo      = tot_vacc / diff_hi,
    nnv_hi      = tot_vacc / diff_low,
    nnv_fatal   = tot_vacc / diff_fatal,
    nnv_fatal_lo= tot_vacc / diff_fatal_hi,
    nnv_fatal_hi= tot_vacc / diff_fatal_low,
    nnv_daly    = tot_vacc / diff_daly,
    nnv_daly_lo = tot_vacc / diff_daly_hi,
    nnv_daly_hi = tot_vacc / diff_daly_low,
    
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

combined_nnv_national_vimkun_dose$setting <- "National"
combined_nnv_setting_sum_vimkun_dose <- rbind(combined_nnv_national_vimkun_dose, 
                                              combined_nnv_setting_sum_vimkun_dose)

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
        VC == "cov1"  ~ "1 %",
        VC == "cov5"  ~ "5 %",
        VC == "cov7" ~ "7 %"),
        levels = c("1 %","5 %","7 %")),
      ## ── VE × Coverage 조합 (fill‧범례용) ────────
      VEVC     = factor(paste0(VE_lab, "_", VC_lab),
                        levels = c("0 %_1 %","0 %_5 %","0 %_7 %",
                                   "97·8 %_1 %","97·8 %_5 %","97·8 %_7 %"))
    )
}

## (2) rbind 후 factor 고정
combined_nnv_setting_sum_vimkun_dose <- combined_nnv_setting_sum_vimkun_dose %>%  fix_factors()

save(combined_nnv_setting_sum_vimkun_vc, file = "00_Data/0_2_Processed/combined_nnv_setting_sum_vimkun_vc.RData")

make_nnv_overlay_dose <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0","VE97.8"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_1 %"    = "#A6CEE3",
    "0 %_5 %"    = "#B2DF8A",
    "0 %_7 %"   = "#FB9A99",
    "97·8 %_1 %" = "#1F78B4",
    "97·8 %_5 %" = "#33A02C",
    "97·8 %_7 %"= "#E31A1C"
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

combined_nnv_dose_vc50 <- combined_nnv_setting_sum_vimkun_dose %>% filter(VC == "cov5" & VE == "VE97.8")

make_nnv_overlay_dose(combined_nnv_dose_vc50,
                    "nnv_fatal",       "NNV to avert a datal case") 
  +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a fatal case",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) + theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y  = element_text(size = 10),
            legend.position = "none") 

