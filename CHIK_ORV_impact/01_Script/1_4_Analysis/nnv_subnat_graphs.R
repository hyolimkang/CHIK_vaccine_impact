make_nnv_overlay_subnat <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0", "VE98.9"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_10 %"     = "#A6CEE3",
    "0 %_50 %"     = "#B2DF8A",
    "0 %_90 %"     = "#FB9A99",
    "98·9 %_10 %"  = "#1F78B4",
    "98·9 %_50 %"  = "#33A02C",
    "98·9 %_90 %"  = "#E31A1C"
  )
  
  dodge_cov <- position_dodge(width = 0.9)
  
  ggplot(df2, aes(x = scenario)) +
    
    # VE 0 %
    geom_col(data = filter(df2, VE_lab == "0 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.8) +
    geom_errorbar(data = filter(df2, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.25) +
    
    # VE 98·9 %
    geom_col(data = filter(df2, VE_lab == "98·9 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "98·9 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.15) +
    
    # 2행 facet
    facet_wrap(~ setting + region, nrow = 2) +
    
    # 색상 및 범례
    scale_fill_manual(
      name   = "VE × Coverage",
      values = cols6,
      breaks = names(cols6)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    
    # y축 log10
    scale_y_log10(
      labels = scales::comma,
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

make_nnv_overlay_subnat_vimkun <- function(df, y_var, y_lab) {
  
  df2 <- df %>% 
    filter(VE %in% c("VE0", "VE97.8"))
  
  y_sym  <- rlang::sym(y_var)
  lo_sym <- rlang::sym(paste0(y_var, "_lo"))
  hi_sym <- rlang::sym(paste0(y_var, "_hi"))
  
  cols6 <- c(
    "0 %_10 %"     = "#A6CEE3",
    "0 %_50 %"     = "#B2DF8A",
    "0 %_90 %"     = "#FB9A99",
    "97·8 %_10 %"  = "#1F78B4",
    "97·8 %_50 %"  = "#33A02C",
    "97·8 %_90 %"  = "#E31A1C"
  )
  
  dodge_cov <- position_dodge(width = 0.9)
  
  ggplot(df2, aes(x = scenario)) +
    
    # VE 0 %
    geom_col(data = filter(df2, VE_lab == "0 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.8) +
    geom_errorbar(data = filter(df2, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.25) +
    
    # VE 98·9 %
    geom_col(data = filter(df2, VE_lab == "97·8 %"),
             aes(y = !!y_sym, fill = VEVC, group = VC_lab),
             position = dodge_cov, width = 0.5) +
    geom_errorbar(data = filter(df2, VE_lab == "97·8 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym, group = VC_lab),
                  position = dodge_cov, width = 0.15) +
    
    # 2행 facet
    facet_wrap(~ setting + region, nrow = 2) +
    
    # 색상 및 범례
    scale_fill_manual(
      name   = "VE × Coverage",
      values = cols6,
      breaks = names(cols6)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    
    # y축 log10
    scale_y_log10(
      labels = scales::comma,
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

combined_nnv_df_region_coverage_model <- combined_nnv_df_region_coverage_model %>% 
  group_by(VE, VC, setting, scenario, region) %>%                   # ★ VE 추가
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
    nnv_daly_hi = tot_vacc / diff_daly_low
  )
    
combined_nnv_df_region_coverage_model <- fix_factors(combined_nnv_df_region_coverage_model)

## vimkun
combined_nnv_df_vimkun_coverage <- combined_nnv_df_vimkun_coverage %>% 
  group_by(VE, VC, setting, scenario, region) %>%                   # ★ VE 추가
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
    nnv_daly_hi = tot_vacc / diff_daly_low
  )

fix_factors_vimkun <- function(df){
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

combined_nnv_df_vimkun_coverage <- fix_factors_vimkun(combined_nnv_df_vimkun_coverage)

write_xlsx(combined_nnv_df_region_coverage_model, path = "02_Outputs/2_2_Tables/table_s5.xlsx")
write_xlsx(combined_nnv_df_vimkun_coverage, path = "02_Outputs/2_2_Tables/table_s6.xlsx")

p1 <- 
make_nnv_overlay_subnat(combined_nnv_df_region_coverage_model,
                    "nnv",       "NNV to avert a symptomatic case") +
  ggtitle(NULL) +                         # <- removes any existing title
  #scale_y_continuous(labels = short_si) +
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
  make_nnv_overlay_subnat(combined_nnv_df_region_coverage_model,
                      "nnv_fatal",       "NNV to avert fatal case") +
  ggtitle(NULL) +                         # <- removes any existing title
  #scale_y_continuous(labels = short_si) +
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
  make_nnv_overlay_subnat(combined_nnv_df_region_coverage_model,
                      "nnv_daly",       "NNV to avert a DALY") +
  ggtitle(NULL) +                         # <- removes any existing title
  #scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a DALY",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  ) 


ggsave(filename = "02_Outputs/2_1_Figures/figs10_p1.jpg", p1, width = 12, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs10_p2.jpg", p2, width = 12, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs10_p3.jpg", p3, width = 12, height = 7, dpi = 1200)


# vimkun
p1 <- 
  make_nnv_overlay_subnat_vimkun(combined_nnv_df_vimkun_coverage ,
                          "nnv",       "NNV to avert a symptomatic case")
p2 <-
  make_nnv_overlay_subnat_vimkun(combined_nnv_df_vimkun_coverage,
                          "nnv_fatal",       "NNV to avert fatal case")

p3 <- 
  make_nnv_overlay_subnat_vimkun(combined_nnv_df_vimkun_coverage,
                          "nnv_daly",       "NNV to avert a DALY") 


ggsave(filename = "02_Outputs/2_1_Figures/figs11_p1.jpg", p1, width = 12, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs11_p2.jpg", p2, width = 12, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs11_p3.jpg", p3, width = 12, height = 7, dpi = 1200)
