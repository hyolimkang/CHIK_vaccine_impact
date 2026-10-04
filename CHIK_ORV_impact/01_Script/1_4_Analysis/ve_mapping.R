
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

make_prevacc <- function(
    posterior, observed, N, region, lhs_sample, sero_vec = NULL
) {

  pre_ui <- simulate_pre_ui_age(
    posterior = posterior,
    bra_foi_state_summ = bra_foi_state_summ,
    age_groups = age_groups,
    N = N,
    region = region,
    observed = observed,
    lhs_sample = lhs_sample,
    sero_vec = sero_vec
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



# 2) all regions (all at once) ------------------------------------------------------------
prevacc_sets <- imap(regions, function(reg_args, region_name) {
  make_prevacc(
    posterior = reg_args$posterior,
    observed  = reg_args$observed,
    N         = reg_args$N,
    region    = region_name,
    lhs_sample = reg_args$lhs_sample
  )
})

sim_results <- imap(regions, function(reg_args, reg_name) {
  map(set_names(ve_inf_set, paste0("VE", ve_inf_set*100)),
      ~ run_simulation_scenarios_ui(
        target_age_list   = target_age_list,
        observed          = reg_args$observed,
        N                 = reg_args$N,
        bra_foi_state_summ = bra_foi_state_summ,
        age_groups        = age_groups,
        region_name       = reg_name,
        hosp              = hosp,
        fatal             = fatal,
        nh_fatal          = nh_fatal,
        lhs_sample_young  = lhs_sample_young,
        lhs_old           = lhs_old,
        le_sample         = le_sample,
        age_gr_levels     = age_gr_levels,
        prevacc_ui        = reg_args$prevacc_ui,
        posterior         = reg_args$posterior,
        ve_inf            = .x        
      )
  )
})

save("sim_results", file = "00_Data/0_2_Processed/sim_results_ve.RData")
save("prevacc_sets", file = "00_Data/0_2_Processed/prevacc_sets_ve.RData")

postsim_all <- imap(regions, function(reg_args, region_name) {
  
  ve_list  <- sim_results[[region_name]]       # 3개 VE 결과
  prevacc  <- prevacc_sets[[region_name]]      # ← ★ 사전-백신 세트 가져오기
  
  imap(ve_list, function(scenario_result, ve_tag) {
    
    postsim_all_ui(
      scenario_result          = scenario_result,
      observed                 = reg_args$observed,
      age_gr_levels            = age_gr_levels,
      pre_summary_cases_age    = prevacc$pre_summary_age,
      pre_summary_cases        = prevacc$prevacc_ui,
      pre_summary_cases_all    = prevacc$prevacc_ui_all,
      region                   = region_name
    )
    
  })                 # imap → VE별 리스트 반환
})  

combined_prepost_case_allve <- imap_dfr(
  postsim_all,                        # 지역 loop
  ~ imap_dfr(                         # VE loop
    .x,                             # .x = ve_list
    ~ .x$summary_week_df |> 
      mutate(VE = .y)             # VE 태그 추가(VE0, VE50, VE98.9)
  )
)

## min, max 
combined_prepost_case_ve_2 <- combined_prepost_case_allve %>% filter(VE != "VE50") 

global_impact_ve <- combined_prepost_case_ve_2 %>%       # ← VE 컬럼 포함된 데이터
  group_by(VE, Scenario, region) %>%                 # ★ VE도 그룹에 추가
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    impact = diff / total_pre_cases * 100
  )

## 2) 주석용 데이터프레임 생성 ---------------------------------------
annotation_ve <- global_impact_ve %>%
  mutate(
    Strategy = dplyr:: recode(Scenario,
                      "Scenario_1" = "Strategy 1",
                      "Scenario_2" = "Strategy 2",
                      "Scenario_3" = "Strategy 3")
  ) %>%
  group_by(region, VE) %>%                          # ★ VE별로 따로 주석
  summarise(
    ann = paste0(Strategy, ": ", round(impact, 1), "%", collapse = "\n"),
    .groups = "drop"
  )

# 2) 관측치도 VE별로 복제
observed_all_full <- tidyr::crossing(
  observed_all,                            # Week, Observed, region
  VE = unique(combined_prepost_case_ve_2$VE)
)

ve_subnat<-
ggplot(combined_prepost_case_ve_2) +
  # ❶ 백신 투입 기간 회색 음영 (facet 전체 공통)
  geom_ribbon(data = combined_prepost_case_ve_2 %>%
                filter(between(Week, vacc_start_week_s1, vacc_end_week_s1)),
              aes(x = Week, ymin = 0, ymax = Inf),
              fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  
  # ❷ 포스트-백신 UI 리본 (시나리오별 색)
  geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = 0.25) +
  
  # ❸ 프리-백신 UI 리본 (연한 회색)
  geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = 0.5) +
  
  # ❹ 포스트-백신 선 (Scenario별 색)
  geom_line(aes(x = Week, y = post_cases, color = Scenario), size = 0.25) +
  
  # ❺ 프리-백신 선 (검정 점선)
  geom_line(aes(x = Week, y = pre_cases, group = Scenario),
            color = "black", linetype = "dashed", size = 0.25) +
  
  # ❻ 관측치 점
  geom_point(data = observed_all_full,
             aes(x = Week, y = Observed),
             size = 0.2, inherit.aes = FALSE) +
  
  # ❼ 지역-별 주석
  geom_text(data = annotation_ve,
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE) +
  
  scale_fill_brewer(name   = "Vaccination strategy",  
                    palette = "Set1",
                    labels = c("Scenario_1" = "1–19 y",
                               "Scenario_2" = "20–59 y",
                               "Scenario_3" = "≥60 y")) +
  scale_color_brewer(name   = "Vaccination strategy",  
                     palette = "Set1",
                     labels = c("Scenario_1" = "1–19 y",
                                "Scenario_2" = "20–59 y",
                                "Scenario_3" = "≥60 y")) +
  scale_y_continuous(labels = comma) +
  labs(x = "Week", y = "Predicted symptomatic cases",
       fill = "Scenario", color = "Scenario") +
  
  # ❽ facet: 행 = 지역, 열 = VE
  facet_grid(region ~ VE, scales = "free_y") +
  #facet_wrap(vars(region, VE), ncol = 6, scales = "free_y")+
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.margin  = margin(5, 20, 5, 5),
    legend.position = "bottom"
  )

## by setting
combined_by_setting <- combined_prepost_case_ve_2 %>% 
  mutate(
    setting = dplyr::recode(region, !!!setting_key)   # add new column
  )

setting_summ <- combined_by_setting %>% 
  group_by(setting, VE, Scenario, Week) %>% 
  summarise(
    post_cases        = sum(post_cases),
    pre_cases         = sum(pre_cases),
    diff              = sum(diff),
    post_fatal        = sum(post_fatal),
    pre_fatal         = sum(pre_fatal),
    diff_fatal        = sum(diff_fatal),
    post_daly         = sum(post_daly),
    pre_daly          = sum(pre_daly),
    diff_daly         = sum(diff_daly),
    post_weekly_median = sum(post_weekly_median),
    post_weekly_low95  = sum(post_weekly_low95),
    post_weekly_hi95   = sum(post_weekly_hi95),
    lo95               = sum(lo95),
    hi95               = sum(hi95),
    .groups = "drop"
  )

combined_nat_setting <- combined_by_setting %>% 
  group_by(VE, Scenario, Week) %>% 
  summarise(
    post_cases        = sum(post_cases),
    pre_cases         = sum(pre_cases),
    diff              = sum(diff),
    post_fatal        = sum(post_fatal),
    pre_fatal         = sum(pre_fatal),
    diff_fatal        = sum(diff_fatal),
    post_daly         = sum(post_daly),
    pre_daly          = sum(pre_daly),
    diff_daly         = sum(diff_daly),
    post_weekly_median = sum(post_weekly_median),
    post_weekly_low95  = sum(post_weekly_low95),
    post_weekly_hi95   = sum(post_weekly_hi95),
    lo95               = sum(lo95),
    hi95               = sum(hi95),
    .groups = "drop"
  )

combined_nat_setting$setting <- "National"

library(forcats)

combined_all <- bind_rows(combined_nat_setting, setting_summ) %>%
  mutate(
    setting = factor(setting,
                     levels = c("Low", "Moderate", "High", "National"))  # bottom → top
  )

## 2 ── one-row data frame for the grey vaccination window
vacc_window_df <- data.frame(
  xmin = vacc_start_week_s1,
  xmax = vacc_end_week_s1
)

global_impact_ve <- combined_all %>%       # ← VE 컬럼 포함된 데이터
  group_by(VE, Scenario, setting) %>%                 # ★ VE도 그룹에 추가
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    impact = diff / total_pre_cases * 100
  )

annotation_setting_ve <- global_impact_ve %>%      # already Low / Moderate / …
  mutate(
    Strategy = dplyr::recode(Scenario,
                      "Scenario_1" = "Strategy 1",
                      "Scenario_2" = "Strategy 2",
                      "Scenario_3" = "Strategy 3",
                      "Scenario_4" = "Strategy 4"),
    setting  = factor(setting,                     # make sure order matches plot
                      levels = c("Low","Moderate","High","National"))
  ) %>% 
  group_by(setting, VE) %>%                        # exactly the facet grid
  summarise(
    ann = paste0(Strategy, ": ", round(impact, 1), "%", collapse = "\n"),
    .groups = "drop"
  )

## 3 ── plot
ve_epi_graph_all <- 
ggplot(combined_all) +
  
  ## ❶ vaccination window (now appears in *all* facets)
  geom_rect(data = vacc_window_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = .20, inherit.aes = FALSE) +
  
  ## ❷ post-vaccine UI ribbon
  geom_ribbon(aes(Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = .25) +
  
  ## ❸ pre-vaccine UI ribbon
  geom_ribbon(aes(Week, ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = .50) +
  
  ## ❹ post-vaccine median
  geom_line(aes(Week, post_cases, colour = Scenario), size = .25) +
  
  ## ❺ pre-vaccine median
  geom_line(aes(Week, pre_cases, group = Scenario),
            colour = "black", linetype = "dashed", size = .25) +
  geom_text(data = annotation_setting_ve,
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE)+
  scale_fill_brewer(palette = "Set1",
                    name   = "Vaccination strategy",
                    labels = c("Scenario_1" = "1–11 y",
                               "Scenario_2" = "12-17 y",
                               "Scenario_3" = "18-59 y",
                               "Scenario_4" = "60+ y")) +
  scale_colour_brewer(palette = "Set1",
                      name   = "Vaccination strategy",
                      labels = c("Scenario_1" = "1–11 y",
                                 "Scenario_2" = "12-17 y",
                                 "Scenario_3" = "18-59 y",
                                 "Scenario_4" = "60+ y")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Predicted symptomatic cases") +
  
  facet_grid(fct_rev(setting) ~ VE, scales = "free_y") +   # National at top
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y     = element_text(angle = 0),
    plot.margin      = margin(5, 20, 5, 5),
    legend.position  = "bottom"
  )

ggsave(filename = "02_Outputs/2_1_Figures/ve_epi_graph_all.jpg", ve_epi_graph_all, width = 8, height = 7, dpi = 1500)

#-------------------------------------------------------------------------------
## 2) 글로벌 임팩트 및 주석 계산 ------------------
combined_nat <- combined_prepost_case_ve_2 %>%          # ← 이미 VE 포함
  group_by(VE, Scenario, Week) %>%                 # ★ VE 추가
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),  # post_cases, pre_cases, hi95 …
    .groups = "drop"
  )

global_impact_nat <- combined_nat %>%          # ← VE, Scenario, Week, post/pre_cases
  group_by(VE, Scenario) %>%                  # ★ VE도 포함
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    impact = diff / total_pre_cases * 100
  )

## 2) 주석 텍스트 만들기 ----------------------------------------------------
annotation_nat <- global_impact_nat %>%
  mutate(
    Strategy = dplyr::recode(Scenario,
                      "Scenario_1" = "Strategy 1",
                      "Scenario_2" = "Strategy 2",
                      "Scenario_3" = "Strategy 3",
                      "Scenario_4" = "Strategy 4")
  ) %>%
  group_by(VE) %>%                            # ★ VE별로 따로 주석
  summarise(
    ann = paste0(Strategy, ": ",
                 round(impact, 1), "%", collapse = "\n"),
    .groups = "drop"
  )

vacc_tbl <- imap_dfr(postsim_all, function(ve_list, reg) {
  imap_dfr(ve_list, function(x, ve) {
    tibble(region = reg, VE = ve,
           start = x$vacc_weeks$scenario1$start,   # (전략 1 기준)
           end   = x$vacc_weeks$scenario1$end)
  })
})%>% 
  filter(VE != "VE50")  

observed_nat <- observed_all %>% 
  group_by(Week) %>%                         # 주차별
  summarise(Observed = sum(Observed, na.rm = TRUE),
            .groups = "drop")

ve_levels <- unique(combined_nat$VE)         # "VE0" "VE50" "VE98.9"

observed_nat_full <- crossing(observed_nat, VE = ve_levels)

ve_nat<- 
ggplot(combined_nat, aes(x = Week)) +
  
  ## (a) 백신 투입 회색 음영 ─ 지역·VE 별로 join 후 geom_rect
  geom_rect(data = vacc_tbl,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "lightgrey", alpha = .1) +
  
  ## (b) 포스트-백신 UI 리본
  geom_ribbon(aes(ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = .18) +
  
  ## (c) 프리-백신 UI 리본
  geom_ribbon(aes(ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = .5) +
  
  ## (d) 포스트-백신 선
  geom_line(aes(y = post_cases, colour = Scenario, linetype = "With"), size = .35) +
  
  ## (e) 프리-백신 선
  geom_line(aes(y = pre_cases, colour = Scenario, linetype = "Without"), color = "black", size = .35) +
  
  ## (f) 관측치
  geom_point(data = observed_nat_full, aes(y = Observed), size = .3) +
  
  scale_colour_brewer(name   = "Vaccination strategy",  
                      palette = "Set1",
                      labels = c("Scenario_1" = "1–11 y",
                                 "Scenario_2" = "12-17 y",
                                 "Scenario_3" = "18-59 y",
                                 "Scenario_4" = "60+ y")) +
  scale_fill_brewer(name   = "Vaccination strategy",  palette = "Set1", guide = "none") +
  scale_linetype_manual(values = c(With = "solid", Without = "dashed"), name = "") +
  scale_y_continuous(labels = comma) +
  labs(x = "Week", y = "Predicted symptomatic cases", colour = "Strategy") +
  
  ## (g) facet : 행 = region, 열 = VE
  facet_wrap(~ VE, scales = "free_y") +
  
  ## (h) 각 패널 오른쪽 위 글로벌 주석
  geom_text(data = annotation_nat,   # region + VE 별 주석 복제본
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE) +
  
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.margin  = margin(5, 25, 5, 5),
    legend.position = "right"
  )

fig1_a <- 
ggarrange(ve_nat ,
          ncol = 1, nrow = 1,
          labels = c("A"),
          common.legend = TRUE,
          legend = "bottom",
          align = "none")

fig1_b <- 
  ggarrange(ve_subnat ,
            ncol = 1, nrow = 1,
            labels = c("B"),
            common.legend = TRUE,
            legend = "bottom",
            align = "none")

ggsave(filename = "02_Outputs/2_1_Figures/fig1_a.jpg", fig1_a, width = 16, height = 6, dpi = 1500)
ggsave(filename = "02_Outputs/2_1_Figures/fig1_b.jpg", fig1_b, width = 12, height = 10, dpi = 1500)

#-------------------------------------------------------------------------------
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

nnv_results <- imap(regions, function(reg_args, region_name) {
  
  ve_list <- postsim_all[[region_name]]
  
  imap(ve_list, function(postsim_ui, ve_tag) {
    
    nnv_list(
      vacc_alloc     = reg_args$vacc_alloc,
      postsim_all_ui = postsim_ui,
      N              = reg_args$N,
      region         = region_name,      
      observed       = reg_args$observed
    )
  })
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

## 1) age_gr 벡터 한 번 정의 -----------------------------------------------
n_scenarios = 4
age_seq <- rep(age_gr[1:length(age_gr_levels)], n_scenarios)          # 길이 18×3 = 54

## 2) nnv_results (지역 → VE → nnv_list 결과) 결합 --------------------------
combined_nnv_df_region <- imap_dfr(nnv_results,      # 지역 loop
                                   function(ve_list, region_name) {
                                     imap_dfr(ve_list,                                # VE loop
                                              function(nnv, ve_tag) {
                                                nnv$final_summ_df %>%                        # (18×3 행)
                                                  mutate(
                                                    region   = region_name,
                                                    VE       = ve_tag,
                                                    setting  = setting_key[region_name],
                                                    age_gr   = factor(age_seq, levels = age_gr_levels)
                                                  )
                                              })
                                   })

combined_nnv_national_ve <- combined_nnv_df_region %>%
  group_by(VE, scenario) %>%                         # ★ VE 추가
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

combined_nnv_setting_ve <- combined_nnv_df_region %>% 
  group_by(VE, setting, scenario, target) %>%           # ★ VE 추가
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

combined_nnv_nat_target_ve <- combined_nnv_df_region %>% 
  group_by(VE, scenario, target) %>%                    # ★ VE 추가
  summarise(
    across(c(tot_vacc, pre_vacc:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(VE, scenario) %>%
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
    nnv_daly_hi = scenario_vacc / diff_daly_low
  )

combined_nnv_setting_summ_ve <- combined_nnv_setting_ve %>% 
  group_by(VE, setting, scenario) %>%                   # ★ VE 추가
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

combined_nnv_national_ve$setting <- "National"
combined_nnv_setting_summ_ve <- rbind(combined_nnv_national_ve, combined_nnv_setting_summ_ve)

combined_nnv_setting_summ_ve2 <- combined_nnv_setting_summ_ve %>% 
  filter(VE != "VE50") 

# pub table ---------------------------------------------------------
library(glue)
metrics <- c("pre_vacc", "post_vacc", "pre_fatal",
             "pre_daly", "post_fatal", "post_daly",
             "diff", "diff_fatal", "diff_daly",
             "nnv", "nnv_fatal", "nnv_daly")

table_nnv_ci_ve <- combined_nnv_setting_summ_ve2 %>%
  rowwise() %>%
  mutate(
    diff_ci = glue(
      "{formatC(diff,        format='f', digits=0, big.mark=',')} ",
      "({formatC(diff_low,  format='f', digits=0, big.mark=',')}–",
      "{formatC(diff_hi,   format='f', digits=0, big.mark=',')})"
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
    scenario, setting, VE,
    ends_with("_ci")
  )


write_xlsx(table_nnv_ci_ve, path = "02_Outputs/2_2_Tables/table_nnv_all.xlsx")


## ── 1.  REWRITE THE NNV PLOT HELPER ────────────────────────────────
short_si <- function(x) {
  dplyr::case_when(
    abs(x) >= 1e6 ~ paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
    abs(x) >= 1e3 ~ paste0(formatC(x / 1e3, format = "f", digits = 1), "K"),
    TRUE          ~ as.character(x)
  )
}

make_nnv_overlay <- function(df, y_var, y_lab) {
  
  df <- df %>% 
    filter(VE %in% c("VE0", "VE98.9")) %>% 
    mutate(
      setting  = factor(setting,  levels = c("National","High","Moderate","Low")),
      scenario = factor(scenario,
                        levels = c("Scenario_1","Scenario_2","Scenario_3", "Scenario_4"),
                        labels = c("1–11 years","12-17 years","18-59 years", "60+ years")),
      VE_lab   = factor(VE, levels = c("VE0","VE98.9"),
                        labels = c("0 %", "98·9 %"))
    )
  
  y_sym  <- sym(y_var)
  lo_sym <- sym(paste0(y_var, "_lo"))
  hi_sym <- sym(paste0(y_var, "_hi"))
  
  ggplot(df, aes(x = scenario, y = !!y_sym, fill = VE_lab)) +
    
    ## rear bar (0 %)
    geom_col(data = subset(df, VE_lab == "0 %"),
             width = .8, alpha = .45) +
    geom_errorbar(data = subset(df, VE_lab == "0 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym),
                  width = .25, alpha = .45) +
    
    ## front bar (98·9 %)
    geom_col(data = subset(df, VE_lab == "98·9 %"),
             width = .5, alpha = .80) +
    geom_errorbar(data = subset(df, VE_lab == "98·9 %"),
                  aes(ymin = !!lo_sym, ymax = !!hi_sym),
                  width = .15, alpha = .80) +
    
    facet_grid(~ setting) +
    
    scale_fill_manual(values = c("0 %" = "#E41A1C", "98·9 %" = "#377EB8"),
                      name   = "Infection-blocking VE") +
    scale_y_continuous(labels = short_si,
                       expand  = expansion(mult = c(0, .05))) +
    
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.x = element_blank()
    )
}

#──────────────────────────────────────────────────────────────────────
#  Build the three panels
#──────────────────────────────────────────────────────────────────────
p_cases <- make_nnv_overlay(combined_nnv_setting_summ_ve2,
                            "nnv",       "NNV to avert a symptomatic case") +
  ggtitle(NULL) +                         # <- removes any existing title
  scale_y_continuous(labels = short_si) +
  labs(
    y        = "NNV to avert a symptomatic case",
    linetype = "Infection-blocking VE",
    colour   = "Infection-blocking VE",
    fill     = "Infection-blocking VE"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 10, colour = "black")
  )

p_death <- make_nnv_overlay(combined_nnv_setting_summ_ve2,
                            "nnv_fatal", "NNV to avert a death")+
  ggtitle(NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank()
  )


p_daly  <- make_nnv_overlay(combined_nnv_setting_summ_ve2,
                            "nnv_daly",  "NNV to avert a DALY")+
  ggtitle(NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    strip.text   = element_blank()
  )

#──────────────────────────────────────────────────────────────────────
#  Combine with common legend & vertical labels
#──────────────────────────────────────────────────────────────────────
combined_nnv_graph <- ggarrange(
  p_cases, p_death, p_daly,
  ncol = 1, labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "bottom",
  align = "v"
)

ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_graph_ve.jpg", combined_nnv_graph, width = 8, height = 9, dpi = 1200)


### heatmpap -------------------------------------------------------------------
ve_grid <- tibble(VE_inf = c(0, 0.989))

regions <- list(
  "Ceará" = list(
    observed   = observed_ce,
    posterior  = posterior_ce,
    N          = N_ceara$Ceará,
    prevacc_ui = prevacc_ui_ce
  ),
  "Bahia" = list(
    observed   = observed_bh,
    posterior  = posterior_bh,
    N          = N_bahia$Bahia,
    prevacc_ui = prevacc_ui_bh
  ),
  "Paraíba" = list(
    observed   = observed_pa,
    posterior  = posterior_pa,
    N          = N_pa$Paraíba,
    prevacc_ui = prevacc_ui_pa
  ),
  "Pernambuco" = list(
    observed   = observed_pn,
    posterior  = posterior_pn,
    N          = N_pemam$Pernambuco,
    prevacc_ui = prevacc_ui_pn
  ),
  "Rio Grande do Norte" = list(
    observed   = observed_rg,
    posterior  = posterior_rg,
    N          = N_rg$`Rio Grande do Norte`,
    prevacc_ui = prevacc_ui_rg
  ),
  "Piauí" = list(
    observed   = observed_pi,
    posterior  = posterior_pi,
    N          = N_pi$Piauí,
    prevacc_ui = prevacc_ui_pi
  ),
  "Tocantins" = list(
    observed   = observed_tc,
    posterior  = posterior_tc,
    N          = N_tc$Tocantins,
    prevacc_ui = prevacc_ui_tc
  ),
  "Alagoas" = list(
    observed   = observed_ag,
    posterior  = posterior_ag,
    N          = N_ag$Alagoas,
    prevacc_ui = prevacc_ui_ag
  ),
  "Minas Gerais" = list(
    observed   = observed_mg,
    posterior  = posterior_mg,
    N          = N_mg$`Minas Gerais`,
    prevacc_ui = prevacc_ui_mg
  ),
  "Sergipe" = list(
    observed   = observed_se,
    posterior  = posterior_se,
    N          = N_se$Sergipe,
    prevacc_ui = prevacc_ui_se
  ),
  "Goiás" = list(
    observed   = observed_go,
    posterior  = posterior_go,
    N          = N_go$Goiás,
    prevacc_ui = prevacc_ui_go
  )
)

regions_tbl <- tibble(
  region = names(regions),        # "Ceará", "Bahia", …
  inputs = regions                # each element is the 4-item list
)


heat_df_VE0 <- map_dfr(
  names(regions),                            # each region name
  \(reg) {
    inp <- regions[[reg]]
    heatmap_delay(
      target_age_list = target_age_list,
      observed        = inp$observed,
      posterior       = inp$posterior,
      N               = inp$N,
      bra_foi_state_summ = bra_foi_state_summ,
      age_groups      = age_groups,
      region_name     = reg,
      hosp            = hosp, fatal = fatal, nh_fatal = nh_fatal,
      lhs_sample_young = lhs_sample_young, lhs_old = lhs_old,
      le_sample       = le_sample,
      pre_summary_cases = inp$prevacc_ui,
      VE_inf          = 0
    ) %>% mutate(VE_inf = 0)
  }
)

heat_df_VE989 <- map_dfr(
  names(regions),
  \(reg) {
    inp <- regions[[reg]]
    heatmap_delay(
      target_age_list = target_age_list,
      observed        = inp$observed,
      posterior       = inp$posterior,
      N               = inp$N,
      bra_foi_state_summ = bra_foi_state_summ,
      age_groups      = age_groups,
      region_name     = reg,
      hosp            = hosp, fatal = fatal, nh_fatal = nh_fatal,
      lhs_sample_young = lhs_sample_young, lhs_old = lhs_old,
      le_sample       = le_sample,
      pre_summary_cases = inp$prevacc_ui,
      VE_inf          = 0.989
    ) %>% mutate(VE_inf = 0.989)
  }
)

# ----- combine & continue -------------------------------------------
heat_df <- bind_rows(heat_df_VE0, heat_df_VE989)
heat_df <- heat_df %>%
  mutate(
    scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
  )

heat_df <- heat_df %>% 
  left_join(vacc_wide, by = "region") %>% 
  mutate(
    # --- 기존: 총효과 (all age groups) ---
    impact_cases_agegr1 =
      (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) /
      totvacc_Scenario_1 * 1e5,
    impact_cases_agegr2 =
      (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) /
      totvacc_Scenario_2 * 1e5,
    impact_cases_agegr3 =
      (diff_cases_agegr1 + diff_cases_agegr2 + diff_cases_agegr3) /
      totvacc_Scenario_3 * 1e5,
    
    # --- Scenario 1 direct / indirect for cases ---
    direct_cases_1    = diff_cases_agegr1,
    indirect_cases_1  = diff_cases_agegr2 + diff_cases_agegr3,
    impact_direct_1   = direct_cases_1   / totvacc_Scenario_1 * 1e5,
    impact_indirect_1 = indirect_cases_1 / totvacc_Scenario_1 * 1e5,
    
    # --- Scenario 2 direct / indirect for cases ---
    direct_cases_2    = diff_cases_agegr2,
    indirect_cases_2  = diff_cases_agegr1 + diff_cases_agegr3,
    impact_direct_2   = direct_cases_2   / totvacc_Scenario_2 * 1e5,
    impact_indirect_2 = indirect_cases_2 / totvacc_Scenario_2 * 1e5,
    
    # --- Scenario 3 direct / indirect for cases ---
    direct_cases_3    = diff_cases_agegr3,
    indirect_cases_3  = diff_cases_agegr1 + diff_cases_agegr2,
    impact_direct_3   = direct_cases_3   / totvacc_Scenario_3 * 1e5,
    impact_indirect_3 = indirect_cases_3 / totvacc_Scenario_3 * 1e5,
    
    # --- 동일 로직을 fatal에 적용 ---
    impact_fatal_agegr1 =
      (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) /
      totvacc_Scenario_1 * 1e5,
    impact_fatal_agegr2 =
      (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) /
      totvacc_Scenario_2 * 1e5,
    impact_fatal_agegr3 =
      (diff_fatal_agegr1 + diff_fatal_agegr2 + diff_fatal_agegr3) /
      totvacc_Scenario_3 * 1e5,
    
    direct_fatal_1    = diff_fatal_agegr1,
    indirect_fatal_1  = diff_fatal_agegr2 + diff_fatal_agegr3,
    impact_fatal_direct_1   = direct_fatal_1   / totvacc_Scenario_1 * 1e5,
    impact_fatal_indirect_1 = indirect_fatal_1 / totvacc_Scenario_1 * 1e5,
    
    direct_fatal_2    = diff_fatal_agegr2,
    indirect_fatal_2  = diff_fatal_agegr1 + diff_fatal_agegr3,
    impact_fatal_direct_2   = direct_fatal_2   / totvacc_Scenario_2 * 1e5,
    impact_fatal_indirect_2 = indirect_fatal_2 / totvacc_Scenario_2 * 1e5,
    
    direct_fatal_3    = diff_fatal_agegr3,
    indirect_fatal_3  = diff_fatal_agegr1 + diff_fatal_agegr2,
    impact_fatal_direct_3   = direct_fatal_3   / totvacc_Scenario_3 * 1e5,
    impact_fatal_indirect_3 = indirect_fatal_3 / totvacc_Scenario_3 * 1e5,
    
    # --- 동일 로직을 DALY에 적용 ---
    impact_daly_agegr1 =
      (diff_daly_agegr1 + diff_daly_agegr2 + diff_daly_agegr3) /
      totvacc_Scenario_1 * 1e5,
    impact_daly_agegr2 =
      (diff_daly_agegr1 + diff_daly_agegr2 + diff_daly_agegr3) /
      totvacc_Scenario_2 * 1e5,
    impact_daly_agegr3 =
      (diff_daly_agegr1 + diff_daly_agegr2 + diff_daly_agegr3) /
      totvacc_Scenario_3 * 1e5,
    
    direct_daly_1    = diff_daly_agegr1,
    indirect_daly_1  = diff_daly_agegr2 + diff_daly_agegr3,
    impact_daly_direct_1   = direct_daly_1   / totvacc_Scenario_1 * 1e5,
    impact_daly_indirect_1 = indirect_daly_1 / totvacc_Scenario_1 * 1e5,
    
    direct_daly_2    = diff_daly_agegr2,
    indirect_daly_2  = diff_daly_agegr1 + diff_daly_agegr3,
    impact_daly_direct_2   = direct_daly_2   / totvacc_Scenario_2 * 1e5,
    impact_daly_indirect_2 = indirect_daly_2 / totvacc_Scenario_2 * 1e5,
    
    direct_daly_3    = diff_daly_agegr3,
    indirect_daly_3  = diff_daly_agegr1 + diff_daly_agegr2,
    impact_daly_direct_3   = direct_daly_3   / totvacc_Scenario_3 * 1e5,
    impact_daly_indirect_3 = indirect_daly_3 / totvacc_Scenario_3 * 1e5
  )
## 0. clean VE labels that will become column suffixes -------------
heat_df <- heat_df %>% 
  mutate(VE_lbl = if_else(VE_inf == 0, "0", "0_989"))   # no dot

## 1. choose the columns that uniquely identify a cell -------------
id_cols <- c("region", "Scenario", "Supply", "Delay", "scenario_num")

## 2. if duplicates exist, collapse them (here: mean) --------------
heat_long <- heat_df %>% 
  group_by(across(all_of(c(id_cols, "VE_lbl")))) %>% 
  summarise(across(starts_with("impact_"), mean), .groups = "drop")

## 3. pivot wider with *those* IDs only ----------------------------
heat_wide <- heat_long %>% 
  pivot_wider(
    id_cols     = all_of(id_cols),
    names_from  = VE_lbl,                      # "0" | "0_989"
    values_from = starts_with("impact_"),
    names_glue  = "{.value}_VE{VE_lbl}"
  )
delta_tbl <- map_dfc(
  str_subset(names(heat_wide), "_VE0$"),     # base columns
  function(col0) {
    col989  <- str_replace(col0, "_VE0$", "_VE0_989")
    delta   <- heat_wide[[col989]] - heat_wide[[col0]]
    setNames(tibble(delta), str_replace(col0, "_VE0$", "_delta"))
  }
)

heat_wide <- bind_cols(heat_wide, delta_tbl)

## 2 ── build the matching delta columns

collapsed_heatmap <- heat_wide %>% 
  group_by(Scenario, Supply, Delay, scenario_num) %>%    # or add `region` if you need it
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


plot_heatmap <- function(data,
                         metric_base,      # e.g. "impact_cases_agegr"
                         suffix,           # "VE0"  or  "delta"
                         legend_title) {
  
  ## ------------------------------------------------------------------
  ## 1. build the regex for pivot_longer()
  ##    column names look like:  impact_cases_agegr1_VE0
  ##                             impact_cases_agegr1_delta
  ## ------------------------------------------------------------------
  pat <- paste0("^", metric_base, "(\\d+)_", suffix, "$")   # two groups
  
  ## ------------------------------------------------------------------
  ## 2. nice axis / legend labeller
  ## ------------------------------------------------------------------
  short_lab <- function(x) {
    dplyr::case_when(
      abs(x) >= 1e6 ~ paste0(round(x / 1e6, 1), "M"),
      abs(x) >= 1e3 ~ formatC(x, format = "d", big.mark = ","),
      TRUE          ~ as.character(x)
    )
  }
  
  data %>% 
    select(Supply, Delay, Scenario, scenario_num, matches(pat)) %>% 
    pivot_longer(
      cols         = matches(pat),
      names_to     = "agegr",
      names_pattern= pat,          # extracts the (\d+) part
      values_to    = "impact"
    ) %>% 
    mutate(agegr = as.integer(agegr)) %>% 
    filter(agegr == scenario_num) %>%            # keep matching bar
    ggplot(aes(Supply, Delay, fill = impact)) +
    geom_tile() +
    scale_fill_distiller(
      palette   = "Spectral",
      direction = -1,
      name      = legend_title,
      labels    = short_lab
    ) +
    scale_x_continuous(
      breaks = seq(0, 0.10, 0.01),
      labels = function(x) sprintf("%.0f", x * 100)
    ) +
    scale_y_continuous(labels = short_lab) +
    facet_wrap(
      ~ Scenario,
      labeller = as_labeller(c(
        Scenario_1 = "1–19 years",
        Scenario_2 = "20–59 years",
        Scenario_3 = "≥60 years"
      ))
    ) +
    coord_cartesian(xlim = c(0, 0.10)) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(x = "Vaccine supply (% of population)")
}
p_cases <- 
  plot_heatmap(
    collapsed_heatmap,
    "impact_cases_agegr",
    "VE0",
    "Cases averted\n(per 100 000 doses)"
  )

p_fatal <-
  plot_heatmap(
    collapsed_heatmap,
    "impact_fatal_agegr",
    "VE0",
    "Deaths averted\n(per 100 000 doses)"
  )

p_daly <-
  plot_heatmap(
    collapsed_heatmap,
    "impact_daly_agegr",
    "VE0",
    "DALYs averted\n(per 100 000 doses)"
  )

p_case_max <-
  plot_heatmap(
    collapsed_heatmap,
    "impact_cases_agegr",
    "VE0_989",
    "Cases averted\n(per 100 000 doses)"
  )


p_cases <- plot_heatmap(
  collapsed_heatmap,
  metric_base = "impact_cases_agegr",
  suffix      = "VE0",            # VE = 0 %
  legend_title= "Cases averted\n(per 100 000 doses)"
) +
  labs(title = "VE = 0 %") +
  theme(legend.position = "none")   # drop its own legend

p_cases_max <- plot_heatmap(
  collapsed_heatmap,
  metric_base = "impact_cases_agegr",
  suffix      = "VE0_989",         # VE = 98·9 %
  legend_title= "Cases averted\n(per 100 000 doses)"
) +
  labs(title = "VE = 98·9 %") +
  theme(legend.position = "none")

p_deaths <- plot_heatmap(
  collapsed_heatmap,
  metric_base = "impact_fatal_agegr",
  suffix      = "VE0",            # VE = 0 %
  legend_title= "Deaths averted\n(per 100 000 doses)"
) +
  labs(title = "VE = 0 %") +
  theme(legend.position = "none")   # drop its own legend

p_deaths_max <- plot_heatmap(
  collapsed_heatmap,
  metric_base = "impact_fatal_agegr",
  suffix      = "VE0_989",         # VE = 98·9 %
  legend_title= "Deaths averted\n(per 100 000 doses)"
) +
  labs(title = "VE = 98·9 %") +
  theme(legend.position = "none")


plot_single <- function(data,
                        metric_base,
                        suffix,
                        legend_title,
                        palette = c("abs", "delta")) {
  
  palette <- match.arg(palette)
  
  pat <- paste0("^", metric_base, "(\\d+)_", suffix, "$")
  
  data %>%                                     # ── same wrangling as before ──
    select(Supply, Delay, Scenario, scenario_num, matches(pat)) %>% 
    tidyr::pivot_longer(
      cols          = matches(pat),
      names_to      = "agegr",
      names_pattern = pat,
      values_to     = "impact"
    ) %>% 
    mutate(agegr = as.integer(agegr)) %>% 
    filter(agegr == scenario_num) %>% 
    ggplot(aes(Supply, Delay, fill = impact)) +
    geom_tile() +
    facet_wrap(
      ~ Scenario,
      labeller = as_labeller(c(
        Scenario_1 = "1–19 years",
        Scenario_2 = "20–59 years",
        Scenario_3 = "≥60 years"
      ))
    ) +
    
    ## ─────────────────────  COLOUR SCALES  ─────────────────────
    {
      if (palette == "abs") {
        ## baseline: your original Spectral (re-versed)
        scale_fill_distiller(
          palette   = "Spectral",
          direction = -1,
          name      = legend_title,
          labels    = short_lab
        )
      } else {
        ## delta: choose a stronger diverging Brewer palette
        ##        (e.g. 'RdBu' or 'PuOr'; both have vivid ends)
        scale_fill_distiller(
          palette   = "RdBu",     # <-- change here
          direction =  1,         # keep red = positive (or set -1 to flip)
          name      = legend_title,
          labels    = short_lab
        )
      }
    } +
    
    scale_x_continuous(
      breaks = seq(0, .10, .02),
      labels = function(x) sprintf("%.0f", x * 100)
    ) +
    scale_y_continuous(labels = short_lab) +
    labs(x = "Vaccine supply (% of total population)", y = "Delay") +
    coord_cartesian(xlim = c(0, .10)) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid   = element_blank(),
      axis.title.y = element_blank(),
      strip.text   = element_text(size = 7)
    )
}


lims_cases <- range(
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_cases_agegr\\d+_VE0$")),
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_cases_agegr\\d+_VE0_989$")),
  na.rm = TRUE
)

## ──────────────────────────────────────────────────────────────
## 2.  VE = 0 %  — baseline tile
## ──────────────────────────────────────────────────────────────
p_abs_cases <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_cases_agegr",
  suffix       = "VE0",
  legend_title = "Cases averted\n(per 100 000 doses)",
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_cases, labels = short_lab,
                       name    = "Cases averted\n(per 100 000 doses)") +
  labs(title = "VE = 0 %") +
  title_style +
  theme(legend.position = "none")

## ──────────────────────────────────────────────────────────────
## 3.  VE = 98·9 %  +  Δ-legend
## ──────────────────────────────────────────────────────────────
p_max_cases <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_cases_agegr",
  suffix       = "VE0_989",
  legend_title = NULL,                # will add our own diverging legend
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_cases, labels = short_lab,
                       name    = "Cases averted\n(per 100 000 doses)") +
  labs(title = "VE = 98·9 %") +
  title_style + hide_y +
  theme(legend.position = "none") +
  
  ## second fill scale for the Δ legend
  ggnewscale::new_scale_fill() +
  scale_fill_distiller(
    palette   = "RdYlBu",
    direction =  1,
    limits    = lims_cases,
    name      = "Δ Cases\n(98·9 % – 0 %)",
    labels    = short_lab
  ) +
  geom_blank()                        # forces ggplot to register the scale

## ──────────────────────────────────────────────────────────────
## 4.  Assemble side-by-side with collected legends
## ──────────────────────────────────────────────────────────────
pair_cases <- (p_abs_cases | p_max_cases) +
  plot_layout(guides = "collect") &
  theme(
    legend.box      = "horizontal",   # legends stack vertically
    legend.position = "right",
    plot.margin     = unit(c(5, 5, 5, 5), "pt")
  )
# -----------------------------------------------------------------------
lims_death <- range(
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_fatal_agegr\\d+_VE0$")),
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_fatal_agegr\\d+_VE0_989$")),
  na.rm = TRUE
)
p_abs_death <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_fatal_agegr",
  suffix       = "VE0",
  legend_title = "Deaths averted\n(per 100 000 doses)",
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_death, labels = short_lab,
                       name    = "Deaths averted\n(per 100 000 doses)") +
  labs(title = "VE = 0 %") +
  title_style +
  theme(legend.position = "none")

## ------------------------------------------------------------------
## 2. VE = 98·9 %  + Δ-scale
## ------------------------------------------------------------------
p_max_death <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_fatal_agegr",
  suffix       = "VE0_989",            # baseline values for the tile
  legend_title = NULL,                 # will be replaced just below
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_death, labels = short_lab,
                       name    = "Deaths averted\n(per 100 000 doses)") +
  labs(title = "VE = 98·9 %") +
  title_style + hide_y +
  theme(legend.position = "none") +
  
  ggnewscale::new_scale_fill() +         # start a *second* fill scale
  scale_fill_distiller(
    palette   = "RdYlBu",
    direction =  1,
    limits    = lims_death,
    name      = "Δ Deaths\n(98·9 % – 0 %)",
    labels    = short_lab
  ) +
  geom_blank()

## ------------------------------------------------------------------
## 3. side-by-side, let patchwork collect both legends
## ------------------------------------------------------------------
pair_death <- (p_abs_death | p_max_death) +
  plot_layout(guides = "collect") &
  theme(legend.box     = "horizontal",
        legend.position = "right",
        plot.margin     = unit(c(5, 5, 5, 5), "pt"))

lims_daly <- range(
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_daly_agegr\\d+_VE0$")),
  dplyr::select(collapsed_heatmap,
                tidyselect::matches("impact_daly_agegr\\d+_VE0_989$")),
  na.rm = TRUE
)

## ──────────────────────────────────────────────────────────────
## 2.  VE = 0 %  ‒ baseline DALYs
## ──────────────────────────────────────────────────────────────
p_abs_daly <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_daly_agegr",
  suffix       = "VE0",
  legend_title = "DALYs averted\n(per 100 000 doses)",
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_daly, labels = short_lab,
                       name    = "DALYs averted\n(per 100 000 doses)") +
  labs(title = "VE = 0 %") +
  title_style +
  theme(legend.position = "none")

## ──────────────────────────────────────────────────────────────
## 3.  VE = 98·9 %  +  Δ-legend
## ──────────────────────────────────────────────────────────────
p_max_daly <- plot_single(
  collapsed_heatmap,
  metric_base  = "impact_daly_agegr",
  suffix       = "VE0_989",
  legend_title = NULL,          # dummy, we'll add our own below
  palette      = "abs"
) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       limits  = lims_daly, labels = short_lab,
                       name    = "DALYs averted\n(per 100 000 doses)") +
  labs(title = "VE = 98·9 %") +
  title_style + hide_y +
  theme(legend.position = "none") +
  
  ggnewscale::new_scale_fill() +          # start 2nd fill scale
  scale_fill_distiller(
    palette   = "RdYlBu",
    direction =  1,
    limits    = lims_daly,
    name      = "Δ DALYs\n(98·9 % – 0 %)",
    labels    = short_lab
  ) +
  geom_blank()

## ──────────────────────────────────────────────────────────────
## 4.  side-by-side pair with 2 legends
## ──────────────────────────────────────────────────────────────
pair_daly <- (p_abs_daly | p_max_daly) +
  plot_layout(guides = "collect") &
  theme(
    legend.box     = "horizontal",
    legend.position = "right",
    plot.margin     = unit(c(5, 5, 5, 5), "pt")
  )




# ─────────────────────────────────────────────────────────────
# extra helper: remove facet-strip label + background
# ─────────────────────────────────────────────────────────────
drop_strip <- function(p) {
  p + theme(strip.text       = element_blank(),
            strip.background = element_blank())
}

# ─────────────────────────────────────────────────────────────
#  rebuild rows
# ─────────────────────────────────────────────────────────────

## --- A • CASES  (keep strips) --------------------------------
pair_cases_row <- {
  p_abs_cases2 <- drop_x_axis(p_abs_cases)          # keep strips
  p_max_cases2 <- drop_x_axis(p_max_cases)
  
  (p_abs_cases2 | p_max_cases2) +
    plot_layout(guides = "collect") &
    theme(legend.box = "horizontal",
          legend.position = "right",
          plot.margin = unit(c(5,5,5,5), "pt"))
}

## --- B • DEATHS  (no strips, no x-axis) -----------------------
pair_death_row <- {
  p_abs_death2 <- p_abs_death  |> drop_title() |> drop_x_axis() |> drop_strip()
  p_max_death2 <- p_max_death  |> drop_title() |> drop_x_axis() |> drop_strip()
  
  (p_abs_death2 | p_max_death2) +
    plot_layout(guides = "collect") &
    theme(legend.box = "horizontal",
          legend.position = "right",
          plot.margin = unit(c(5,5,5,5), "pt"))
}

## --- C • DALYs  (no strips, keep x-axis) ----------------------
pair_daly_row <- {
  p_abs_daly2 <- p_abs_daly  |> drop_title() |> drop_strip()        # keep x-axis
  p_max_daly2 <- p_max_daly  |> drop_title() |> drop_strip()
  
  (p_abs_daly2 | p_max_daly2) +
    plot_layout(guides = "collect") &
    theme(legend.box = "horizontal",
          legend.position = "right",
          plot.margin = unit(c(5,5,5,5), "pt"))
}

# ─────────────────────────────────────────────────────────────
# final stacked figure
# ─────────────────────────────────────────────────────────────
final_heatmap <- ggpubr::ggarrange(
  pair_cases_row,
  pair_death_row,
  pair_daly_row,
  ncol   = 1, nrow = 3,
  labels = c("A", "B", "C"),
  align  = "v",
  common.legend = FALSE,
  legend = "right",
  heights  = c(1, 0.8, 1)
)


ggsave(filename = "02_Outputs/2_1_Figures/pair_cases_heatmap.jpg", pair_cases, width = 8, height = 2, dpi = 1200)

ggsave(filename = "02_Outputs/2_1_Figures/final_heatmap.jpg", final_heatmap, width = 9, height = 6, dpi = 1200)

write_xlsx(heat_wide, path = "02_Outputs/2_2_Tables/heat_wide.xlsx")
write_xlsx(collapsed_heatmap, path = "02_Outputs/2_2_Tables/collapsed_heatmap_3.xlsx")




