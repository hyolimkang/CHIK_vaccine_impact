load("00_Data/0_2_Processed/combined_heatmap_samedose.RData")
target_age_list <- list(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 1-11 years
                        c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 12-17 years
                        c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0), # 18-59
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1) # 60+
                        
) 
simulate_pre_ui_age_med <- function(posterior, bra_foi_state_summ, age_groups, N, region, observed,
                                    A = 20,
                                    delay = 53, VE_block = 0, VE_inf = 0,
                                    target_age = rep(0, A),
                                    coverage_threshold = 0, total_coverage = 0,
                                    weekly_delivery_speed = 0) {
  # 1) 시간 길이 설정
  T <- nrow(observed)
  
  # 2) posterior에서 median 파라미터 추출
  base_beta_med <- apply(posterior$base_beta, 2, median)
  I0_med        <- apply(posterior$I0,      2, median)
  gamma_med     <- median(posterior$gamma)
  rho_med       <- median(posterior$rho)
  sigma_med     <- median(posterior$sigma)
  
  # 3) 단일 median 기반 시뮬레이션
  sim_out <- sirv_sim_coverageSwitch(
    T                    = T,
    A                    = A,
    N                    = N,
    r                    = rep(0, A),
    base_beta            = base_beta_med,
    I0_draw              = I0_med,
    R0                   = 1 - exp(- bra_foi_state_summ$avg_foi[
      bra_foi_state_summ$NAME_1 == region
    ] * age_groups),
    rho                  = rho_med,
    gamma                = gamma_med,
    sigma                = sigma_med,
    delay                = delay,
    VE_block             = VE_block,
    VE_inf               = VE_inf,
    target_age           = target_age,
    coverage_threshold   = coverage_threshold,
    total_coverage       = total_coverage,
    weekly_delivery_speed= weekly_delivery_speed
  )
  
  # 4) age-stratified symptomatic median (A×T)
  median_by_age_rawsymp <- sim_out$true_symptomatic
  # no uncertainty at median-run
  low95_by_age_rawsymp  <- median_by_age_rawsymp
  hi95_by_age_rawsymp   <- median_by_age_rawsymp
  
  # 5) 주별 합계
  weekly_cases_rawsymp <- colSums(median_by_age_rawsymp)
  
  # 6) weekly UI = median only
  weekly_low95_rawsymp <- weekly_cases_rawsymp
  weekly_hi95_rawsymp  <- weekly_cases_rawsymp
  
  # 7) Data frames
  df_rawsymp <- data.frame(
    week   = 1:T,
    median = weekly_cases_rawsymp,
    low95  = weekly_low95_rawsymp,
    hi95   = weekly_hi95_rawsymp
  )
  
  # 반환
  list(
    sim_out                     = sim_out,
    median_by_age_rawsymp       = median_by_age_rawsymp,
    low95_by_age_rawsymp        = low95_by_age_rawsymp,
    hi95_by_age_rawsymp         = hi95_by_age_rawsymp,
    weekly_cases_median_rawsymp = df_rawsymp$median,
    weekly_cases_low95_rawsymp  = df_rawsymp$low95,
    weekly_cases_hi95_rawsymp   = df_rawsymp$hi95,
    df_rawsymp                  = df_rawsymp
  )
}

heatmap_delay <- function(target_age_list,
                          observed,
                          posterior, 
                          N, 
                          bra_foi_state_summ, 
                          age_groups, 
                          region_name,
                          hosp, 
                          fatal, 
                          nh_fatal,
                          lhs_sample_young, 
                          lhs_old, 
                          le_sample, 
                          pre_summary_cases,
                          VE_inf = 0 ) {
  
  # no vacc scenario (수정 없음)
  pre_ui <- simulate_pre_ui_age_med(
    posterior, bra_foi_state_summ, age_groups, N, region_name, observed,
    delay       = Inf,     
    VE_block    = 0,   
    VE_inf      = 0,
    target_age = rep(0, length(N)),
    coverage_threshold = 0,
    total_coverage      = 0,
    weekly_delivery_speed = 0
  )
  
  pre_summary_cases <- summarise_presim_ui(
    sim_result       = pre_ui,
    observed         = observed,
    age_gr_levels    =  c(
      "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
      "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
      "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
    ),
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = hosp,
    fatal            = fatal,
    nh_fatal         = nh_fatal,
    region           = region_name
  )$summary_cases_pre 
  
  # Set T dynamically (수정 없음)
  T = nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <-  c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  
  age_gr_levels <- default_age_vector # 변수 이름 일관성을 위해 추가
  
  # Create an age group vector for the age-by-week summary:
  # Each age group label is repeated for every week.
  age_gr <- rep(default_age_vector, T)
  
  # supply and delay (수정 없음)
  pct <- region_coverage[region_name]
  supply <- seq(0, 1, by = 0.1)
  delay_steps  <- seq(1, T, by = 1)
  
  # extract median params (수정 없음)
  base_beta        <- apply(posterior$base_beta, 2, median)    
  I0               <- apply(posterior$I0, 2, median)           
  gamma            <- median(posterior$gamma)
  rho              <- median(posterior$rho)
  sigma            <- median(posterior$sigma)
  
  # ########################################################################## #
  # ### 아래부터 새로운 로직 적용 ###
  # ########################################################################## #
  
  # 결과를 저장할 비어있는 데이터프레임 초기화
  impact_summ <- data.frame()
  
  # Loop over scenarios, supply rates, and delay steps
  for (scenario_index in seq_along(target_age_list)) {
    
    target <- target_age_list[[scenario_index]]
    
    # ★ 현재 시나리오에 맞춰 직/간접 효과 그룹을 동적으로 정의
    direct_groups   <- which(target == 1)
    indirect_groups <- which(target == 0)
    
    for (supply_rate in supply) {
      for (delay_step in delay_steps) {
        
        # 시뮬레이션 실행 (기존과 동일)
        sim_result <- sirv_sim_coverageSwitch(
          T = T, 
          A = length(age_gr_levels),
          N = N,
          r = rep(0, 20),
          base_beta = base_beta,    
          I0 = I0,               
          R0  = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
          rho = rho,                
          gamma = gamma,        
          sigma = sigma,
          delay = delay_step,
          VE_block = 0.978,
          VE_inf = 0.978,
          coverage_threshold = 1,
          target_age = target,
          total_coverage = supply_rate,
          weekly_delivery_speed = 0.1
        )
        
        # 시뮬레이션 결과를 DALY까지 계산된 데이터프레임으로 변환 (기존 로직 유지)
        sim_df <- as.data.frame.table(sim_result$true_symptomatic, responseName = "Cases") %>%
          mutate(
            Scenario = scenario_index,
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          ) %>%
          mutate(
            hosp_rate = rep(hosp, T),
            hospitalised = Cases * hosp_rate,
            non_hospitalised = Cases - hospitalised,
            fatality = rep(fatal, T),
            nh_fatality = rep(nh_fatal, T),
            fatal = (hospitalised * fatality + non_hospitalised * nh_fatality)
          ) %>%
          mutate(
            age_numeric = case_when(
              age_gr == "<1" ~ 0, age_gr == "1-4" ~ 2.5, age_gr == "5–9" ~ 7,
              age_gr == "10-11" ~ 10.5, age_gr == "12-17" ~ 14.5, age_gr == "18–19" ~ 18.5,
              age_gr == "20–24" ~ 22, age_gr == "25–29" ~ 27, age_gr == "30–34" ~ 32,
              age_gr == "35–39" ~ 37, age_gr == "40–44" ~ 42, age_gr == "45–49" ~ 47,
              age_gr == "50–54" ~ 52, age_gr == "55–59" ~ 57, age_gr == "60–64" ~ 62,
              age_gr == "65–69" ~ 67, age_gr == "70–74" ~ 72, age_gr == "75–79" ~ 77,
              age_gr == "80–84" ~ 82, age_gr == "85+" ~ 87.5,
              TRUE ~ NA_real_
            )
          ) %>%
          mutate(
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
            dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
            dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
          ) %>%
          mutate(
            subac_prop = case_when(
              age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),
              age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)
            ),
            chr_6m = case_when(
              age_numeric < 40 ~ quantile(lhs_sample_young$chr6m, 0.5),
              age_numeric >= 40 ~ quantile(lhs_old$chr6m, 0.5)
            ),
            chr_12m = case_when(
              age_numeric < 40 ~ quantile(lhs_sample_young$chr12m, 0.5),
              age_numeric >= 40 ~ quantile(lhs_old$chr12m, 0.5)
            ),
            chr_30m = case_when(
              age_numeric < 40 ~ quantile(lhs_sample_young$chr30m, 0.5),
              age_numeric >= 40 ~ quantile(lhs_old$chr30m, 0.5)
            ),
            chr_prop = chr_6m + chr_12m + chr_30m
          ) %>%
          mutate(
            yld_acute = (hospitalised * dw_hosp * dur_acute) + (non_hospitalised * dw_nonhosp * dur_acute),
            yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
            yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) + (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
            yld_total = yld_acute + yld_subacute + yld_chronic
          ) %>%
          mutate(
            le_left = case_when(
              age_numeric <= 1 ~ quantile(le_sample$le_1, 0.5),
              age_numeric > 1 & age_numeric < 10 ~ quantile(le_sample$le_1, 0.5),
              age_numeric >= 10 & age_numeric < 20 ~ quantile(le_sample$le_2, 0.5),
              age_numeric >= 20 & age_numeric < 30 ~ quantile(le_sample$le_2, 0.5),
              age_numeric >= 30 & age_numeric < 40 ~ quantile(le_sample$le_3, 0.5),
              age_numeric >= 40 & age_numeric < 50 ~ quantile(le_sample$le_4, 0.5),
              age_numeric >= 50 & age_numeric < 60 ~ quantile(le_sample$le_5, 0.5),
              age_numeric >= 60 & age_numeric < 70 ~ quantile(le_sample$le_6, 0.5),
              age_numeric >= 70 & age_numeric < 80 ~ quantile(le_sample$le_7, 0.5),
              age_numeric >= 80 ~ quantile(le_sample$le_8, 0.5),
              TRUE ~ NA_real_
            )
          ) %>%
          mutate(
            yll = fatal * le_left,
            daly_tot = yld_total + yll,
            cum_daly = cumsum(daly_tot)
          )
        
        # ★★★ 직접/간접 효과 분리 집계 (새로운 로직) ★★★
        
        # 1. 접종 후(Post-vaccination) 결과 집계
        post_summary <- sim_df %>%
          group_by(AgeGroup) %>%
          summarise(across(c(Cases, fatal, daly_tot), sum, na.rm = TRUE), .groups = "drop") %>%
          summarise(
            post_direct_cases = sum(Cases[AgeGroup %in% direct_groups], na.rm = TRUE),
            post_indirect_cases = sum(Cases[AgeGroup %in% indirect_groups], na.rm = TRUE),
            post_direct_fatal = sum(fatal[AgeGroup %in% direct_groups], na.rm = TRUE),
            post_indirect_fatal = sum(fatal[AgeGroup %in% indirect_groups], na.rm = TRUE),
            post_direct_daly = sum(daly_tot[AgeGroup %in% direct_groups], na.rm = TRUE),
            post_indirect_daly = sum(daly_tot[AgeGroup %in% indirect_groups], na.rm = TRUE)
          )
        
        # 2. 접종 전(Pre-vaccination, Baseline) 결과 집계
        pre_summary <- pre_summary_cases %>%
          group_by(AgeGroup) %>%
          summarise(across(c(Median, fatal, daly_tot), sum, na.rm = TRUE), .groups = "drop") %>%
          summarise(
            pre_direct_cases = sum(Median[AgeGroup %in% direct_groups], na.rm = TRUE),
            pre_indirect_cases = sum(Median[AgeGroup %in% indirect_groups], na.rm = TRUE),
            pre_direct_fatal = sum(fatal[AgeGroup %in% direct_groups], na.rm = TRUE),
            pre_indirect_fatal = sum(fatal[AgeGroup %in% indirect_groups], na.rm = TRUE),
            pre_direct_daly = sum(daly_tot[AgeGroup %in% direct_groups], na.rm = TRUE),
            pre_indirect_daly = sum(daly_tot[AgeGroup %in% indirect_groups], na.rm = TRUE)
          )
        
        # 3. 효과(Averted) 계산 및 결과 행 생성
        total_doses <- sum(sim_result$raw_allocation_age, na.rm = TRUE)
        
        summary_row <- data.frame(
          Scenario = scenario_index,
          Supply = supply_rate,
          Delay = delay_step,
          Region = region_name,
          TotalDoses = total_doses,
          
          DirectCasesAverted = pre_summary$pre_direct_cases - post_summary$post_direct_cases,
          IndirectCasesAverted = pre_summary$pre_indirect_cases - post_summary$post_indirect_cases,
          
          DirectFatalAverted = pre_summary$pre_direct_fatal - post_summary$post_direct_fatal,
          IndirectFatalAverted = pre_summary$pre_indirect_fatal - post_summary$post_indirect_fatal,
          
          DirectDalyAverted = pre_summary$pre_direct_daly - post_summary$post_direct_daly,
          IndirectDalyAverted = pre_summary$pre_indirect_daly - post_summary$post_indirect_daly,
          
          Pre_direct_cases = pre_summary$pre_direct_cases,
          Pre_indirect_cases = pre_summary$pre_indirect_cases,
          Pre_cases  = pre_summary$pre_direct_cases + pre_summary$pre_indirect_cases,
          Pre_direct_deaths = pre_summary$pre_direct_fatal,
          Pre_indirect_deaths = pre_summary$pre_indirect_fatal,
          Pre_deaths = pre_summary$pre_direct_fatal + pre_summary$pre_indirect_fatal,
          Pre_direct_dalys = pre_summary$pre_direct_daly,
          Pre_indirect_dalys = pre_summary$pre_indirect_daly,
          Pre_dalys  = pre_summary$pre_direct_daly + pre_summary$pre_indirect_daly
        )
        
        # 4. 최종 결과 데이터프레임에 추가
        impact_summ <- rbind(impact_summ, summary_row)
      }
    }
  }
  
  return(impact_summ)
}

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

heatmap_summ <- combined_heatmap_ixchiq %>%
  mutate(
    scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
  )

# 'heatmap_summ'은 이전에 수정한 heatmap_delay 함수의 결과물이라고 가정합니다.

heatmap_summ_final <- heatmap_summ %>%
  mutate(
    
    # --- 1. 총 예방 효과 계산 (10만 도즈 당) ---
    # 분모(TotalDoses)가 0이 되어 발생하는 오류를 방지하기 위해 ifelse 사용
    
    CasesAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectCasesAverted + IndirectCasesAverted) / TotalDoses * 1e5,
      0
    ),
    
    FatalAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectFatalAverted + IndirectFatalAverted) / TotalDoses * 1e5,
      0
    ),
    
    DalyAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectDalyAverted + IndirectDalyAverted) / TotalDoses * 1e5,
      0
    ),

    # --- 2. (선택사항) 직접/간접 효과의 '기여도(비중)' 계산 ---
    
    # 먼저 전체 예방 케이스 수를 합산
    TotalCasesAverted = DirectCasesAverted + IndirectCasesAverted,
    
    # 전체 예방 케이스 중 직접/간접 효과의 비율 계산
    DirectContribution_Cases = ifelse(
      TotalCasesAverted > 0,
      DirectCasesAverted / TotalCasesAverted,
      0
    ),
    IndirectContribution_Cases = ifelse(
      TotalCasesAverted > 0,
      IndirectCasesAverted / TotalCasesAverted,
      0
    ),
    
    Case_averted_prop = ifelse(
      TotalDoses > 0,
      TotalCasesAverted / Pre_cases,
      0
    )
  )

collapsed_heatmap_ixchiq <- heatmap_summ_final %>%
  group_by(Scenario, Supply, Delay) %>%
  summarise(
    across(
      .cols = where(is.numeric),  
      .fns = mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

collapsed_heatmap_ixchiq <- collapsed_heatmap_ixchiq %>% filter(TotalDoses!=0)

plot_heatmap <- function(data, metric_to_plot, legend_title) {
  
  custom_labeller <- function(x) {
    sapply(x, function(v) {
      if (abs(v) >= 1e6) {
        paste0(round(v / 1e6, 1), "M")
      } else if (abs(v) >= 1e3) {
        formatC(v, format = "d", big.mark = ",")
      } else {
        as.character(v)
      }
    })
  }
    data_to_plot <- data %>%
    mutate(
      ScenarioLabel = factor(
        Scenario,
        levels = 1:4, 
        labels = c(
          "Target: 1-11 years", 
          "Target: 12-17 years", 
          "Target: 18-59 years", 
          "Target: 60+ years"
        )
      )
    )
  
  ggplot(data_to_plot, aes(x = Supply, y = Delay, fill = .data[[metric_to_plot]])) +
    geom_tile() +
    scale_fill_distiller(
      palette   = "Spectral",
      direction = -1,
      name      = legend_title,
      labels    = function(x) sprintf("%.0f", x * 100)
    ) + 
    scale_x_continuous(
      breaks = seq(0.1, 1, by = 0.1),  
      labels = function(x) sprintf("%.0f", x * 100) # 가독성을 위해 % 기호 추가
    ) +
    scale_y_continuous(labels = custom_labeller) +
    coord_cartesian(xlim = c(0.1, 1)) +
    facet_wrap(~ ScenarioLabel) + # ★ 새로 만든 라벨을 사용
    theme_minimal() +
    theme(
      axis.text.x = element_text(),
      panel.grid  = element_blank()
    ) +
    labs(x = "Vaccine supply (% of total population)", y = "Start of vaccination campaign (Week)")
}

plot_heatmap(
  data = collapsed_heatmap_ixchiq,
  metric_to_plot = "FatalAverted_per_100k_doses",
  legend_title = "FatalAverted_per_100k_doses"
)

plot_heatmap(
  data = collapsed_heatmap_ixchiq,
  metric_to_plot = "DalyAverted_per_100k_doses",
  legend_title = "DalyAverted_per_100k_doses"
)

p2<-
plot_heatmap(
  data = collapsed_heatmap_ixchiq,
  metric_to_plot = "IndirectContribution_Cases",
  legend_title = "% of Symptomatic cases averted \n among unvaccinated"
)

p1<-
plot_heatmap(
  data = collapsed_heatmap_ixchiq,
  metric_to_plot = "Case_averted_prop",
  legend_title = "% of Total symptomatic  \n cases averted"
)

combined_heatmap <- 
p1 / p2 +
  plot_layout(guides = "collect") +       # 범례 한 번만
  labs(
    x = "Vaccine supply (% of total population)",
    y = "Start of vaccination (Week)"
  )  +
  plot_annotation(
    tag_levels = "A",    # 첫 번째 패치웍 = "A", 두 번째 = "B"
    tag_suffix = ""      # "A." 대신 "A" 로 붙이려면 비워두기
  )

ggsave(filename = "02_Outputs/2_1_Figures/combined_heatmap_coverage.jpg", combined_heatmap, width = 7, height = 8, dpi = 1200)


heatmap_filter <- collapsed_heatmap
heatmap_filter <- heatmap_filter %>% filter(Delay == 2)
heatmap_filter <- heatmap_filter %>% filter(scenario_num == 2)



write_xlsx(collapsed_heatmap_ixchiq, path = "02_Outputs/2_2_Tables/collapsed_heatmap_ixchiq.xlsx")
write_xlsx(collapsed_heatmap_vimkun, path = "02_Outputs/2_2_Tables/collapsed_heatmap_vimkun.xlsx")

# -------------------------------------------------------------------------------
load("00_Data/0_2_Processed/combined_heatmap_coverage_ixchiq.RData")
combined_heatmap_ixchiq <- combined_heatmap
load("00_Data/0_2_Processed/combined_heatmap_coverage_vimkun.RData")
combined_heatmap_vimkun <- combined_heatmap

heatmap_summ_vimkun <- combined_heatmap_vimkun %>%
  mutate(
    scenario_num = as.integer(sub("Scenario_(\\d+)", "\\1", Scenario))
  )

heatmap_summ_vimkun <- heatmap_summ_vimkun %>%
  mutate(
    
    # --- 1. 총 예방 효과 계산 (10만 도즈 당) ---
    # 분모(TotalDoses)가 0이 되어 발생하는 오류를 방지하기 위해 ifelse 사용
    
    CasesAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectCasesAverted + IndirectCasesAverted) / TotalDoses * 1e5,
      0
    ),
    
    FatalAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectFatalAverted + IndirectFatalAverted) / TotalDoses * 1e5,
      0
    ),
    
    DalyAverted_per_100k_doses = ifelse(
      TotalDoses > 0,
      (DirectDalyAverted + IndirectDalyAverted) / TotalDoses * 1e5,
      0
    ),
    
    # --- 2. (선택사항) 직접/간접 효과의 '기여도(비중)' 계산 ---
    
    # 먼저 전체 예방 케이스 수를 합산
    TotalCasesAverted = DirectCasesAverted + IndirectCasesAverted,
    
    # 전체 예방 케이스 중 직접/간접 효과의 비율 계산
    DirectContribution_Cases = ifelse(
      TotalCasesAverted > 0,
      DirectCasesAverted / TotalCasesAverted,
      0
    ),
    IndirectContribution_Cases = ifelse(
      TotalCasesAverted > 0,
      IndirectCasesAverted / TotalCasesAverted,
      0
    ),
    
    Case_averted_prop = ifelse(
      TotalDoses > 0,
      TotalCasesAverted / Pre_cases,
      0
    )
  )

collapsed_heatmap_vimkun <- heatmap_summ_vimkun %>%
  group_by(Scenario, Supply, Delay) %>%
  summarise(
    across(
      .cols = where(is.numeric),  
      .fns = mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

collapsed_heatmap_vimkun <- collapsed_heatmap_vimkun %>% filter(TotalDoses!=0)

library(dplyr)
library(ggplot2)
library(scales)
plot_df <- collapsed_heatmap_vimkun %>% 
  ## ── 축 순서 고정용 factor 생성 ──────────────────
  mutate(
    Supply_f = factor(Supply, levels = seq(0.1, 1, 0.1)),     # x축
    Delay_f  = factor(Delay,  levels = seq(52, 1, -1)),       # y축 (위→아래 52→1)
    
    ## 보기 좋은 시나리오 라벨
    Scenario_f = factor(
      Scenario,
      levels = c(1, 2, 3, 4),
      labels = c("1–11 yrs", "12–17 yrs", "18–59 yrs", "60 + yrs")
    )
  )%>% 
  mutate(
    ## ① 퍼센트 값(0 – 100) 따로 계산
    Case_averted_pct = Case_averted_prop * 100,
    ## ② 라벨 문자열: 정수/소수 한 자리 등 원하는 형식으로
    pct_label        = sprintf("%.0f", Case_averted_pct),
    Indirect_ct = IndirectContribution_Cases * 100
  )


delay_thresh <- 18      # 예: 40주 초과는 너무 늦음
# --------------------------------------------------
plot_df2 <- collapsed_heatmap_vimkun %>% 
  filter(Delay <= delay_thresh) %>%                # ★ 1) Delay 기준 필터
  mutate(
    ## 축‧라벨용 factor 지정 -------------------------
    Supply_f = factor(Supply, levels = seq(0.1, 1, 0.1)),
    # Delay_f  = factor(Delay,  levels = seq(delay_thresh, 1, -1)), # ★ 2) 아래→위 1→20
    Delay_f  = factor(Delay,  levels = 1:delay_thresh), 
    Scenario_f = factor(
      Scenario,
      levels = 1:4,
      labels = c("1–11 yrs", "12–17 yrs", "18–59 yrs", "60+ yrs")
    ),
    Case_averted_pct        = Case_averted_prop * 100,
    pct_label_case          = sprintf("%.0f", Case_averted_pct),
    Indirect_pct            = IndirectContribution_Cases * 100,
    pct_label_indirect      = sprintf("%.0f", Indirect_pct)
    
  )

h1 <- 
ggplot(plot_df2,
       aes(Supply_f, Delay_f, fill = Case_averted_pct)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  #geom_text(aes(label = pct_label_case), size = 2.2) +
  facet_wrap(~ Scenario_f, nrow = 1) +
  scale_fill_gradient(low = "#fff7ec", high = "#fb6a4a",
                      limits = c(0,100),
                      name   = "% Total symptomatic \n cases averted",
                      labels = function(x) paste0(x, "%")) +
  labs(x = "Vaccine coverage", y = "Delay (weeks)") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold"))+ theme(axis.title.x = element_blank())
h2 <- 
ggplot(plot_df2,
       aes(Supply_f, Delay_f, fill = Indirect_pct     )) +
  geom_tile(colour = "white", linewidth = 0.25) +
  #geom_text(aes(label = pct_label_indirect), size = 2.2) +
  facet_wrap(~ Scenario_f, nrow = 1) +
  scale_fill_gradient(low = "#fff7ec", high = "#fb6a4a",
                      limits = c(0,100),
                      name   = "% Indirectly averted \n cases",
                      labels = function(x) paste0(x, "%")) +
  labs(x = "Vaccine coverage", y = "Delay (weeks)") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold"))
heat <- 
(h1 / h2) + 
  plot_annotation(tag_levels = list(c("A","B")))+
  plot_layout(guides = "collect") 

ggsave(filename = "02_Outputs/2_1_Figures/fig3_heat.jpg", heat, width = 10, height = 7, dpi = 1200)
