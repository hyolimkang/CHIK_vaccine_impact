df_long <- df_nnv_ce %>% 
  pivot_longer(
    cols      = c(direct_per1M, indirect_per1M),
    names_to  = "effect_type",
    values_to = "per1M_value"
  )
df_long <- df_long %>%
  mutate(
    # 백신맞은 그룹 여부
    vaccinated_age = case_when(
      scenario == "Scenario_1" & AgeGroup %in% c("1","2","3","4")   ~ TRUE,
      scenario == "Scenario_2" & AgeGroup %in% c("5","6","7","8","9","10","11") ~ TRUE,
      scenario == "Scenario_3" & AgeGroup %in% c("12","13","14","15","16","17","18") ~ TRUE,
      TRUE ~ FALSE
    )
  )

ggplot(df_long, aes(x = AgeGroup, y = per1M_value, fill = effect_type)) +
  geom_col(position = "stack", width = 0.7, color = "gray20") +
  facet_wrap(~ scenario) +
  scale_fill_manual(
    values = c(
      "direct_per1M"   = "#1b9e77",
      "indirect_per1M" = "#d95f02"
    ),
    labels = c("Direct", "Indirect")
  ) +
  labs(
    x = "Age Group",
    y = "Per-dose Impact (per 1 FVP)",
    fill = "Effect Type"
  ) +
  #theme_minimal(base_size = 12) +
  theme(
    legend.position      = "bottom",
    panel.grid.major.x   = element_blank()
  )

ggplot(df_nnv_ce, aes(x=AgeGroup, y=pre_vacc)) + 
  geom_bar(stat = "identity")

ggplot(df_nnv_ce, aes(x=AgeGroup, y=indirect_per1M)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~scenario)

ggplot(df_nnv_ce, aes(x=AgeGroup, y=overall_per1M)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~scenario)

ggplot(df_nnv_ce, aes(x=AgeGroup, y=direct_per1M)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~scenario)


ggplot(df_nnv_ce, aes(x=AgeGroup, y=tot_pop)) + 
  geom_bar(stat = "identity")

ggplot(df_nnv_ce, aes(x=AgeGroup, y=tot_vacc_prop)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~scenario)


prevacc_prop <- df_nnv_go$pre_vacc / df_nnv_go$tot_pop

df_nnv_go$prevacc_prop <- prevacc_prop

ggplot(df_nnv_go, aes(x=AgeGroup, y=prevacc_prop)) + 
  geom_bar(stat = "identity")

vacc_ce <- df_nnv_ce %>% group_by(scenario) %>% summarise(
  tot_vacc = sum(tot_vacc)
)

vacc_impact_ce <- df_nnv_ce %>% group_by(scenario) %>% summarise(
  pre_vacc_tot  = sum(pre_vacc),
  post_vacc_tot = sum(Median),
  RR = 1 - (post_vacc_tot / pre_vacc_tot),
  overall_per_fvp = sum(overall_per1M),
  overall_daly_fvp = sum(overall_daly_per1M),
  overall_fatal_fvp = sum(overall_fatal_per1M),
  diff_fatal = sum(diff_fatal),
  diff_daly = sum(diff_daly),
  diff = sum(diff),
  tot_vacc= sum(tot_vacc)
)

View(preui_ce)

R_ce <- as.data.frame(preui_ce$sim_out$R)
S_ce <- as.data.frame(preui_ce$sim_out$S)
I_ce <- as.data.frame(preui_ce$sim_out$I)

phi_prevacc    <- as.data.frame(preui_ce[[1]]$sim_out$phi)
phi_vec_s1 <- as.data.frame(postsim_ce_ui[[1]]$sim_out$phi)
phi_vec_s2 <- as.data.frame(postsim_ce_ui[[2]]$sim_out$phi)
phi_vec_s3 <- as.data.frame(postsim_ce_ui[[3]]$sim_out$phi)
phi_prevacc <- as.data.frame(preui_ce$sim_out$phi)

phi_vec_s1$scenario <- "1-19"
phi_vec_s2$scenario <- "20-59"
phi_vec_s3$scenario <- ">60"
phi_prevacc$scenario <- "prevacc"

colnames(phi_vec_s1) <- c("FOI", "scenario")
colnames(phi_vec_s2) <- c("FOI", "scenario")
colnames(phi_vec_s3) <- c("FOI", "scenario")
colnames(phi_prevacc) <- c("FOI", "scenario")

phi_vec <- rbind(phi_vec_s1, phi_vec_s2, phi_vec_s3, phi_prevacc)
phi_vec$time <- rep(1:52, 4)

ggplot(phi_vec)+
  geom_line(aes(x = time, y = FOI, color = scenario))


# scenario 3 post
R_post_ce <- as.data.frame(postsim_all_ce_ui$scenario_result[[1]]$sim_out$R)
S_post_ce <- as.data.frame(postsim_all_ce_ui$scenario_result[[1]]$sim_out$S)
I_post_ce <- as.data.frame(postsim_all_ce_ui$scenario_result[[1]]$sim_out$I)

T = 52
df_S <- as.data.frame(S_ce)
colnames(df_S) <- paste0("S", seq_len(T))
df_S <- df_S %>%
  mutate(AgeGroup = age_gr_levels) %>%  # 행에 age_labels 할당
  select(AgeGroup, everything())

df_S_post <- as.data.frame(S_post_ce)
colnames(df_S_post) <- paste0("S", seq_len(T))
df_S_post <- df_S_post %>%
  mutate(AgeGroup = age_gr_levels) %>%  # 행에 age_labels 할당
  select(AgeGroup, everything())

df_I <- as.data.frame(I_ce)
colnames(df_I) <- paste0("I", seq_len(T))
df_I <- df_I %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything())

df_I_post <- as.data.frame(I_post_ce)
colnames(df_I_post) <- paste0("I", seq_len(T))
df_I_post <- df_I_post %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything())

df_R <- as.data.frame(R_ce)
colnames(df_R) <- paste0("R", seq_len(T))
df_R <- df_R %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything())

df_R_post <- as.data.frame(R_post_ce)
colnames(df_R_post) <- paste0("R", seq_len(T))
df_R_post <- df_R_post %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything())

df_all <- df_S %>%
  left_join(df_I, by = "AgeGroup") %>%
  left_join(df_R, by = "AgeGroup")

df_all_post <- df_S_post %>%
  left_join(df_I_post, by = "AgeGroup") %>%
  left_join(df_R_post, by = "AgeGroup")

df_all_pre  <- df_all  %>% mutate(type = "pre")
df_all_post <- df_all_post %>% mutate(type = "post")
df_combined <- bind_rows(df_all_pre, df_all_post)

df_long <- df_combined %>%
  pivot_longer(
    cols      = -c(AgeGroup, type),   
    names_to  = "metric",
    values_to = "count"
  ) %>%
  mutate(
    compartment = case_when(
      startsWith(metric, "S") ~ "S",
      startsWith(metric, "I") ~ "I",
      startsWith(metric, "R") ~ "R"
    ),
    week = as.integer(sub("^[SIR]", "", metric))  # "S12" → 12, "I34" → 34 ...
  ) %>%
  select(AgeGroup, week, compartment, type, count)

pop_ce <- as.data.frame(N_ceara$Ceará)
pop_ce$AgeGroup <- age_gr_levels
colnames(pop_ce) <- c("tot_pop", "AgeGroup")

df_long <- df_long %>%
  left_join(pop_ce, by = "AgeGroup") %>%
  mutate(
    prop = count / tot_pop 
  )%>% 
  mutate(AgeGroup = factor(AgeGroup, levels = age_gr_levels))

df_long %>%
  filter(compartment == "S") %>%
  ggplot(aes(x = week, y = prop, color = type)) +
  geom_line(size = 0.8) +
  facet_wrap(~AgeGroup)

df_long %>%
  filter(compartment == "I") %>%
  ggplot(aes(x = week, y = count, color = type)) +
  geom_line(size = 0.8) +
  facet_wrap(~AgeGroup)

df_long %>%
  filter(compartment == "R") %>%
  ggplot(aes(x = week, y = prop, color = type)) +
  geom_line(size = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01))+
  facet_wrap(~AgeGroup)

df_long %>%
  filter(compartment == "R") %>%
  ggplot(aes(x = week, y = count, color = type)) +
  geom_line(size = 0.8) +
  #scale_y_continuous(labels = scales::percent_format(accuracy = 0.01))+
  facet_wrap(~AgeGroup)



##
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) 시나리오별 반복 처리를 위해 인덱스 벡터 정의
scenarios <- 1:3

# 2) 시나리오별로 pre/post 데이터를 뽑아서 하나의 데이터프레임으로 합치는 함수
build_df_long_per_scenario <- function(sc_idx) {
  df_S_pre <- as.data.frame(S_ce)
  colnames(df_S_pre) <- paste0("S", seq_len(ncol(df_S_pre)))
  df_S_pre <- df_S_pre %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  df_I_pre <- as.data.frame(I_ce)
  colnames(df_I_pre) <- paste0("I", seq_len(ncol(df_I_pre)))
  df_I_pre <- df_I_pre %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  df_R_pre <- as.data.frame(R_ce)
  colnames(df_R_pre) <- paste0("R", seq_len(ncol(df_R_pre)))
  df_R_pre <- df_R_pre %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  df_all_pre <- df_S_pre %>%
    left_join(df_I_pre, by = "AgeGroup") %>%
    left_join(df_R_pre, by = "AgeGroup") %>%
    mutate(type = "pre", scenario = paste0("Scenario_", sc_idx))
  
  # --- post (시나리오 적용 후) 데이터 추출 ---
  S_post <- as.data.frame(postsim_all_ce_ui$scenario_result[[sc_idx]]$sim_out$S)
  I_post <- as.data.frame(postsim_all_ce_ui$scenario_result[[sc_idx]]$sim_out$I)
  R_post <- as.data.frame(postsim_all_ce_ui$scenario_result[[sc_idx]]$sim_out$R)
  
  colnames(S_post) <- paste0("S", seq_len(ncol(S_post)))
  df_S_post <- S_post %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  colnames(I_post) <- paste0("I", seq_len(ncol(I_post)))
  df_I_post <- I_post %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  colnames(R_post) <- paste0("R", seq_len(ncol(R_post)))
  df_R_post <- R_post %>%
    mutate(AgeGroup = age_gr_levels) %>%
    select(AgeGroup, everything())
  
  df_all_post <- df_S_post %>%
    left_join(df_I_post, by = "AgeGroup") %>%
    left_join(df_R_post, by = "AgeGroup") %>%
    mutate(type = "post", scenario = paste0("Scenario_", sc_idx))
  
  # 3) pre & post 합치기
  df_combined <- bind_rows(df_all_pre, df_all_post)
  
  # 4) wide-to-long, compartment/week 열 생성
  df_long <- df_combined %>%
    pivot_longer(
      cols      = -c(AgeGroup, type, scenario),
      names_to  = "metric",
      values_to = "count"
    ) %>%
    mutate(
      compartment = case_when(
        startsWith(metric, "S") ~ "S",
        startsWith(metric, "I") ~ "I",
        startsWith(metric, "R") ~ "R"
      ),
      week = as.integer(sub("^[SIR]", "", metric))
    ) %>%
    select(AgeGroup, week, compartment, type, scenario, count)
  
  # 5) 인구수 정보 병합 → prop 계산 (원한다면)
  pop_ce_df <- as.data.frame(N_ceara$Ceará)
  pop_ce_df$AgeGroup <- age_gr_levels
  colnames(pop_ce_df) <- c("tot_pop", "AgeGroup")
  
  df_long <- df_long %>%
    left_join(pop_ce_df, by = "AgeGroup") %>%
    mutate(prop = count / tot_pop)
  
  # 6) AgeGroup 순서 보장
  df_long <- df_long %>%
    mutate(AgeGroup = factor(AgeGroup, levels = age_gr_levels))
  
  return(df_long)
}

# 3) 세 시나리오 모두 하나로 묶기
df_long_all <- bind_rows(lapply(scenarios, build_df_long_per_scenario))


# 4) ggplot 예시: compartment = "S", 선 색깔 = 시나리오, 선종류 = pre/post
ggplot(
  df_long_all %>% filter(compartment == "S"),
  aes(x = week, y = count,
      color    = scenario,
      linetype = type       # pre는 실선(“solid”), post는 점선(“dashed”) 등으로 자동 지정됨
  )
) +
  geom_line(size = 0.8) +
  facet_wrap(~AgeGroup) +
  labs(
    title    = "각 연령별 Susceptible(S) 인원수 변화",
    subtitle = "색깔: 시나리오 1/2/3, 선종류: pre vs post",
    x        = "Week",
    y        = "Count (log10 scale)",
    color    = "시나리오",
    linetype = "Type"
  ) 

ggplot(
  df_long_all %>% 
    filter(compartment == "I"),
  aes(x = week, y = count,
      color    = scenario,
      linetype = type
  )
) +
  geom_line(size = 0.8) +
  facet_wrap(~AgeGroup) +
  labs(
    title    = "각 연령별 Infectious(I) 인원수 변화",
    subtitle = "색깔: 시나리오 1/2/3, 선종류: pre vs post",
    x        = "Week",
    y        = "Count",
    color    = "시나리오",
    linetype = "Type"
  ) 


ggplot(
  df_long_all %>% filter(scenario == "Scenario_3") %>% filter(compartment == "R"),
  aes(x = week, y = prop,
      color    = scenario,
      linetype = type
  )
) +
  geom_line(size = 0.8) +
  facet_wrap(~AgeGroup) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title    = "각 연령별 Recovered(R) 비율 변화",
    subtitle = "색깔: 시나리오 1/2/3, 선종류: pre vs post",
    x        = "Week",
    y        = "Proportion",
    color    = "시나리오",
    linetype = "Type"
  ) 


S_init <- S_ce[, 1]      # length 18인 벡터
target_list <- list(
  Scenario1 = 1:4,        # <20세 풀
  Scenario2 = 5:12,       # 20–59세 풀
  Scenario3 = 13:18       # ≥60세 풀
)

# 시나리오별로 “week 1 S 합계” 계산
sus_sum <- sapply(target_list, function(idx) sum(S_init[idx]))

total_supply <- floor(min(sus_sum)) 

loss_frac  <- 1 - exp(-sum(phi_vec[1:2]))

eff1 <- sus_sum[1] * (1 - loss_frac)
eff2 <- sus_sum[2] * (1 - loss_frac)
eff3 <- sus_sum[3] * (1 - loss_frac)


