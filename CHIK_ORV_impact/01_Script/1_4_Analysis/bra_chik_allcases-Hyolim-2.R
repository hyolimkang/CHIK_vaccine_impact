library(lubridate)
chikv_vac_cases_fix <- readRDS("00_Data/0_1_Raw/chikv_vac_cases_fix.RDS")
chikv_vac_cases_fix <- na.omit(chikv_vac_cases_fix)
#chikv_vac_cases_fix <- chikv_vac_cases_fix %>% filter(age != 0)


chikv_vac_cases_sum <- chikv_vac_cases_fix %>%
  mutate(
    symptom_onset = ymd(symptom_onset),              
    epi_year = epiweek(symptom_onset) %/% 100 + 2000, 
    epi_week = epiweek(symptom_onset),
    epi_year = year(symptom_onset + days(3 - wday(symptom_onset))) 
  ) %>%
  group_by(epi_year, epi_week, code_state, age) %>%
  summarise(cases = n(), .groups = "drop")

chikv_2022 <- chikv_vac_cases_sum %>% filter(epi_year == 2022)

age_gr_levels <- c(
  "<1", "1-4", "5–9", "10-11", "12–17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
)

age_breaks <- c(0, 1, 5, 10, 12, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf)

chikv_2022 <- chikv_2022 %>%
  mutate(
    age_group = cut(
      age,
      breaks = age_breaks,         
      labels = age_gr_levels,
      right = FALSE,               
      include.lowest = TRUE        
    )
  )

chikv_2022_sum <- chikv_2022 %>%
  group_by(epi_year, epi_week, code_state, age_group) %>%
  summarise(cases = sum(cases), .groups = "drop")

state_map <- data.frame(
  code_state = sprintf("%02d", c(
    11, 12, 13, 14, 15, 16, 17,  # North
    21, 22, 23, 24, 25, 26, 27, 28, 29,  # Northeast
    31, 32, 33, 35,              # Southeast
    41, 42, 43,                  # South
    50, 51, 52, 53               # Central-West
  )),
  state_abbr = c(
    "RO", "AC", "AM", "RR", "PA", "AP", "TO",
    "MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA",
    "MG", "ES", "RJ", "SP",
    "PR", "SC", "RS",
    "MS", "MT", "GO", "DF"
  ),
  state_full = c(
    "Rondônia", "Acre", "Amazonas", "Roraima", "Pará", "Amapá", "Tocantins",
    "Maranhão", "Piauí", "Ceará", "Rio Grande do Norte", "Paraíba", "Pernambuco", "Alagoas", "Sergipe", "Bahia",
    "Minas Gerais", "Espírito Santo", "Rio de Janeiro", "São Paulo",
    "Paraná", "Santa Catarina", "Rio Grande do Sul",
    "Mato Grosso do Sul", "Mato Grosso", "Goiás", "Distrito Federal"
  )
)

chikv_2022_sum <- chikv_2022_sum %>%
  mutate(code_state = sprintf("%02d", as.numeric(code_state))) %>%
  left_join(state_map, by = "code_state")

chikv_2022_sum_collapsed <- chikv_2022_sum %>% group_by(epi_week, state_full) %>% summarise(
  cases = sum(cases)
)

ggplot(chikv_2022_sum_collapsed) +
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~state_full)

peak_cases_by_state <- chikv_2022_sum_collapsed %>%
  group_by(state_full) %>%
  summarise(
    peak_cases = max(cases, na.rm = TRUE),
    peak_week  = epi_week[which.max(cases)]
  ) %>%
  arrange(desc(peak_cases))

print(peak_cases_by_state)

peak_cases_by_state <- left_join(peak_cases_by_state, tot_pop_df, by = "state_full")
peak_cases_by_state <- peak_cases_by_state %>% mutate(
  peak_per_M  = peak_cases / tot_pop * 1e6 
)


target_states <- c("17", "22", "23", "24", "25", "26", "27", "29", "31", "28", "52")
chikv_2022_filter <- chikv_2022_sum %>% mutate(code_state = sprintf("%02d", as.numeric(code_state))) %>%
  filter(code_state %in% target_states)

weeks <- 1:52
age_levels <- levels(factor(chikv_2022_filter$age_group))  

chikv_2022_full <- chikv_2022_filter %>%
  complete(
    nesting(epi_year, code_state, state_full),  
    epi_week = 1:52,
    age_group = age_levels,
    fill = list(cases = 0)
  )%>%
  mutate(
    target = case_when(
      age_group %in% c("0-1 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("1-4 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("5-9 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("10-11 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("12-17 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("18-19 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("20-24 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("25-29 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("30-34 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("35-39 years")                  ~ "Scenario_1",  # infant
      age_group %in% c("40-44 years")                 ~ "Scenario_2",  # children
      age_group %in% c("45-49 years")                ~ "Scenario_3",  # adolescents
      age_group %in% c("50-54 years")                ~ "Scenario_4",  # adults,
      age_group %in% c("55-59 years")                ~ "Scenario_4", # adults,
      age_group %in% c("60-64 years")                  ~ "Scenario_5",  # elderlies
      age_group %in% c("65-69 years")                  ~ "Scenario_5",  # elderlies
      age_group %in% c("70-74 years")                  ~ "Scenario_5",  # elderlies
      age_group %in% c("75-79 years")                  ~ "Scenario_5",  # elderlies
      age_group %in% c("80-84 years")                  ~ "Scenario_5",  # elderlies
      age_group %in% c("85+ years")                  ~ "Scenario_5",  # elderlies
      TRUE                                           ~ NA_character_
    )
  )

chikv_2022_full_other <- chikv_2022_sum %>%
  complete(
    nesting(epi_year, code_state, state_full),  
    epi_week = 1:52,
    age_group = age_levels,
    fill = list(cases = 0)
  )

## filter by state
bra_ag_22 <- chikv_2022_full %>% filter(state_full == "Alagoas")
bra_bh_22 <- chikv_2022_full %>% filter(state_full == "Bahia")
bra_ce_22 <- chikv_2022_full %>% filter(state_full == "Ceará")
bra_mg_22 <- chikv_2022_full %>% filter(state_full == "Minas Gerais")
bra_pn_22 <- chikv_2022_full %>% filter(state_full == "Pernambuco")
bra_pa_22 <- chikv_2022_full %>% filter(state_full == "Paraíba")
bra_rg_22 <- chikv_2022_full %>% filter(state_full == "Rio Grande do Norte")
bra_pi_22 <- chikv_2022_full %>% filter(state_full == "Piauí")
bra_tc_22 <- chikv_2022_full %>% filter(state_full == "Tocantins")
bra_se_22 <- chikv_2022_full %>% filter(state_full == "Sergipe")
bra_go_22 <- chikv_2022_full %>% filter(state_full == "Goiás")


# aggregated case
chikv_all <- chikv_2022_full %>% group_by(epi_year, code_state, epi_week, state_full) %>% summarise(
  cases = sum(cases)
)

chikv_all_strategy <- chikv_2022_full %>% group_by(epi_year, code_state, epi_week, state_full, target) %>% summarise(
  cases = sum(cases)
)

  
chikv_all_other <- chikv_2022_full_other %>% group_by(epi_year, code_state, epi_week, state_full) %>% summarise(
  cases = sum(cases)
)


bra_all_sum <- chikv_all
colnames(bra_all_sum)[3] <- c("Week")
colnames(bra_all_sum)[5] <- c("Observed")

bra_all_sum_target <- chikv_all_strategy
colnames(bra_all_sum_target)[3] <- c("Week")
colnames(bra_all_sum_target)[5] <- c("Scenario")
colnames(bra_all_sum_target)[6] <- c("Observed")

bra_all_sum <- bra_all_sum %>% group_by(Week) %>% summarise(
  Observed = sum(Observed)
)

bra_all_sum_target <- bra_all_sum_target %>% group_by(Week, Scenario) %>% summarise(
  Observed = sum(Observed)
)

bra_sum_ag <- chikv_all %>% filter(state_full == "Alagoas")
bra_sum_bh <- chikv_all %>% filter(state_full == "Bahia")
bra_sum_ce <- chikv_all %>% filter(state_full == "Ceará")
bra_sum_mg <- chikv_all %>% filter(state_full == "Minas Gerais")
bra_sum_pn <- chikv_all %>% filter(state_full == "Pernambuco")
bra_sum_pa <- chikv_all %>% filter(state_full == "Paraíba")
bra_sum_rg <- chikv_all %>% filter(state_full == "Rio Grande do Norte")
bra_sum_pi <- chikv_all %>% filter(state_full == "Piauí")
bra_sum_tc <- chikv_all %>% filter(state_full == "Tocantins")
bra_sum_se <- chikv_all %>% filter(state_full == "Sergipe")
bra_sum_go <- chikv_all %>% filter(state_full == "Goiás")

# sanity check
ggplot(bra_bh_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_ce_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_mg_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_pn_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_pa_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_rg_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_pi_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_ag_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_tc_22)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

save(bra_ag_22, bra_bh_22, bra_ce_22,
     bra_mg_22, bra_pa_22, bra_pi_22,
     bra_pn_22, bra_rg_22, bra_tc_22,
     bra_se_22, bra_go_22,
     bra_sum_ag, bra_sum_bh, bra_sum_ce,
     bra_sum_mg, bra_sum_pa, bra_sum_pi,
     bra_sum_pn, bra_sum_rg, bra_sum_tc,
     bra_sum_se, bra_sum_go, bra_all_sum,
     file = "00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")

save(bra_all_sum_target, file = "00_Data/0_2_Processed/bra_all_sum_target.RData")


#### per capita 
chikv_all <- left_join(chikv_all, tot_pop_df, by = "state_full")
chikv_all_other <- left_join(chikv_all_other, tot_pop_df, by = "state_full")

chikv_all <- chikv_all %>% mutate(
  per_capita_inc = cases / tot_pop * 1e6
)

chikv_all_other <- chikv_all_other %>% mutate(
  per_capita_inc = cases / tot_pop * 1e6
)

ggplot(chikv_all)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~state_full)+
  theme_pubclean()

ggplot(chikv_all_other)+
  geom_line(aes(x = epi_week, y = per_capita_inc))+
  facet_wrap(~state_full)+
  theme_minimal()

summary_all <- chikv_all %>%
  group_by(state_full) %>%
  summarise(
    max_cases            = max(cases, na.rm = TRUE),
    max_per_capita_inc   = max(per_capita_inc, na.rm = TRUE)
  )

summary_other <- chikv_all_other %>%
  group_by(state_full) %>%
  summarise(
    max_cases            = max(cases, na.rm = TRUE),
    max_per_capita_inc   = max(per_capita_inc, na.rm = TRUE)
  )

chikv_all <- chikv_all_other %>% group_by(state_full) %>% summarise(
  tot_cases = sum(cases),
  tot_pop   = first(tot_pop),
  .groups = "drop" 
) %>% mutate(
  per_capita_inc = tot_cases / tot_pop * 1e5
)

colnames(chikv_all)[1] <- "region"

chikv_all <- left_join(chikv_all, setting_key_df, by = "region") 
chikv_all <- chikv_all %>%  mutate(setting = replace_na(setting, "Very low"))

chikv_all <- chikv_all %>% mutate(
  inc_bin = cut(
    per_capita_inc,
    breaks = quantile(per_capita_inc, probs = c(0, .25, .5, .75, 1), na.rm=TRUE),
    include.lowest = TRUE,
    labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)")
  )
)

chikv_all_tot_cases <- chikv_all[,1:2]
names(chikv_all_tot_cases) <- c("state_full", "tot_cases")
peak_cases_by_state <- left_join(peak_cases_by_state, chikv_all_tot_cases, by = "state_full")

save(chikv_all, file = "00_Data/0_2_Processed/chikv_all.RData")

chikv_long <- chikv_all %>%
  select(region, setting, tot_cases, per_capita_inc) %>%
  pivot_longer(
    cols = c(tot_cases, per_capita_inc),
    names_to = "indicator",
    values_to = "value"
  )

# 2. 막대 정렬용 factor 지정 (state_full은 그대로 유지)
chikv_long <- chikv_long %>%
  mutate(region = factor(region, levels = unique(region[order(-value)])))

y_labels <- c(
  tot_cases = "Total reported cases",
  per_capita_inc = "Per 100,000 symptomatic cases"
)

# 3. 시각화: setting별 색상 지정
p <- 
ggplot(chikv_long, aes(x = region, y = value, fill = setting)) +
  geom_col() +
  facet_grid(rows = vars(indicator), scales = "free_y", switch = "y", labeller = labeller(indicator = y_labels)) +
  scale_y_continuous(labels = comma) +
  labs(x = "State", y = NULL, fill = "Per-capita incidence") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.text.y = element_text(face = "bold")
  ) 

ggsave(filename = "02_Outputs/2_1_Figures/figs3.jpg", p, width = 10, height = 6, dpi = 1200)

chikv_plotdata <- chikv_2022_full %>%
  filter(!is.na(cases)) %>%
  group_by(state_abbr, age_group) %>%
  summarise(total_cases = sum(cases), .groups = "drop")

state_lookup <- chikv_2022_full %>%
  select(state_abbr, state_full) %>%
  distinct()

# state_full 추가
chikv_plotdata_full <- chikv_plotdata %>%
  left_join(state_lookup, by = "state_abbr")
chikv_with_pop <- chikv_plotdata_full %>%
  left_join(ibge_long, by = c("state_full", "age_group"))%>%
  filter(!is.na(population), !is.na(total_cases)) %>%
  mutate(incidence_per_100k = (total_cases / population) * 100000)

chikv_with_pop <- chikv_with_pop %>%
  mutate(age_group = factor(age_group, levels = age_gr_levels))

age_colors <- rep(brewer.pal(12, "Set3"), length.out = length(unique(chikv_with_pop$age_group)))

ggplot(chikv_with_pop, aes(x = age_group, y = incidence_per_100k)) +
  geom_col(fill = "#2c7fb8") +  # 한 가지 색상
  facet_wrap(~ state_full, scales = "free_y") +
  labs(
    x = "Age group",
    y = "Incidence per 100,000"
  ) +
  theme_pubclean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  )

scale_factor <- max(chikv_with_pop$population, na.rm = TRUE) / 
  max(chikv_with_pop$total_cases, na.rm = TRUE)

scale_factor <- max(chikv_with_pop$population, na.rm = TRUE) / 
  max(chikv_with_pop$total_cases, na.rm = TRUE)

ggplot(chikv_with_pop, aes(x = age_group)) +
  # 발생 건수 (케이스): 막대
  geom_col(aes(y = total_cases), fill = "#3182bd", alpha = 0.7) +
  
  # 인구 수: 선
  geom_line(aes(y = population / scale_factor, group = 1),
            color = "#e6550d", linewidth = 1.2) +
  
  # 이중 y축
  scale_y_continuous(
    name = "Total Cases",
    labels = comma,
    sec.axis = sec_axis(~ . * scale_factor, name = "Population", labels = comma)
  ) +
  
  # Facet by state
  facet_wrap(~ state_full, scales = "free_y") +
  
  # 라벨 및 테마
  labs(x = "Age Group") +
  theme_pubclean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  )
#-------------------------------------------------------------------------------
bra_ac_22 <- chikv_2022_full_other %>% filter(state_full == "Acre")
bra_ap_22 <- chikv_2022_full_other %>% filter(state_full == "Amapá")
bra_am_22 <- chikv_2022_full_other %>% filter(state_full == "Amazonas")
bra_df_22 <- chikv_2022_full_other %>% filter(state_full == "Distrito Federal")
bra_es_22 <- chikv_2022_full_other %>% filter(state_full == "Espírito Santo")
bra_go_22 <- chikv_2022_full_other %>% filter(state_full == "Goiás")
bra_ma_22 <- chikv_2022_full_other %>% filter(state_full == "Maranhão")
bra_mt_22 <- chikv_2022_full_other %>% filter(state_full == "Mato Grosso")
bra_ms_22 <- chikv_2022_full_other %>% filter(state_full == "Mato Grosso do Sul")
bra_pr_22 <- chikv_2022_full_other %>% filter(state_full == "Paraná")
bra_para_22 <- chikv_2022_full_other %>% filter(state_full == "Pará")
bra_rs_22 <- chikv_2022_full_other %>% filter(state_full == "Rio Grande do Sul")
bra_rj_22 <- chikv_2022_full_other %>% filter(state_full == "Rio de Janeiro")
bra_ro_22 <- chikv_2022_full_other %>% filter(state_full == "Rondônia")
bra_rr_22 <- chikv_2022_full_other %>% filter(state_full == "Roraima")
bra_sc_22 <- chikv_2022_full_other %>% filter(state_full == "Santa Catarina")
bra_se_22 <- chikv_2022_full_other %>% filter(state_full == "Sergipe")
bra_sp_22 <- chikv_2022_full_other %>% filter(state_full == "São Paulo")

bra_sum_ac <- chikv_all_other %>% filter(state_full == "Acre")
bra_sum_ap <- chikv_all_other %>% filter(state_full == "Amapá")
bra_sum_am <- chikv_all_other %>% filter(state_full == "Amazonas")
bra_sum_df <- chikv_all_other %>% filter(state_full == "Distrito Federal")
bra_sum_es <- chikv_all_other %>% filter(state_full == "Espírito Santo")
bra_sum_go <- chikv_all_other %>% filter(state_full == "Goiás") ## peak 200 
bra_sum_ma <- chikv_all_other %>% filter(state_full == "Maranhão")
bra_sum_mt <- chikv_all_other %>% filter(state_full == "Mato Grosso")
bra_sum_ms <- chikv_all_other %>% filter(state_full == "Mato Grosso do Sul")
bra_sum_pr <- chikv_all_other %>% filter(state_full == "Paraná")
bra_sum_para <- chikv_all_other %>% filter(state_full == "Pará")
bra_sum_rs <- chikv_all_other %>% filter(state_full == "Rio Grande do Sul")
bra_sum_rj <- chikv_all_other %>% filter(state_full == "Rio de Janeiro")
bra_sum_ro <- chikv_all_other %>% filter(state_full == "Rondônia")
bra_sum_rr <- chikv_all_other %>% filter(state_full == "Roraima")
bra_sum_sc <- chikv_all_other %>% filter(state_full == "Santa Catarina")
bra_sum_se <- chikv_all_other %>% filter(state_full == "Sergipe") # peak 170
bra_sum_sp <- chikv_all_other %>% filter(state_full == "São Paulo")

ggplot(bra_sum_sp)+
  geom_line(aes(x = epi_week, y = cases))+
  theme_minimal()


save(bra_ac_22, bra_ap_22, bra_am_22,
     bra_df_22, bra_es_22, bra_go_22,
     bra_ma_22, bra_mt_22, bra_ms_22,
     bra_pr_22, bra_para_22, bra_rs_22,
     bra_rj_22, bra_ro_22, bra_rr_22, 
     bra_sc_22, bra_se_22, bra_sp_22,
     bra_sum_ac, bra_sum_ap, bra_sum_am,
     bra_sum_df, bra_sum_es, bra_sum_go,
     bra_sum_ma, bra_sum_mt, bra_sum_ms,
     bra_sum_pr, bra_sum_para, bra_sum_rs,
     bra_sum_rj, bra_sum_ro, bra_sum_rr,
     bra_sum_sc, bra_sum_se, bra_sum_sp, 
     file = "00_Data/0_2_Processed/bra_cases_2022_cleaned_otherstates.RData")
