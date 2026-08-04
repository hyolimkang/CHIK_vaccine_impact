library(lubridate)
load("00_Data/0_2_Processed/tot_pop_df.RData")
chikv_vac_cases_fix <- readRDS("00_Data/0_1_Raw/chikv_vac_cases_fix.RDS")
chikv_vac_cases_fix <- na.omit(chikv_vac_cases_fix)
chikv_vac_cases_fix <- chikv_vac_cases_fix %>% filter(age != 0)
chikv_vac_2022 <- chikv_vac_cases_fix %>% 
  filter(year(symptom_onset) == 2022)

chikv_vac_cases_sum <- chikv_vac_cases_fix %>%
  mutate(
    symptom_onset = ymd(symptom_onset),              # 날짜 형식 변환
    epi_year = epiweek(symptom_onset) %/% 100 + 2000, # lubridate의 epiweek는 epiyear 미포함
    epi_week = epiweek(symptom_onset),
    epi_year = year(symptom_onset + days(3 - wday(symptom_onset))) # 정확한 epiyear 보정
  ) %>%
  group_by(epi_year, epi_week, code_state, age) %>%
  summarise(cases = n(), .groups = "drop")

chikv_2022 <- chikv_vac_cases_sum %>% filter(epi_year == 2022)


age_gr_levels <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                   "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                   "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                   "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                   "80-84 years", "85+ years")

chikv_2022 <- chikv_2022 %>%
  mutate(
    age_group = cut(
      age,
      breaks = c(seq(0, 85, by = 5), Inf),  # 90세 이상도 마지막 구간에 포함
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

target_states <- c("17", "22", "23", "24", "25", "26", "27", "28", "29", "31", "52")
chikv_2022_filter <- chikv_2022_sum %>% mutate(code_state = sprintf("%02d", as.numeric(code_state))) %>%
  filter(code_state %in% target_states)

weeks <- 1:52
age_levels <- levels(factor(chikv_2022_filter$age_group))  # 18개 age group 자동 추출

chikv_2022_full <- chikv_2022_filter %>%
  complete(
    nesting(epi_year, code_state, state_full),  # 각 주 단위로
    epi_week = 1:52,
    age_group = age_levels,
    fill = list(cases = 0)
  )

state_pop <- tibble(
  state_full = c("Tocantins", "Bahia", "Alagoas",
                 "Ceará", "Minas Gerais",
                 "Pernambuco", "Paraíba", "Rio Grande do Norte",
                 "Piauí", "Sergipe", "Goiás"),
  total_pop = c(sum(N_tc$Tocantins), sum(N_bahia$Bahia), sum(N_ag$Alagoas),
                sum(N_ceara$Ceará), sum(N_mg$`Minas Gerais`), sum(N_pemam$Pernambuco),
                sum(N_pa$Paraíba), sum(N_rg$`Rio Grande do Norte`), sum(N_pi$Piauí), sum(N_se$Sergipe), sum(N_go$Goiás))
)

chikv_2022_full <- chikv_2022_full %>%
  left_join(tot_pop_df, by = "state_full")

chik_2022_full_summ <- chikv_2022_sum %>% group_by(state_full, epi_week) %>% summarise(
  tot_cases = sum(cases)
)

chik_2022_full_summ <- chik_2022_full_summ %>% 
  left_join(tot_pop_df, by = "state_full")

chik_2022_full_summ <- chik_2022_full_summ %>% mutate(
  per_capita_case = tot_cases / tot_pop * 1000000
)

chik_2022_full_agg <- chik_2022_full_summ %>% group_by(state_full, tot_pop) %>% summarise(
  tot_cases = sum(tot_cases)
) %>% mutate(
  per_capita_case = tot_cases / tot_pop * 1000000
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
chikv_all <- chikv_2022_full %>% group_by(epi_year, code_state, epi_week, state_full, tot_pop) %>% summarise(
  cases = sum(cases)
) %>% group_by(code_state) %>% mutate(
  peak_per_million = max(cases) / tot_pop * 1000000
)


chikv_all_summ <- chikv_all %>% group_by(state_full, tot_pop) %>% summarise(
  cum_case = sum(cases)
) %>% mutate(
  case_per_1m = cum_case / tot_pop * 1000000
)

bra_all_sum <- chikv_all
colnames(bra_all_sum)[3] <- c("Week")
colnames(bra_all_sum)[6] <- c("Observed")

bra_all_sum <- bra_all_sum %>% group_by(Week) %>% summarise(
  Observed = sum(Observed)
)

p <- 
ggplot(chik_2022_full_summ)+
  geom_line(aes(x = epi_week, y= tot_cases))+
  facet_wrap(~state_full, scales = "free_y")+
  theme_pubclean()
ggsave(filename = "02_Outputs/2_1_Figures/full_raw_cases.jpg", p, width = 12, height = 10, dpi = 1200)

p1 <- 
  chik_2022_full_agg %>%
  mutate(state_full = fct_reorder(state_full, per_capita_case, .desc = TRUE)) %>%
  ggplot(aes(x = state_full, y = per_capita_case)) +
  geom_bar(stat = "identity", fill = state_full) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Attack rate") +
  scale_y_continuous(labels = comma)
p2 <- 
  chik_2022_full_agg %>%
  mutate(state_full = fct_reorder(state_full, tot_cases, .desc = TRUE)) %>%  # Reorder by total_cases in descending order
  ggplot() +
  geom_bar(aes(x = state_full, y = tot_cases, fill = state_full), stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Total cumulative reported cases")+
  scale_y_continuous(labels = comma)

region_order <- chik_2022_full_agg %>%
  arrange(desc(tot_cases)) %>%
  pull(state_full)

combined_metric <- chik_2022_full_agg %>%
  select(state_full, tot_cases, per_capita_case) %>%
  pivot_longer(cols = c(per_capita_case, tot_cases),
               names_to = "metric",
               values_to = "value") %>%
  mutate(region = factor(state_full, levels = region_order),
         metric = factor(metric, levels = c("tot_cases", "per_capita_case")))


p <- combined_metric %>%
  ggplot(aes(
    x = region,  # ← no reordering here!
    y = value,
    fill = region
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left",
             labeller = labeller(metric = c(
               tot_cases = "Total reported cases",
               per_capita_case = "Per capita case (per million)"
             ))) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.placement = "outside",
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(filename = "02_Outputs/2_1_Figures/outbreak_chr_state.jpg", p, width = 12, height = 6, dpi = 1200)


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


# peak case per million


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
     bra_sum_ag, bra_sum_bh, bra_sum_ce,
     bra_sum_mg, bra_sum_pa, bra_sum_pi,
     bra_sum_pn, bra_sum_rg, bra_sum_tc,
     file = "00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")

#-------------------------------------------------------------------------------
bra_ac_22 <- chikv_2022_full %>% filter(state_full == "Acre")
bra_ap_22 <- chikv_2022_full %>% filter(state_full == "Amapa")
bra_am_22 <- chikv_2022_full %>% filter(state_full == "Amazonas")
bra_df_22 <- chikv_2022_full %>% filter(state_full == "Distrito Federal")
bra_es_22 <- chikv_2022_full %>% filter(state_full == "Espírito Santo")
bra_go_22 <- chikv_2022_full %>% filter(state_full == "Goiás")
bra_ma_22 <- chikv_2022_full %>% filter(state_full == "Maranhão")
bra_mt_22 <- chikv_2022_full %>% filter(state_full == "Mato Grosso")
bra_ms_22 <- chikv_2022_full %>% filter(state_full == "Mato Grosso do Sul")
bra_pr_22 <- chikv_2022_full %>% filter(state_full == "Paraná")
bra_para_22 <- chikv_2022_full %>% filter(state_full == "Pará")
bra_rs_22 <- chikv_2022_full %>% filter(state_full == "Rio Grande do Sul")
bra_rj_22 <- chikv_2022_full %>% filter(state_full == "Rio de Janeiro")
bra_ro_22 <- chikv_2022_full %>% filter(state_full == "Rondônia")
bra_rr_22 <- chikv_2022_full %>% filter(state_full == "Roraima")
bra_sc_22 <- chikv_2022_full %>% filter(state_full == "Santa Catarina")
bra_se_22 <- chikv_2022_full %>% filter(state_full == "Sergipe")
bra_sp_22 <- chikv_2022_full %>% filter(state_full == "São Paulo")

bra_sum_ac <- chikv_all %>% filter(state_full == "Acre")
bra_sum_ap <- chikv_all %>% filter(state_full == "Amapa")
bra_sum_am <- chikv_all %>% filter(state_full == "Amazonas")
bra_sum_df <- chikv_all %>% filter(state_full == "Distrito Federal")
bra_sum_es <- chikv_all %>% filter(state_full == "Espírito Santo")
bra_sum_go <- chikv_all %>% filter(state_full == "Goiás")
bra_sum_ma <- chikv_all %>% filter(state_full == "Maranhão")
bra_sum_mt <- chikv_all %>% filter(state_full == "Mato Grosso")
bra_sum_ms <- chikv_all %>% filter(state_full == "Mato Grosso do Sul")
bra_sum_pr <- chikv_all %>% filter(state_full == "Paraná")
bra_sum_para <- chikv_all %>% filter(state_full == "Pará")
bra_sum_rs <- chikv_all %>% filter(state_full == "Rio Grande do Sul")
bra_sum_rj <- chikv_all %>% filter(state_full == "Rio de Janeiro")
bra_sum_ro <- chikv_all %>% filter(state_full == "Rondônia")
bra_sum_rr <- chikv_all %>% filter(state_full == "Roraima")
bra_sum_sc <- chikv_all %>% filter(state_full == "Santa Catarina")
bra_sum_se <- chikv_all %>% filter(state_full == "Sergipe")
bra_sum_sp <- chikv_all %>% filter(state_full == "São Paulo")

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
