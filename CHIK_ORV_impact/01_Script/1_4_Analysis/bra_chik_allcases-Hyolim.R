library(lubridate)
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

target_states <- c("17", "22", "23", "24", "25", "26", "27", "29", "31")
chikv_2022_filter <- chikv_2022_sum %>% mutate(code_state = sprintf("%02d", as.numeric(code_state))) %>%
  filter(code_state %in% target_states)

## filter by state
bra_ag <- chikv_2022_filter %>% filter(state_full == "Alagoas")
bra_ce <- chikv_2022_filter %>% filter(state_full == "Ceará")
ggplot(bra_ag)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()

ggplot(bra_ce)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_group)+
  theme_minimal()


