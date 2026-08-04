ibge_bra_pop <- read_excel("00_Data/0_1_Raw/ibge_bra_pop.xlsx", 
                           sheet = "Sheet1")

ibge_long <- ibge_bra_pop %>% pivot_longer(
  cols = -c(code, state),
  names_to = "age",
  values_to = "pop"
)

ibge_long <- ibge_long %>% mutate(
  age_group = case_when(
    age == "<1" ~ "<1",
    age %in% as.character(1:4) ~ "1-4",
    age %in% as.character(5:9) ~ "5–9",
    age %in% as.character(10:11) ~ "10-11",
    age %in% as.character(12:17) ~ "12-17",
    age %in% as.character(18:19) ~ "18–19",
    age %in% as.character(20:24) ~ "20–24",
    age %in% as.character(25:29) ~ "25–29",
    age %in% as.character(30:34) ~ "30–34",
    age %in% as.character(35:39) ~ "35–39",
    age %in% as.character(40:44) ~ "40–44",
    age %in% as.character(45:49) ~ "45–49",
    age %in% as.character(50:54) ~ "50–54",
    age %in% as.character(55:59) ~ "55–59",
    age %in% as.character(60:64) ~ "60–64",
    age %in% as.character(65:69) ~ "65–69",
    age %in% as.character(70:74) ~ "70–74",
    age %in% as.character(75:79) ~ "75–79",
    age %in% as.character(80:84) ~ "80–84",
    TRUE ~ "85+"
  )
)

ibge_agegr <- ibge_long %>% group_by(code, state, age_group) %>%
  summarise(pop = sum(pop, na.rm = T), .groups = "drop")

age_group_levels <- c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
)

ibge_agegr <- ibge_agegr %>% mutate(
  age_group = factor(age_group, levels = age_group_levels, ordered = TRUE)
)

ibge_agegr_clean <- ibge_agegr %>%
  group_by(age_group, state) %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

ibge_agegr_clean <- ibge_agegr_clean %>%
  mutate(age_group = factor(age_group, levels = age_group_levels, ordered = TRUE))

ibge_agegr_wide <- ibge_agegr_clean %>%
  pivot_wider(
    names_from = state,
    values_from = pop
  ) %>%
  arrange(age_group)

save(ibge_agegr_wide, file = "00_Data/0_2_Processed/ibge_agegr.RData")


ibge_long <- ibge_agegr_wide %>%
  pivot_longer(
    cols = -age_group,
    names_to = "state_full",
    values_to = "population"
  )

pop_by_state <- ibge_long %>% group_by(state_full) %>% summarise(
  tot_pop = sum(population)
)%>%
  dplyr::rename(region = state_full)

save(pop_by_state, file = "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_benefit_risk/01_Data/pop_by_state.RData")

