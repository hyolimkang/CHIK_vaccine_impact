ibge_bra_pop <- read_excel("00_Data/0_1_Raw/ibge_bra_pop.xlsx", 
                           sheet = "Sheet1")

ibge_long <- ibge_bra_pop %>% pivot_longer(
  cols = -c(code, state),
  names_to = "age",
  values_to = "pop"
)

ibge_long <- ibge_long %>% mutate(
  age_group = case_when(
    age == "<1" ~ "0-1",  
    age %in% as.character(1:11)   ~ "1-11",
    age %in% as.character(12:17)   ~ "12-17",
    age %in% as.character(18:39)  ~ "18-39",
    age %in% as.character(40:59) ~ "40-59",
    TRUE                          ~ "60+"
  )
)

ibge_agegr <- ibge_long %>% group_by(code, state, age_group) %>%
  summarise(pop = sum(pop, na.rm = T), .groups = "drop")

ibge_agegr_collapsed <- ibge_agegr %>%
  mutate(
    age_group = case_when(
      age_group == "0-1"                ~ "0-1 years",
      age_group == "1-11"              ~ "1-11 years",
      age_group == "12–17"             ~ "12-17 years",
      age_group == "18-39"             ~ "18-39 years",
      age_group == "40–59"             ~ "40-59 years",
      age_group == "60+"               ~ "60+ years",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  group_by(state, age_group) %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
  mutate(age_group = factor(
    age_group,
    levels = c("0-1 years", "1-11 years", "12-17 years", "18-39 years", "40-59 years", "60+ years"),
    ordered = TRUE
  ))

# ④ wide format으로 변환
ibge_agegr_wide <- ibge_agegr_collapsed %>%
  pivot_wider(
    names_from = state,
    values_from = pop
  ) %>%
  arrange(age_group)

save(ibge_agegr_wide, file = "00_Data/0_2_Processed/ibge_agegr.RData")

