bra_all_ages_21 <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                               sheet = "bra_all_ages_2021")

bra_all_ages_22 <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                              sheet = "bra_all_ages_2022")
bra_all_ages_region <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                              sheet = "bra_region_age_week_2022")


bra_all_sum_21 <- bra_all_ages_21 %>% group_by(epi_week) %>% # use epi_week when only using 2022 data
  summarise(
  epi_week  = first(epi_week),
  tot_cases = sum(cases),
  region    = first(region)
)

bra_all_sum_22 <- bra_all_ages_22 %>% group_by(week) %>% # use epi_week when only using 2022 data
  summarise(
    epi_week  = first(week),
    tot_cases = sum(cases),
    region    = first(region)
  )


age_gr_levels <- c("<1",
                    "1-4",
                    "5-9",
                    "10-14",
                    "15-19",
                    "20-39",
                    "40-59",
                    "60-64",
                    "65-69",
                    "70-79",
                    ">80")
  
bra_all_ages$age_gr <- factor(bra_all_ages$age_gr, levels = age_gr_levels)

model_age_bins <- c("<5 years",
                   "5-9 years",
                   "10-14 years",
                   "15-19 years",
                   "20-24 years",
                   "25-29 years",
                   "30-34 years",
                   "35-39 years",
                   "40-44 years",
                   "45-49 years",
                   "50-54 years",
                   "55-59 years",
                   "60-64 years",
                   "65-69 years",
                   "70-74 years",
                   "75-79 years",
                   "80-84 years",
                   "85-89 years")

bra_expanded_21 <- bra_all_ages_21 %>%
  mutate(age_gr = if_else(age_gr %in% c("<1", "1-4"), "<5", age_gr)) %>%
  group_by(epi_week, age_gr, region) %>%
  summarise(cases = sum(cases), .groups = "drop")

bra_expanded_22 <- bra_all_ages_22 %>%
  mutate(age_gr = if_else(age_gr %in% c("<1", "1-4"), "<5", age_gr)) %>%
  group_by(week, age_gr, region) %>%
  summarise(cases = sum(cases), .groups = "drop")

mapping_list <- list(
  "<5"     = 1,         # Aggregated "<5" goes to model bin 1 ("<5 years")
  "5-9"    = 2,         # "5-9" goes to model bin 2 ("5-9 years")
  "10-14"  = 3,         # "10-14" goes to model bin 3 ("10-14 years")
  "15-19"  = 4,         # "15-19" goes to model bin 4 ("15-19 years")
  "20-39"  = 5:8,       # "20-39" spans bins 5-8 ("20-24", "25-29", "30-34", "35-39")
  "40-59"  = 9:12,      # "40-59" spans bins 9-12 ("40-44", "45-49", "50-54", "55-59")
  "60-64"  = 13,        # "60-64" goes to bin 13 ("60-64 years")
  "65-69"  = 14,        # "65-69" goes to bin 14 ("65-69 years")
  "70-79"  = 15:16,     # "70-79" spans bins 15-16 ("70-74", "75-79")
  ">80"    = 17:18      # ">80" spans bins 17-18 ("80-84", "85-89")
)

bra_expanded_21 <- bra_expanded_21 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_expanded_21 <- bra_expanded_21 %>%
  group_by(epi_week, age_gr, region) %>%
  mutate(pop_total = sum(N_2023[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_2023[age_bin] / pop_total))

bra_expanded_21 <- bra_expanded_21 %>%
  group_by(epi_week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(epi_week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_expanded_21$age_bin_label <- factor(bra_expanded_21$age_bin_label, levels = model_age_bins)


## 2022
bra_expanded_22 <- bra_expanded_22 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_expanded_22 <- bra_expanded_22 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_2023[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_2023[age_bin] / pop_total))

bra_expanded_22 <- bra_expanded_22 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_expanded_22$age_bin_label <- factor(bra_expanded_22$age_bin_label, levels = model_age_bins)


ggplot(bra_all_ages)+
  geom_line(aes(x = epi_week, y = cases, color = factor(age_gr)))+
  theme_minimal()

ggplot(bra_expanded_21)+
  geom_line(aes(x = epi_week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_expanded_22)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_all_sum_21)+
  geom_line(aes(x = epi_week, y = tot_cases))+
  theme_minimal()

ggplot(bra_all_sum_22)+
  geom_line(aes(x = week, y = tot_cases))+
  theme_minimal()
