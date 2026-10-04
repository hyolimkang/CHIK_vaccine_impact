bra_all_ages_region_23 <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                                     sheet = "bra_region_age_week_2023")

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

brazil_pop_transposed <- brazil_pop_transposed %>%
  tibble::rownames_to_column(var = "age_index")

brazil_pop_long <- brazil_pop_transposed %>%
  mutate(age_index = as.integer(age_index)) %>%
  mutate(age_gr = model_age_bins[age_index]) %>%  # Map each numeric index to your age group labels
  select(-age_index) %>%                          # Remove the temporary index column
  pivot_longer(
    cols = -age_gr,                              # Pivot all columns except age_gr
    names_to = "region",                         # The column names become the region
    values_to = "pop"                            # The values become population counts
  )


bra_pop_all <- brazil_pop_long %>% group_by(region) %>% summarise(
  tot_pop = sum(pop)
)

bra_pop_all <- bra_pop_all[-28,]

bra_all_ages_region_23_sum <- bra_all_ages_region_23 %>% group_by(region, week) %>% # use epi_week when only using 2022 data
  summarise(
    epi_week  = first(week),
    tot_cases = sum(cases)
  )

bra_all_sum_23 <- bra_all_ages_region_23_sum %>% group_by(week) %>% # use epi_week when only using 2022 data
  summarise(
    epi_week  = first(week),
    tot_cases = sum(tot_cases)
  )

colnames(bra_all_sum_23)[1] <- "Week"
colnames(bra_all_sum_23)[3] <- "Observed"

brazil_pop_long %>%
  filter(region != "tot_pop") %>%
  mutate(age_gr = factor(age_gr, levels = unique(age_gr))) %>%
  ggplot() +
  geom_bar(aes(x = age_gr, y = pop), stat = "identity") +
  facet_wrap(~region) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# split into regions 
bra_sum_bh_23 <- bra_all_ages_region_23_sum %>% filter(region == "Bahia"& week <=41)
bra_sum_ce_23 <- bra_all_ages_region_23_sum %>% filter(region == "Ceará")
bra_sum_mg_23 <- bra_all_ages_region_23_sum %>% filter(region == "Minas Gerais"& week <=36)
bra_sum_ms_23 <- bra_all_ages_region_23_sum %>% filter(region == "Mato Grosso do Sul"& week <=41)
bra_sum_sp_23 <- bra_all_ages_region_23_sum %>% filter(region == "São Paulo"& week <=31)
bra_sum_tc_23 <- bra_all_ages_region_23_sum %>% filter(region == "Tocantins"& week <=30)


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

bra_all_ages_region_23$age_gr <- factor(bra_all_ages_region_23$age_gr, levels = age_gr_levels)

bra_expanded_region_23 <- bra_all_ages_region_23 %>%
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

bra_expanded_region_23 <- bra_expanded_region_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

# by region
bra_bh_23 <- bra_expanded_region_23 %>% filter(region == "Bahia", week <= 41)
bra_ce_23 <- bra_expanded_region_23 %>% filter(region == "Ceará")
bra_mg_23 <- bra_expanded_region_23 %>% filter(region == "Minas Gerais", week <= 36) 
bra_ms_23 <- bra_expanded_region_23 %>% filter(region == "Mato Grosso do Sul", week <= 41)
bra_sp_23 <- bra_expanded_region_23 %>% filter(region == "São Paulo", week <= 31)
bra_tc_23 <- bra_expanded_region_23 %>% filter(region == "Tocantins", week <= 30)

N_bahia <- brazil_pop_transposed[,6]
N_ceara <- brazil_pop_transposed[,7]
N_mg    <- brazil_pop_transposed[,14]
N_pemam <- brazil_pop_transposed[,18]
N_go    <- brazil_pop_transposed[,10]
N_mt    <- brazil_pop_transposed[,12]
N_ms    <- brazil_pop_transposed[,13]
N_sp    <- brazil_pop_transposed[,27]
N_tc    <- brazil_pop_transposed[,28]


# bahia
bra_bh_23 <- bra_bh_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_bh_23 <- bra_bh_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_bahia[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_bahia[age_bin] / pop_total))

bra_bh_23 <- bra_bh_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_bh_23$age_bin_label <- factor(bra_bh_23$age_bin_label, levels = model_age_bins)

# ceara
bra_ce_23 <- bra_ce_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_ce_23 <- bra_ce_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_ceara[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_ceara[age_bin] / pop_total))

bra_ce_23 <- bra_ce_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_ce_23$age_bin_label <- factor(bra_ce_23$age_bin_label, levels = model_age_bins)

## mg
bra_mg_23 <- bra_mg_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_mg_23 <- bra_mg_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_mg[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_mg[age_bin] / pop_total))

bra_mg_23 <- bra_mg_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_mg_23$age_bin_label <- factor(bra_mg_23$age_bin_label, levels = model_age_bins)

# Mato Grosso do Sul
bra_ms_23 <- bra_ms_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_ms_23 <- bra_ms_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_ms[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_ms[age_bin] / pop_total))

bra_ms_23 <- bra_ms_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_ms_23$age_bin_label <- factor(bra_ms_23$age_bin_label, levels = model_age_bins)

# São Paulo
bra_sp_23 <- bra_sp_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_sp_23 <- bra_sp_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_sp[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_sp[age_bin] / pop_total))

bra_sp_23 <- bra_sp_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_sp_23$age_bin_label <- factor(bra_sp_23$age_bin_label, levels = model_age_bins)

## tocantins
bra_tc_23 <- bra_tc_23 %>%
  rowwise() %>%
  mutate(model_bins = list(mapping_list[[age_gr]])) %>%  # Assign target bins from mapping_list
  ungroup() %>%
  unnest(cols = c(model_bins)) %>%                      # Expand: one row per model bin
  mutate(age_bin = model_bins) %>%                      # Rename for clarity
  select(-model_bins)

bra_tc_23 <- bra_tc_23 %>%
  group_by(week, age_gr, region) %>%
  mutate(pop_total = sum(N_tc[age_bin])) %>%  # Sum only over the bins that the record covers
  ungroup() %>%
  mutate(cases_bin = cases * (N_tc[age_bin] / pop_total))

bra_tc_23 <- bra_tc_23 %>%
  group_by(week, region, age_bin) %>%
  summarise(cases = sum(cases_bin), .groups = "drop") %>%
  arrange(week, age_bin) %>%
  mutate(age_bin_label = model_age_bins[age_bin])

bra_tc_23$age_bin_label <- factor(bra_tc_23$age_bin_label, levels = model_age_bins)

## graphs

ggplot(bra_bh_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_ce_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_mg_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_ms_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_sp_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()

ggplot(bra_tc_23)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~age_bin_label)+
  theme_minimal()


ggplot(bra_sum_tc_23)+
  geom_line(aes(x= week, y = tot_cases))+
  theme_minimal()

combined_23_all <- bind_rows(
  bra_sum_bh_23, bra_sum_ce_23, bra_sum_mg_23,
  bra_sum_tc_23, bra_sum_ms_23, bra_sum_sp_23
)

ggplot(combined_23_all)+
  geom_line(aes(x= week, y = tot_cases))+
  theme_minimal()+
  facet_wrap(~region)

ggplot(combined_23_all)+
  geom_line(aes(x= week, y = tot_cases))+
  theme_light()+
  facet_wrap(~region, scales = "free_y")


combined_2023 <- combined_23_all %>% group_by(week) %>% summarise(
  tot_cases = sum(tot_cases)
)

ggplot(combined_2023)+
  geom_line(aes(x= week, y = tot_cases))+
  theme_light()+
  scale_y_continuous(label = comma)
