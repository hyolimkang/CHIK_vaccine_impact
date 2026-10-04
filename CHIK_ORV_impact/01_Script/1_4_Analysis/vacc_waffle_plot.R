full_palette <- c(
  "<5 years"   = "#1f78b4", "5-9 years"   = "#33a02c",
  "10-14 years"= "#e31a1c", "15-19 years" = "#ff7f00",
  "20-24 years"= "#6a3d9a", "25-29 years" = "#b15928",
  "30-34 years"= "#a6cee3", "35-39 years" = "#b2df8a",
  "40-44 years"= "#fb9a99", "45-49 years" = "#fdbf6f",
  "50-54 years"= "#cab2d6", "55-59 years" = "#ffff99",
  "60-64 years"= "#8dd3c7", "65-69 years" = "#ffffb3",
  "70-74 years"= "#bebada", "75-79 years" = "#fb8072",
  "80-84 years"= "#80b1d3", "85-89 years" = "#fdb462",
  "90+ years"  = "#b3de69"
)

age_levels <- c(
  "<5 years", "5-9 years", "10-14 years", "15-19 years", "20-24 years",
  "25-29 years", "30-34 years", "35-39 years", "40-44 years", "45-49 years",
  "50-54 years", "55-59 years", "60-64 years", "65-69 years", "70-74 years",
  "75-79 years", "80-84 years", "85-89 years"
)

cb_palette_fixed <- c(
  "#E69F00",  # <5 years     (orange)
  "#56B4E9",  # 5-9 years    (blue)
  "#009E73",  # 10-14 years  (green)
  "#F0E442",  # 15-19 years  (yellow)
  "#0072B2",  # 20-24 years  (dark blue)
  "#D55E00",  # 25-29 years  (vermillion)
  "#CC79A7",  # 30-34 years  (pink)
  "#800000",  # 35-39 years  (maroon) ✅ NEW
  "#6A3D9A",  # 40-44 years  (purple)
  "#1F78B4",  # 45-49 years  (steel blue)
  "#B2DF8A",  # 50-54 years  (light green)
  "#FB9A99",  # 55-59 years  (salmon pink) ✅ NEW
  "#33A02C",  # 60-64 years  (forest green)
  "#A6CEE3",  # 65-69 years  (light blue)
  "#FF7F00",  # 70-74 years  (orange)
  "#B15928",  # 75-79 years  (brown)
  "#CAB2D6",  # 80-84 years  (lavender)
  "#6B6BD6"   # 85-89 years  (vivid indigo) ✅ NEW
)

full_palette <- setNames(cb_palette_fixed, age_levels)

scenario_labels <- c(
  Scenario1 = "<20 years only",
  Scenario2 = "20-59 years only",
  Scenario3 = ">60 years only"
)

waffle_data18 <- combined_vacc_alloc_summ %>%
  # 1) 시나리오×연령대별 총합
  group_by(Scenario, age_gr) %>%
  summarise(Doses = sum(Vaccinated, na.rm = TRUE), .groups = "drop") %>%
  
  # 2) 시나리오별 100칸 보정
  group_by(Scenario) %>%
  mutate(
    raw_pct     = Doses / sum(Doses) * 100,
    pct_floor   = floor(raw_pct),
    remainder   = 100 - sum(pct_floor),
    ranker      = rank(desc(raw_pct - pct_floor), ties.method = "first"),
    pct         = pct_floor + if_else(ranker <= remainder, 1, 0),
    total_million = sum(Doses) / 1e6
  ) %>%
  ungroup() %>%
  
  # 3) 시나리오 라벨 붙이기
  mutate(
    Scenario_label = scenario_labels[Scenario],
    strip_label = paste0(Scenario_label, "\nTotal: ", round(total_million, 1), "M doses"),
    
    # 팩터 레벨 고정
    Scenario    = factor(Scenario, levels = names(scenario_labels)),
    strip_label = factor(strip_label, levels = paste0(
      scenario_labels,
      "\nTotal: ",
      round(
        c(
          sum(filter(., Scenario == "Scenario1")$Doses) / 1e6,
          sum(filter(., Scenario == "Scenario2")$Doses) / 1e6,
          sum(filter(., Scenario == "Scenario3")$Doses) / 1e6
        ),
        1
      ),
      "M doses"
    )),
    age_gr = factor(age_gr, levels = names(full_palette))
  )

p <- 
ggplot(waffle_data18, aes(values = pct, fill = age_gr)) +
  geom_waffle(
    n_rows = 10,
    size   = 0.25,
    colour = "white"
  ) +
  facet_wrap(~ strip_label, nrow = 1) +
  scale_fill_manual(
    values = full_palette,
    name   = "Age group"
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    #strip.text      = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm")
  )

ggsave(filename = "02_Outputs/2_1_Figures/fig_vacc_waffle.jpg", p, width = 12, height = 7, dpi = 1200)

###
selected_weeks <- c(1, 3, 6, 10, 15, 20)

waffle_weekly_data_subset <- waffle_weekly_data %>%
  filter(Week %in% selected_weeks) %>%
  mutate(
    week = factor(Week, levels = selected_weeks),
    Scenario_label = paste0(scenario_labels[Scenario], "\nWeek ", Week)
  )

ggplot(waffle_weekly_data_subset, aes(values = pct, fill = age_gr)) +
  geom_waffle(n_rows = 10, size = 0.25, colour = "white") +
  facet_wrap(~ Scenario_label, nrow = 3) +
  scale_fill_manual(values = full_palette, name = "Age group") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 11)
  )
