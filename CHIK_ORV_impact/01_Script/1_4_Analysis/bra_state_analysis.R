bra_region <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                               sheet = "bra_region_2022")
bra_region_22 <- bra_region %>% filter(year == 2022)
bra_region_23 <- bra_region %>% filter(year == 2023)
bra_region_24 <- bra_region %>% filter(year == 2024)

colnames(bra_region_22)[2] <- "case_22"
colnames(bra_region_22)[4] <- "case_per_100k_22"
colnames(bra_region_23)[2] <- "case_23"
colnames(bra_region_23)[4] <- "case_per_100k_23"
colnames(bra_region_24)[2] <- "case_24"
colnames(bra_region_24)[4] <- "case_per_100k_24"

bra_foi <- allfoi %>% filter(country == "Brazil")
# state
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))


bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)


bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = quantile(foi_mid[foi_mid > 0], 0.5, na.rm = TRUE),
            foi_lo  = quantile(foi_mid[foi_mid > 0], 0.025, na.rm = TRUE),
            foi_hi  = quantile(foi_mid[foi_mid > 0], 0.975, na.rm = TRUE)
            )%>%
  filter(!is.na(NAME_1))

save("bra_foi_state_summ", file = "00_Data/0_2_Processed/bra_foi_state_summ.RData")
save("bra_foi_states", file = "00_Data/0_2_Processed/bra_foi_states.RData")

names(bra_foi_state_summ) <- c("region", "avg_foi", "foi_lo", "foi_hi", "geometry")

bra_foi_state_summ <- merge(bra_foi_state_summ, bra_region, by = "region")

br_states <- dplyr::rename(br_states, region = NAME_1)

library(dplyr)

br_states <- br_states %>%
  left_join(bra_region_22, by = "region") %>%
  left_join(bra_region_23, by = "region") %>%
  left_join(bra_region_24, by = "region")

bra_foi_state_summ_df <- st_set_geometry(bra_foi_state_summ, NULL)
br_states <- br_states %>% 
  left_join(bra_foi_state_summ_df, by = "region")

# number of cases in 2022
p <- 
ggplot(br_states) +
  geom_sf(aes(fill = case_per_100k_22)) +
  scale_fill_distiller(
    palette = "Reds", 
    direction = 1, 
    name = "per 100k cases"
    ) +
  theme_minimal()+
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text labels (lat/lon)
    axis.ticks = element_blank(),        # Remove axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )  +
  geom_sf_text(aes(label = region), 
               size = 2,     # Text size
               color = "black")
ggsave(filename = "02_Outputs/2_1_Figures/bra_ob_24_per100k.jpg", plot = p,
       width = 8, height = 3)

p <- 
ggplot(br_states) +
  geom_sf(aes(fill = case_22)) +
  scale_fill_distiller(
    palette = "Reds", 
    direction = 1, 
    name = "reported cases"
  ) +
  theme_minimal()+
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text labels (lat/lon)
    axis.ticks = element_blank(),        # Remove axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )  +
  geom_sf_text(aes(label = region), 
               size = 2,     # Text size
               color = "black")

ggsave(filename = "02_Outputs/2_1_Figures/bra_ob_23.jpg", plot = p,
       width = 8, height = 3)

## foi
ggplot(br_states) +
  geom_sf(aes(fill = avg_foi)) +
  scale_fill_distiller(
    palette = "Reds", 
    direction = 1, 
    name = "reported cases"
  ) +
  theme_minimal()+
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text labels (lat/lon)
    axis.ticks = element_blank(),        # Remove axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )  +
  geom_sf_text(aes(label = region), 
               size = 2,     # Text size
               color = "black")




## epidemic diversity by subnational region
bra_region_week <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                         sheet = "bra_region_week_2022")

bra_region_week <- left_join(bra_region_week, bra_pop_all, by = "region")
bra_region_week <- bra_region_week %>% mutate(
  cases_per_million = cases / tot_pop * 1000000
)

p <- 
ggplot(bra_region_week)+
  geom_line(aes(x = week, y = cases))+
  facet_wrap(~region, scales = "free_y")+
  theme_minimal()

ggplot(bra_region_week)+
  geom_line(aes(x = week, y = cases_per_million))+
  facet_wrap(~region, scales = "free_y")+
  theme_light()


ggsave(filename = "02_Outputs/2_1_Figures/fig_allstate_raw.jpg", p, width = 10, height = 4, dpi = 1200)


threshold <- 20

metrics_df <- bra_region_week %>%
  group_by(region) %>%
  summarize(
    total_cases    = sum(cases, na.rm = TRUE),
    peak_incidence = max(cases, na.rm = TRUE),
    week_of_peak   = week[which.max(cases)],            # first week of max
    outbreak_duration = sum(cases > threshold, na.rm = TRUE)  # # of wks > threshold
  ) %>%
  ungroup()

metrics_df <- left_join(metrics_df, bra_pop_all, by = "region")
metrics_df <- metrics_df %>% mutate(
  cases_per_100k = total_cases / tot_pop * 100000
)

bra_foi_state_summ_22 <- bra_foi_state_summ %>% filter(year == 2022)

metrics_df <- left_join(metrics_df, bra_foi_state_summ_22[,c(1:2)], by = "region")
metrics_df <- metrics_df %>% mutate(
  s_pop_med   = exp(-avg_foi * 50),
  s_pop       = s_pop_med * tot_pop,
  attack_rate = total_cases / s_pop 
)

metrics_df <- metrics_df %>% mutate(
  rel_cases = total_cases / tot_pop,
  rel_cases_norm = rel_cases / max(rel_cases),
  abs_cases_norm = total_cases / max(total_cases),
  comp_score = 0.5 * rel_cases_norm + 0.5 * abs_cases_norm
) %>% mutate(
  case_per_million = (peak_incidence / tot_pop) * 1000000
) 

metrics_df <- metrics_df %>%
  mutate(
    outbreak_burden = case_when(
      total_cases   > 25000  |
        peak_incidence > 1500  |
        case_per_million > 500  ~ "High",
      
      total_cases   >= 10000 & 
        (peak_incidence >= 500 | case_per_million >= 200) ~ "Moderate",
      
      total_cases   < 10000 &
        (peak_incidence >=  500 &
        case_per_million >= 100)       ~ "Low",
      
      total_cases < 10000 &
        peak_incidence < 500 &
        case_per_million < 150 ~ "Very low"
    )
  ) 

library(forcats)

bra_region_week_threshold <- bra_region_week %>% filter(cases_per_million > 100)

p1 <- 
metrics_df %>%
  mutate(region = fct_reorder(region, attack_rate, .desc = TRUE)) %>%  # Reorder by total_cases in descending order
  ggplot() +
  geom_bar(aes(x = region, y = attack_rate, fill = region), stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Attack rate")+
  scale_y_continuous(labels = comma)

p2 <- 
  metrics_df %>%
  mutate(region = fct_reorder(region, total_cases, .desc = TRUE)) %>%  # Reorder by total_cases in descending order
  ggplot() +
  geom_bar(aes(x = region, y = total_cases, fill = region), stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Total cumulative reported cases")+
  scale_y_continuous(labels = comma)

region_order <- metrics_df %>%
  arrange(desc(total_cases)) %>%
  pull(region)

combined_metric <- metrics_df %>%
  select(region, total_cases, attack_rate) %>%
  pivot_longer(cols = c(attack_rate, total_cases),
               names_to = "metric",
               values_to = "value") %>%
  mutate(region = factor(region, levels = region_order),
         metric = factor(metric, levels = c("total_cases", "attack_rate")))


p <- combined_metric %>%
  ggplot(aes(
    x = region,  # ← no reordering here!
    y = value,
    fill = region
  )) +
  geom_bar(stat = "identity") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left",
             labeller = labeller(metric = c(
               total_cases = "Total reported cases",
               attack_rate = "Attack rate (%)"
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


tot_10_comp <- metrics_df %>% filter(comp_score > 0.15)
tot_5_attack <- metrics_df %>% filter(attack_rate > 0.0068)
tot_10_totcase <- metrics_df %>% filter(total_cases > 9000)

ggarrange(p1, p2, 
          nrow = 1, ncol = 2)

ggplot(metrics_df)+
  geom_histogram(aes(x = week_of_peak))

cluster_data <- metrics_df %>%
  select(-region)

cluster_data_scaled <- scale(cluster_data)
rownames(cluster_data_scaled) <- metrics_df$region

## looking at the relation between indexP and epidemic pattern
bra_admin <- do.call(rbind, lapply(bra_combined_output, function(x){
   x$admin[1]
  }))

colnames(bra_admin) <- "state"
 
state_mapping <- data.frame(
  state = c("BRRS", "BRRR", "BRPA", "BRAC", "BRAP", "BRMS", "BRPR", "BRSC", 
            "BRAM", "BRRO", "BRMT", "BRMA", "BRPI", "BRCE", "BRRN", "BRPB", 
            "BRPE", "BRAL", "BRSE", "BRBA", "BRES", "BRRJ", "BRSP", "BRGO", 
            "BRDF", "BRMG", "BRTO"),
  full_name = c("Rio Grande do Sul", "Roraima", "Pará", "Acre", "Amapá", 
                "Mato Grosso do Sul", "Paraná", "Santa Catarina", "Amazonas", 
                "Rondônia", "Mato Grosso", "Maranhão", "Piauí", "Ceará", 
                "Rio Grande do Norte", "Paraíba", "Pernambuco", "Alagoas", 
                "Sergipe", "Bahia", "Espírito Santo", "Rio de Janeiro", 
                "São Paulo", "Goiás", "Distrito Federal", "Minas Gerais", 
                "Tocantins"),
  stringsAsFactors = FALSE
)

colnames(indexP_matrix) <- state_mapping$full_name
indexP_matrix <- as.data.frame(indexP_matrix)

indexP_matrix <- indexP_matrix %>% 
  mutate(week = row_number())

indexP_matrix_long <- indexP_matrix %>%
  pivot_longer(
    cols = -week,
    names_to = "state",
    values_to = "indexP"
  )

indexP_matrix_long <- arrange(indexP_matrix_long, state)

ggplot(indexP_matrix_long, aes(x = week, y = indexP)) +
  geom_line() +
  facet_wrap(~ state) +
  labs(title = "Weekly indexP",
       x = "Week",
       y = "Cases")+
  theme_light()

## correlation analysis
corr <- cor(bra_region_week$cases[1:52],
            indexP_matrix_long$indexP[1:52])

colnames(bra_region_week)[3] <- "state"
combined_df <- merge(bra_region_week, indexP_matrix_long, by = c("state", "week"))

plot_df <- combined_df %>%
  pivot_longer(cols = c(cases, indexP), names_to = "variable", values_to = "value")

ggplot(plot_df, aes(x = week, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ state, scales = "free_y") +  # scales free_y allows each state to adjust to its range
  labs(x = "Week",
       y = "Value",
       color = "Metric")+
  scale_y_log10()+
  theme_minimal()


correlations <- combined_df %>%
  group_by(state) %>%
  summarise(correlation = cor(cases, indexP, use = "complete.obs"))

metrics_df$correlations <- correlations$correlation

ggplot(correlations, aes(x = reorder(state, correlation), y = correlation)) +
  geom_col(fill = "skyblue") +
  coord_flip()

## classification
# what should be the threshold for high outbreak? 
threshold_cases <- median(metrics_df$total_cases) 
threshold_corr <- 0.5  

metrics_df <- metrics_df %>%
  mutate(
    group = case_when(
      total_cases >= threshold_cases & abs(correlations) >= threshold_corr ~ "A",  # high outbreak, strong correlation
      total_cases >= threshold_cases & abs(correlations) <  threshold_corr ~ "B",  # high outbreak, weak correlation
      total_cases <  threshold_cases & abs(correlations)  >= threshold_corr ~ "C",  # small outbreak, strong corrleation                                  ~ "C"   # low outbreak
      total_cases <  threshold_cases & abs(correlations) <  threshold_corr ~ "D"   # small outbreak, weak correlation
    )
  )
colnames(metrics_df)[1] <- "state"
plot_df <- left_join(plot_df, metrics_df %>% select(state, group), by = "state")

colnames(bra_region_week)[3] <- "state"
bra_region_week <- bra_region_week %>% 
  left_join(metrics_df %>% select(state, peak_incidence, total_cases), by = "state")%>%
mutate(state = reorder(state, -peak_incidence))

bra_region_week <- bra_region_week  %>%
  mutate(
    group = case_when(
      total_cases >= quantile(total_cases, probs = 2/3, na.rm = TRUE) ~ "High",
      total_cases >= quantile(total_cases, probs = 1/3, na.rm = TRUE) & 
      total_cases < quantile(total_cases, probs = 2/3, na.rm = TRUE) ~ "Medium",
      TRUE ~ "Low"
    )
  )

bra_region_week$group <- factor(bra_region_week$group, levels = c("High", "Medium", "Low"))

ggplot(bra_region_week, aes(x = week, y = cases, color = group)) +
  geom_line() +
  facet_wrap(~ group + state, scales = "free_y") +
  labs(x = "Week",
       y = "Value",
       color = "Metric",
       caption = "Group A: Large outbreak, strong correlation\nGroup B: Large outbreak, weak correlation\nGroup C: Small outbreak, strong correlation\nGroup D: Small outbreak, weak correlation") +
  theme_minimal()

bra_mg <- bra_region_week %>% filter(state == "Minas Gerais")

ggplot(bra_mg)+
  geom_line(aes(x = week, y = cases))

