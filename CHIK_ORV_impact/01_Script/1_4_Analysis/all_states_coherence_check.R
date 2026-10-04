## -----------------------------------------------------------------------
## all_states_coherence_check.R
##
## For all 11 states: decompose the municipalities covering ~90% of 2022
## case burden into spatially-contiguous clusters, compute each cluster's
## own weekly case curve and peak week, and classify whether state-level
## aggregate SEIR fitting is defensible or whether the state should be
## split into sub-regions for fitting -- and if so, how many.
##
## Decision rule (stated explicitly so it's auditable / adjustable):
##   - Build spatial clusters (border-touching connected components) among
##     municipalities covering the top ~90% of state 2022 case burden.
##   - Clusters with < 5% of state case burden are treated as minor and
##     folded into "residual" (not separately fit).
##   - Among clusters with >= 5% share ("significant" clusters), compute
##     the range of case-weighted peak weeks.
##   - If there is only 1 significant cluster, OR peak weeks across
##     significant clusters fall within a 6-week window: AGGREGATE OK.
##   - Otherwise: SPLIT, with recommended sub-regions = number of
##     significant clusters (after merging any that ARE within 6 weeks of
##     each other into one fitting group).
## -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

sf::sf_use_s2(FALSE)  # planar geometry engine; the IBGE adm2 shapefile has a
                       # few topologically invalid polygons that trip up s2

SHARE_THRESHOLD    <- 0.05   # min case share for a cluster to count as "significant"
SYNC_WEEK_TOLERANCE <- 6     # weeks; peak-week spread within this = "synchronized"
COVERAGE_TARGET     <- 0.90  # cumulative case coverage used to define the analysis set

target_states <- c(
  "Bahia" = 29, "Ceará" = 23, "Minas Gerais" = 31, "Pernambuco" = 26,
  "Paraíba" = 25, "Rio Grande do Norte" = 24, "Piauí" = 22, "Alagoas" = 27,
  "Tocantins" = 17, "Sergipe" = 28, "Goiás" = 52
)

## ---- 1. Load municipality case + population data (reuse prior pipeline) --

load("00_Data/0_2_Processed/muni_hotspot_2022.RData")  # muni_2022_pop

adm2 <- st_read(
  "00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm2_ibge_2020.shp",
  quiet = TRUE
)
st_crs(adm2) <- 4674
adm2 <- adm2 %>%
  mutate(muni6 = as.character(as.numeric(substr(ADM2_PCODE, 3, 9)) %/% 10L))

muni_week <- readRDS("00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds") %>%
  filter(epi_year == 2022) %>%
  dplyr::select(muni6, epi_week, cases_confirmed)

## ---- 2. Per-state: build coverage set, spatial clusters, cluster time series -

analyze_state <- function(state_name) {

  df <- muni_2022_pop %>%
    filter(state_full == state_name, cases_2022 > 0) %>%
    arrange(desc(cases_2022)) %>%
    mutate(cum_share = cumsum(cases_2022) / sum(cases_2022))

  # Coverage set: smallest top-N covering >= COVERAGE_TARGET of state cases.
  cutoff_row <- min(which(df$cum_share >= COVERAGE_TARGET))
  cov_set <- df %>% slice(1:cutoff_row)

  sf_set <- adm2 %>%
    inner_join(cov_set %>% dplyr::select(muni6, muni_name, cases_2022), by = "muni6")

  # Spatial clusters via border-touching connected components (BFS).
  touch_list <- st_touches(sf_set)
  n_muni <- nrow(sf_set)
  comp_id <- rep(0L, n_muni)
  cur <- 0L
  for (i in seq_len(n_muni)) {
    if (comp_id[i] != 0L) next
    cur <- cur + 1L
    queue <- i
    while (length(queue) > 0) {
      node <- queue[1]; queue <- queue[-1]
      if (comp_id[node] != 0L) next
      comp_id[node] <- cur
      queue <- c(queue, touch_list[[node]])
    }
  }
  sf_set$cluster <- comp_id

  cluster_share <- as.data.frame(sf_set) %>%
    group_by(cluster) %>%
    summarise(cluster_cases = sum(cases_2022), .groups = "drop") %>%
    mutate(cluster_share = cluster_cases / sum(df$cases_2022)) %>%
    arrange(desc(cluster_share))

  sig_clusters <- cluster_share %>% filter(cluster_share >= SHARE_THRESHOLD)

  # Weekly time series per significant cluster.
  muni6_by_cluster <- as.data.frame(sf_set) %>% dplyr::select(muni6, cluster)

  cluster_weekly <- muni_week %>%
    inner_join(muni6_by_cluster, by = "muni6") %>%
    filter(cluster %in% sig_clusters$cluster) %>%
    group_by(cluster, epi_week) %>%
    summarise(cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")

  peak_by_cluster <- cluster_weekly %>%
    group_by(cluster) %>%
    slice_max(cases, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(cluster, peak_week = epi_week)

  peak_week_range <- if (nrow(peak_by_cluster) > 1) {
    max(peak_by_cluster$peak_week) - min(peak_by_cluster$peak_week)
  } else 0L

  verdict <- if (nrow(sig_clusters) <= 1 || peak_week_range <= SYNC_WEEK_TOLERANCE) {
    "Aggregate OK"
  } else {
    "Split recommended"
  }

  n_subregions <- if (verdict == "Aggregate OK") 1L else nrow(sig_clusters)

  list(
    state_full            = state_name,
    n_muni_in_covset       = nrow(cov_set),
    coverage_pct           = 100 * cov_set$cum_share[nrow(cov_set)],
    n_spatial_clusters     = nrow(cluster_share),
    n_significant_clusters = nrow(sig_clusters),
    peak_week_range        = peak_week_range,
    verdict                = verdict,
    n_subregions            = n_subregions,
    cluster_weekly          = cluster_weekly %>% mutate(state_full = state_name),
    cluster_share           = sig_clusters %>% mutate(state_full = state_name),
    sf_set                  = sf_set %>% mutate(state_full = state_name)
  )
}

results <- lapply(names(target_states), analyze_state)
names(results) <- names(target_states)

## ---- 3. Summary table ------------------------------------------------------

summary_table <- bind_rows(lapply(results, function(r) {
  tibble::tibble(
    state_full              = r$state_full,
    n_municipalities_90pct  = r$n_muni_in_covset,
    coverage_pct             = round(r$coverage_pct, 1),
    n_spatial_clusters       = r$n_spatial_clusters,
    n_significant_clusters   = r$n_significant_clusters,
    peak_week_range          = r$peak_week_range,
    verdict                  = r$verdict,
    recommended_subregions   = r$n_subregions
  )
})) %>%
  arrange(desc(peak_week_range))

cat("===== State-level aggregation feasibility =====\n")
print(as.data.frame(summary_table))

dir.create("00_Data/0_2_Processed", showWarnings = FALSE, recursive = TRUE)
save(results, summary_table, file = "00_Data/0_2_Processed/all_states_coherence.RData")
readr::write_csv(summary_table, "00_Data/0_2_Processed/all_states_coherence_summary.csv")

## ---- 4. Figure A: normalized weekly curves per cluster, all 11 states -----

all_cluster_weekly <- bind_rows(lapply(results, function(r) r$cluster_weekly)) %>%
  left_join(bind_rows(lapply(results, function(r) r$cluster_share)),
            by = c("state_full", "cluster")) %>%
  group_by(state_full, cluster) %>%
  mutate(share_of_cluster_total = cases / sum(cases)) %>%
  ungroup() %>%
  mutate(state_full = factor(state_full, levels = summary_table$state_full)) %>%
  # Rank clusters 1..K within each state (by share, descending) so the
  # color scale only ever needs as many colors as the largest per-state
  # cluster count, regardless of the raw (state-specific) cluster IDs.
  group_by(state_full) %>%
  mutate(cluster_rank = dense_rank(desc(cluster_share))) %>%
  ungroup() %>%
  mutate(cluster_lab = paste0("cluster ", cluster, " (", scales::percent(cluster_share, accuracy = 1), ")"))

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
grid_col      <- "#e1e0d9"
axis_col      <- "#c3c2b7"

palette_clusters <- c("#2a78d6", "#eb6834", "#3f9f5a", "#a24fb0", "#c9a227")

p_all_synchrony <- ggplot(
  all_cluster_weekly,
  aes(x = epi_week, y = share_of_cluster_total, color = factor(cluster_rank), group = cluster)
) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~state_full, ncol = 3, scales = "free_y") +
  scale_color_manual(values = palette_clusters, guide = "none") +
  scale_x_continuous(breaks = seq(0, 52, by = 13)) +
  labs(
    title    = "Timing of major case clusters within each state, 2022",
    subtitle = sprintf(
      "Each line = one spatially-contiguous cluster covering >= %.0f%% of state cases (weekly share of that cluster's own 2022 total)",
      100 * SHARE_THRESHOLD
    ),
    x = "Epidemiological week", y = "Share of cluster's 2022 total"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(color = ink_primary),
    plot.title          = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle       = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    strip.text          = element_text(face = "bold", size = 9.5, color = ink_primary),
    axis.line.x         = element_line(color = axis_col, linewidth = 0.3),
    panel.grid.major    = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor    = element_blank(),
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(12, 16, 10, 12)
  )

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("02_Outputs/2_1_Figures/all_states_cluster_synchrony.pdf",
       plot = p_all_synchrony, width = 12, height = 12, device = cairo_pdf)
ggsave("02_Outputs/2_1_Figures/all_states_cluster_synchrony.png",
       plot = p_all_synchrony, width = 12, height = 12, dpi = 300, bg = "white")

## ---- 5. Figure B: summary dumbbell of peak-week range per state -----------

peak_range_plot_data <- bind_rows(lapply(results, function(r) {
  wk <- r$cluster_weekly %>%
    group_by(cluster) %>%
    slice_max(cases, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(r$cluster_share, by = c("state_full", "cluster"))
  wk
})) %>%
  mutate(state_full = factor(state_full, levels = rev(summary_table$state_full)))

p_peak_range <- ggplot(peak_range_plot_data, aes(x = epi_week, y = state_full)) +
  geom_line(aes(group = state_full), color = grid_col, linewidth = 3) +
  geom_point(aes(size = cluster_share), color = "#8a1f11", alpha = 0.85) +
  scale_size_area(max_size = 7, labels = scales::percent, name = "Cluster share\nof state cases") +
  labs(
    title = "Peak week of each significant case cluster, by state (2022)",
    subtitle = "Wide spread = clusters likely reflect separate, staggered local outbreaks",
    x = "Epidemiological week of cluster peak", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(color = ink_primary),
    plot.title          = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle       = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    panel.grid.major.x  = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(12, 16, 10, 12)
  )

ggsave("02_Outputs/2_1_Figures/all_states_peak_week_range.pdf",
       plot = p_peak_range, width = 8, height = 6, device = cairo_pdf)
ggsave("02_Outputs/2_1_Figures/all_states_peak_week_range.png",
       plot = p_peak_range, width = 8, height = 6, dpi = 400, bg = "white")

message("[save] 00_Data/0_2_Processed/all_states_coherence_summary.csv")
message("[save] 02_Outputs/2_1_Figures/all_states_cluster_synchrony.pdf/.png")
message("[save] 02_Outputs/2_1_Figures/all_states_peak_week_range.pdf/.png")
