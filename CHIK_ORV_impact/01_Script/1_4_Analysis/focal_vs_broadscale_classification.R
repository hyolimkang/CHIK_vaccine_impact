## -----------------------------------------------------------------------
## focal_vs_broadscale_classification.R
##
## Reclassify the 25 analysis units (from transmission_tier_classification.R)
## into two epidemiologically distinct archetypes, justified by a real,
## tested relationship rather than an arbitrary split of attack rate:
##
##   Focal      = small population, few constituent municipalities, near-
##                complete local burnout (very high attack rate)
##   Broad-scale = large population spanning many municipalities, modest
##                attack rate (diluted by areas the outbreak didn't reach)
##
## Evidence: Spearman correlation between log(population) and estimated
## infection attack rate across the 25 units = -0.57 (p = 0.0035) -- see
## check_pombal_and_popsize.R. This script formalizes that into a
## population-threshold classification, reports it as a table, and maps
## which municipalities fall into each category (so the "why focal" is
## visually obvious, not just a number).
## -----------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(ggpubr)

sf::sf_use_s2(FALSE)

load("00_Data/0_2_Processed/transmission_tier_units.RData")  # units_all
load("00_Data/0_2_Processed/all_states_coherence.RData")     # results, summary_table
load("00_Data/0_2_Processed/muni_hotspot_2022.RData")        # muni_2022_pop

POP_THRESHOLD <- 200000  # sits in the visual gap between the small-town
                          # cluster (<130k) and the next unit up (~285k)
AR_THRESHOLD  <- 3.63    # the Jenks natural-break point on infection AR (%)
                          # found in transmission_tier_classification.R

## ---- 1. Count constituent municipalities per unit -------------------------

aggregate_ok_states <- summary_table %>% filter(verdict == "Aggregate OK") %>% pull(state_full)
split_states        <- summary_table %>% filter(verdict == "Split recommended") %>% pull(state_full)

muni_counts <- bind_rows(
  # Aggregate-OK states: n = municipalities across ALL significant clusters
  bind_rows(lapply(aggregate_ok_states, function(st) {
    r <- results[[st]]
    n <- sf::st_drop_geometry(r$sf_set) %>% filter(cluster %in% r$cluster_share$cluster) %>% nrow()
    tibble::tibble(unit_id = st, n_municipalities = n)
  })),
  # Split states: n per significant cluster
  bind_rows(lapply(split_states, function(st) {
    r <- results[[st]]
    sf_df <- sf::st_drop_geometry(r$sf_set)
    sf_df %>%
      filter(cluster %in% r$cluster_share$cluster) %>%
      count(cluster, name = "n_municipalities") %>%
      mutate(unit_id = paste0(st, " - cluster ", cluster)) %>%
      dplyr::select(unit_id, n_municipalities)
  }))
)

## ---- 2. Classify and build the summary table -------------------------------

## "Focal" requires BOTH small population AND high attack rate -- either
## criterion alone misclassifies a real case: Bahia: Bom Jesus da Lapa is
## small (141k) but low-AR (2.1%, an ordinary modest outbreak in a small
## town, not an explosive burnout); Ceará: Brejo Santo (the Cariri
## cluster) is large (925k, 17 municipalities) but high-AR (9.2%, a
## genuinely intense REGIONAL outbreak, not a single-town phenomenon).
## Requiring both conditions keeps the label tied to the actual mechanism
## (near-complete local burnout in a small, tightly-mixed population)
## rather than either variable alone.
units_classified <- units_all %>%
  left_join(muni_counts, by = "unit_id") %>%
  mutate(
    setting_type = ifelse(
      population < POP_THRESHOLD & infection_ar_pct > AR_THRESHOLD,
      "Focal", "Broad-scale"
    ),
    setting_type = factor(setting_type, levels = c("Focal", "Broad-scale"))
  ) %>%
  arrange(population)

cat("===== Focal vs. broad-scale classification =====\n")
print(
  as.data.frame(
    units_classified %>%
      dplyr::select(unit_label, setting_type, n_municipalities, population,
             cases_2022, attack_rate_pct, infection_ar_pct)
  ),
  digits = 4
)

cat(sprintf("\nPopulation threshold: %s\n", format(POP_THRESHOLD, big.mark = ",")))
cat("\nFocal units:", sum(units_classified$setting_type == "Focal"),
    " | Broad-scale units:", sum(units_classified$setting_type == "Broad-scale"), "\n")

cat("\nException check -- any Broad-scale unit with high infection AR (>8%)?\n")
print(
  as.data.frame(
    units_classified %>% filter(setting_type == "Broad-scale", infection_ar_pct > 8) %>%
      dplyr::select(unit_label, n_municipalities, population, infection_ar_pct)
  )
)

dir.create("00_Data/0_2_Processed", showWarnings = FALSE, recursive = TRUE)
save(units_classified, POP_THRESHOLD, file = "00_Data/0_2_Processed/focal_broadscale_units.RData")
readr::write_csv(units_classified, "00_Data/0_2_Processed/focal_broadscale_units.csv")

## ---- 3. Scatter: population vs infection AR, colored by classification ----

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
grid_col      <- "#e1e0d9"
col_focal     <- "#a3271f"
col_broad     <- "#2a78d6"

p_scatter <- ggplot(units_classified, aes(x = population, y = infection_ar_pct, color = setting_type)) +
  geom_vline(xintercept = POP_THRESHOLD, linetype = "22", color = grid_col, linewidth = 0.6) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_number(suffix = "%")) +
  scale_color_manual(values = c(Focal = col_focal, "Broad-scale" = col_broad), name = NULL) +
  labs(
    title    = "Focal vs. broad-scale transmission settings",
    subtitle = sprintf("Population threshold = %s (dashed line); Spearman r = -0.57, p = 0.0035", format(POP_THRESHOLD, big.mark = ",")),
    x = "Population (log scale)", y = "Estimated infection attack rate (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(color = ink_primary),
    plot.title          = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle       = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    panel.grid.minor    = element_blank(),
    legend.position     = "top",
    legend.justification = "left",
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(12, 16, 10, 12)
  )

ggsave("02_Outputs/2_1_Figures/focal_broadscale_scatter.png", plot = p_scatter, width = 8, height = 6.5, dpi = 400, bg = "white")
ggsave("02_Outputs/2_1_Figures/focal_broadscale_scatter.pdf", plot = p_scatter, width = 8, height = 6.5, device = cairo_pdf)

## ---- 4. Map: ONE whole-Brazil map, clusters dissolved into single shapes --
## Previous version faceted into 11 tiny per-state panels with every
## municipality border still visible -- gave no sense of where states sit
## relative to each other, and made it look like each municipality was
## being treated as its own separate thing. Fixed here: (a) one map, full
## Brazil extent, with ALL state borders as background context; (b) each
## unit's municipalities are DISSOLVED (st_union) into one shape, so a
## 17-municipality cluster reads as one region, matching how it's
## actually treated in the fitting.

adm1 <- st_read(
  "00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm1_ibge_2020.shp",
  quiet = TRUE
)
st_crs(adm1) <- 4674

adm2 <- st_read(
  "00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm2_ibge_2020.shp",
  quiet = TRUE
)
st_crs(adm2) <- 4674
adm2 <- adm2 %>% mutate(muni6 = as.character(as.numeric(substr(ADM2_PCODE, 3, 9)) %/% 10L))

all_states <- names(results)

muni_setting <- bind_rows(lapply(all_states, function(st) {

  r <- results[[st]]
  sf_df <- sf::st_drop_geometry(r$sf_set) %>%
    filter(cluster %in% r$cluster_share$cluster)

  sf_df <- if (st %in% aggregate_ok_states) {
    sf_df %>% mutate(unit_id = st)
  } else {
    sf_df %>% mutate(unit_id = paste0(st, " - cluster ", cluster))
  }

  sf_df %>%
    left_join(units_classified %>% dplyr::select(unit_id, setting_type), by = "unit_id") %>%
    dplyr::select(muni6, state_full, unit_id, setting_type)
}))

map_data_muni <- adm2 %>% inner_join(muni_setting, by = "muni6")

# Dissolve to one polygon per unit (cluster), dropping internal municipality
# borders.
map_data_dissolved <- map_data_muni %>%
  group_by(unit_id, state_full, setting_type) %>%
  summarise(.groups = "drop")

cat("\nUnits dissolved to single shapes:", nrow(map_data_dissolved), "\n")

# Label points: centroid of each state's LARGEST (by population) colored
# unit, not the state's own polygon centroid -- a state-centroid label can
# land far from every colored shape (Minas Gerais) or sit on top of one
# (Tocantins). Anchoring to the biggest actual cluster keeps the label
# near real content.
target_state_names <- names(results)
state_labels <- map_data_dissolved %>%
  left_join(units_classified %>% dplyr::select(unit_id, population), by = "unit_id") %>%
  group_by(state_full) %>%
  slice_max(population, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label_geom = st_centroid(geometry)) %>%
  st_drop_geometry() %>%
  mutate(geometry = label_geom) %>%
  st_as_sf()

p_map <- ggplot() +
  # Background: every Brazilian state, for geographic context.
  geom_sf(data = adm1, fill = "#f4f1ea", color = "#c3c2b7", linewidth = 0.2) +
  # Outline the 11 target states a bit more visibly.
  geom_sf(
    data = adm1 %>% filter(ADM1_PT %in% target_state_names),
    fill = NA, color = "#898781", linewidth = 0.5
  ) +
  # The actual classified units, dissolved to one shape each.
  geom_sf(data = map_data_dissolved, aes(fill = setting_type), color = "white", linewidth = 0.15) +
  geom_sf_text(
    data = state_labels, aes(label = state_full),
    size = 2.9, color = ink_secondary, fontface = "bold",
    nudge_y = 0.35
  ) +
  scale_fill_manual(values = c(Focal = col_focal, "Broad-scale" = col_broad), name = NULL, na.translate = FALSE) +
  coord_sf(xlim = c(-50, -34), ylim = c(-21, 0), expand = FALSE) +
  labs(
    title    = "Focal vs. broad-scale transmission settings",
    subtitle = stringr::str_wrap(
      "Each colored shape = the population/case footprint used for THAT unit's attack rate (a cluster of 1-27 municipalities, NOT the full state).",
      width = 78
    ),
    caption  = stringr::str_wrap(
      sprintf(
        paste(
          "25 colored units total (5 Focal, 20 Broad-scale). Beige/uncolored = zero 2022 cases OR part of a cluster below the 5%% state-case-share cutoff (excluded from this analysis).",
          "Focal: population < %s AND estimated infection attack rate > %.1f%%. Ceará's Cariri cluster (Crato/Juazeiro do Norte area) is Broad-scale despite a high attack rate -- see scatter plot.",
          "Population here is CLUSTER-restricted (where the 2022 outbreak actually was), not full state population -- a different denominator than the state-level SEIR fits used elsewhere in this analysis.",
          sep = " "
        ),
        format(POP_THRESHOLD, big.mark = ","), AR_THRESHOLD
      ),
      width = 105
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                   = element_text(color = ink_primary),
    plot.title             = element_text(face = "bold", size = 14, color = ink_primary),
    plot.subtitle          = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 8)),
    plot.caption           = element_text(size = 7.5, color = ink_muted, hjust = 0, margin = ggplot2::margin(t = 8)),
    axis.text              = element_blank(),
    axis.title             = element_blank(),
    panel.grid             = element_blank(),
    legend.position        = "top",
    legend.justification   = "left",
    plot.title.position    = "plot",
    plot.caption.position  = "plot",
    plot.margin            = ggplot2::margin(12, 16, 10, 12)
  )

ggsave("02_Outputs/2_1_Figures/focal_broadscale_map.pdf", plot = p_map, width = 9, height = 11.6, device = cairo_pdf)
ggsave("02_Outputs/2_1_Figures/focal_broadscale_map.png", plot = p_map, width = 9, height = 11.6, dpi = 400, bg = "white")

message("[save] 00_Data/0_2_Processed/focal_broadscale_units.csv")
message("[save] 02_Outputs/2_1_Figures/focal_broadscale_scatter.png/.pdf")
message("[save] 02_Outputs/2_1_Figures/focal_broadscale_map.png/.pdf")
