## -----------------------------------------------------------------------
## ceara_hotspot_coherence_check.R
##
## Before treating N_eff (population of the top effective-N hotspot
## municipalities) as a valid substitute for full state population in the
## state-aggregated SEIR fit, check whether that substitution is coherent:
##
##   1) Coverage: what share of the state's 2022 case burden do the hotspot
##      municipalities actually account for?
##   2) Spatial coherence: do the hotspot municipalities form one
##      contiguous cluster, or several geographically separate ones?
##   3) Temporal coherence: did the hotspot municipalities' outbreaks rise
##      and peak around the same weeks, or were they staggered in time?
##
## If cases are concentrated (1), geographically contiguous (2), and
## temporally synchronized (3), summing them into one state-level curve
## and fitting one SEIR with N_eff is a reasonable approximation. If not,
## that's the same "several overlapping local outbreaks fit as one curve"
## problem raised in peer review, just at a smaller spatial scale.
##
## Pilot state: Ceará.
## -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

## ---- 0. Load hotspot definitions from the previous pipeline --------------

load("00_Data/0_2_Processed/muni_hotspot_2022.RData")     # muni_2022_pop, concentration_summary
load("00_Data/0_2_Processed/n_eff_table_2022.RData")        # n_eff_table

ce_hotspot_row <- n_eff_table %>% filter(state_full == "Ceará")
ce_hotspot_names <- strsplit(ce_hotspot_row$top_muni_names, "; ")[[1]]

ce_all <- muni_2022_pop %>% filter(state_full == "Ceará", cases_2022 > 0)

ce_hotspot <- ce_all %>% filter(muni_name %in% ce_hotspot_names)

## ---- 1. Coverage: hotspot municipalities vs. full state -------------------

coverage_pct <- 100 * sum(ce_hotspot$cases_2022) / sum(ce_all$cases_2022)

cat("===== Ceará hotspot coverage =====\n")
cat("Hotspot municipalities:", paste(ce_hotspot_names, collapse = ", "), "\n")
cat("N municipalities:", nrow(ce_hotspot), "of", nrow(ce_all), "with any 2022 cases\n")
cat("Cases in hotspot set:", sum(ce_hotspot$cases_2022), "\n")
cat("Total Ceará 2022 confirmed cases:", sum(ce_all$cases_2022), "\n")
cat("Coverage:", round(coverage_pct, 1), "% of state case burden\n")
cat("Hotspot population:", format(round(sum(ce_hotspot$population)), big.mark = ","), "\n")
cat("Full Ceará population:", format(round(sum(ce_all$population)), big.mark = ","), "\n")

## Also report full coverage-vs-N curve (how coverage % grows as more
## municipalities are added, ranked by case count) for context.
ce_ranked <- ce_all %>%
  arrange(desc(cases_2022)) %>%
  mutate(
    rank            = row_number(),
    cum_cases       = cumsum(cases_2022),
    cum_coverage_pct = 100 * cum_cases / sum(cases_2022)
  )

cat("\nCumulative coverage by rank (top 15 municipalities):\n")
print(as.data.frame(ce_ranked %>% dplyr::select(rank, muni_name, cases_2022, cum_coverage_pct) %>% head(15)))

## ---- 2. Spatial coherence: are the hotspot municipalities contiguous? ----

adm2 <- st_read(
  "00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm2_ibge_2020.shp",
  quiet = TRUE
)
st_crs(adm2) <- 4674

adm2 <- adm2 %>%
  mutate(muni6 = as.character(as.numeric(substr(ADM2_PCODE, 3, 9)) %/% 10L))

ce_hotspot_sf <- adm2 %>% inner_join(ce_hotspot %>% dplyr::select(muni6, muni_name, cases_2022), by = "muni6")

stopifnot(nrow(ce_hotspot_sf) == nrow(ce_hotspot))

# Build adjacency (shares a border) among hotspot municipalities and find
# connected components via a simple breadth-first search (avoids adding an
# igraph dependency).
touch_list <- st_touches(ce_hotspot_sf)

n_muni <- nrow(ce_hotspot_sf)
component_id <- rep(0L, n_muni)
current_component <- 0L

for (i in seq_len(n_muni)) {
  if (component_id[i] != 0L) next
  current_component <- current_component + 1L
  queue <- i
  while (length(queue) > 0) {
    node <- queue[1]; queue <- queue[-1]
    if (component_id[node] != 0L) next
    component_id[node] <- current_component
    queue <- c(queue, touch_list[[node]])
  }
}

ce_hotspot_sf$component <- component_id

cat("\n===== Spatial coherence: connected components among hotspot municipalities =====\n")
print(
  as.data.frame(ce_hotspot_sf) %>%
    dplyr::select(muni_name, cases_2022, component) %>%
    arrange(component, desc(cases_2022))
)
cat("Number of spatially separate clusters:", length(unique(component_id)), "\n")

# Pairwise centroid distances (km) for context, regardless of clustering.
centroids <- st_centroid(st_geometry(ce_hotspot_sf))
dist_km <- st_distance(centroids) / 1000
rownames(dist_km) <- ce_hotspot_sf$muni_name
colnames(dist_km) <- ce_hotspot_sf$muni_name
cat("\nPairwise centroid distances (km):\n")
print(round(dist_km, 0))

## ---- 3. Temporal coherence: weekly timing of each hotspot municipality ----

muni_week <- readRDS("00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds")

ce_hotspot_weekly <- muni_week %>%
  filter(epi_year == 2022, muni6 %in% ce_hotspot$muni6) %>%
  left_join(ce_hotspot %>% dplyr::select(muni6, muni_name), by = "muni6")

rm(muni_week); gc()

peak_weeks <- ce_hotspot_weekly %>%
  group_by(muni_name) %>%
  slice_max(cases_confirmed, n = 1, with_ties = FALSE) %>%
  dplyr::select(muni_name, peak_week = epi_week, peak_cases = cases_confirmed)

cat("\n===== Peak week per hotspot municipality (2022) =====\n")
print(as.data.frame(peak_weeks %>% arrange(peak_week)))

save(
  ce_hotspot, ce_ranked, ce_hotspot_sf, coverage_pct, ce_hotspot_weekly, peak_weeks,
  file = "00_Data/0_2_Processed/ceara_hotspot_coherence.RData"
)

## ---- 4. Plot: weekly case curves for each hotspot municipality -----------
## Normalized to each municipality's own 2022 total, so the plot compares
## TIMING (when each place's outbreak rose/peaked/fell), not magnitude.

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
grid_col      <- "#e1e0d9"
axis_col      <- "#c3c2b7"
col_recent    <- "#2a78d6"

ce_hotspot_weekly_norm <- ce_hotspot_weekly %>%
  group_by(muni_name) %>%
  mutate(share_of_annual_total = cases_confirmed / sum(cases_confirmed)) %>%
  ungroup()

p_synchrony <- ggplot(
  ce_hotspot_weekly_norm,
  aes(x = epi_week, y = share_of_annual_total)
) +
  geom_line(color = col_recent, linewidth = 0.8) +
  facet_wrap(~muni_name, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 52, by = 13)) +
  labs(
    title    = "Ceará hotspot municipalities: timing of the 2022 outbreak",
    subtitle = "Weekly confirmed cases as a share of each municipality's own 2022 total (compares timing, not size)",
    x        = "Epidemiological week",
    y        = "Share of municipality's 2022 total"
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
ggsave(
  "02_Outputs/2_1_Figures/ceara_hotspot_synchrony.pdf",
  plot = p_synchrony, width = 8, height = 8, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/ceara_hotspot_synchrony.png",
  plot = p_synchrony, width = 8, height = 8, dpi = 400, bg = "white"
)

message("[save] 02_Outputs/2_1_Figures/ceara_hotspot_synchrony.pdf")
message("[save] 02_Outputs/2_1_Figures/ceara_hotspot_synchrony.png")
