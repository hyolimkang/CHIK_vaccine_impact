## -----------------------------------------------------------------------
## muni_hotspot_map_2022.R
##
## Municipality-level spatial concentration of the 2022 chikungunya outbreak
## within each of the 11 states used for age-structured fitting.
##
## Purpose: visualize where within each state the 2022 outbreak was
## concentrated, and quantify that concentration (effective number of
## municipalities, via the Simpson/participation-ratio index) as a
## non-arbitrary basis for defining an "effective population" (N_eff) --
## an alternative to the full state population currently used as the
## SEIR mixing denominator / susceptible pool base.
##
## Inputs:
##   - 00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds
##     municipality x week case counts, 2015-2024
##   - 00_Data/0_1_Raw/estimativa_dou_2024.xls (sheet "MUNICÍPIOS")
##     IBGE municipality population estimates (2024, closest available)
##   - 00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm2_ibge_2020.shp
##     official IBGE municipality boundaries
## -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readxl)
library(sf)
library(ggplot2)
library(ggpubr)

## ---- 0. Config ---------------------------------------------------------

target_states <- c(
  "Bahia" = 29, "Ceará" = 23, "Minas Gerais" = 31, "Pernambuco" = 26,
  "Paraíba" = 25, "Rio Grande do Norte" = 24, "Piauí" = 22, "Alagoas" = 27,
  "Tocantins" = 17, "Sergipe" = 28, "Goiás" = 52
)

## ---- 1. Municipality-level 2022 case counts -----------------------------

muni_week <- readRDS("00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds")

muni_2022 <- muni_week %>%
  filter(epi_year == 2022) %>%
  mutate(uf_code = as.integer(substr(muni6, 1, 2))) %>%
  filter(uf_code %in% target_states) %>%
  group_by(muni6, uf_code) %>%
  summarise(cases_2022 = sum(cases_confirmed, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    state_full = names(target_states)[match(uf_code, target_states)],
    # SINAN convention: muniXX0000 marks "municipality unknown/ignored"
    is_unknown_muni = substr(muni6, 3, 6) == "0000"
  )

rm(muni_week); gc()

cat("Share of state 2022 cases with unknown/unassigned municipality:\n")
print(
  muni_2022 %>%
    group_by(state_full) %>%
    summarise(
      total_cases   = sum(cases_2022),
      unknown_cases = sum(cases_2022[is_unknown_muni]),
      pct_unknown   = 100 * unknown_cases / total_cases
    )
)

## ---- 2. Municipality population (IBGE 2024 estimate) --------------------

pop_raw <- readxl::read_excel(
  "00_Data/0_1_Raw/estimativa_dou_2024.xls",
  sheet = "MUNICÍPIOS", skip = 1
)

colnames(pop_raw) <- c("uf_abbr", "cod_uf", "cod_munic", "muni_name", "population")

pop_muni <- pop_raw %>%
  filter(!is.na(cod_uf), !is.na(cod_munic)) %>%
  mutate(
    cod_uf     = as.integer(cod_uf),
    cod_munic  = as.integer(cod_munic),
    population = suppressWarnings(as.numeric(gsub("[^0-9]", "", population))),
    muni7      = cod_uf * 100000L + cod_munic,
    muni6      = as.character(muni7 %/% 10L)
  ) %>%
  filter(!is.na(population)) %>%
  dplyr::select(muni6, muni_name, population)

cat("\nMunicipality population rows parsed:", nrow(pop_muni), "\n")
cat("Total Brazil population (sanity check, should be ~212M):",
    format(sum(pop_muni$population), big.mark = ","), "\n")

## ---- 3. Join cases + population, compute attack rate ---------------------

muni_2022_pop <- muni_2022 %>%
  filter(!is_unknown_muni) %>%
  left_join(pop_muni, by = "muni6") %>%
  filter(!is.na(population), population > 0) %>%
  mutate(attack_rate_per100k = 1e5 * cases_2022 / population)

cat("\nMunicipalities matched to population (of non-unknown case rows):",
    nrow(muni_2022_pop), "\n")
cat("Unmatched (population join failed):",
    sum(!is.na(muni_2022$muni6) & !muni_2022$is_unknown_muni) - nrow(muni_2022_pop), "\n")

## ---- 4. Concentration index per state ------------------------------------
## Effective number of municipalities (Simpson reciprocal / participation
## ratio): N_eff_muni = (sum cases)^2 / sum(cases^2). Equals the number of
## municipalities if cases were spread perfectly evenly; approaches 1 if
## concentrated in a single municipality. No arbitrary threshold required.

concentration_summary <- muni_2022_pop %>%
  filter(cases_2022 > 0) %>%
  group_by(state_full) %>%
  summarise(
    n_municipalities_with_cases = n(),
    total_cases                 = sum(cases_2022),
    n_eff_municipalities        = (sum(cases_2022)^2) / sum(cases_2022^2),
    pct_muni_effective          = 100 * n_eff_municipalities / n_municipalities_with_cases
  ) %>%
  arrange(pct_muni_effective)

cat("\n===== Municipality-level concentration of the 2022 outbreak =====\n")
print(as.data.frame(concentration_summary), digits = 3)

dir.create("00_Data/0_2_Processed", showWarnings = FALSE, recursive = TRUE)
save(muni_2022_pop, concentration_summary, file = "00_Data/0_2_Processed/muni_hotspot_2022.RData")
readr::write_csv(concentration_summary, "00_Data/0_2_Processed/muni_concentration_summary_2022.csv")

## ---- 5. Load municipality geometries, join, map ---------------------------

adm2 <- st_read(
  "00_Data/0_1_Raw/rasterfiles/bra_adm_ibge_2020/bra_admbnda_adm2_ibge_2020.shp",
  quiet = TRUE
)
st_crs(adm2) <- 4674  # SIRGAS 2000, standard IBGE datum

adm2 <- adm2 %>%
  mutate(muni6 = as.character(as.numeric(substr(ADM2_PCODE, 3, 9)) %/% 10L))

map_data <- adm2 %>%
  inner_join(muni_2022_pop, by = "muni6") %>%
  filter(!is.na(state_full))

cat("\nMunicipalities matched to geometry:", nrow(map_data), "\n")

## ---- 6. Choropleth: ABSOLUTE case count, one small map per state ----------
## Colored by raw case count (not per-capita rate) so the map visually
## matches what the effective-N concentration index actually measures --
## case BURDEN concentration, not per-capita risk. A rate map makes small,
## low-population municipalities with a handful of cases look "hot" even
## though they contribute almost nothing to the state's case burden.
##
## geom_sf's default coordinate system does not support facet_wrap(scales =
## "free"), which is required here since each state has a different extent.
## Build one ggplot per state (so each gets its own auto-zoomed coord_sf)
## with a SHARED fill scale (same limits/breaks) for cross-state
## comparability, then arrange them into one composite figure.

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"

fill_limits <- range(map_data$cases_2022, na.rm = TRUE)
fill_breaks <- c(0, 10, 100, 1000, 5000)

make_state_map <- function(state_name) {
  ggplot(map_data %>% filter(state_full == state_name)) +
    geom_sf(aes(fill = cases_2022), color = "white", linewidth = 0.05) +
    scale_fill_gradient(
      low = "#f4f1ea", high = "#8a1f11",
      trans = "log1p",
      limits = fill_limits,
      breaks = fill_breaks,
      labels = scales::label_number(accuracy = 1),
      name = "Confirmed cases\n(count, 2022)"
    ) +
    labs(title = state_name) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title  = element_text(face = "bold", size = 10, color = ink_primary, hjust = 0.5),
      axis.text   = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 8.5, color = ink_secondary),
      legend.text     = element_text(size = 8, color = ink_secondary)
    )
}

state_maps <- lapply(names(target_states), make_state_map)

p_muni_map <- ggpubr::ggarrange(
  plotlist      = state_maps,
  ncol          = 4, nrow = 3,
  common.legend = TRUE, legend = "right"
)

p_muni_map <- ggpubr::annotate_figure(
  p_muni_map,
  top = grid::textGrob(
    "Municipality-level concentration of the 2022 chikungunya outbreak",
    gp = grid::gpar(fontface = "bold", fontsize = 14, col = ink_primary), hjust = 0.5
  ),
  bottom = grid::textGrob(
    "Confirmed case count (not per-capita rate), by municipality of residence, 2022. Off-white: zero confirmed cases.",
    gp = grid::gpar(fontsize = 8.5, col = ink_muted), hjust = 0.5
  )
)

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "02_Outputs/2_1_Figures/muni_hotspot_map_2022.pdf",
  plot = p_muni_map, width = 14, height = 11, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/muni_hotspot_map_2022.png",
  plot = p_muni_map, width = 14, height = 11, dpi = 300, bg = "white"
)

message("[save] 02_Outputs/2_1_Figures/muni_hotspot_map_2022.pdf")
message("[save] 02_Outputs/2_1_Figures/muni_hotspot_map_2022.png")
