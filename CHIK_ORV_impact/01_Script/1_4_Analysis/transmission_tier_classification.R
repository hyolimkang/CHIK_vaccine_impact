## -----------------------------------------------------------------------
## transmission_tier_classification.R
##
## Redefine the analysis unit as "state" for the 5 states where state-level
## aggregation checks out (see all_states_coherence_check.R), and as
## "significant cluster within state" for the 6 states that need splitting.
## Compute a simple attack rate (2022 confirmed cases / population) for
## each resulting unit and classify into high / moderate / low transmission
## tiers (tertiles) -- a purely data-driven classification that does not
## depend on SEIR fitting, so it sidesteps the S0-beta identifiability
## problem entirely.
## -----------------------------------------------------------------------

library(dplyr)
library(ggplot2)

load("00_Data/0_2_Processed/all_states_coherence.RData")  # results, summary_table
load("00_Data/0_2_Processed/muni_hotspot_2022.RData")     # muni_2022_pop
load("00_Data/0_2_Processed/rho_by_state.RData")           # rho_table (state-level fitted reporting probability)

aggregate_ok_states <- summary_table %>% filter(verdict == "Aggregate OK") %>% pull(state_full)
split_states        <- summary_table %>% filter(verdict == "Split recommended") %>% pull(state_full)

## ---- 1. Units for "Aggregate OK" states: significant clusters combined ----
## IMPORTANT: use the SAME denominator logic as the split states below --
## population of the significant (>=5% share) clusters only, not the full
## state population. "Aggregate OK" means those clusters are synchronized
## enough to fit as ONE region; it does not mean the whole state (including
## municipalities with negligible case counts) should dilute the
## denominator. Using full state population here while using
## cluster-restricted population for split states would make attack rates
## incomparable across the two groups.

units_aggregate_ok <- bind_rows(lapply(aggregate_ok_states, function(st) {

  r <- results[[st]]
  sf_df <- sf::st_drop_geometry(r$sf_set) %>%
    left_join(muni_2022_pop %>% dplyr::select(muni6, population), by = "muni6")

  sig <- sf_df %>% filter(cluster %in% r$cluster_share$cluster)

  tibble::tibble(
    state_full = st,
    unit_id    = st,
    unit_label = st,
    unit_type  = "State (aggregate OK)",
    cases_2022 = sum(sig$cases_2022),
    population = sum(sig$population)
  )
}))

## ---- 2. Units for "Split" states: one row per significant cluster ----------

units_split <- bind_rows(lapply(split_states, function(st) {

  r <- results[[st]]
  sf_df <- sf::st_drop_geometry(r$sf_set) %>%
    left_join(muni_2022_pop %>% dplyr::select(muni6, population), by = "muni6")

  cluster_pop <- sf_df %>%
    filter(cluster %in% r$cluster_share$cluster) %>%
    group_by(cluster) %>%
    summarise(
      # top_muni MUST be computed before cases_2022 is overwritten by
      # sum() below, or which.max() sees the already-summed scalar
      # instead of the per-municipality vector.
      top_muni     = muni_name[which.max(cases_2022)],
      population   = sum(population),
      cases_2022    = sum(cases_2022),
      .groups = "drop"
    )

  cluster_pop %>%
    mutate(
      state_full = st,
      unit_id    = paste0(st, " - cluster ", cluster),
      unit_label = paste0(st, ": ", top_muni, " cluster"),
      unit_type  = "Sub-region (state split)"
    )
}))

## ---- 3. Combine, compute attack rate, classify into tiers -------------------

units_all <- bind_rows(
  units_aggregate_ok %>% dplyr::select(state_full, unit_id, unit_label, unit_type, cases_2022, population),
  units_split %>% dplyr::select(state_full, unit_id, unit_label, unit_type, cases_2022, population)
) %>%
  left_join(rho_table %>% dplyr::select(state_full, rho_median), by = "state_full") %>%
  mutate(
    attack_rate_per100k = 1e5 * cases_2022 / population,
    attack_rate_pct      = attack_rate_per100k / 1000,
    # Estimated INFECTION attack rate: confirmed-case AR divided by the
    # state's own fitted reporting probability (rho), which in this model
    # already combines symptomatic probability and case-ascertainment
    # probability. Sub-region clusters inherit their parent state's rho
    # (no cluster-level rho is fitted).
    infection_ar_pct     = attack_rate_pct / rho_median
  ) %>%
  arrange(infection_ar_pct)

## ---- Tier classification: Jenks natural breaks, not forced tertiles ------
## Forced tertiles (equal group SIZE) produced a "High" group spanning
## 2-21% infection AR -- internally incoherent, and defensible only by
## "these happened to be the top third of OUR 25 units," not by any
## external or structural criterion. Jenks natural-breaks instead finds
## the cutpoints that minimize within-group variance / maximize between-
## group separation in the data actually observed -- groups may be
## unequal in size, but each is more internally homogeneous, and the
## break itself falls in a genuine gap in the distribution rather than at
## an arbitrary percentile.
jenks <- classInt::classIntervals(units_all$infection_ar_pct, n = 3, style = "jenks")
tier_cuts_infection_ar <- jenks$brks[2:3]

units_all <- units_all %>%
  mutate(
    transmission_tier = case_when(
      infection_ar_pct <= tier_cuts_infection_ar[1] ~ "Low",
      infection_ar_pct <= tier_cuts_infection_ar[2] ~ "Moderate",
      TRUE                                            ~ "High"
    ),
    transmission_tier = factor(transmission_tier, levels = c("Low", "Moderate", "High"))
  )

cat("===== Units and transmission tiers (Jenks natural breaks on estimated infection AR) =====\n")
print(
  as.data.frame(
    units_all %>%
      dplyr::select(unit_label, unit_type, cases_2022, population, rho_median,
             attack_rate_pct, infection_ar_pct, transmission_tier)
  ),
  digits = 4
)

cat("\nJenks break points (estimated infection attack rate, %):\n")
print(round(tier_cuts_infection_ar, 2))

cat("\nUnits per tier:\n")
print(table(units_all$transmission_tier))

cat("\nWithin-tier range (min-max infection AR %), to check homogeneity:\n")
print(
  units_all %>%
    group_by(transmission_tier) %>%
    summarise(n = n(), min_pct = min(infection_ar_pct), max_pct = max(infection_ar_pct))
)

dir.create("00_Data/0_2_Processed", showWarnings = FALSE, recursive = TRUE)
save(units_all, tier_cuts_infection_ar, file = "00_Data/0_2_Processed/transmission_tier_units.RData")
readr::write_csv(units_all, "00_Data/0_2_Processed/transmission_tier_units.csv")

## ---- 4. Plot: ranked attack rate, colored by tier ---------------------------

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
grid_col      <- "#e1e0d9"
axis_col      <- "#c3c2b7"

col_low  <- "#2a78d6"
col_mod  <- "#c9a227"
col_high <- "#a3271f"

units_all$unit_label <- factor(units_all$unit_label, levels = units_all$unit_label)

p_tiers <- ggplot(units_all, aes(x = infection_ar_pct, y = unit_label, color = transmission_tier)) +
  geom_vline(xintercept = tier_cuts_infection_ar, linetype = "22", color = grid_col, linewidth = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Low = col_low, Moderate = col_mod, High = col_high), name = "Transmission\ntier") +
  scale_x_continuous(labels = scales::label_number(accuracy = 1, suffix = "%")) +
  labs(
    title    = "Classification units by estimated infection attack rate, into transmission tiers",
    subtitle = "Confirmed-case AR / state rho (reporting probability); tiers = Jenks natural breaks (not forced-equal tertiles)",
    x        = "Estimated infection attack rate, 2022 (%)",
    y        = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    text                   = element_text(color = ink_primary),
    plot.title             = element_text(face = "bold", size = 12.5, color = ink_primary),
    plot.subtitle          = element_text(size = 8.5, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    axis.text.y            = element_text(size = 8, color = ink_primary),
    axis.line.x            = element_line(color = axis_col, linewidth = 0.3),
    panel.grid.major.x     = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.major.y     = element_blank(),
    panel.grid.minor       = element_blank(),
    legend.position        = "top",
    legend.justification   = "left",
    plot.title.position    = "plot",
    plot.margin            = ggplot2::margin(12, 16, 10, 12)
  )

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("02_Outputs/2_1_Figures/transmission_tier_classification.pdf",
       plot = p_tiers, width = 9, height = 8, device = cairo_pdf)
ggsave("02_Outputs/2_1_Figures/transmission_tier_classification.png",
       plot = p_tiers, width = 9, height = 8, dpi = 400, bg = "white")

message("[save] 00_Data/0_2_Processed/transmission_tier_units.csv")
message("[save] 02_Outputs/2_1_Figures/transmission_tier_classification.pdf/.png")
