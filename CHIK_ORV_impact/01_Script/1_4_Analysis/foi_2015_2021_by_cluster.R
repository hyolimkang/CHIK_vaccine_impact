## -----------------------------------------------------------------------
## foi_2015_2021_by_cluster.R
##
## Cluster-level analogue of foi_2015_2021_from_reported_cases.R.
##
## Back-calculates FOI (2015-2021) for each of the 25 classification units
## (state-level units where aggregation checks out, or significant
## sub-region clusters where it doesn't -- see all_states_coherence.R /
## focal_vs_broadscale_classification.R), using the SAME method as the
## state-level version: cumulative reported confirmed cases -> infections
## (via Monte Carlo over reporting probability rho and symptomatic
## probability) -> cumulative attack rate -> FOI.
##
##   Infection_cum   = Cases_reported / (rho * P_symp)
##   q                = Infection_cum / N_unit
##   FOI(2015-2021)   = -log(1 - q) / 7
##
## Population and case counts are restricted to each unit's constituent
## municipalities (the same cluster used for the Focal/Broad-scale
## classification), NOT the full state -- consistent with how the
## classification's attack rates were computed.
##
## No cluster-level rho has been fitted; each unit inherits its PARENT
## STATE's fitted rho (same approximation already used for the infection
## attack rate in transmission_tier_classification.R).
##
## Inputs:
##   - 00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds
##   - 00_Data/0_2_Processed/rho_df.RData          (state-level rho quantiles)
##   - 00_Data/0_2_Processed/lhs_sample_young.RDS   (symptomatic probability draws)
##   - 00_Data/0_2_Processed/all_states_coherence.RData   (results, summary_table)
##   - 00_Data/0_2_Processed/focal_broadscale_units.RData (units_classified)
## -----------------------------------------------------------------------

library(dplyr)
library(sf)

set.seed(123)

N_YEARS <- 7   # 2015-2021 inclusive
N_DRAWS <- 1000

OUT_DIR   <- "00_Data/0_2_Processed"
OUT_RDATA <- file.path(OUT_DIR, "foi_2015_2021_by_cluster.RData")
OUT_CSV   <- file.path(OUT_DIR, "foi_2015_2021_by_cluster.csv")

## ---- 1. Rebuild the municipality -> unit_id mapping -----------------------
## Same logic as focal_vs_broadscale_classification.R section 4, kept in
## sync deliberately rather than re-sourcing that whole script.

load(file.path(OUT_DIR, "all_states_coherence.RData"))     # results, summary_table
load(file.path(OUT_DIR, "focal_broadscale_units.RData"))    # units_classified, POP_THRESHOLD

aggregate_ok_states <- summary_table %>% dplyr::filter(verdict == "Aggregate OK") %>% dplyr::pull(state_full)

muni_to_unit <- bind_rows(lapply(names(results), function(st) {

  r <- results[[st]]
  sf_df <- sf::st_drop_geometry(r$sf_set) %>%
    dplyr::filter(cluster %in% r$cluster_share$cluster)

  sf_df <- if (st %in% aggregate_ok_states) {
    sf_df %>% dplyr::mutate(unit_id = st)
  } else {
    sf_df %>% dplyr::mutate(unit_id = paste0(st, " - cluster ", cluster))
  }

  sf_df %>% dplyr::select(muni6, state_full, unit_id)
}))

cat("Municipalities mapped to units:", nrow(muni_to_unit), "\n")
cat("Units covered:", dplyr::n_distinct(muni_to_unit$unit_id), "of", nrow(units_classified), "\n")

## ---- 2. Cumulative 2015-2021 confirmed cases, by unit ----------------------

muni_week <- readRDS("00_Data/0_1_Raw/chik_brazil_muni_week_2015_2024.rds")

cases_2015_2021_by_unit <- muni_week %>%
  dplyr::filter(epi_year >= 2015, epi_year <= 2021) %>%
  dplyr::inner_join(muni_to_unit, by = "muni6") %>%
  dplyr::group_by(unit_id) %>%
  dplyr::summarise(cases_reported = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")

rm(muni_week); gc()

## ---- 3. Reporting probability (rho) and symptomatic probability draws ----

load(file.path(OUT_DIR, "rho_df.RData"))          # rho_df (region, rho_p50, rho_p2.5, rho_p97.5)
lhs_sample_young <- readRDS(file.path(OUT_DIR, "lhs_sample_young.RDS"))
stopifnot(N_DRAWS == nrow(lhs_sample_young))

sample_rho_logitnormal <- function(p50, p2.5, p97.5, n) {
  mu_logit <- qlogis(p50)
  sd_logit <- (qlogis(p97.5) - qlogis(p2.5)) / (2 * 1.959964)
  plogis(rnorm(n, mean = mu_logit, sd = sd_logit))
}

rho_state <- rho_df %>% dplyr::rename(state_full = region)

symp_draws <- lhs_sample_young$symp_overall

## ---- 4. Monte Carlo FOI per unit --------------------------------------------

units_info <- units_classified %>%
  dplyr::select(unit_id, unit_label, state_full, setting_type, population) %>%
  dplyr::inner_join(cases_2015_2021_by_unit, by = "unit_id") %>%
  dplyr::left_join(rho_state %>% dplyr::select(state_full, rho_p50, rho_p2.5, rho_p97.5), by = "state_full")

missing_rho <- units_info %>% dplyr::filter(is.na(rho_p50))
if (nrow(missing_rho) > 0) {
  cat("\n[WARNING] Units with no parent-state rho match:\n")
  print(as.data.frame(missing_rho %>% dplyr::select(unit_id, state_full)))
}

foi_results <- lapply(seq_len(nrow(units_info)), function(i) {

  row <- units_info[i, ]

  rho_draws <- sample_rho_logitnormal(row$rho_p50, row$rho_p2.5, row$rho_p97.5, N_DRAWS)

  infection_cum_draws <- row$cases_reported / (rho_draws * symp_draws)
  q_draws              <- infection_cum_draws / row$population

  n_clipped       <- sum(q_draws >= 1)
  q_draws_clipped <- pmin(q_draws, 0.999999)

  foi_draws <- -log(1 - q_draws_clipped) / N_YEARS

  tibble::tibble(
    unit_id            = row$unit_id,
    unit_label         = row$unit_label,
    state_full         = row$state_full,
    setting_type       = row$setting_type,
    population         = row$population,
    cases_reported     = row$cases_reported,
    rho_used_p50       = row$rho_p50,
    q_median           = median(q_draws_clipped),
    q_low95            = quantile(q_draws_clipped, 0.025),
    q_hi95             = quantile(q_draws_clipped, 0.975),
    foi_median         = median(foi_draws),
    foi_low95          = quantile(foi_draws, 0.025),
    foi_hi95           = quantile(foi_draws, 0.975),
    pct_draws_clipped  = 100 * n_clipped / N_DRAWS
  )
})

foi_2015_2021_by_cluster <- dplyr::bind_rows(foi_results) %>%
  dplyr::arrange(dplyr::desc(foi_median))

## ---- 5. Save + report ---------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save(foi_2015_2021_by_cluster, file = OUT_RDATA)
readr::write_csv(foi_2015_2021_by_cluster, OUT_CSV)

cat("\n===== FOI (2015-2021) by cluster =====\n")
print(as.data.frame(foi_2015_2021_by_cluster), digits = 3)

cat(sprintf("\n[save] %s\n[save] %s\n", OUT_RDATA, OUT_CSV))

if (any(foi_2015_2021_by_cluster$pct_draws_clipped > 0)) {
  message(
    "\n[WARNING] Some units had Monte Carlo draws with cumulative infections >= unit population ",
    "(q >= 1), clipped before computing FOI -- check 'pct_draws_clipped'. This is most likely for ",
    "small-population Focal units where low rho/symptomatic draws push implausibly high implied ",
    "attack rates."
  )
  print(
    as.data.frame(
      foi_2015_2021_by_cluster %>%
        dplyr::filter(pct_draws_clipped > 0) %>%
        dplyr::select(unit_label, population, cases_reported, pct_draws_clipped)
    )
  )
}
