## -----------------------------------------------------------------------
## sinan_state_age_cumulative_2015_2021.R
##
## Goal: cumulative reported CONFIRMED chikungunya cases by state x age
##       group (2015-2021) for the 11 states used in the vaccine-impact
##       model, built from the individual-level SINAN dataset already
##       cleaned in the CHIK_CLIM repo.
##
## Source data (external repo, not part of CHIK_ORV_impact):
##   CHIK_CLIM/01_Data/chik_sinan_individual_2015_2024.rds
##   (produced by CHIK_CLIM/02_Script/clean_chik_sinan_brazil.R)
##
## NOTE ON YEAR RANGE:
##   The raw SINAN extracts in CHIK_CLIM start at CHIKBR15.csv.zip
##   (2015) - there is no 2014 file, so "2014-2021" as requested is not
##   available from this source. This script covers 2015-2021 instead.
##
## Age groups match the 20-bin scheme used throughout the SEIR fitting
## pipeline (e.g. age_struc_bra_seir_m2.stan / age_struc_fitting_region_func.R):
##   <1, 1-4, 5-9, 10-11, 12-17, 18-19, 20-24, 25-29, 30-34, 35-39,
##   40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85+
## -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)

## ---- 0. Config -----------------------------------------------------------

YEAR_START <- 2015
YEAR_END   <- 2021

# External repo path (CHIK_CLIM) - update here if the repo moves.
SINAN_INDIVIDUAL_RDS <- file.path(
  "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine",
  "Documents/GitHub/CHIK_CLIM/CHIK_CLIM/01_Data",
  "chik_sinan_individual_2015_2024.rds"
)

OUT_DIR <- "00_Data/0_2_Processed"
OUT_RDATA <- file.path(OUT_DIR, "sinan_confirmed_cases_age_state_2015_2021.RData")
OUT_CSV   <- file.path(OUT_DIR, "sinan_confirmed_cases_age_state_2015_2021.csv")

## ---- 1. State map (same convention as bra_chik_allcases.R) --------------

state_map <- data.frame(
  code_state = sprintf("%02d", c(
    11, 12, 13, 14, 15, 16, 17,               # North
    21, 22, 23, 24, 25, 26, 27, 28, 29,       # Northeast
    31, 32, 33, 35,                           # Southeast
    41, 42, 43,                               # South
    50, 51, 52, 53                            # Central-West
  )),
  state_full = c(
    "Rondônia", "Acre", "Amazonas", "Roraima", "Pará", "Amapá", "Tocantins",
    "Maranhão", "Piauí", "Ceará", "Rio Grande do Norte", "Paraíba",
    "Pernambuco", "Alagoas", "Sergipe", "Bahia",
    "Minas Gerais", "Espírito Santo", "Rio de Janeiro", "São Paulo",
    "Paraná", "Santa Catarina", "Rio Grande do Sul",
    "Mato Grosso do Sul", "Mato Grosso", "Goiás", "Distrito Federal"
  ),
  stringsAsFactors = FALSE
)

# The 11 states used in the vaccine-impact model.
target_states <- c("17", "22", "23", "24", "25", "26", "27", "28", "29", "31", "52")

## ---- 2. Load individual-level SINAN data ---------------------------------

if (!file.exists(SINAN_INDIVIDUAL_RDS)) {
  stop(
    "Cannot find cleaned SINAN individual-level data at:\n  ", SINAN_INDIVIDUAL_RDS,
    "\nRun CHIK_CLIM/02_Script/clean_chik_sinan_brazil.R first (in the CHIK_CLIM repo)."
  )
}

sinan <- readRDS(SINAN_INDIVIDUAL_RDS)
message(sprintf("[load] %s rows read from %s", format(nrow(sinan), big.mark = ","), basename(SINAN_INDIVIDUAL_RDS)))

## ---- 3. Filter: 11 states, confirmed cases, valid age, 2015-2021 --------

sinan_filt <- sinan %>%
  dplyr::mutate(
    code_state = sprintf("%02d", suppressWarnings(as.numeric(dplyr::coalesce(uf_residence, uf_notif)))),
    event_year = as.integer(format(event_date, "%Y"))
  ) %>%
  dplyr::filter(
    code_state %in% target_states,
    is_confirmed_chik,
    !is.na(age_years), age_years >= 0, age_years <= 105,
    !is.na(event_year), event_year >= YEAR_START, event_year <= YEAR_END
  ) %>%
  dplyr::left_join(state_map, by = "code_state")

message(sprintf(
  "[filter] %s confirmed cases kept (11 states, %d-%d, valid age)",
  format(nrow(sinan_filt), big.mark = ","), YEAR_START, YEAR_END
))

## ---- 4. Assign the 20-group age scheme -----------------------------------

age_gr_levels <- c(
  "<1", "1-4", "5-9", "10-11", "12-17", "18-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
  "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
)
age_breaks <- c(0, 1, 5, 10, 12, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf)

sinan_filt <- sinan_filt %>%
  dplyr::mutate(
    age_group = cut(
      age_years,
      breaks = age_breaks,
      labels = age_gr_levels,
      right = FALSE,
      include.lowest = TRUE
    )
  )

## ---- 5. Aggregate ----------------------------------------------------------

# Main output: cumulative confirmed cases, state x age group, 2015-2021 combined.
cumulative_state_age <- sinan_filt %>%
  dplyr::count(state_full, age_group, name = "confirmed_cases") %>%
  tidyr::complete(
    state_full = state_map$state_full[state_map$code_state %in% target_states],
    age_group  = age_gr_levels,
    fill = list(confirmed_cases = 0L)
  ) %>%
  dplyr::mutate(age_group = factor(age_group, levels = age_gr_levels)) %>%
  dplyr::arrange(state_full, age_group)

# Secondary output: same, but broken down by year (for QC / trend checks).
cumulative_state_age_year <- sinan_filt %>%
  dplyr::count(state_full, event_year, age_group, name = "confirmed_cases") %>%
  tidyr::complete(
    state_full = state_map$state_full[state_map$code_state %in% target_states],
    event_year = YEAR_START:YEAR_END,
    age_group  = age_gr_levels,
    fill = list(confirmed_cases = 0L)
  ) %>%
  dplyr::mutate(age_group = factor(age_group, levels = age_gr_levels)) %>%
  dplyr::arrange(state_full, event_year, age_group)

## ---- 6. Save ---------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

save(
  cumulative_state_age, cumulative_state_age_year,
  file = OUT_RDATA
)
readr::write_csv(cumulative_state_age, OUT_CSV)

message(sprintf("[save] %s", OUT_RDATA))
message(sprintf("[save] %s", OUT_CSV))

## ---- 7. Quick sanity summary ------------------------------------------------

message("\n[summary] Cumulative confirmed cases by state (2015-2021):")
print(
  cumulative_state_age %>%
    dplyr::group_by(state_full) %>%
    dplyr::summarise(total_confirmed = sum(confirmed_cases), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_confirmed))
)

message(sprintf(
  "\n[NOTE] 2014 not available in the SINAN source (earliest raw extract is CHIKBR15 = 2015); this table covers %d-%d only.",
  YEAR_START, YEAR_END
))
