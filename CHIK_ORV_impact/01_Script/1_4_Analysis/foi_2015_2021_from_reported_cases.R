## -----------------------------------------------------------------------
## foi_2015_2021_from_reported_cases.R
##
## Back-calculate state-level force of infection (FOI) for 2015-2021 from
## cumulative reported confirmed cases, propagating full uncertainty in
## the reporting probability (rho) and the symptomatic probability
## (P_symp) via Monte Carlo.
##
##   Infection_cum   = Cases_reported / (rho * P_symp)
##   q                = Infection_cum / N_state                 (cumulative attack rate)
##   FOI(2015-2021)   = -log(1 - q) / 7                          (7 = number of years)
##
## Inputs:
##   - cumulative_state_age  (00_Data/0_2_Processed/sinan_confirmed_cases_age_state_2015_2021.RData)
##     -> built by sinan_state_age_cumulative_2015_2021.R
##   - rho_df                (00_Data/0_2_Processed/rho_df.RData)
##     -> per-state reporting probability: rho_p50 / rho_p2.5 / rho_p97.5
##        (summary quantiles only, not raw posterior draws)
##   - lhs_sample_young$symp_overall (00_Data/0_2_Processed/lhs_sample_young.RDS)
##     -> 1000 LHS draws of the overall symptomatic-probability-given-infection
##   - tot_pop_df            (00_Data/0_2_Processed/tot_pop_df.RData)
##     -> all-age population by state
##
## Uncertainty propagation:
##   rho_df only gives 3 quantiles per state (no raw draws), so a full
##   distribution is reconstructed on the logit scale (rho is a bounded
##   probability): rho_draw = plogis(rnorm(n, logit(rho_p50), sd_logit)),
##   with sd_logit chosen so the reconstructed 2.5/97.5 percentiles match
##   rho_p2.5/rho_p97.5. This is then paired draw-by-draw with
##   lhs_sample_young$symp_overall (n = 1000) for the Monte Carlo.
## -----------------------------------------------------------------------

library(dplyr)
library(readr)

set.seed(123)

## ---- 0. Config -------------------------------------------------------------

N_YEARS <- 7   # 2015-2021 inclusive
N_DRAWS <- 1000

OUT_DIR   <- "00_Data/0_2_Processed"
OUT_RDATA <- file.path(OUT_DIR, "foi_2015_2021_reported.RData")
OUT_CSV   <- file.path(OUT_DIR, "foi_2015_2021_reported.csv")

target_states <- c("Ceará", "Bahia", "Paraíba", "Pernambuco",
                   "Rio Grande do Norte", "Piauí", "Tocantins", "Alagoas",
                   "Minas Gerais", "Sergipe", "Goiás")

## ---- 1. Load inputs ---------------------------------------------------------

load(file.path(OUT_DIR, "sinan_confirmed_cases_age_state_2015_2021.RData"))  # cumulative_state_age
load(file.path(OUT_DIR, "rho_df.RData"))                                     # rho_df
load(file.path(OUT_DIR, "tot_pop_df.RData"))                                 # tot_pop_df
load(file.path(OUT_DIR, "bra_foi_state_summ.RData"))                         # bra_foi_state_summ (sf)
lhs_sample_young <- readRDS(file.path(OUT_DIR, "lhs_sample_young.RDS"))

stopifnot(N_DRAWS == nrow(lhs_sample_young))

# Long-term average FOI (catalytic model, whole lifetime exposure), for
# comparison against the 2015-2021-only FOI computed below.
long_term_foi_state <- bra_foi_state_summ %>%
  sf::st_drop_geometry() %>%
  dplyr::transmute(
    state_full         = NAME_1,
    long_term_avg_foi  = unname(avg_foi),
    long_term_foi_lo   = unname(foi_lo),
    long_term_foi_hi   = unname(foi_hi)
  ) %>%
  dplyr::filter(state_full %in% target_states)

stopifnot(nrow(long_term_foi_state) == length(target_states))

# State-level cumulative reported confirmed cases (sum across age groups).
cases_reported_state <- cumulative_state_age %>%
  dplyr::group_by(state_full) %>%
  dplyr::summarise(cases_reported = sum(confirmed_cases), .groups = "drop") %>%
  dplyr::filter(state_full %in% target_states)

pop_state <- tot_pop_df %>%
  dplyr::filter(state_full %in% target_states) %>%
  dplyr::select(state_full, tot_pop)

rho_state <- rho_df %>%
  dplyr::rename(state_full = region) %>%
  dplyr::filter(state_full %in% target_states)

stopifnot(
  nrow(cases_reported_state) == length(target_states),
  nrow(pop_state)            == length(target_states),
  nrow(rho_state)            == length(target_states)
)

## ---- 2. Reconstruct full rho draws from (p2.5, p50, p97.5) -----------------

# rho is a bounded (0,1) probability -> sample on the logit scale so draws
# stay in range, with sd chosen to reproduce the reported 95% UI.
sample_rho_logitnormal <- function(p50, p2.5, p97.5, n) {
  mu_logit <- qlogis(p50)
  sd_logit <- (qlogis(p97.5) - qlogis(p2.5)) / (2 * 1.959964)
  plogis(rnorm(n, mean = mu_logit, sd = sd_logit))
}

symp_draws <- lhs_sample_young$symp_overall

## ---- 3. Monte Carlo propagation, per state ----------------------------------

foi_results <- lapply(target_states, function(st) {

  cases_st <- cases_reported_state$cases_reported[cases_reported_state$state_full == st]
  n_st     <- pop_state$tot_pop[pop_state$state_full == st]
  rho_row  <- rho_state[rho_state$state_full == st, ]

  rho_draws <- sample_rho_logitnormal(
    rho_row$rho_p50, rho_row$rho_p2.5, rho_row$rho_p97.5, N_DRAWS
  )

  infection_cum_draws <- cases_st / (rho_draws * symp_draws)
  q_draws              <- infection_cum_draws / n_st

  n_clipped <- sum(q_draws >= 1)
  q_draws_clipped <- pmin(q_draws, 0.999999)

  foi_draws <- -log(1 - q_draws_clipped) / N_YEARS

  tibble::tibble(
    state_full        = st,
    cases_reported    = cases_st,
    tot_pop           = n_st,
    q_median           = median(q_draws_clipped),
    q_low95            = quantile(q_draws_clipped, 0.025),
    q_hi95             = quantile(q_draws_clipped, 0.975),
    foi_median          = median(foi_draws),
    foi_low95           = quantile(foi_draws, 0.025),
    foi_hi95            = quantile(foi_draws, 0.975),
    pct_draws_clipped  = 100 * n_clipped / N_DRAWS
  )
})

foi_2015_2021_reported <- dplyr::bind_rows(foi_results) %>%
  dplyr::left_join(long_term_foi_state, by = "state_full") %>%
  dplyr::arrange(state_full)

## ---- 4. Save -----------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save(foi_2015_2021_reported, file = OUT_RDATA)
readr::write_csv(foi_2015_2021_reported, OUT_CSV)

message(sprintf("[save] %s", OUT_RDATA))
message(sprintf("[save] %s", OUT_CSV))

## ---- 5. Print summary ----------------------------------------------------------

message("\n[summary] FOI (2015-2021), reporting- and symptomatic-probability-adjusted:")
print(foi_2015_2021_reported, width = Inf)

if (any(foi_2015_2021_reported$pct_draws_clipped > 0)) {
  message(
    "\n[WARNING] Some Monte Carlo draws implied cumulative infections >= state population ",
    "(q >= 1) and were clipped to 0.999999 before computing FOI. See 'pct_draws_clipped' ",
    "column above -- a large percentage means rho/P_symp uncertainty is pushing implausibly ",
    "low reporting rates for that state."
  )
}
