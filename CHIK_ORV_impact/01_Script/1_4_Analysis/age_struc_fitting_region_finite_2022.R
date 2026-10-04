load("00_Data/0_2_Processed/bra_cases_2022_cleaned.RData")
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")

observed_cases_bahia <- as.matrix(bra_bh_22[,6])
observed_cases_bahia <- matrix(observed_cases_bahia, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_ce <- as.matrix(bra_ce_22[,6])
observed_cases_ce <- matrix(observed_cases_ce, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_mg <- as.matrix(bra_mg_22[,6])
observed_cases_mg <- matrix(observed_cases_mg, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pn <- as.matrix(bra_pn_22[,6])
observed_cases_pn <- matrix(observed_cases_pn, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pa <- as.matrix(bra_pa_22[,6])
observed_cases_pa <- matrix(observed_cases_pa, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_rg <- as.matrix(bra_rg_22[,6])
observed_cases_rg <- matrix(observed_cases_rg, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_pi <- as.matrix(bra_pi_22[,6])
observed_cases_pi <- matrix(observed_cases_pi, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_ag <- as.matrix(bra_ag_22[,6])
observed_cases_ag <- matrix(observed_cases_ag, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_tc <- as.matrix(bra_tc_22[,6])
observed_cases_tc <- matrix(observed_cases_tc, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_se <- as.matrix(bra_se_22[,6])
observed_cases_se <- matrix(observed_cases_se, nrow = 20, ncol = 52, byrow = FALSE)

observed_cases_go <- as.matrix(bra_go_22[,6])
observed_cases_go <- matrix(observed_cases_go, nrow = 20, ncol = 52, byrow = FALSE)

stan_model_age <- stan_model("01_Script/1_2_SIR_models/age_struc_bra_seir_m3.stan")

# state specific foi (long-term average, kept for reference/comparison)
bra_foi <- allfoi %>% filter(country == "Brazil")
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))
bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)

bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = mean(foi_mid, na.rm = TRUE))%>%
  filter(!is.na(NAME_1))

# state specific short-term FOI (H_median column) -- this is what's used for
# the finite-history baseline seroprevalence below
#load("00_Data/0_2_Processed/bra_state_short_term_foi.RData")  # results_df
#bra_state_short_term_foi <- results_df

# age groups
age_groups <- c(mean(0:1),
                mean(1:4),
                mean(5:9),
                mean(10:11),
                mean(12:17),
                mean(18:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)


### function
make_finite_prevacc_data <- function(
    old_data,
    state_name,
    fitting_year,
    introduction_year = 2014L,
    age_midpoints = age_groups,
    B = 2L
) {
  
  # Observed age-by-week case matrix
  observed_cases <- round(
    as.matrix(old_data$observed_cases_by_age)
  )
  
  storage.mode(observed_cases) <- "integer"
  
  A <- nrow(observed_cases)
  T <- ncol(observed_cases)
  
  stopifnot(
    length(age_midpoints) == A
  )
  
  # Extract state-specific short-term FOI (bra_state_short_term_foi$H_median)
  #state_foi <- bra_state_short_term_foi$H_median[
  #  bra_state_short_term_foi$state == state_name
  #]
  #stopifnot(
  #  "state_foi not found for state_name in bra_state_short_term_foi" =
  #    length(state_foi) == 1
  #)

  # Extract state-specific long-term average FOI (bra_foi_state_summ$avg_foi)
  # -- same 8-year finite-history exposure cap as before, just a different
  # FOI level feeding it.
  state_foi <- bra_foi_state_summ$avg_foi[
    bra_foi_state_summ$NAME_1 == state_name
  ]

  stopifnot(
    "state_foi not found for state_name in bra_foi_state_summ" =
      length(state_foi) == 1
  )
  
  # Maximum duration of possible exposure before fitting year
  years_since_introduction <-
    fitting_year - introduction_year
  
  # Age-specific exposure duration capped by circulation history
  exposure_years_finite <- pmin(
    age_midpoints,
    years_since_introduction
  )
  
  # Finite-history baseline seroprevalence
  sero_finite <- 1 - exp(
    -state_foi * exposure_years_finite
  )
  
  # Initial infectious prior
  prior_I0 <- pmax(
    as.numeric(old_data$prior_I0),
    1
  )
  
  finite_data <- list(
    T = as.integer(T),
    B = as.integer(B),
    A = as.integer(A),
    
    observed_cases_by_age =
      observed_cases,
    
    N =
      as.numeric(old_data$N),
    
    r =
      rep(0, A),
    
    prior_I0 =
      prior_I0,
    
    prior_sd_I0 =
      as.numeric(old_data$prior_sd_I0),
    
    sero =
      as.numeric(sero_finite),

    # Fixed, externally-supplied symptomatic probability (age_struc_bra_seir_m3.stan
    # separates rho = p_sym * rho_sym; p_sym is data, rho_sym is estimated).
    # Matches symp_const_fixed used downstream in 02b_setup_ar_by_state.R.
    p_sym =
      0.524
  )

  return(finite_data)
}

fitting_years <- c(
  "Bahia" = 2022L,
  "Ceará" = 2022L,
  "Minas Gerais" = 2022L,
  "Pernambuco" = 2022L,
  "Paraíba" = 2022L,
  "Rio Grande do Norte" = 2022L,
  "Piauí" = 2022L,
  "Alagoas" = 2022L,
  "Tocantins" = 2022L,
  "Sergipe" = 2022L,
  "Goiás" = 2022L
)

# Built self-contained from this file's own observed_cases_* / N_* objects
# (loaded/constructed at the top of this file) -- no dependency on
# age_struc_fitting_region_2022.R's stan_data_prevacc_xx_longterm objects.
# make_finite_prevacc_data() only reads observed_cases_by_age/N/prior_I0/
# prior_sd_I0 from old_data (prior_I0 gets pmax(., 1)'d internally, so the
# raw week-1 case count can be passed as-is).
prevacc_longterm_data <- list(
  "Bahia" = list(
    observed_cases_by_age = observed_cases_bahia,
    N            = as.numeric(N_bahia$Bahia),
    prior_I0     = observed_cases_bahia[, 1],
    prior_sd_I0  = 100
  ),

  "Ceará" = list(
    observed_cases_by_age = observed_cases_ce,
    N            = as.numeric(N_ceara$Ceará),
    prior_I0     = observed_cases_ce[, 1],
    prior_sd_I0  = 100
  ),

  "Minas Gerais" = list(
    observed_cases_by_age = observed_cases_mg,
    N            = as.numeric(N_mg$`Minas Gerais`),
    prior_I0     = observed_cases_mg[, 1],
    prior_sd_I0  = 100
  ),

  "Pernambuco" = list(
    observed_cases_by_age = observed_cases_pn,
    N            = as.numeric(N_pemam$Pernambuco),
    prior_I0     = observed_cases_pn[, 1],
    prior_sd_I0  = 100
  ),

  "Paraíba" = list(
    observed_cases_by_age = observed_cases_pa,
    N            = as.numeric(N_pa$Paraíba),
    prior_I0     = observed_cases_pa[, 1],
    prior_sd_I0  = 100
  ),

  "Rio Grande do Norte" = list(
    observed_cases_by_age = observed_cases_rg,
    N            = as.numeric(N_rg$`Rio Grande do Norte`),
    prior_I0     = observed_cases_rg[, 1],
    prior_sd_I0  = 100
  ),

  "Piauí" = list(
    observed_cases_by_age = observed_cases_pi,
    N            = as.numeric(N_pi$Piauí),
    prior_I0     = observed_cases_pi[, 1],
    prior_sd_I0  = 100
  ),

  "Alagoas" = list(
    observed_cases_by_age = observed_cases_ag,
    N            = as.numeric(N_ag$Alagoas),
    prior_I0     = observed_cases_ag[, 1],
    prior_sd_I0  = 100
  ),

  "Tocantins" = list(
    observed_cases_by_age = observed_cases_tc,
    N            = as.numeric(N_tc$Tocantins),
    prior_I0     = observed_cases_tc[, 1],
    prior_sd_I0  = 100
  ),

  "Sergipe" = list(
    observed_cases_by_age = observed_cases_se,
    N            = as.numeric(N_se$Sergipe),
    prior_I0     = observed_cases_se[, 1],
    prior_sd_I0  = 100
  ),

  "Goiás" = list(
    observed_cases_by_age = observed_cases_go,
    N            = as.numeric(N_go$Goiás),
    prior_I0     = observed_cases_go[, 1],
    prior_sd_I0  = 100
  )
)

prevacc_finite_data <- lapply(
  names(prevacc_longterm_data),
  function(state_name) {
    
    make_finite_prevacc_data(
      old_data =
        prevacc_longterm_data[[state_name]],
      
      state_name =
        state_name,
      
      fitting_year =
        fitting_years[[state_name]],
      
      introduction_year =
        2014L,
      
      age_midpoints =
        age_groups,
      
      B =
        2L
    )
  }
)

names(prevacc_finite_data) <-
  names(prevacc_longterm_data)

stan_data_prevacc_bh_finite <-
  prevacc_finite_data[["Bahia"]]

stan_data_prevacc_ce_finite <-
  prevacc_finite_data[["Ceará"]]

stan_data_prevacc_mg_finite <-
  prevacc_finite_data[["Minas Gerais"]]

stan_data_prevacc_pn_finite <-
  prevacc_finite_data[["Pernambuco"]]

stan_data_prevacc_pa_finite <-
  prevacc_finite_data[["Paraíba"]]

stan_data_prevacc_rg_finite <-
  prevacc_finite_data[["Rio Grande do Norte"]]

stan_data_prevacc_pi_finite <-
  prevacc_finite_data[["Piauí"]]

stan_data_prevacc_ag_finite <-
  prevacc_finite_data[["Alagoas"]]

stan_data_prevacc_tc_finite <-
  prevacc_finite_data[["Tocantins"]]

stan_data_prevacc_se_finite <-
  prevacc_finite_data[["Sergipe"]]

stan_data_prevacc_go_finite <-
  prevacc_finite_data[["Goiás"]]


### -----------------------------------------------------------------------
### Sanity check: finite-history baseline seroprevalence-by-age curve,
### all states -- inspect before committing to the expensive Stan fits below
### -----------------------------------------------------------------------

sero_by_age_finite <- purrr::imap_dfr(
  prevacc_finite_data,
  function(dat, state_name) {
    tibble::tibble(
      state = state_name,
      age   = age_groups,
      sero  = as.numeric(dat[["sero"]])
    )
  }
)

# Legend/line order follows seroprevalence at the oldest age group, so the
# color key reads top-to-bottom in the same order the curves stack on the
# right edge of the plot.
sero_state_order <- sero_by_age_finite %>%
  filter(age == max(age)) %>%
  arrange(desc(sero)) %>%
  pull(state)

sero_by_age_finite <- sero_by_age_finite %>%
  mutate(state = factor(state, levels = sero_state_order))

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
grid_col      <- "#e1e0d9"
axis_col      <- "#c3c2b7"

sero_state_cols <- colorRampPalette(ggsci::pal_lancet("lanonc")(9))(
  length(sero_state_order)
)
names(sero_state_cols) <- sero_state_order

# Recompute sero = 1 - exp(-foi * pmin(age, years_since_introduction)) over a
# fine, continuous age grid (rather than just the 20 sparse age_group
# midpoints), so the pre-cap portion renders as a smooth exponential rise
# instead of a couple of straight segments between distant points.
age_grid_fine <- seq(0, max(age_groups), length.out = 300)

sero_by_age_finite_smooth <- purrr::imap_dfr(
  fitting_years,
  function(fyear, state_name) {

    state_foi <- bra_foi_state_summ$avg_foi[
      bra_foi_state_summ$NAME_1 == state_name
    ]

    years_since_introduction <- fyear - 2014L

    tibble::tibble(
      state = state_name,
      age   = age_grid_fine,
      sero  = 1 - exp(-state_foi * pmin(age_grid_fine, years_since_introduction))
    )
  }
) %>%
  mutate(state = factor(state, levels = sero_state_order))

p_sero_finite <- ggplot(
  sero_by_age_finite_smooth,
  aes(x = age, y = sero * 100, color = state)
) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = sero_state_cols, breaks = sero_state_order) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title    = "Finite-history baseline seroprevalence by age, before model fitting",
    subtitle = "sero = 1 - exp(-FOI x exposure years capped at years since 2014 introduction), by state",
    x        = "Age (years)",
    y        = "Seropositive (%)",
    color    = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                  = element_text(color = ink_primary),
    plot.title            = element_text(face = "bold", size = 12.5, color = ink_primary),
    plot.subtitle         = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    axis.title            = element_text(size = 9.5, color = ink_secondary),
    axis.text             = element_text(color = ink_primary, size = 8.5),
    axis.line.x           = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks            = element_blank(),
    panel.grid.major      = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor      = element_blank(),
    legend.position       = "right",
    legend.text           = element_text(size = 8.5, color = ink_secondary),
    legend.key.height     = unit(14, "pt"),
    plot.title.position   = "plot",
    plot.margin           = ggplot2::margin(12, 16, 10, 12)
  )

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "02_Outputs/2_1_Figures/sero_finite_by_age.pdf",
  plot = p_sero_finite, width = 8.5, height = 5.5, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/sero_finite_by_age.png",
  plot = p_sero_finite, width = 8.5, height = 5.5, dpi = 400, bg = "white"
)

print(p_sero_finite)


### -----------------------------------------------------------------------
### Overlay: finite-history (capped) vs long-term average (uncapped) curve,
### faceted by state
### -----------------------------------------------------------------------

# Reuses the same bra_foi_state_summ (mean of all foi_mid, computed above at
# L47-50) that feeds sero_finite/make_finite_prevacc_data(), so this curve
# differs from the capped one ONLY by the exposure cap -- not also by FOI
# estimator. (Previously this loaded 00_Data/0_2_Processed/bra_foi_state_summ.RData,
# which is built in bra_state_analysis.R from a different statistic --
# median of foi_mid[foi_mid > 0] -- and differed from the in-script mean by
# up to ~10% for some states, confounding the cap-only comparison.)
sero_longterm_by_age <- purrr::map_dfr(
  sero_state_order,
  function(state_name) {

    foi_longterm <- bra_foi_state_summ$avg_foi[
      bra_foi_state_summ$NAME_1 == state_name
    ]

    stopifnot(
      "state_name not found in bra_foi_state_summ" =
        length(foi_longterm) == 1
    )

    tibble::tibble(
      state = state_name,
      age   = age_groups,
      sero  = 1 - exp(-foi_longterm * age_groups)
    )
  }
)

col_capped   <- "#2a78d6"   # finite-history (capped at introduction year)
col_uncapped <- "#eb6834"   # long-term average (uncapped, lifetime exposure)

sero_compare <- bind_rows(
  sero_by_age_finite %>%
    mutate(assumption = "Finite-history (capped)"),
  sero_longterm_by_age %>%
    mutate(
      state      = factor(state, levels = sero_state_order),
      assumption = "Long-term average (uncapped)"
    )
) %>%
  mutate(
    assumption = factor(
      assumption,
      levels = c("Finite-history (capped)", "Long-term average (uncapped)")
    )
  )

p_sero_compare <- ggplot(
  sero_compare,
  aes(x = age, y = sero * 100, color = assumption)
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state, ncol = 4) +
  scale_color_manual(values = c(col_capped, col_uncapped)) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title    = "Finite-history vs long-term average baseline seroprevalence, by state",
    subtitle = "Capped: exposure begins at 2014 introduction. Uncapped: lifetime exposure at the long-term average FOI (bra_foi_state_summ$avg_foi)",
    x        = "Age (years)",
    y        = "Seropositive (%)",
    color    = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    text                  = element_text(color = ink_primary),
    plot.title            = element_text(face = "bold", size = 12.5, color = ink_primary),
    plot.subtitle         = element_text(size = 8.5, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    axis.title            = element_text(size = 9, color = ink_secondary),
    axis.text             = element_text(color = ink_primary, size = 7.5),
    axis.line.x           = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks            = element_blank(),
    panel.grid.major      = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor      = element_blank(),
    panel.spacing         = unit(10, "pt"),
    strip.text            = element_text(face = "bold", size = 8.5, color = ink_primary),
    strip.background      = element_blank(),
    legend.position       = "top",
    legend.justification  = "left",
    legend.text           = element_text(size = 9, color = ink_secondary),
    plot.title.position   = "plot",
    plot.margin           = ggplot2::margin(12, 16, 10, 12)
  )

ggsave(
  "02_Outputs/2_1_Figures/sero_finite_vs_longterm_by_state.pdf",
  plot = p_sero_compare, width = 11, height = 8, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/sero_finite_vs_longterm_by_state.png",
  plot = p_sero_compare, width = 11, height = 8, dpi = 400, bg = "white"
)

print(p_sero_compare)


### -----------------------------------------------------------------------
### Post-hoc sensitivity: +/-10% uncertainty on the introduction-year
### assumption (2014), WITHOUT refitting. make_finite_prevacc_data() is
### reused as-is with a shifted introduction_year, so years_since_introduction
### (and thus the exposure cap) moves by +/-10% around the base 8 years used
### for fitting (7.2y / 8y / 8.8y). This only recomputes the sero_finite
### curve for display -- it does NOT touch fit_prevacc_*_finite, which stays
### fit under the fixed introduction_year = 2014L assumption.
### -----------------------------------------------------------------------

introduction_year_variants <- c(
  "+10% exposure (8.8y)" = 2022 - 8.8,
  "Base (8y, as fitted)" = 2014,
  "-10% exposure (7.2y)" = 2022 - 7.2
)

sero_introyear_sensitivity <- purrr::map_dfr(
  names(introduction_year_variants),
  function(variant_name) {
    intro_year <- introduction_year_variants[[variant_name]]

    purrr::map_dfr(
      sero_state_order,
      function(state_name) {

        finite_data <- make_finite_prevacc_data(
          old_data          = prevacc_longterm_data[[state_name]],
          state_name        = state_name,
          fitting_year      = fitting_years[[state_name]],
          introduction_year = intro_year,
          age_midpoints     = age_groups,
          B                 = 2L
        )

        tibble::tibble(
          state     = state_name,
          age       = age_groups,
          sero      = finite_data$sero,
          variant   = variant_name
        )
      }
    )
  }
) %>%
  mutate(
    state   = factor(state, levels = sero_state_order),
    variant = factor(variant, levels = names(introduction_year_variants))
  )

col_introyear <- c(
  "+10% exposure (8.8y)" = "#eb6834",
  "Base (8y, as fitted)" = "#2a78d6",
  "-10% exposure (7.2y)" = "#3aa76d"
)

p_sero_introyear_sensitivity <- ggplot(
  sero_introyear_sensitivity,
  aes(x = age, y = sero * 100, color = variant)
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~state, ncol = 4) +
  scale_color_manual(values = col_introyear) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title    = "Sensitivity of finite-history baseline seroprevalence to the introduction-year assumption",
    subtitle = "+/-10% on the 8-year exposure cap (2014 introduction, 2022 fitting year) -- sero_finite recomputed post-hoc, fit_prevacc_*_finite is not refit under these variants",
    x        = "Age (years)",
    y        = "Seropositive (%)",
    color    = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    text                  = element_text(color = ink_primary),
    plot.title            = element_text(face = "bold", size = 12.5, color = ink_primary),
    plot.subtitle         = element_text(size = 8.5, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    axis.title            = element_text(size = 9, color = ink_secondary),
    axis.text             = element_text(color = ink_primary, size = 7.5),
    axis.line.x           = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks            = element_blank(),
    panel.grid.major      = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor      = element_blank(),
    panel.spacing         = unit(10, "pt"),
    strip.text            = element_text(face = "bold", size = 8.5, color = ink_primary),
    strip.background      = element_blank(),
    legend.position       = "top",
    legend.justification  = "left",
    legend.text           = element_text(size = 9, color = ink_secondary),
    plot.title.position   = "plot",
    plot.margin           = ggplot2::margin(12, 16, 10, 12)
  )

ggsave(
  "02_Outputs/2_1_Figures/sero_finite_introyear_sensitivity_by_state.pdf",
  plot = p_sero_introyear_sensitivity, width = 11, height = 8, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/sero_finite_introyear_sensitivity_by_state.png",
  plot = p_sero_introyear_sensitivity, width = 11, height = 8, dpi = 400, bg = "white"
)

print(p_sero_introyear_sensitivity)


### fitting: finite-history baseline immunity

fit_prevacc_bh_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_bh_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

fit_prevacc_ce_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ce_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

fit_prevacc_mg_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_mg_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_pn_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pn_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_pa_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pa_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.8,
    max_treedepth = 15
  )
)

fit_prevacc_rg_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_rg_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_pi_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_pi_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_ag_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_ag_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_tc_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_tc_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_se_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_se_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

fit_prevacc_go_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_go_finite,
  iter = 2000,
  chains = 1,
  warmup = 1000,
  thin = 1,
  seed = 123,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  )
)

save(
  "fit_prevacc_bh_finite", "fit_prevacc_ce_finite",
  "fit_prevacc_mg_finite", "fit_prevacc_pn_finite",
  "fit_prevacc_pa_finite", "fit_prevacc_rg_finite",
  "fit_prevacc_pi_finite", "fit_prevacc_ag_finite",
  "fit_prevacc_tc_finite", "fit_prevacc_se_finite",
  "fit_prevacc_go_finite",
  
  file = "00_Data/0_2_Processed/fits_prevacc_finite.RData"

)


## post processing
list_bh <- create_summary_df(fit_prevacc_bh_finite ,
                             bra_sum_bh,
                             region = "Bahia")

list_ce <- create_summary_df(fit_prevacc_ce_finite ,
                             bra_sum_ce,
                             region = "Ceará")

list_mg <- create_summary_df(fit_prevacc_mg_finite ,
                             bra_sum_mg,
                             region = "Minas Gerais")

list_pn <- create_summary_df(fit_prevacc_pn_finite ,
                             bra_sum_pn,
                             region = "Pernambuco")

list_pa <- create_summary_df(fit_prevacc_pa_finite ,
                             bra_sum_pa,
                             region = "Paraíba")

list_rg <- create_summary_df(fit_prevacc_rg_finite ,
                             bra_sum_rg,
                             region = "Rio Grande do Norte")

list_pi <- create_summary_df(fit_prevacc_pi_finite ,
                             bra_sum_pi,
                             region = "Piauí")

list_ag <- create_summary_df(fit_prevacc_ag_finite ,
                             bra_sum_ag,
                             region = "Alagoas")

list_tc <- create_summary_df(fit_prevacc_tc_finite ,
                             bra_sum_tc,
                             region = "Tocantins")

list_se <- create_summary_df(fit_prevacc_se_finite ,
                             bra_sum_se,
                             region = "Sergipe")

list_go <- create_summary_df(fit_prevacc_go_finite ,
                             bra_sum_go,
                             region = "Goiás")

## bahia
df_bh <- list_bh$df_out

df_bh_summ <- list_bh$df_summ

observed_bh <- list_bh$observed

overall_fit_gg(observed_bh, df_bh_summ)


## ceara
df_ce <- list_ce$df_out

df_ce_summ <- list_ce$df_summ

observed_ce <- list_ce$observed


overall_fit_gg(observed_ce, df_ce_summ)


## minas gerais
df_mg <- list_mg$df_out

df_mg_summ <- list_mg$df_summ

observed_mg <- list_mg$observed

age_strat_gg(df_mg)

overall_fit_gg(observed_mg, df_mg_summ)


## pernambuco
df_pn <- list_pn$df_out

df_pn_summ <- list_pn$df_summ

observed_pn <- list_pn$observed

age_strat_gg(df_pn)

overall_fit_gg(observed_pn, df_pn_summ)


## paraiba
df_pa <- list_pa$df_out

df_pa_summ <- list_pa$df_summ

observed_pa <- list_pa$observed


overall_fit_gg(observed_pa, df_pa_summ)


## rio grande norte
df_rg <- list_rg$df_out

df_rg_summ <- list_rg$df_summ

observed_rg <- list_rg$observed

age_strat_gg(df_rg)

overall_fit_gg(observed_rg, df_rg_summ)

## Piuai
df_pi <- list_pi$df_out

df_pi_summ <- list_pi$df_summ

observed_pi <- list_pi$observed

age_strat_gg(df_pi)

overall_fit_gg(observed_pi, df_pi_summ)

## Alagoas
df_ag <- list_ag$df_out

df_ag_summ <- list_ag$df_summ

observed_ag <- list_ag$observed

age_strat_gg(df_ag)

overall_fit_gg(observed_ag, df_ag_summ)

## tocantins
df_tc <- list_tc$df_out

df_tc_summ <- list_tc$df_summ

observed_tc <- list_tc$observed

age_strat_gg(df_tc)

overall_fit_gg(observed_tc, df_tc_summ)

## Sergipe
df_se_22 <- list_se$df_out

df_se_summ_22 <- list_se$df_summ

observed_se <- list_se$observed


overall_fit_gg(observed_se, df_se_summ_22)

## goias
df_go_22 <- list_go$df_out

df_go_summ_22 <- list_go$df_summ

observed_go <- list_go$observed

age_strat_gg(df_go_22)

overall_fit_gg(observed_go, df_go_summ_22)

## overall graph 
observed_ce$region <- "Ceará" 
observed_bh$region <- "Bahia" 
observed_ag$region <- "Alagoas" 
observed_mg$region <- "Minas Gerais" 
observed_tc$region <- "Tocantins"
observed_pa$region <- "Paraíba" 
observed_pi$region <- "Piauí"
observed_pn$region <- "Pernambuco"
observed_rg$region <- "Rio Grande do Norte"
observed_se$region <- "Sergipe"
observed_go$region <- "Goiás"

observed_all <- bind_rows(observed_ce, observed_bh, observed_ag,
                          observed_mg, observed_tc, observed_pa,
                          observed_pi, observed_pn, observed_rg,
                          observed_se, observed_go)

pred_all <- bind_rows(
  df_ag_summ, df_bh_summ, df_ce_summ,
  df_mg_summ, df_pa_summ, df_pi_summ,
  df_rg_summ, df_pn_summ, df_tc_summ,
  df_go_summ_22, df_se_summ_22
)

model_fit <- 
  
  ggplot()+
  geom_point(data = observed_all, aes(x = Week, y = Observed), size = 0.8)+
  facet_wrap(~region, scales = "free_y")+
  # Predicted line
  geom_line(data = pred_all, aes(x = Week, y = Median, color = Type), size = 1) +
  
  # Prediction interval (ribbon)
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), 
              alpha = 0.2)+
  theme_pubclean()+
  theme(legend.position = "right")+
  ylab("Predicted and observed reported symptomatic cases")+
  scale_y_continuous(labels = comma)

model_fit_2 <- 
  ggarrange(model_fit, 
            ncol = 1,
            labels = c("D"),
            common.legend = TRUE,
            legend = "bottom",
            align = "none")

ggsave(filename = "02_Outputs/2_1_Figures/figs5.jpg", model_fit, width = 12, height = 8, dpi = 1200)


save("observed_ce", "observed_bh", "observed_ag", 
     "observed_mg", "observed_tc", "observed_pa",
     "observed_se", "observed_go",
     "observed_pi", "observed_pn", "observed_rg", "observed_all", "bra_all_sum",
     file = "00_Data/0_2_Processed/observed_2022.RData")


