## -----------------------------------------------------------------------
## plot_foi_2015_2021_vs_longterm.R
##
## Comparison figure: state-level force of infection (FOI) estimated from
## 2015-2021 reported confirmed cases (with 95% uncertainty interval, via
## Monte Carlo propagation of reporting- and symptomatic-probability
## uncertainty) vs. the long-term average FOI (catalytic model).
##
## Input : 00_Data/0_2_Processed/foi_2015_2021_reported.RData
##         (built by foi_2015_2021_from_reported_cases.R)
## Output: 02_Outputs/2_1_Figures/foi_2015_2021_vs_longterm.pdf / .png
## -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

load("00_Data/0_2_Processed/foi_2015_2021_reported.RData")  # foi_2015_2021_reported

## ---- 1. Reshape to long format for the two comparison series --------------

plot_df <- dplyr::bind_rows(
  foi_2015_2021_reported %>%
    dplyr::transmute(
      state_full, period = "2015–2021 (reported cases)",
      value = foi_median, low = foi_low95, high = foi_hi95
    ),
  foi_2015_2021_reported %>%
    dplyr::transmute(
      state_full, period = "Long-term average (catalytic)",
      value = long_term_avg_foi, low = long_term_foi_lo, high = long_term_foi_hi
    )
) %>%
  dplyr::mutate(
    period = factor(
      period,
      levels = c("2015–2021 (reported cases)", "Long-term average (catalytic)")
    )
  )

# Order states by the 2015-2021 FOI median (ascending, so the largest
# estimate ends up at the top after coord_flip()).
order_levels <- foi_2015_2021_reported %>%
  dplyr::arrange(foi_median) %>%
  dplyr::pull(state_full)
plot_df$state_full <- factor(plot_df$state_full, levels = order_levels)

## ---- 2. Palette (validated categorical slots 1-2; see dataviz skill) ------

col_recent    <- "#2a78d6"   # slot 1 - blue
col_longterm  <- "#eb6834"   # slot 2 - orange

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
grid_col      <- "#e1e0d9"
axis_col      <- "#c3c2b7"

## ---- 3. Build the figure ----------------------------------------------------

p <- ggplot(plot_df, aes(x = state_full, y = value, color = period,
                        shape = period, linetype = period)) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    width = 0, linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.6, stroke = 0.8
  ) +
  coord_flip() +
  scale_color_manual(values = c(col_recent, col_longterm)) +
  scale_shape_manual(values = c(16, 18)) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    title    = "FOI by state: 2015–2021 vs. long-term average",
    #subtitle = "2015–2021 estimated from cumulative reported confirmed cases (95% uncertainty interval);\nlong-term average from the catalytic (whole-lifetime exposure) model",
    x        = NULL,
    y        = "average FOI/year",
    color    = NULL,
    shape    = NULL,
    linetype = NULL,
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(color = ink_primary),
    plot.title          = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle       = element_text(size = 9.5, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    plot.caption        = element_text(size = 7.5, color = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)),
    axis.text           = element_text(color = ink_primary, size = 10),
    axis.title.x        = element_text(color = ink_secondary, size = 10, margin = ggplot2::margin(t = 8)),
    axis.line.x         = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks          = element_blank(),
    panel.grid.major.x  = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.position     = "top",
    legend.justification = "left",
    legend.text         = element_text(size = 9.5, color = ink_secondary),
    legend.key.spacing.x = unit(12, "pt"),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin         = ggplot2::margin(12, 16, 10, 12)
  )

## ---- 4. Save -----------------------------------------------------------------

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "02_Outputs/2_1_Figures/foi_2015_2021_vs_longterm.pdf",
  plot = p, width = 7.2, height = 5.5, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/foi_2015_2021_vs_longterm.png",
  plot = p, width = 7.2, height = 5.5, dpi = 400, bg = "white"
)

message("[save] 02_Outputs/2_1_Figures/foi_2015_2021_vs_longterm.pdf")
message("[save] 02_Outputs/2_1_Figures/foi_2015_2021_vs_longterm.png")

## -----------------------------------------------------------------------
## Immunity profile by state: % seropositive vs. age, from the 2015-2021
## reported-case FOI (median + 95% UI).
##
##   sero(age) = 1 - exp( -FOI * min(age, assessment_year - introduction_year) )
##
## i.e. exposure is capped at the number of years CHIKV has been circulating
## before the assessment year (chikungunya was introduced to Brazil in 2014;
## the 2022 fits cap exposure at 2022 - 2014 = 8 years).
## -----------------------------------------------------------------------

introduction_year <- 2014L
assessment_year   <- 2022L
years_exposure_cap <- assessment_year - introduction_year

age_seq <- seq(0, 90, by = 0.5)

immunity_df <- foi_2015_2021_reported %>%
  dplyr::select(state_full, foi_median, foi_low95, foi_hi95) %>%
  tidyr::crossing(age = age_seq) %>%
  dplyr::mutate(
    exposure_years = pmin(age, years_exposure_cap),
    sero_median     = 1 - exp(-foi_median * exposure_years),
    sero_low        = 1 - exp(-foi_low95  * exposure_years),
    sero_high       = 1 - exp(-foi_hi95   * exposure_years)
  )

# Facet order: highest median FOI first (reuse the ordering already computed
# above for the comparison plot).
immunity_df$state_full <- factor(immunity_df$state_full, levels = rev(order_levels))

p_immunity <- ggplot(immunity_df, aes(x = age, y = sero_median)) +
  geom_ribbon(
    aes(ymin = sero_low, ymax = sero_high),
    fill = col_recent, alpha = 0.16
  ) +
  geom_line(
    color = col_recent, linewidth = 0.7
  ) +
  facet_wrap(~state_full, ncol = 4) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 30),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title    = "Predicted baseline immunity profile by state",
    subtitle = sprintf(
      "1 − exp(−FOI × min(age, %d)); FOI from 2015–2021 reported cases (median, 95%% UI)",
      years_exposure_cap
    ),
    x        = "Age (years)",
    y        = "% seropositive",
    caption  = sprintf(
      "Exposure duration capped at %d years (chikungunya introduction in %d to the %d assessment year).",
      years_exposure_cap, introduction_year, assessment_year
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                   = element_text(color = ink_primary),
    plot.title             = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle          = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    plot.caption           = element_text(size = 7.5, color = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)),
    axis.text              = element_text(color = ink_primary, size = 8.5),
    axis.title.x           = element_text(color = ink_secondary, size = 10, margin = ggplot2::margin(t = 8)),
    axis.title.y           = element_text(color = ink_secondary, size = 10, margin = ggplot2::margin(r = 8)),
    axis.line.x            = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks             = element_blank(),
    panel.grid.major       = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor       = element_blank(),
    panel.spacing          = unit(14, "pt"),
    strip.text             = element_text(face = "bold", size = 9, color = ink_primary, hjust = 0),
    plot.title.position    = "plot",
    plot.caption.position  = "plot",
    plot.margin            = ggplot2::margin(12, 16, 10, 12)
  )

ggsave(
  "02_Outputs/2_1_Figures/immunity_profile_by_state.pdf",
  plot = p_immunity, width = 9, height = 8, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/immunity_profile_by_state.png",
  plot = p_immunity, width = 9, height = 8, dpi = 400, bg = "white"
)

message("[save] 02_Outputs/2_1_Figures/immunity_profile_by_state.pdf")
message("[save] 02_Outputs/2_1_Figures/immunity_profile_by_state.png")

## -----------------------------------------------------------------------
## Immunity profile by state, LONG-TERM AVERAGE FOI: finite-history
## (8-year post-introduction cap, flat thereafter) vs. the simple
## whole-lifetime catalytic model (uncapped, 1 - exp(-FOI * age)).
##
## Below the cap the two models are IDENTICAL (exposure = age either way),
## so only one line is visible; above the cap they diverge -- the capped
## model plateaus (solid), while the uncapped model is shown as a dashed
## counterfactual extension of what the same FOI would imply if exposure
## had accumulated over the person's whole life rather than just since
## the 2014 introduction.
## -----------------------------------------------------------------------

longterm_immunity_df <- foi_2015_2021_reported %>%
  dplyr::select(state_full, long_term_avg_foi) %>%
  tidyr::crossing(age = age_seq) %>%
  dplyr::mutate(
    exposure_years_capped = pmin(age, years_exposure_cap),
    sero_finite            = 1 - exp(-long_term_avg_foi * exposure_years_capped),
    sero_simple            = 1 - exp(-long_term_avg_foi * age)
  )

longterm_immunity_df$state_full <- factor(longterm_immunity_df$state_full, levels = rev(order_levels))

longterm_immunity_long <- dplyr::bind_rows(
  longterm_immunity_df %>%
    dplyr::transmute(state_full, age, model = "Finite-history (8-yr cap)", sero = sero_finite),
  longterm_immunity_df %>%
    dplyr::filter(age >= years_exposure_cap) %>%
    dplyr::transmute(state_full, age, model = "Simple lifetime catalytic (uncapped)", sero = sero_simple)
) %>%
  dplyr::mutate(
    model = factor(model, levels = c("Finite-history (8-yr cap)", "Simple lifetime catalytic (uncapped)"))
  )

p_immunity_longterm <- ggplot(longterm_immunity_long, aes(x = age, y = sero, linetype = model)) +
  geom_line(color = col_longterm, linewidth = 0.7) +
  facet_wrap(~state_full, ncol = 4) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 30),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title    = "Predicted baseline immunity profile by state: long-term average FOI",
    subtitle = sprintf(
      "Solid: 1 − exp(−FOI × min(age, %d)), capped at introduction. Dashed: 1 − exp(−FOI × age), uncapped whole-lifetime catalytic model.",
      years_exposure_cap
    ),
    x        = "Age (years)",
    y        = "% seropositive",
    caption  = sprintf(
      "Long-term average FOI (catalytic model). Exposure capped at %d years (chikungunya introduction in %d to the %d assessment year); the two models are identical below the cap.",
      years_exposure_cap, introduction_year, assessment_year
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                   = element_text(color = ink_primary),
    plot.title             = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle          = element_text(size = 8.5, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    plot.caption           = element_text(size = 7.5, color = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)),
    axis.text              = element_text(color = ink_primary, size = 8.5),
    axis.title.x           = element_text(color = ink_secondary, size = 10, margin = ggplot2::margin(t = 8)),
    axis.title.y           = element_text(color = ink_secondary, size = 10, margin = ggplot2::margin(r = 8)),
    axis.line.x            = element_line(color = axis_col, linewidth = 0.3),
    axis.ticks             = element_blank(),
    panel.grid.major       = element_line(color = grid_col, linewidth = 0.3),
    panel.grid.minor       = element_blank(),
    panel.spacing          = unit(14, "pt"),
    strip.text             = element_text(face = "bold", size = 9, color = ink_primary, hjust = 0),
    legend.position         = "top",
    legend.justification    = "left",
    legend.text             = element_text(size = 9, color = ink_secondary),
    plot.title.position    = "plot",
    plot.caption.position  = "plot",
    plot.margin            = ggplot2::margin(12, 16, 10, 12)
  )

ggsave(
  "02_Outputs/2_1_Figures/immunity_profile_longterm_finite_vs_simple.pdf",
  plot = p_immunity_longterm, width = 9, height = 8.5, device = cairo_pdf
)
ggsave(
  "02_Outputs/2_1_Figures/immunity_profile_longterm_finite_vs_simple.png",
  plot = p_immunity_longterm, width = 9, height = 8.5, dpi = 400, bg = "white"
)

message("[save] 02_Outputs/2_1_Figures/immunity_profile_longterm_finite_vs_simple.pdf")
message("[save] 02_Outputs/2_1_Figures/immunity_profile_longterm_finite_vs_simple.png")
