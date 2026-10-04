# =============================================================================
# posterior_diagnostics_finite.R
#
# Posterior-diagnostic figures/table for the 11-state finite-history Stan fit,
# supporting the manuscript Methods text (chat record 2026-09-02):
#   "Posterior convergence was assessed using trace plots, effective sample
#    sizes... model fit was evaluated using posterior predictive checks
#    comparing aggregated observed weekly case counts with the corresponding
#    model-fitted reported symptomatic case counts."
#
# Built from whichever finite-history stanfit list is on disk under FIT_FILE
# below -- rstan::summary()/as.array()/extract() are chain-count-generic, so
# this runs unchanged against either the current 1-chain interim fit or the
# 4-chain/5000-iter refit (age_struc_fitting_region_finite_2022_4chains.R,
# launched in background 2026-09-02, saving to
# 00_Data/0_2_Processed/fits_prevacc_finite_4chains.RData) once that
# completes -- just flip FIT_FILE/FIT_LABEL below and re-run.
#
# Diagnostic parameters: rho (reporting rate), gamma (recovery rate), sigma
# (E->I rate), alpha_log_beta / sigma_beta_rw (transmission-rate random-walk
# hyperparameters), lp__ (log posterior) -- the model's scalar hyperparameters;
# age/week-indexed vectors (I0, beta, S/E/I/R, ...) are excluded from
# trace/ESS (not meaningfully traceable/tabulated one-by-one at this scale).
#
# Output (project convention: 02_Outputs/2_1_Figures for figures,
# 02_Outputs/2_2_Tables for tables -- matches posterior_analysis.R's
# table_s3.xlsx):
#   02_Outputs/2_1_Figures/figure_trace_plots_finite.png
#   02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite.png
#   02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.xlsx
#   02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(bayesplot)
  library(rstan)
  library(patchwork)
  library(scales)
  library(flextable)
  library(officer)
  library(writexl)
})

# ---- Config: point at the fit to diagnose -----------------------------------
FIT_FILE  <- "00_Data/0_2_Processed/fits_prevacc_finite.RData"
FIT_LABEL <- "1 chain, 2,000 iterations (1,000 warm-up) -- interim, pre-refit"
# Once the background 4-chain refit finishes, switch to:
# FIT_FILE  <- "00_Data/0_2_Processed/fits_prevacc_finite_4chains.RData"
# FIT_LABEL <- "4 chains, 5,000 iterations each (1,000 warm-up), 4,000 post-warm-up draws/chain"

message("Loading ", FIT_FILE, " ...")
load(FIT_FILE)

state_fits <- list(
  "Ceará" = fit_prevacc_ce_finite, "Bahia" = fit_prevacc_bh_finite,
  "Paraíba" = fit_prevacc_pa_finite, "Pernambuco" = fit_prevacc_pn_finite,
  "Rio Grande do Norte" = fit_prevacc_rg_finite, "Piauí" = fit_prevacc_pi_finite,
  "Alagoas" = fit_prevacc_ag_finite, "Tocantins" = fit_prevacc_tc_finite,
  "Minas Gerais" = fit_prevacc_mg_finite, "Sergipe" = fit_prevacc_se_finite,
  "Goiás" = fit_prevacc_go_finite
)

DIAG_PARS <- c("rho", "gamma", "sigma", "alpha_log_beta", "sigma_beta_rw", "lp__")
par_labels <- c(
  rho = "Reporting rate (\u03c1)", gamma = "Recovery rate (\u03b3)",
  sigma = "Incubation rate (\u03c3)", alpha_log_beta = "log(\u03b2) intercept",
  sigma_beta_rw = "\u03b2 random-walk SD", lp__ = "Log posterior"
)

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("02_Outputs/2_2_Tables", showWarnings = FALSE, recursive = TRUE)

# ---- 1) Trace plots (bayesplot::mcmc_trace, chain-count-generic) -----------
trace_plots <- lapply(names(state_fits), function(region) {
  arr <- as.array(state_fits[[region]], pars = DIAG_PARS)
  dimnames(arr)$parameters <- unname(par_labels[dimnames(arr)$parameters])
  bayesplot::mcmc_trace(arr, facet_args = list(nrow = 2)) +
    ggtitle(region) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none",
          strip.text = element_text(size = 8),
          axis.text = element_text(size = 6))
})
names(trace_plots) <- names(state_fits)

p_trace <- patchwork::wrap_plots(trace_plots, ncol = 3) +
  patchwork::plot_annotation(
    title = "Posterior trace plots by state",
    subtitle = paste0("Scalar hyperparameters; fit: ", FIT_LABEL),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, colour = "grey30"))
  )

ggsave("02_Outputs/2_1_Figures/figure_trace_plots_finite.png", p_trace,
       width = 20, height = 24, dpi = 250, bg = "white", limitsize = FALSE)
message("Saved: 02_Outputs/2_1_Figures/figure_trace_plots_finite.png")

# ---- 2) ESS (n_eff) + Rhat table (rstan::summary, chain-count-generic) ------
diag_table <- bind_rows(lapply(names(state_fits), function(region) {
  summ <- rstan::summary(state_fits[[region]], pars = DIAG_PARS)$summary
  tibble::tibble(
    State = region,
    Parameter = par_labels[rownames(summ)],
    Median = summ[, "50%"],
    Bulk_ESS = summ[, "n_eff"],
    Rhat = summ[, "Rhat"]
  )
}))

diag_table_fmt <- diag_table %>%
  mutate(
    Median = sprintf("%.4f", Median),
    Bulk_ESS = sprintf("%.0f", Bulk_ESS),
    Rhat = sprintf("%.3f", Rhat)
  )

ft <- flextable::flextable(diag_table_fmt) %>%
  flextable::set_header_labels(Bulk_ESS = "Effective sample size", Rhat = "R\u0302") %>%
  flextable::merge_v(j = "State") %>%
  flextable::valign(j = "State", valign = "top") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("State", "Parameter"), align = "left", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Posterior convergence diagnostics, by state", style = "heading 2") %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste0(
      "Fit: ", FIT_LABEL, ". Effective sample size (bulk n_eff) and R\u0302 (potential scale ",
      "reduction factor) for each state's key scalar hyperparameters, from rstan::summary(). ",
      "R\u0302 close to 1 (< 1.01) and effective sample size in the thousands indicate good mixing ",
      "and convergence; see accompanying trace plots for a visual check."
    ),
    style = "Normal"
  )

print(doc, target = "02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.docx")
message("Saved: 02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.docx")

writexl::write_xlsx(diag_table, "02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.xlsx")
message("Saved: 02_Outputs/2_2_Tables/table_posterior_diagnostics_finite.xlsx")

# ---- 3) Posterior predictive check: aggregated observed vs fitted cases ----
load("00_Data/0_2_Processed/observed_2022.RData")  # -> observed_all

pred_all <- bind_rows(lapply(names(state_fits), function(region) {
  post <- rstan::extract(state_fits[[region]], pars = "pred_cases")$pred_cases  # draws x 52 weeks
  tibble::tibble(
    Week = seq_len(ncol(post)),
    Median = apply(post, 2, median, na.rm = TRUE),
    Lower = apply(post, 2, quantile, probs = 0.025, na.rm = TRUE),
    Upper = apply(post, 2, quantile, probs = 0.975, na.rm = TRUE),
    Type = "Predicted",
    region = region
  )
}))

ink_primary <- "#0b0b0b"; ink_secondary <- "#52514e"
grid_col <- "#e7e6e0"; axis_col <- "#c3c2b7"; col_predicted <- "#C0392B"

p_ppc <- ggplot() +
  geom_ribbon(data = pred_all, aes(x = Week, ymin = Lower, ymax = Upper, fill = Type), alpha = 0.18, colour = NA) +
  geom_line(data = pred_all, aes(x = Week, y = Median, colour = Type), linewidth = 0.8) +
  geom_point(data = observed_all, aes(x = Week, y = Observed, shape = Type), size = 1.5, colour = ink_primary, alpha = 0.85) +
  facet_wrap(~region, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c(Predicted = col_predicted), name = NULL) +
  scale_fill_manual(values = c(Predicted = col_predicted), name = NULL) +
  scale_shape_manual(values = c(Observed = 16), name = NULL) +
  scale_x_continuous(name = "Epidemiological week (2022)", breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(name = "Reported symptomatic cases", labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Posterior predictive check: fitted vs. observed weekly reported cases",
    subtitle = paste0("Points = observed; line/ribbon = posterior median and 95% credible interval. Fit: ", FIT_LABEL)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(colour = ink_primary),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10.5, colour = ink_secondary, margin = margin(b = 10), lineheight = 1.15),
    plot.title.position = "plot",
    axis.title = element_text(size = 12), axis.text = element_text(size = 10),
    axis.line.x = element_line(colour = axis_col, linewidth = 0.3), axis.ticks = element_blank(),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3), panel.grid.minor = element_blank(),
    panel.spacing = unit(12, "pt"), strip.text = element_text(face = "bold", size = 12),
    strip.background = element_blank(), legend.position = "bottom", legend.title = element_blank(),
    plot.margin = margin(14, 18, 12, 14)
  )

ggsave("02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite.png", p_ppc,
       width = 14, height = 10, dpi = 350, bg = "white")
message("Saved: 02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite.png")
