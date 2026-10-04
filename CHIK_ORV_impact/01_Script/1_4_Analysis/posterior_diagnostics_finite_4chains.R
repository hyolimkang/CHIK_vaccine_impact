# =============================================================================
# posterior_diagnostics_finite_4chains.R
#
# Posterior-diagnostic figures/table for the 11-state finite-history 4-chain
# Stan refit -- same 3 outputs as posterior_diagnostics_finite.R (trace
# plots, ESS+Rhat table, posterior predictive check), but built from the
# LIGHTWEIGHT diagnostics bundle (diagnostics_finite_4chains.RData) that
# age_struc_fitting_region_finite_2022_4chains.R computes per-state right
# after fitting, instead of from a saved raw stanfit list -- the 4-chain run
# deliberately never writes fits_prevacc_finite_4chains.RData to disk (see
# that script's header for why: ~20GB+ of unused per-age-per-week matrices
# at 16,000 draws/state). trace_arr/diag_summary/pred_cases were already
# computed via as.array()/rstan::summary()/rstan::extract() while each
# state's stanfit was still in memory, so this script only needs to re-plot/
# re-tabulate them -- no Stan object, no rstan call, needed here.
#
# Called automatically at the end of age_struc_fitting_region_finite_2022_
# 4chains.R; can also be re-run standalone once diagnostics_finite_4chains.
# RData exists on disk.
#
# Output (*_4chains suffix -- does not overwrite the 1-chain interim
# diagnostics from posterior_diagnostics_finite.R):
#   02_Outputs/2_1_Figures/trace_plots_finite_4chains/figure_trace_plots_finite_4chains_<State>.png
#     -- one file per state (chat record 2026-09-05), not a single combined
#        image -- 11 states x 6 parameters in one plot was unreadable
#   02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite_4chains.png
#   02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.xlsx
#   02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(bayesplot)
  library(scales)
  library(flextable)
  library(officer)
  library(writexl)
})

FIT_LABEL <- "4 chains, 5,000 iterations each (1,000 warm-up), 4,000 post-warm-up draws/chain"

message("Loading 00_Data/0_2_Processed/diagnostics_finite_4chains.RData ...")
load("00_Data/0_2_Processed/diagnostics_finite_4chains.RData")
# -> trace_arr_list, diag_summary_list, pred_cases_list, DIAG_PARS

par_labels <- c(
  rho = "Reporting rate (\u03c1)", gamma = "Recovery rate (\u03b3)",
  sigma = "Incubation rate (\u03c3)", alpha_log_beta = "log(\u03b2) intercept",
  sigma_beta_rw = "\u03b2 random-walk SD", lp__ = "Log posterior"
)

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("02_Outputs/2_2_Tables", showWarnings = FALSE, recursive = TRUE)

# ---- 1) Trace plots (bayesplot::mcmc_trace on the precomputed arrays) ------
# One file PER STATE (chat record 2026-09-05) -- the previous version packed
# all 11 states x 6 parameters (66 tiny facets) into a single patchwork
# image, which was unreadable (chain colours indistinguishable, axis text
# near-illegible at that scale). Splitting one-state-per-file lets each
# state's 6 panels render at a size where the 4 chains and mixing/divergence
# pattern are actually visible, with its own chain-colour legend.
trace_dir <- "02_Outputs/2_1_Figures/trace_plots_finite_4chains"
dir.create(trace_dir, showWarnings = FALSE, recursive = TRUE)

for (region in names(trace_arr_list)) {
  arr <- trace_arr_list[[region]]
  dimnames(arr)$parameters <- unname(par_labels[dimnames(arr)$parameters])
  p_state <- bayesplot::mcmc_trace(arr, facet_args = list(nrow = 2)) +
    labs(title = region, subtitle = paste0("Fit: ", FIT_LABEL)) +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 9))

  fname <- file.path(trace_dir, sprintf(
    "figure_trace_plots_finite_4chains_%s.png",
    gsub("[^A-Za-z0-9]+", "_", iconv(region, from = "UTF-8", to = "ASCII//TRANSLIT"))
  ))
  ggsave(fname, p_state, width = 12, height = 6.5, dpi = 300, bg = "white")
  message("Saved: ", fname)
}

# ---- 2) ESS (n_eff) + Rhat table (from the precomputed summaries) ----------
diag_table <- bind_rows(lapply(names(diag_summary_list), function(region) {
  summ <- diag_summary_list[[region]]
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

print(doc, target = "02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.docx")
message("Saved: 02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.docx")

writexl::write_xlsx(diag_table, "02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.xlsx")
message("Saved: 02_Outputs/2_2_Tables/table_posterior_diagnostics_finite_4chains.xlsx")

# ---- 3) Posterior predictive check (from the precomputed pred_cases draws) -
load("00_Data/0_2_Processed/observed_2022.RData")  # -> observed_all

pred_all <- bind_rows(lapply(names(pred_cases_list), function(region) {
  post <- pred_cases_list[[region]]  # draws x 52 weeks
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
    plot.subtitle = element_text(size = 10.5, colour = ink_secondary, margin = ggplot2::margin(b = 10), lineheight = 1.15),
    plot.title.position = "plot",
    axis.title = element_text(size = 12), axis.text = element_text(size = 10),
    axis.line.x = element_line(colour = axis_col, linewidth = 0.3), axis.ticks = element_blank(),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3), panel.grid.minor = element_blank(),
    panel.spacing = unit(12, "pt"), strip.text = element_text(face = "bold", size = 12),
    strip.background = element_blank(), legend.position = "bottom", legend.title = element_blank(),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

ggsave("02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite_4chains.png", p_ppc,
       width = 14, height = 10, dpi = 350, bg = "white")
message("Saved: 02_Outputs/2_1_Figures/figure_posterior_predictive_check_finite_4chains.png")
