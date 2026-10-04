## -----------------------------------------------------------------------
## foi_vs_2022_ar_plot.R
##
## Visualize the relationship between each cluster's 2015-2021 cumulative
## FOI (prior exposure, back-calculated from reported cases) and its 2022
## outbreak attack rate -- the inverse relationship that explains why the
## Focal (small, explosive) clusters are exactly the ones with the LOWEST
## prior exposure: low FOI 2015-2021 -> large naive susceptible pool ->
## explosive local burnout when the 2022 outbreak reached them.
## -----------------------------------------------------------------------

library(dplyr)
library(ggplot2)

load("00_Data/0_2_Processed/foi_2015_2021_by_cluster.RData")  # foi_2015_2021_by_cluster
load("00_Data/0_2_Processed/focal_broadscale_units.RData")     # units_classified

plot_data <- foi_2015_2021_by_cluster %>%
  dplyr::select(unit_id, foi_median, foi_low95, foi_hi95, pct_draws_clipped) %>%
  dplyr::inner_join(
    units_classified %>% dplyr::select(unit_id, unit_label, setting_type, infection_ar_pct, population),
    by = "unit_id"
  ) %>%
  dplyr::mutate(
    data_flag = ifelse(pct_draws_clipped > 5, "Flagged (unreliable FOI)", as.character(setting_type))
  )

cat("===== Correlation: FOI (2015-2021) vs 2022 infection attack rate =====\n")
clean <- plot_data %>% dplyr::filter(pct_draws_clipped <= 5)
print(cor.test(log(clean$foi_median), log(clean$infection_ar_pct), method = "spearman"))

## ---- Plot ------------------------------------------------------------------

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
grid_col      <- "#e1e0d9"

col_focal   <- "#a3271f"
col_broad   <- "#2a78d6"
col_flagged <- "#c3c2b7"

plot_data$data_flag <- factor(
  plot_data$data_flag,
  levels = c("Focal", "Broad-scale", "Flagged (unreliable FOI)")
)

p_foi_vs_ar <- ggplot(
  plot_data,
  aes(x = foi_median, y = infection_ar_pct, color = data_flag)
) +
  geom_errorbarh(
    aes(xmin = foi_low95, xmax = foi_hi95),
    height = 0, alpha = 0.35, linewidth = 0.5
  ) +
  geom_point(size = 3.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = plot_data %>% dplyr::filter(data_flag != "Broad-scale" | infection_ar_pct > 3),
    aes(label = unit_label),
    size = 2.6, color = ink_secondary, seed = 1, max.overlaps = 20,
    segment.color = grid_col, segment.size = 0.3
  ) +
  scale_x_log10(labels = scales::label_percent(accuracy = 0.001)) +
  scale_y_log10(labels = scales::label_number(suffix = "%")) +
  scale_color_manual(
    values = c(Focal = col_focal, "Broad-scale" = col_broad, "Flagged (unreliable FOI)" = col_flagged),
    name = NULL
  ) +
  labs(
    title    = "Prior exposure (2015-2021 FOI) vs. the 2022 outbreak attack rate",
    subtitle = "Low prior FOI -> large naive susceptible pool -> explosive local 2022 outbreak (the Focal clusters)",
    x        = "Cumulative FOI, 2015-2021 (median, 95% UI; log scale)",
    y        = "Estimated infection attack rate, 2022 (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    text                = element_text(color = ink_primary),
    plot.title          = element_text(face = "bold", size = 13, color = ink_primary),
    plot.subtitle       = element_text(size = 9, color = ink_secondary, margin = ggplot2::margin(b = 10)),
    panel.grid.minor    = element_blank(),
    legend.position     = "top",
    legend.justification = "left",
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(12, 16, 10, 12)
  )

dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("02_Outputs/2_1_Figures/foi_vs_2022_ar.png", plot = p_foi_vs_ar, width = 9, height = 7.5, dpi = 400, bg = "white")
ggsave("02_Outputs/2_1_Figures/foi_vs_2022_ar.pdf", plot = p_foi_vs_ar, width = 9, height = 7.5, device = cairo_pdf)

message("[save] 02_Outputs/2_1_Figures/foi_vs_2022_ar.png/.pdf")
