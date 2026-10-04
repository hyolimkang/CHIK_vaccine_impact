## model calibration check
library(tidyr)
library(dplyr)
library(ggplot2)

sero_comparison_long <- sero_comparison_bahia %>%
  pivot_longer(
    cols = c(LongTerm, FiniteHistory),
    names_to = "Assumption",
    values_to = "Seroprevalence"
  )

ggplot(
  sero_comparison_long,
  aes(
    x = AgeMidpoint,
    y = Seroprevalence,
    linetype = Assumption
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  labs(
    x = "Age midpoint",
    y = "Baseline seroprevalence"
  ) +
  theme_bw()

population_bahia <- as.numeric(
  N_bahia$Bahia
)

overall_sero_longterm <- weighted.mean(
  sero_bahia_longterm,
  w = population_bahia
)

overall_sero_finite <- weighted.mean(
  sero_bahia_finite,
  w = population_bahia
)

c(
  LongTerm = overall_sero_longterm,
  FiniteHistory = overall_sero_finite
)

stan_data_prevacc_bh_finite <- stan_data_prevacc_bh

stan_data_prevacc_bh_finite$sero <-
  as.numeric(sero_bahia_finite)

fit_prevacc_bh_finite <- sampling(
  object = stan_model_age,
  data = stan_data_prevacc_bh_finite,
  iter = 1200,            # Reduced from 4000
  chains = 1,             # Reduced from 4
  warmup = 600,           # Specify warmup period
  thin = 2,               # Add thinning
  seed = 123,
  control = list(
    adapt_delta = 0.95,    # Start with lower adaptation
    max_treedepth = 12    # Reduced tree depth
  )
)


################################################################################
# downstream ggplots
################################################################################
list_bh_finite <- create_summary_df(fit_prevacc_bh_finite,
                             bra_sum_bh,
                             region = "Bahia")

df_bh_finite <- list_bh_finite$df_out

df_bh_summ_finite <- list_bh_finite$df_summ

observed_bh <- list_bh_finite$observed

overall_fit_gg(observed_bh, df_bh_summ_finite)

post_bh <- rstan::extract(
  fit_prevacc_bh_finite,
  pars = c(
    "beta_observed",
    "sigma_beta_rw"
  )
)

beta_draws <- post_bh$beta_observed

beta_summary <- data.frame(
  Week = seq_len(ncol(beta_draws)),
  Median = apply(
    beta_draws,
    2,
    median
  ),
  Lower = apply(
    beta_draws,
    2,
    quantile,
    probs = 0.025
  ),
  Upper = apply(
    beta_draws,
    2,
    quantile,
    probs = 0.975
  )
)

ggplot(
  beta_summary,
  aes(
    x = Week,
    y = Median
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    alpha = 0.2
  ) +
  geom_line() +
  labs(
    x = "Week",
    y = "Transmission rate"
  ) +
  theme_bw()
