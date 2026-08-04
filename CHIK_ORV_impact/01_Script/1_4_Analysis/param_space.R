# all possible parameter space search
n_ages <- 18

all_strategies <- list()
idx <- 1

# 1 block (with as many continuous age bins as possible)
for(start1 in 1:n_ages){
  for(end1 in start1:n_ages){
    if ((end1 - start1 + 1) == 5) {
    vac_vec <- rep(0, n_ages)
    vac_vec[start1: end1] <- 1
    all_strategies[[idx]] <- vac_vec
    idx <- idx + 1
    }
  }
}

for(start1 in 1:n_ages){
  for(end1 in start1:n_ages){
      vac_vec <- rep(0, n_ages)
      vac_vec[start1: end1] <- 1
      all_strategies[[idx]] <- vac_vec
      idx <- idx + 1
    }
}

# reverse
for (start1 in seq(n_ages, 1, by = -1)) {
  for (end1 in seq(start1, 1, by = -1)) {
    vac_vec <- rep(0, n_ages)
    # Because start1 >= end1 here, we fill vac_vec[end1 : start1].
    vac_vec[end1:start1] <- 1
    all_strategies[[idx]] <- vac_vec
    idx <- idx + 1
  }
}


# 2 contiguous blocks 
for (start1 in 1:n_ages) {
  for (end1 in start1:n_ages) {
    if(end1 < n_ages){
      for (start2 in (end1+1):n_ages) {
        for (end2 in start2:n_ages) {
          vac_vec <- rep(0, n_ages)
          vac_vec[start1:end1] <- 1
          vac_vec[start2:end2] <- 1
          all_strategies[[idx]] <- vac_vec
          idx <- idx + 1
        }
      }
    }
  }
}

all_strat_mat <- do.call(rbind, all_strategies)
all_strat_mat <- unique(all_strat_mat)
all_strat_df  <- as.data.frame(all_strat_mat)

all_strat_df$strat <- seq_len(nrow(all_strat_df))

all_strat_long <- all_strat_df %>%
  pivot_longer(
    # The first 18 columns are your age groups.
    # If your data frame has exactly 18 columns before adding 'Strategy', you can do:
    cols = -strat,  # i.e., pivot all columns except 'Strategy'
    names_to = "age_group",
    values_to = "targeted"
  ) %>%
  # 'AgeIndex' might be "V1", "V2", etc. Convert that to a numeric 1..18
  mutate(age_group = rep(age_gr[1:18], nrow(all_strat_df)))

all_strat_long$age_group <- factor(all_strat_long$age_group, levels = age_gr_levels)

ggplot(all_strat_long, aes(x = strat, y = age_group, fill = factor(targeted))) +
  geom_tile() +
  # Choose colors: "white" for 0, "grey50" for 1, for example
  scale_fill_manual(values = c("0" = "white", "1" = "grey50")) +
  labs(x = "Strategy", y = "Age group", fill = "Targeted") +
  theme_minimal()


##
prevacc_cases <- sum(summary_cases_pre_all$Median)
prevacc_fatal <- sum(summary_cases_pre_all$fatal)

supplies <- seq(0.1, 1.0, by = 0.05)  

# We'll store all results here:
df_results <- data.frame()
scenario_index <- 1
for (s in supplies) {
  for (i in seq_len(nrow(all_strat_mat))) {
    
    # Extract the i-th strategy
    strategy_vec <- all_strat_mat[i, ]
    
    # Run the model at supply = s (converted to decimal if needed)
    sim_out <- sirv_sim_coverageSwitch(
      T = 52,
      A = 18,
      N = N_2023,
      r = rep(0, 18),
      base_beta = base_beta,
      I0_draw = I0,
      R0  = 1 - exp(-lambda * age_groups),
      rho = rho,
      gamma = gamma,
      delay = 6,
      VE_block = 0.75,
      coverage_threshold = 1,
      target_age = strategy_vec,
      total_coverage = s,
      weekly_delivery_speed = 0.1
    )
    
    sim_df <- as.data.frame.table(sim_out$age_stratified_cases, responseName = "Cases") %>%
      # First, create the basic data frame
      mutate(
        Scenario  = scenario_index,
        AgeGroup  = as.numeric(Var1),
        Week      = as.numeric(Var2)
      ) %>%
      mutate(
        hosp_rate         = rep(hosp, 52),
        hospitalised      = Cases * hosp_rate,
        non_hospitalised  = Cases - hospitalised,
        fatality          = rep(fatal, 52),   # Ensure fatal is correctly replicated
        nh_fatality       = rep(nh_fatal, 52),
        fatal             = (hospitalised * fatality + non_hospitalised * nh_fatality)
      ) %>%
      mutate(
        age_numeric = case_when(
          str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
          TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
        )
      ) %>%
      # Add YLD parameters that are independent of age
      mutate(
        dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Add age-dependent parameters for subacute and chronic outcomes
      mutate(
        subac_prop = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)
        ),
        chr_6m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr6m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr6m, 0.5)
        ),
        chr_12m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr12m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr12m, 0.5)
        ),
        chr_30m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr30m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr30m, 0.5)
        ),
        chr_prop = chr_6m + chr_12m + chr_30m
      ) %>%
      # Calculate YLD estimates
      mutate(
        yld_acute    = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_total    = yld_acute + yld_subacute + yld_chronic
      ) %>%
      # Assign life expectancy based on age_numeric and calculate YLL and DALY metrics
      mutate(
        le_left = case_when(
          age_numeric %in% c(0, 5)    ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10, 15)  ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20, 25)  ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30, 35)  ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40, 45)  ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50, 55)  ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60, 65)  ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70, 75)  ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80, 85)  ~ quantile(le_sample$le_9, 0.5)
        )
      ) %>%
      mutate(
        yll      = fatal * le_left,
        daly_tot = yld_total + yll,
        cum_daly = cumsum(daly_tot)
      )
    
    post_vacc_cases <- sum(sim_df$Cases)
    cases_averted   <- prevacc_cases - post_vacc_cases / sum(N_2023) * 100000
    
    post_vacc_fatal <- sum(sim_df$fatal)
    fatal_averted   <- prevacc_fatal - post_vacc_fatal / sum(N_2023) * 100000

    # Add a row to df_results
    df_results <- rbind(
      df_results,
      data.frame(
        Scenario            = scenario_index,
        Supply              = s * 100,  # store as % if you like
        StrategyID          = i,
        cases_averted       = cases_averted,
        fatal_averted       = fatal_averted
      )
    )
    scenario_index <- scenario_index + 1
  }
}


df_frontier <- df_results %>%
  group_by(Supply) %>%
  summarise(
    max_case_averted  = max(cases_averted),
    max_fatal_averted = max(fatal_averted)
  )

# 2) Plot
ggplot(df_results, aes(x = Supply, y = fatal_averted, color = factor(StrategyID))) +
  #geom_point(position = position_jitter(width = 0.5, height = 0), alpha = 0.4) +  # black by default
  geom_point(alpha = 0.3)+
  geom_line(data = df_frontier, aes(x = Supply, y = max_fatal_averted), color = "red", size = 1) +
  geom_point(data = df_frontier, aes(x = Supply, y = max_fatal_averted), color = "red", size = 2) +
  labs(
    x = "Supply (%)",
    y = "Cases Averted per 1000",
    title = "Frontier of Strategies"
  ) +
  theme_minimal()

ggplot(df_results, aes(x = Supply, y = cases_averted, color = factor(StrategyID))) + 
  geom_point(alpha = 0.3) +
  geom_line(data = df_frontier, aes(x = Supply, y = max_case_averted), 
            color = "red", size = 1) +
  geom_point(data = df_frontier, aes(x = Supply, y = max_case_averted), 
             color = "red", size = 2) +
  labs(
    x = "Supply (%)",
    y = "Cases Averted per 1000",
    title = "Frontier of Strategies"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(71900, 71920))


best_strategies <- df_results %>%
  group_by(Supply) %>%
  filter(fatal_averted == max(fatal_averted)) %>%
  ungroup()

ggplot(df_results, aes(x = Supply, y = cases_averted, color = factor(StrategyID))) +
  geom_point()
