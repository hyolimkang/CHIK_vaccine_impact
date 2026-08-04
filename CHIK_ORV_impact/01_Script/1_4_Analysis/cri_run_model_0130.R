# prevacc data------------------------------------------------------------------
posterior_prevacc <- rstan::extract(fit_prevacc)

base_beta       <- apply(posterior_prevacc$beta, 2, median)       # length T
I0              <- apply(posterior_prevacc$I0, 2, median)           
gamma           <- median(posterior_prevacc$gamma)
rho             <- median(posterior_prevacc$rho)

months <- seq(1, 12, by = 1)
weeks <- seq(1, 12, length.out = 52)
loess_model <- loess(base_beta ~ months, span = 0.5)  # Adjust 'span' for smoothness

# Predict weekly beta values
weekly_base_beta <- predict(loess_model, newdata = weeks)
plot(weekly_base_beta)

N_cri <- as.matrix(cri_pop_transposed$tot_pop)

# costa rica FOI ? (mean FOI)
lambda <- 0.003843 # 0.00279100 0.005134
# age groups
age_groups <- c(mean(0:4),
                mean(5:9),
                mean(10:14),
                mean(15:19),
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

sim_result_novacc <- sirv_sim_coverageSwitch(
  T = 52, 
  A = 18,
  N = N_cri,
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = weekly_base_beta,
  I0_draw = I0,
  R0  = 1 - exp(-lambda * age_groups),
  rho = rho,
  gamma = gamma,
  
  delay = 52,                   # Vaccination starts from Week x
  VE_block = 0,             # Vaccine efficacy
  target_age = c(rep(0,18)),  # Targeting specific age groups
  coverage_threshold = 0,
  total_coverage = 0,
  weekly_delivery_speed = 0
)

age_gr <- rep(c("<5 years",
                "5-9 years",
                "10-14 years",
                "15-19 years",
                "20-24 years",
                "25-29 years",
                "30-34 years",
                "35-39 years",
                "40-44 years",
                "45-49 years",
                "50-54 years",
                "55-59 years",
                "60-64 years",
                "65-69 years",
                "70-74 years",
                "75-79 years",
                "80-84 years",
                "85-89 years"
), 52)
age_gr_levels <- c("<5 years",
                   "5-9 years",
                   "10-14 years",
                   "15-19 years",
                   "20-24 years",
                   "25-29 years",
                   "30-34 years",
                   "35-39 years",
                   "40-44 years",
                   "45-49 years",
                   "50-54 years",
                   "55-59 years",
                   "60-64 years",
                   "65-69 years",
                   "70-74 years",
                   "75-79 years",
                   "80-84 years",
                   "85-89 years")

# pre-vacc
summary_cases_pre <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Median")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )
summary_cases_pre$age_gr <- age_gr
summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)

summary_cases_pre$Scenario <- "Pre-Vaccination"

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_R_pre <- as.data.frame.table(sim_result_novacc$R, responseName = "R")
colnames(summary_R_pre) <- c("AgeGroup", "Week","Median")

summary_R_pre <- summary_R_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week),
    Median_rho = Median * rho
  )

summary_R_pre$Scenario <- "Pre-Vaccination"

summary_R_pre_all <- summary_R_pre %>% group_by(Week) %>% 
  summarise(
    Median_raw = sum(Median),
    Median   = Median_raw * rho
  )

ggplot()+
  geom_line(data = summary_cases_pre, aes(x = Week, y = Median, color = Scenario))+
  facet_wrap(~age_gr)+
  theme_bw()

summary_cases_pre <- summary_cases_pre %>%
  mutate(
    hospitalised = Median * 0.04,
    fatality = rep(fatal, 52),          # Ensure fatal is correctly replicated
    fatal = hospitalised * fatality    # Calculate fatal cases
  ) %>%
  arrange(Week, AgeGroup) %>%          # Order data by Week first, then AgeGroup
  group_by(Week, AgeGroup) %>%         # Group by Week and AgeGroup
  mutate(
    cum_fatal = cumsum(fatal)          # Calculate cumulative fatal cases
  ) %>%
  ungroup()  

summary_cases_pre_age <- summary_cases_pre %>% group_by(AgeGroup) %>%
  summarise(
    Median       = sum(Median),
    hospitalised = sum(hospitalised),
    fatalilty    = first(fatality),
    fatal        = sum(fatal)
  ) %>% mutate(
    Scenario = "Pre-vaccination"
  )

summary_cases_pre_age$age_gr <- age_gr[1:18]
summary_cases_pre_age$age_gr <- factor(summary_cases_pre$age_gr[1:18], levels = age_gr_levels)

summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median),
    hospitalised = sum(hospitalised),
    fatalilty    = first(fatality),
    fatal        = sum(fatal)
  ) %>% mutate(
    Scenario = "Pre-vaccination",
    cum_fatal    = cumsum(fatal),
    cum_hosp     = cumsum(hospitalised)
  )

target_age_list <- list(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1))

scenario_result <- list()

for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  sim_result <- sirv_sim_coverageSwitch(
    T = 52,
    A = 18,
    N = N_cri,
    r = rep(0, 18),
    base_beta = weekly_base_beta,
    I0_draw = I0,
    R0  = 1 - exp(-lambda * age_groups),
    rho = rho,
    gamma = gamma,
    delay = 2,
    VE_block = 1,
    coverage_threshold = 1,
    target_age = target,
    total_coverage = 0.5,
    weekly_delivery_speed = 0.1
  )
  
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
    mutate(
      Scenario            = scenario_index,
      AgeGroup            = as.numeric(Var1),
      Week                = as.numeric(Var2)
    )
  
  scenario_result[[scenario_index]] <- list(sim_result = sim_result, 
                                            sim_df     = sim_df)
}

scenario_data <- lapply(seq_along(scenario_result), function(idx){
  list <- scenario_result[[idx]]
  data <- list$sim_result
  return(data)
})

summary_list <- lapply(seq_along(scenario_result), function(idx){
  list <- scenario_result[[idx]]
  df <- list$sim_df
  return(df)
})

final_summ <- lapply(summary_list, function(df) {
  df %>%
    group_by(AgeGroup) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(Cases, na.rm = TRUE),  # Summarize MedianCases
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_cases_pre_age$Median,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100
    ) 
})

summary_week <- lapply(seq_along(summary_list), function(i) {
  summary_list[[i]] %>% 
    mutate(Scenario = paste0("Scenario_", i),
           pre_vacc = summary_cases_pre$Median) %>% group_by(Week, Scenario) %>%
    summarise(post_cases = sum(Cases),
              pre_cases  = sum(pre_vacc),
              diff       = pre_cases - post_cases,
              Scenario   = first(Scenario)
    ) 
  
})

summary_week_df <- do.call(rbind, summary_week)

global_impact_week <- summary_week_df %>%
  group_by(Scenario) %>%
  summarise(
    total_post_cases = sum(post_cases),  # Sum cumulative fatal cases for each scenario and type
    total_pre_case   = sum(pre_cases),
    .groups = "drop"
  ) %>%
  mutate(
    diff = `total_pre_case` - `total_post_cases`,
    impact = diff / `total_pre_case`  * 100
  )

annotation_text <- paste0(
  global_impact_week$Scenario, ": ", 
  round(global_impact_week$impact, 1), "%",
  collapse = "\n"
)

total_duration <- scenario_data[[1]]$total_vaccine_duration_age[1]

ggplot(summary_week_df) +
  geom_ribbon(data = summary_week_df[summary_week_df$Week <= total_duration, ],
              aes(x = Week, ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.3) +
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases")) +
  geom_line(aes(x = Week, y = pre_cases, color = Scenario, linetype = "Pre Cases"), color = "black") +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  #scale_color_brewer(type = "qual", palette = 1)+
  labs(color = "Scenario", linetype = "Type") +
  labs(
    title = "Coverage:40%, delivery speed: 10%",  # Add title
    color = "Scenario", 
    linetype = "Type", 
    x = "Week", 
    y = "Cases"
  )+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 10)
  ) + 
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1,
           label = annotation_text,  # from the paste step
           color = "black", size = 3)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5))


# by coverage and delay (heatmap generation) -----------------------------------------------------


supply <- seq(0.1, 1, by = 0.05)
delay_steps  <- seq(1, 52, by = 1)

debug_impact_data <- list()
# Iterate over scenarios, supply rates, and delays
for (scenario_index in seq_along(target_age_list)) {
  
  target <- target_age_list[[scenario_index]]
  
  for (supply_rate in supply) { 
    for (delay_step in delay_steps) {
      # Run the simulation
      sim_result <- sirv_sim_coverageSwitch(
        T = 52,
        A = 18,
        N = N_cri,
        r = rep(0, 18),
        base_beta = weekly_base_beta,
        I0_draw = I0,
        R0  = 1 - exp(-lambda * age_groups),
        rho = rho,
        gamma = gamma,
        delay = delay_step,
        VE_block = 1,
        coverage_threshold = 1,
        target_age = target,
        total_coverage = supply_rate,
        weekly_delivery_speed = 0.05
      )
      
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
        mutate(
          Scenario = scenario_index,
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      pre_vacc_cases <- summary_cases_pre %>%
        group_by(Week) %>%
        summarise(tot_pre_cases = sum(Median, na.rm = TRUE), .groups = "drop")
      
      post_vacc_cases <- sim_df %>%
        group_by(Week) %>%
        summarise(tot_post_cases = sum(Cases, na.rm = TRUE), .groups = "drop")
      
      impact_data <- pre_vacc_cases %>%
        left_join(post_vacc_cases, by = "Week") %>%
        summarise(
          tot_pre_cases = sum(tot_pre_cases, na.rm = TRUE),
          tot_post_cases = sum(tot_post_cases, na.rm = TRUE),
          tot_diff = tot_pre_cases - tot_post_cases,
          impact = tot_diff / tot_pre_cases * 100,
          .groups = "drop"
        )
      
      debug_impact_data[[paste0("Scenario_", scenario_index, "_Supply_", supply_rate, "_Delay_", delay_step)]] <- list(
        pre_vacc_cases = pre_vacc_cases,
        post_vacc_cases = post_vacc_cases,
        impact_data = impact_data
      )
      
    }
  }
}

impact_summ <- do.call(rbind, lapply(names(debug_impact_data), function(idx) {
  list <- debug_impact_data[[idx]]  # Access the element using its name
  df <- list$impact_data           
  df$scenario_id <- idx          
  return(df)                        
}))

impact_summ <- impact_summ %>%
  mutate(
    Scenario = sub("_Supply.*", "", scenario_id),  # Extract only the Scenario part
    Supply = as.numeric(sub("Supply_", "", stringr::str_extract(scenario_id, "Supply_\\d+(\\.\\d+)?"))),
    Delay = as.numeric(sub("Delay_", "", stringr::str_extract(scenario_id, "Delay_\\d+(\\.\\d+)?")))
  )

p <- ggplot(impact_summ, aes(x = Supply, y = Delay, fill = impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% cases averted") +  # Add color legend breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~Scenario)+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig11.jpg", plot = p,
       width = 8, height = 3)
