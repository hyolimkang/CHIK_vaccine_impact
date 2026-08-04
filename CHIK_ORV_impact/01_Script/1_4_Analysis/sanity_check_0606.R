preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups, 
                                N = N_ceara$Ceará, 
                                region = "Ceará",
                                observed = observed_ce
)

pre_results_ce_ui <- summarise_presim_ui(sim_result     = preui_ce, 
                                         observed         = observed_ce,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Ceará"
)

prevacc_ui_ce <- pre_results_ce_ui$summary_cases_pre
prevacc_ui_ce_all <- pre_results_ce_ui$summary_cases_pre_all
pre_summary_age_ce <- pre_results_ce_ui$summary_cases_pre_age

postsim_ce_ui <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_ce,
  observed          = observed_ce,
  N                 = N_ceara$Ceará,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Ceará",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_ce,
  posterior         = posterior_ce
)

postsim_all_ce_ui <- postsim_all_ui(
  scenario_result      = postsim_ce_ui,
  observed             = observed_ce,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_ce,
  pre_summary_cases     = prevacc_ui_ce,
  pre_summary_cases_all     = prevacc_ui_ce_all,
  region                = "Ceará"
)

epi_graph(postsim_all_ce_ui, observed_ce)


vacc_alloc_ce <- vacc_allocation(postsim_all_ce_ui, observed_ce, "Ceará")


vacc_alloc_graph(vacc_alloc_ce)

nnv_ce <- nnv_list(vacc_alloc_ce, postsim_all_ce_ui, N_ceara$Ceará, "Ceará", observed_ce)


df_nnv_ce <- nnv_ce$final_summ_df

ggplot(df_nnv_ce, aes(x=AgeGroup, y=tot_vacc_prop)) + 
  geom_bar(stat = "identity")+
  facet_wrap(~scenario)

vacc_ce <- df_nnv_ce %>% group_by(scenario) %>% summarise(
  tot_vacc = sum(tot_vacc)
)

vacc_impact_ce <- df_nnv_ce %>% group_by(scenario, target) %>% summarise(
  diff  = sum(diff),
  diff_fatal = sum(diff_fatal),
  diff_daly = sum(diff_daly),
  tot_vacc= sum(tot_vacc)
)

vacc_impact_ce <- vacc_impact_ce %>% mutate(
  overall_case_1m  = (diff / tot_vacc) * 1e6,
  overall_death_1m = (diff_fatal / tot_vacc) * 1e6,
  overall_daly_1m  = (diff_daly / tot_vacc) * 1e6
)

vacc_impact_ce <- vacc_impact_ce %>%
  filter(is.finite(overall_case_1m)) %>%        
  mutate(
    target = factor(
      target,
      levels = c("<20 years", "20-59 years", ">60 years")    
    )
  )

ggplot(vacc_impact_ce, aes(x = target, y = overall_death_1m, fill = scenario)) + 
  geom_bar(stat = "identity")+
  ylab("vaccine impact ratio to avert death (per 1 FVP)")+
  facet_wrap(~scenario)

ggplot(vacc_impact_ce, aes(x = target, y = overall_case_1m, fill = scenario)) + 
  geom_bar(stat = "identity")+
  ylab("vaccine impact ratio to avert case (per 1 FVP)")+
  facet_wrap(~scenario)

ggplot(vacc_impact_ce, aes(x = target, y = overall_daly_1m, fill = scenario)) + 
  geom_bar(stat = "identity")+
  ylab("vaccine impact ratio to avert daly (per 1 FVP)")+
  facet_wrap(~scenario)


# wasted dose track 

# 시나리오 1
T = 52
raw_s1   <- as.data.frame(postsim_all_ce_ui$scenario_result[[1]]$sim_out$raw_allocation_age)
wasted_s1 <- as.data.frame(postsim_all_ce_ui$scenario_result[[1]]$sim_out$wasted_dose)

colnames(raw_s1)   <- paste0("W", seq_len(T))
colnames(wasted_s1)<- paste0("W", seq_len(T))

raw_s1 <- raw_s1 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S1")

wasted_s1 <- wasted_s1 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S1")

# ───────────────────────────────────────────────────────────────────
# 시나리오 2
raw_s2   <- as.data.frame(postsim_all_ce_ui$scenario_result[[2]]$sim_out$raw_allocation_age)
wasted_s2 <- as.data.frame(postsim_all_ce_ui$scenario_result[[2]]$sim_out$wasted_dose)

colnames(raw_s2)   <- paste0("W", seq_len(T))
colnames(wasted_s2)<- paste0("W", seq_len(T))

raw_s2 <- raw_s2 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S2")

wasted_s2 <- wasted_s2 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S2")

# ───────────────────────────────────────────────────────────────────
# 시나리오 3
raw_s3   <- as.data.frame(postsim_all_ce_ui$scenario_result[[3]]$sim_out$raw_allocation_age)
wasted_s3 <- as.data.frame(postsim_all_ce_ui$scenario_result[[3]]$sim_out$wasted_dose)

colnames(raw_s3)   <- paste0("W", seq_len(T))
colnames(wasted_s3)<- paste0("W", seq_len(T))

raw_s3 <- raw_s3 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S3")

wasted_s3 <- wasted_s3 %>%
  mutate(AgeGroup = age_gr_levels) %>%
  select(AgeGroup, everything()) %>%
  mutate(type = "S3")

# ───────────────────────────────────────────────────────────────────
# 합치기
raw_combined   <- bind_rows(raw_s1, raw_s2, raw_s3)
wasted_combined<- bind_rows(wasted_s1, wasted_s2, wasted_s3)


# 1) raw_allocation_age long 포맷으로 변환
raw_long <- raw_combined %>%
  pivot_longer(
    cols      = starts_with("W"),  # W1, W2, ... W<T>
    names_to  = "metric",
    values_to = "raw_alloc"
  ) %>%
  mutate(
    week = as.integer(str_remove(metric, "^W"))
  ) %>%
  select(AgeGroup, type, week, raw_alloc)

# 2) wasted_dose long 포맷으로 변환
wasted_long2 <- wasted_combined %>%
  pivot_longer(
    cols      = starts_with("W"),
    names_to  = "metric",
    values_to = "wasted"
  ) %>%
  mutate(
    week = as.integer(str_remove(metric, "^W"))
  ) %>%
  select(AgeGroup, type, week, wasted)

# 3) raw와 wasted를 합쳐서 “effective” 계산
alloc_waste <- raw_long %>%
  left_join(wasted_long2, by = c("AgeGroup", "type", "week")) %>%
  mutate(
    effective = raw_alloc - wasted  # “할당된 raw 도즈 중 낭비되지 않고 면역 후보가 된 양”
  )

# 1) raw_allocation_age long 포맷으로 변환
raw_long <- raw_combined %>%
  pivot_longer(
    cols      = starts_with("W"),  # W1, W2, ... W<T>
    names_to  = "metric",
    values_to = "raw_alloc"
  ) %>%
  mutate(
    week = as.integer(str_remove(metric, "^W"))
  ) %>%
  select(AgeGroup, type, week, raw_alloc)

# 2) wasted_dose long 포맷으로 변환
wasted_long2 <- wasted_combined %>%
  pivot_longer(
    cols      = starts_with("W"),
    names_to  = "metric",
    values_to = "wasted"
  ) %>%
  mutate(
    week = as.integer(str_remove(metric, "^W"))
  ) %>%
  select(AgeGroup, type, week, wasted)

# 3) raw와 wasted를 합쳐서 “effective” 계산
alloc_waste <- raw_long %>%
  left_join(wasted_long2, by = c("AgeGroup", "type", "week")) %>%
  mutate(
    effective = raw_alloc - wasted  # “할당된 raw 도즈 중 낭비되지 않고 면역 후보가 된 양”
  )

df_summary2 <- alloc_waste %>%
  group_by(type, week) %>%
  summarise(
    total_raw_alloc  = sum(raw_alloc, na.rm = TRUE),
    total_wasted     = sum(wasted,    na.rm = TRUE),
    total_effective  = sum(effective, na.rm = TRUE),
    .groups = "drop"
  )


# 1) 스택형 그래프용 long 포맷으로 변환
df_stack <- df_summary2 %>%
  pivot_longer(
    cols      = c(total_wasted, total_effective),
    names_to  = "status",            # “ineffective” vs “effective”
    values_to = "count"
  ) %>%
  mutate(
    status = dplyr::recode(status,
                    "total_wasted"    = "ineffective",
                    "total_effective" = "effective")
  )
  #%>%
  #mutate(
  #  AgeGroup = factor(AgeGroup, levels = age_gr_levels)
  #)

# 2) 시각화: 시나리오(type)별·주차(week)별 스택형 막대
ggplot(df_stack, aes(x = week, y = count, fill = status)) +
  geom_col(width = 0.8) +
  facet_wrap(~type) +
  scale_fill_manual(
    values = c("ineffective" = "#D55E00",  # 짙은 주황 (낭비)
               "effective"   = "#0072B2")  # 짙은 파랑 (유효)
  ) +
  scale_x_continuous(
    breaks = seq(1, max(df_stack$week), by = 4)
  ) + ylab("Total doses given")+
  scale_x_continuous(labels = comma)


### age group (scenario specific) epi curves -----------------------------------
summary_list_ce <- postsim_all_ce_ui$summary_list

# assume each summary_list[[i]] has columns: Week, AgeGroup, Cases, fatal, daly_tot
summary_age_week <- lapply(seq_along(summary_list_ce), function(i) {
  summary_list_ce[[i]] %>%
    mutate(
      Scenario = paste0("Scenario_", i),
      pre_vacc_cases = prevacc_ui_ce$Median,      # if pre_summary_cases also age-specific
      pre_vacc_fatal = prevacc_ui_ce$fatal,      # same length/order as AgeGroup
      pre_vacc_daly  = prevacc_ui_ce$daly_tot
    ) %>%
    group_by(Week, AgeGroup, Scenario) %>%
    summarise(
      post_cases   = sum(Cases,          na.rm = TRUE),
      pre_cases    = sum(pre_vacc_cases, na.rm = TRUE),
      diff_cases   = pre_cases - post_cases,
      post_fatal   = sum(fatal,          na.rm = TRUE),
      pre_fatal    = sum(pre_vacc_fatal, na.rm = TRUE),
      diff_fatal   = pre_fatal - post_fatal,
      post_daly    = sum(daly_tot,       na.rm = TRUE),
      pre_daly     = sum(pre_vacc_daly,  na.rm = TRUE),
      diff_daly    = pre_daly - post_daly,
      .groups = "drop"
    )
})
summary_age_week_df <- bind_rows(summary_age_week)

age_targets <- list(
  Scenario_1 = 1:4,
  Scenario_2 = 5:12,
  Scenario_3 = 13:18
)

# 2) target 플래그 추가
summary_age_week_df <- summary_age_week_df %>%
  mutate(
    target = map2_lgl(Scenario, AgeGroup,
                      ~ .y %in% age_targets[[.x]])
  )
summary_target_week <- summary_age_week_df %>%
  filter(target) %>%
  group_by(Week, Scenario) %>%
  summarise(
    post_cases = sum(post_cases, na.rm = TRUE),
    pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  )


summary_target_week %>%
  pivot_longer(
    cols = c(post_cases, pre_cases),
    names_to  = "Type",
    values_to = "Cases"
  ) %>%
  ggplot(aes(x = Week, y = Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario) +
  labs(
    x = "Week",
    y = "Cases"
  ) +
  theme_pubclean()+
  theme(legend.position = "bottom")


summary_age_week_df %>% 
  filter(target == FALSE) %>% 
  group_by(Week, Scenario) %>% 
  summarise(
    post_cases = sum(post_cases),
    pre_cases  = sum(pre_cases),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = c(post_cases, pre_cases),
    names_to  = "Type",
    values_to = "Cases"
  ) %>% 
  ggplot(aes(Week, Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario) +
  theme_pubclean()+
  theme(legend.position = "bottom")

summary_age_week_df %>% 
  filter(target == TRUE) %>% 
  group_by(Week, Scenario) %>% 
  summarise(
    post_cases = sum(post_cases),
    pre_cases  = sum(pre_cases),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = c(post_cases, pre_cases),
    names_to  = "Type",
    values_to = "Cases"
  ) %>% 
  ggplot(aes(Week, Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario) +
  theme_pubclean()+
  theme(legend.position = "bottom")

summary_age_week_df %>% 
  filter(target == TRUE) %>%
  group_by(Week, Scenario) %>% 
  summarise(
    post_fatal = sum(post_fatal),
    pre_fatal  = sum(pre_fatal),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = c(post_fatal, pre_fatal),
    names_to  = "Type",
    values_to = "Cases"
  ) %>% 
  ggplot(aes(Week, Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario) +
  theme_pubclean()+
  theme(legend.position = "bottom")


post_fatal <- summary_age_week_df %>% 
  filter(target == TRUE) %>%
  group_by(Scenario) %>% 
  summarise(
    post_fatal = sum(post_fatal),
    pre_fatal  = sum(pre_fatal),
    diff = pre_fatal - post_fatal,
    reduction = (diff) /pre_fatal * 100,
    .groups = "drop"
  )
### only infections 
summarise_presim_ui_inf <- function(
    sim_result,        # output from simulate_pre_ui_age
    observed,
    age_gr_levels = c("<5 years",
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
                      "85-89 years"),
    lhs_sample_young,  # data for younger ages (for YLD calculation)
    lhs_old,           # data for older ages (for YLD calculation)
    le_sample,         # life-expectancy sample
    hosp,              # hospitalization rate vector or single numeric
    fatal,             # fatality rate for hospitalized
    nh_fatal,           # fatality rate for non-hospitalized
    region
) {
  
  # Set T dynamically
  T =  nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  age_gr <- rep(default_age_vector, each = T)
  
  # --- Create Age-Stratified Summary Data Frame ---
  median_df <- as.data.frame.table(sim_result$median_by_age_rawinf, responseName = "Median")
  colnames(median_df) <- c("AgeGroup", "Week", "Median")
  
  median_df <- median_df %>%
    mutate(
      AgeGroup = as.numeric(AgeGroup),
      Week = as.numeric(Week)
    )
  
  low95_df <- as.data.frame.table(sim_result$low95_by_age_rawinf, responseName = "low95")
  hi95_df  <- as.data.frame.table(sim_result$hi95_by_age_rawinf, responseName = "hi95")
  
  summary_cases_pre <- median_df %>%
    mutate(
      low95 = low95_df$low95,
      hi95  = hi95_df$hi95,
      age_gr = default_age_vector[AgeGroup]  # assign provided age group labels
    )
  summary_cases_pre$age_gr <- factor(summary_cases_pre$age_gr, levels = age_gr_levels)
  
  # --- Weekly Summary ---
  weekly_df <- data.frame(
    Week = 1:length(sim_result$weekly_cases_median_rawinf),
    weekly_median = sim_result$weekly_cases_median_rawinf,
    weekly_low95  = sim_result$weekly_cases_low95_rawinf,
    weekly_hi95   = sim_result$weekly_cases_hi95_rawinf
  )
  
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(Median = sum(Median)) %>%
    mutate(Scenario = "Pre-vaccination") %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      lo95  = weekly_low95,
      hi95   = weekly_hi95
    )
  
  
  # --- Calculate Hospitalizations, Fatalities, and Cumulative Totals ---
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      hosp_rate = rep(hosp, T),
      hospitalised = Median * hosp_rate,
      hospitalised_lo = low95 * hosp_rate,
      hospitalised_hi = hi95 * hosp_rate,
      non_hospitalised = Median - hospitalised,
      non_hospitalised_lo = low95 - hospitalised_hi,  # conservative: lower bound for non-hosp
      non_hospitalised_hi = hi95 - hospitalised_lo,
      fatality = rep(fatal, T),          # Ensure fatal is correctly replicated
      nh_fatality = rep(nh_fatal, T),
      fatal = (hospitalised * fatality +   
                 non_hospitalised * nh_fatality),
      fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
      fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality)
    )%>%
    arrange(Week, AgeGroup) %>%          # Order data by Week first, then AgeGroup
    group_by(Week, AgeGroup) %>%         # Group by Week and AgeGroup
    mutate(
      cum_fatal = cumsum(fatal),          # Calculate cumulative fatal cases
      cum_hosp  = cumsum(hospitalised)
    ) %>%
    ungroup()
  
  # --- Summarize by AgeGroup ---
  summary_cases_pre_age <- summary_cases_pre %>% group_by(AgeGroup) %>%
    summarise(
      Median       = sum(Median),
      hospitalised = sum(hospitalised),
      hospitalised_lo = sum(hospitalised_lo),
      hospitalised_hi = sum(hospitalised_hi),
      fatalilty    = first(fatality),
      nh_fatality  = first(nh_fatality),
      fatal        = sum(fatal)
    ) %>% mutate(
      Scenario = "Pre-vaccination"
    )
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  # --- Summarize by Week for Entire Population ---
  summary_cases_pre_all <- summary_cases_pre %>%
    group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      hospitalised = sum(hospitalised),
      fatality     = first(fatality),
      fatal        = sum(fatal),
      fatal_lo     = sum(fatal_lo),
      fatal_hi     = sum(fatal_hi)
    ) %>%
    mutate(
      Scenario     = "Pre-vaccination",
      cum_fatal    = cumsum(fatal),
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      cum_hosp     = cumsum(hospitalised)
    ) %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      pre_weekly_low95  = weekly_low95,
      pre_weekly_hi95   = weekly_hi95
    )
  
  # --- DALY Calculations ---
  summary_cases_pre <- summary_cases_pre %>% mutate(
    # YLD parameters
    dw_hosp      = quantile(lhs_sample_young$dw_hosp, 0.5),
    dur_acute    = quantile(lhs_sample_young$dur_acute, 0.5),
    dw_nonhosp   = quantile(lhs_sample_young$dw_nonhosp, 0.5),
    dw_chronic   = quantile(lhs_sample_young$dw_chronic, 0.5),
    dur_chronic  = quantile(lhs_sample_young$dur_chronic, 0.5),
    dw_subacute  = quantile(lhs_sample_young$dw_subac, 0.5),
    dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5),
    subac_prop   = quantile(lhs_sample_young$subac, 0.5),
    chr_6m       = quantile(lhs_sample_young$chr6m, 0.5),
    chr_12m      = quantile(lhs_sample_young$chr12m, 0.5),
    chr_30m      = quantile(lhs_sample_young$chr30m, 0.5),
    chr_prop     = chr_6m + chr_12m + chr_30m 
  )
  summary_cases_pre <- summary_cases_pre %>%
    mutate(
      age_numeric = case_when(
        str_detect(age_gr, "<5") ~ 0,  # Set age_numeric to 0 for "<5 years"
        TRUE ~ as.numeric(str_extract(age_gr, "^\\d+")) # Extract first number otherwise
      ),  # Extract numbers from age_gr
      subac_prop  = case_when(
        age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),   # Replace with the desired value for <39
        age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)   # Replace with the desired value for ≥40
      ),
      chr_6m  = case_when(
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
      chr_prop = chr_6m + chr_12m + chr_30m  # Recalculate chr_prop based on updated values
    ) %>% mutate(
      # YLD estimates
      yld_acute    = (hospitalised * dw_hosp * dur_acute) +
        (non_hospitalised * dw_nonhosp * dur_acute),
      yld_acute_lo = (hospitalised_lo * dw_hosp * dur_acute) +
        ((low95 - hospitalised_hi) * dw_nonhosp * dur_acute),
      yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
        ((hi95 - hospitalised_lo) * dw_nonhosp * dur_acute),
      
      yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) + 
        (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
        ((low95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
      yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
      
      yld_chronic  = (hospitalised * chr_prop * dw_chronic * dur_chronic) + 
        (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
        ((low95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
      yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
        ((hi95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
      
      yld_total    = yld_acute + yld_subacute + yld_chronic,
      yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
      yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
    ) %>% mutate(
      # le left 
      le_left = case_when(
        age_numeric %in% c(0, 5)   ~ quantile(le_sample$le_1, 0.5),
        age_numeric %in% c(10, 15) ~ quantile(le_sample$le_2, 0.5),
        age_numeric %in% c(20, 25) ~ quantile(le_sample$le_3, 0.5),
        age_numeric %in% c(30, 35) ~ quantile(le_sample$le_4, 0.5),
        age_numeric %in% c(40, 45) ~ quantile(le_sample$le_5, 0.5),
        age_numeric %in% c(50, 55) ~ quantile(le_sample$le_6, 0.5),
        age_numeric %in% c(60, 65) ~ quantile(le_sample$le_7, 0.5),
        age_numeric %in% c(70, 75) ~ quantile(le_sample$le_8, 0.5),
        age_numeric %in% c(80, 85) ~ quantile(le_sample$le_9, 0.5)
      ) 
    ) %>% mutate(
      # YLL  
      yll = fatal * le_left,
      yll_lo = fatal_lo * le_left,
      yll_hi = fatal_hi * le_left,
      # total DALY
      daly_tot = yld_total + yll,
      daly_tot_lo = yld_total_lo + yll_lo,
      daly_tot_hi = yld_total_hi + yll_hi,
      cum_daly = cumsum(daly_tot)
    )
  
  # --- Final Summaries ---
  summary_cases_pre_age <- summary_cases_pre %>% 
    group_by(AgeGroup) %>%
    summarise(
      Median           = sum(Median, na.rm = TRUE),
      low95            = sum(low95, na.rm = TRUE),
      hi95             = sum(hi95, na.rm = TRUE),
      hospitalised     = sum(hospitalised, na.rm = TRUE),
      hospitalised_lo  = sum(hospitalised_lo, na.rm = TRUE),
      hospitalised_hi  = sum(hospitalised_hi, na.rm = TRUE),
      yld_acute        = sum(yld_acute, na.rm = TRUE),
      yld_acute_lo     = sum(yld_acute_lo, na.rm = TRUE),
      yld_acute_hi     = sum(yld_acute_hi, na.rm = TRUE),
      yld_subacute     = sum(yld_subacute, na.rm = TRUE),
      yld_subacute_lo  = sum(yld_subacute_lo, na.rm = TRUE),
      yld_subacute_hi  = sum(yld_subacute_hi, na.rm = TRUE),
      yld_chronic      = sum(yld_chronic, na.rm = TRUE),
      yld_chronic_lo   = sum(yld_chronic_lo, na.rm = TRUE),
      yld_chronic_hi   = sum(yld_chronic_hi, na.rm = TRUE),
      yld_total        = sum(yld_total, na.rm = TRUE),
      yld_total_lo     = sum(yld_total_lo, na.rm = TRUE),
      yld_total_hi     = sum(yld_total_hi, na.rm = TRUE),
      yll              = sum(yll, na.rm = TRUE),
      yll_lo           = sum(yll_lo, na.rm = TRUE),
      yll_hi           = sum(yll_hi, na.rm = TRUE),
      daly_tot         = sum(daly_tot, na.rm = TRUE),
      daly_tot_lo      = sum(daly_tot_lo, na.rm = TRUE),
      daly_tot_hi      = sum(daly_tot_hi, na.rm = TRUE),
      fatal            = sum(fatal, na.rm = TRUE),
      fatal_lo         = sum(fatal_lo, na.rm = TRUE),
      fatal_hi         = sum(fatal_hi, na.rm = TRUE)
    ) %>% 
    mutate(
      Scenario = "Pre-vaccination"
    ) %>%
    ungroup()
  
  summary_cases_pre_age$age_gr <- age_gr[1:18]
  summary_cases_pre_age$age_gr <- factor(summary_cases_pre_age$age_gr, levels = age_gr_levels)
  
  summary_cases_pre_all <- summary_cases_pre %>% group_by(Week) %>%
    summarise(
      Median       = sum(Median),
      yld_acute    = sum(yld_acute),
      yld_subacute = sum(yld_subacute),
      yld_chronic  = sum(yld_chronic),
      yld_total    = sum(yld_total),
      yll          = sum(yll),
      daly_tot     = sum(daly_tot),
      hospitalised = sum(hospitalised),
      fatalilty    = first(fatality),
      fatal        = sum(fatal),
      fatal_lo     = sum(fatal_lo),
      fatal_hi     = sum(fatal_hi)
    ) %>% mutate(
      Scenario   = "Pre-vaccination",
      cum_fatal    = cumsum(fatal),
      cum_fatal_lo = cumsum(fatal_lo),
      cum_fatal_hi = cumsum(fatal_hi),
      cum_hosp   = cumsum(hospitalised),
      cum_yld_acute = cumsum(yld_acute),
      cum_yld_subacute = cumsum(yld_subacute),
      cum_yld_chronic = cumsum(yld_chronic),
      cum_yld_total   = cumsum(yld_total),
      cum_yll         = cumsum(yll),
      cum_daly_tot    = cumsum(daly_tot)
    ) %>%
    ungroup() %>%
    left_join(weekly_df, by = "Week") %>%
    dplyr::rename(
      pre_weekly_median = weekly_median,
      lo95  = weekly_low95,
      hi95   = weekly_hi95
    )
  
  summary_cases_pre$region     <- region
  summary_cases_pre_all$region <- region
  summary_cases_pre_age$region <- region
  
  return(list(
    summary_cases_pre = summary_cases_pre,
    summary_cases_pre_all = summary_cases_pre_all,
    summary_cases_pre_age = summary_cases_pre_age
  ))
}

run_simulation_scenarios_ui_inf <- function(target_age_list, 
                                        observed,
                                        N, bra_foi_state_summ, age_groups, region_name,
                                        hosp, fatal, nh_fatal,
                                        lhs_sample_young, lhs_old, le_sample,
                                        age_gr_levels,
                                        prevacc_ui = NULL,
                                        posterior
) {
  # target_age_list: list of 0/1 vectors for each scenario
  # observed: observed data (used to set T)
  # N: population vector (length = number of age groups)
  # bra_foi_state_summ: data frame with avg_foi and NAME_1 for FOI calculation
  # age_groups: numeric vector of age group midpoints used in FOI calculation
  # hosp, fatal, nh_fatal: hospitalization, fatality rates (numeric or vector)
  # lhs_sample_young, lhs_old: samples for YLD parameters for young/old ages
  # le_sample: life-expectancy sample (for YLL calculation)
  # age_gr_levels: factor levels for age group labels
  # prevacc_ui: (optional) pre-vaccination UI data frame with columns: AgeGroup, Week, Median, low95, hi95
  # posterior: list with posterior draws (e.g., base_beta, I0, gamma, rho, etc.)
  
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  n_draws <- length(posterior$gamma)
  
  T <- nrow(observed)
  
  # Define default age groups (assumed 18 groups)
  default_age_vector <- c("<5 years", "5-9 years", "10-14 years", "15-19 years",
                          "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                          "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                          "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                          "80-84 years", "85-89 years")
  
  # Create an age group vector for the age-by-week summary:
  # Here, each default age label is repeated T times.
  age_gr <- rep(default_age_vector, T)
  
  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    
    for (d in 1:n_draws) {
      # Extract d-th draw parameters
      base_beta_draw <- posterior$base_beta[d, ]   # vector length = T
      I0_draw        <- posterior$I0[d, ]          # vector length = 18 (A)
      gamma_draw     <- posterior$gamma[d]
      rho_draw       <- posterior$rho[d]
      
      sim_out <- sirv_sim_coverageSwitch(
        T = T,
        A = 18,
        N = N,
        r = rep(0, 18),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        R0 = 1 - exp(- bra_foi_state_summ$avg_foi[bra_foi_state_summ$NAME_1 == region_name] * age_groups),
        #R0 = 1 - exp(- 0.1 * age_groups),
        rho = rho_draw,
        gamma = gamma_draw,
        delay = 1,
        VE_block = 0.989,
        #coverage_threshold = 1,
        target_age = target,
        #total_coverage = 1e6 * 0.2 / sum(N),
        total_coverage = 0.12,
        #cap_supply =  1e6,
        weekly_delivery_speed = 0.1
      )
      
      draw_results_raw_inf[[d]]  <- sim_out$age_stratified_cases_raw  # A x T matrix
      draw_results_raw_symp[[d]] <- sim_out$true_symptomatic  # A x T matrix
      
    }
    
    # Combine draw results into a 3D array [18, T, n_draws]
    age_array_raw_inf  <- array(unlist(draw_results_raw_inf), dim = c(18, T, n_draws))
    age_array_raw_symp <- array(unlist(draw_results_raw_symp), dim = c(18, T, n_draws))
    
    # Compute quantiles for each age and week over draws:
    median_by_age_rawinf <- apply(age_array_raw_inf, c(1, 2), median)
    low95_by_age_rawinf  <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawinf   <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.975)
    
    median_by_age_rawsymp <- apply(age_array_raw_symp, c(1, 2), median)
    low95_by_age_rawsymp  <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawsymp   <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.975)
    
    # Weekly totals (sum over ages) for each draw:
    weekly_totals_rawinf       <- apply(age_array_raw_inf, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawinf <- apply(weekly_totals_rawinf, 1, median)
    weekly_cases_low95_rawinf  <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawinf   <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.975)
    
    weekly_totals_rawsymp       <- apply(age_array_raw_symp, c(2, 3), sum)  # [T, n_draws]
    weekly_cases_median_rawsymp <- apply(weekly_totals_rawsymp, 1, median)
    weekly_cases_low95_rawsymp  <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawsymp   <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.975)
    
    # Create simulation data frame from median_by_age: (raw symptomatic cases)
    sim_df <- as.data.frame.table(median_by_age_rawinf, responseName = "Cases")
    
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      mutate(
        Scenario = s,
        # Convert AgeGroup from factor to numeric via as.character
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    
    sim_df <- sim_df %>%
      mutate(region_full = region_name)
    
    # Attach UI for cases:
    low95_df <- as.data.frame.table(low95_by_age_rawinf, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age_rawinf, responseName = "post_hi95")
    
    sim_df <- sim_df %>%
      mutate(
        post_low95        = low95_df$post_low95,
        post_hi95         = hi95_df$post_hi95
      )
    
    # Weekly summary data frame:
    weekly_df <- data.frame(
      Week = 1:T,
      weekly_median = weekly_cases_median_rawinf,
      weekly_low95 = weekly_cases_low95_rawinf,
      weekly_hi95 = weekly_cases_hi95_rawinf
    )
    
    # Attach age group labels using our default vector:
    sim_df <- sim_df %>%
      mutate(
        age_gr = rep(default_age_vector, T)
      )
    sim_df$age_gr <- factor(sim_df$age_gr, levels = age_gr_levels)
    
    # Calculate hospitalisations, fatalities, and uncertainty bounds:
    sim_df <- sim_df %>%
      mutate(
        hosp_rate = rep(hosp, T),
        hospitalised = Cases * hosp_rate,
        hospitalised_lo = post_low95 * hosp_rate,
        hospitalised_hi = post_hi95 * hosp_rate,
        non_hospitalised = Cases - hospitalised,
        non_hospitalised_lo = post_low95 - hospitalised_hi,   # conservative
        non_hospitalised_hi = post_hi95 - hospitalised_lo,
        fatality = rep(fatal, T),
        nh_fatality = rep(nh_fatal, T),
        fatal = (hospitalised * fatality + non_hospitalised * nh_fatality),
        fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
        fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality),
        
      ) %>%
      arrange(Week, AgeGroup) %>%
      group_by(Week, AgeGroup) %>%
      mutate(
        cum_fatal = cumsum(fatal),
        cum_hosp = cumsum(hospitalised)
      ) %>%
      ungroup()
    
    # Calculate age_numeric from age_gr for age-dependent parameters.
    sim_df <- sim_df %>%
      mutate(
        age_numeric = case_when(
          str_detect(as.character(age_gr), "<5") ~ 0,
          TRUE ~ as.numeric(str_extract(as.character(age_gr), "^\\d+"))
        )
      )
    
    # Add YLD parameters (age independent)
    sim_df <- sim_df %>%
      mutate(
        dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      # Age-dependent parameters:
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
      # Calculate YLD with uncertainty bounds
      mutate(
        yld_acute = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_acute_lo = (hospitalised_lo * dw_hosp * dur_acute) +
          ((post_low95 - hospitalised_hi) * dw_nonhosp * dur_acute),
        yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
          ((post_hi95 - hospitalised_lo) * dw_nonhosp * dur_acute),
        
        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
          ((post_low95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
          ((post_hi95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
        
        yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
          ((post_low95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
          ((post_hi95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
        
        yld_total = yld_acute + yld_subacute + yld_chronic,
        yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
        yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
      ) %>%
      # Determine remaining life expectancy by age and calculate YLL and DALY
      mutate(
        le_left = case_when(
          age_numeric %in% c(0,5) ~ quantile(le_sample$le_1, 0.5),
          age_numeric %in% c(10,15) ~ quantile(le_sample$le_2, 0.5),
          age_numeric %in% c(20,25) ~ quantile(le_sample$le_3, 0.5),
          age_numeric %in% c(30,35) ~ quantile(le_sample$le_4, 0.5),
          age_numeric %in% c(40,45) ~ quantile(le_sample$le_5, 0.5),
          age_numeric %in% c(50,55) ~ quantile(le_sample$le_6, 0.5),
          age_numeric %in% c(60,65) ~ quantile(le_sample$le_7, 0.5),
          age_numeric %in% c(70,75) ~ quantile(le_sample$le_8, 0.5),
          age_numeric %in% c(80,85) ~ quantile(le_sample$le_9, 0.5)
        )
      ) %>%
      mutate(
        yll = fatal * le_left,
        yll_lo = fatal_lo * le_left,
        yll_hi = fatal_hi * le_left,
        daly_tot = yld_total + yll,
        daly_tot_lo = yld_total_lo + yll_lo,
        daly_tot_hi = yld_total_hi + yll_hi,
        cum_daly = cumsum(daly_tot)
      )
    
    if(!is.null(prevacc_ui)) {
      sim_df <- left_join(sim_df, 
                          dplyr::select(prevacc_ui, AgeGroup, Week, 
                                        pre_Median = Median, 
                                        pre_low95 = low95, 
                                        pre_hi95 = hi95), 
                          by = c("AgeGroup", "Week"))
    }
    
    scenario_result[[s]] <- list(sim_result = list(age_array_raw_symp           = age_array_raw_symp,
                                                   age_array_raw_inf            = age_array_raw_inf,
                                                   weekly_cases_median_rawsymp  = weekly_cases_median_rawsymp,
                                                   weekly_cases_median_rawinf   = weekly_cases_median_rawinf,
                                                   weekly_cases_low95_rawsymp   = weekly_cases_low95_rawsymp,
                                                   weekly_cases_low95_rawinf    = weekly_cases_low95_rawinf,
                                                   weekly_cases_hi95_rawsymp    = weekly_cases_hi95_rawsymp),
                                 sim_out                      = sim_out,
                                 sim_df                       = sim_df,
                                 weekly_df                    = weekly_df
    )
  }
  
  return(scenario_result)
}



preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups, 
                                N = N_ceara$Ceará, 
                                region = "Ceará",
                                observed = observed_ce
)

pre_results_ce_ui <- summarise_presim_ui_inf(sim_result     = preui_ce, 
                                         observed         = observed_ce,
                                         age_gr_levels    = age_gr_levels,
                                         lhs_sample_young = lhs_sample_young,
                                         lhs_old          = lhs_old,
                                         le_sample        = le_sample,
                                         hosp             = hosp,
                                         fatal            = fatal,
                                         nh_fatal         = nh_fatal,
                                         region           = "Ceará"
)

prevacc_ui_ce <- pre_results_ce_ui$summary_cases_pre
prevacc_ui_ce_all <- pre_results_ce_ui$summary_cases_pre_all
pre_summary_age_ce <- pre_results_ce_ui$summary_cases_pre_age

postsim_ce_ui <- run_simulation_scenarios_ui_inf(
  target_age_list   = target_age_list,
  #supply            = supply_ce,
  observed          = observed_ce,
  N                 = N_ceara$Ceará,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Ceará",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_ce,
  posterior         = posterior_ce
)

postsim_all_ce_ui <- postsim_all_ui(
  scenario_result      = postsim_ce_ui,
  observed             = observed_ce,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_ce,
  pre_summary_cases     = prevacc_ui_ce,
  pre_summary_cases_all     = prevacc_ui_ce_all,
  region                = "Ceará"
)

summary_list_ce <- postsim_all_ce_ui$summary_list
# assume each summary_list[[i]] has columns: Week, AgeGroup, Cases, fatal, daly_tot
summary_age_week <- lapply(seq_along(summary_list_ce), function(i) {
  summary_list_ce[[i]] %>%
    mutate(
      Scenario = paste0("Scenario_", i),
      pre_vacc_cases = prevacc_ui_ce$Median,      # if pre_summary_cases also age-specific
      pre_vacc_fatal = prevacc_ui_ce$fatal,      # same length/order as AgeGroup
      pre_vacc_daly  = prevacc_ui_ce$daly_tot
    ) %>%
    group_by(Week, AgeGroup, Scenario) %>%
    summarise(
      post_cases   = sum(Cases,          na.rm = TRUE),
      pre_cases    = sum(pre_vacc_cases, na.rm = TRUE),
      diff_cases   = pre_cases - post_cases,
      post_fatal   = sum(fatal,          na.rm = TRUE),
      pre_fatal    = sum(pre_vacc_fatal, na.rm = TRUE),
      diff_fatal   = pre_fatal - post_fatal,
      post_daly    = sum(daly_tot,       na.rm = TRUE),
      pre_daly     = sum(pre_vacc_daly,  na.rm = TRUE),
      diff_daly    = pre_daly - post_daly,
      .groups = "drop"
    )
})
summary_age_week_df <- bind_rows(summary_age_week)

age_targets <- list(
  Scenario_1 = 1:4,
  Scenario_2 = 5:12,
  Scenario_3 = 13:18
)

# 2) target 플래그 추가
summary_age_week_df <- summary_age_week_df %>%
  mutate(
    target = map2_lgl(Scenario, AgeGroup,
                      ~ .y %in% age_targets[[.x]])
  )
summary_target_week <- summary_age_week_df %>%
  filter(target) %>%
  group_by(Week, Scenario) %>%
  summarise(
    post_cases = sum(post_cases, na.rm = TRUE),
    pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  )


summary_target_week %>%
  pivot_longer(
    cols = c(post_cases, pre_cases),
    names_to  = "Type",
    values_to = "Cases"
  ) %>%
  ggplot(aes(x = Week, y = Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ Scenario) +
  labs(
    x = "Week",
    y = "Cases"
  ) +
  theme_pubclean()+
  theme(legend.position = "bottom")


summary_age_week_df %>% 
  filter(target == FALSE) %>% 
  group_by(Week, Scenario) %>% 
  summarise(
    post_cases = sum(post_cases),
    pre_cases  = sum(pre_cases),
    .groups = "drop"
  ) %>% 
  pivot_longer(
    cols = c(post_cases, pre_cases),
    names_to  = "Type",
    values_to = "Cases"
  ) %>% 
  ggplot(aes(Week, Cases, color = Type)) +
  geom_line() +
  facet_wrap(~ Scenario) +
  theme_minimal()
