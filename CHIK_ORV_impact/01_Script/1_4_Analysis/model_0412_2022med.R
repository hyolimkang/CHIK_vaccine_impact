load("00_Data/0_2_Processed/fit_prevacc_sg_22.RData")
load("00_Data/0_2_Processed/fit_prevacc_go_22.RData")
load("00_Data/0_2_Processed/fit_prevacc_sp_22.RData")
load("00_Data/0_2_Processed/fit_prevacc_ms_22.RData")
load("00_Data/0_2_Processed/fit_prevacc_mh_22.RData")
load("00_Data/0_2_Processed/brazil_pop_transposed.RData")
load("00_Data/0_2_Processed/observed_2022_med.RData")
load("00_Data/0_2_Processed/bra_foi_state_summ.RData")
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")
lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")
le_sample <- readRDS("00_Data/0_2_Processed/le_sample.RDS")
source("01_Script/1_1_Functions/sim_function_final.R")
source("01_Script/1_1_Functions/age_struc_fitting_region_func.R")
source("01_Script/1_1_Functions/library.R")

# extract params-----------------------------------------------------------------
param_sg_22 <- extract_params(fit_prevacc_sg_22)
param_go_22 <- extract_params(fit_prevacc_go_22)
param_sp_22 <- extract_params(fit_prevacc_sp_22)
param_ms_22 <- extract_params(fit_prevacc_ms_22)
param_mh_22 <- extract_params(fit_prevacc_mh_22)

posterior_go_22 <- param_go_22$posterior_prevacc
posterior_sg_22 <- param_sg_22$posterior_prevacc
posterior_sp_22 <- param_sp_22$posterior_prevacc
posterior_ms_22 <- param_ms_22$posterior_prevacc
posterior_mh_22 <- param_mh_22$posterior_prevacc

# scenarios
target_age_list <- list(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 20-59 yrs old
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)) # >60 yrs old) # <20 yrs old

# prevacc 95%UI simulation ----------------------------------------------------

## pre-processing ----------------------------------------------------
preui_go_22 <- simulate_pre_ui_age(posterior = posterior_go_22, bra_foi_state_summ, age_groups, 
                                   N = N_go, 
                                   region = "Goiás",
                                   observed = observed_go_22
)

preui_sp_22 <- simulate_pre_ui_age(posterior = posterior_sp_22, bra_foi_state_summ, age_groups, 
                                   N = N_sp, 
                                   region = "São Paulo",
                                   observed = observed_sp_22
)

preui_sg_22 <- simulate_pre_ui_age(posterior = posterior_sg_22, bra_foi_state_summ, age_groups, 
                                   N = N_sg, 
                                   region = "Sergipe",
                                   observed = observed_sg_22
)

preui_ms_22 <- simulate_pre_ui_age(posterior = posterior_ms_22, bra_foi_state_summ, age_groups, 
                                   N = N_ms, 
                                   region = "Mato Grosso do Sul",
                                   observed = observed_ms_22
)

preui_mh_22 <- simulate_pre_ui_age(posterior = posterior_mh_22, bra_foi_state_summ, age_groups, 
                                   N = N_mh, 
                                   region = "Maranhão",
                                   observed = observed_mh_22
)

## summarise ----------------------------------------------------
pre_results_go_ui_22 <- summarise_presim_ui(sim_result       = preui_go_22, 
                                            observed         = observed_go_22,
                                            age_gr_levels    = age_gr_levels,
                                            lhs_sample_young = lhs_sample_young,
                                            lhs_old          = lhs_old,
                                            le_sample        = le_sample,
                                            hosp             = hosp,
                                            fatal            = fatal,
                                            nh_fatal         = nh_fatal,
                                            region           = "Goiás"
)

pre_results_sp_ui_22 <- summarise_presim_ui(sim_result       = preui_sp_22, 
                                            observed         = observed_sp_22,
                                            age_gr_levels    = age_gr_levels,
                                            lhs_sample_young = lhs_sample_young,
                                            lhs_old          = lhs_old,
                                            le_sample        = le_sample,
                                            hosp             = hosp,
                                            fatal            = fatal,
                                            nh_fatal         = nh_fatal,
                                            region           = "São Paulo"
)

pre_results_sg_ui_22 <- summarise_presim_ui(sim_result       = preui_sg_22, 
                                            observed         = observed_sg_22,
                                            age_gr_levels    = age_gr_levels,
                                            lhs_sample_young = lhs_sample_young,
                                            lhs_old          = lhs_old,
                                            le_sample        = le_sample,
                                            hosp             = hosp,
                                            fatal            = fatal,
                                            nh_fatal         = nh_fatal,
                                            region           = "Sergipe"
)

pre_results_ms_ui_22 <- summarise_presim_ui(sim_result       = preui_ms_22, 
                                            observed         = observed_ms_22,
                                            age_gr_levels    = age_gr_levels,
                                            lhs_sample_young = lhs_sample_young,
                                            lhs_old          = lhs_old,
                                            le_sample        = le_sample,
                                            hosp             = hosp,
                                            fatal            = fatal,
                                            nh_fatal         = nh_fatal,
                                            region           = "Mato Grosso do Sul"
)

pre_results_mh_ui_22 <- summarise_presim_ui(sim_result       = preui_mh_22, 
                                            observed         = observed_mh_22,
                                            age_gr_levels    = age_gr_levels,
                                            lhs_sample_young = lhs_sample_young,
                                            lhs_old          = lhs_old,
                                            le_sample        = le_sample,
                                            hosp             = hosp,
                                            fatal            = fatal,
                                            nh_fatal         = nh_fatal,
                                            region           = "Maranhão"
)

## prevacc outputs ----------------------------------------------------

prevacc_ui_sp_22 <- pre_results_sp_ui_22$summary_cases_pre
prevacc_ui_sp_all_22 <- pre_results_sp_ui_22$summary_cases_pre_all
pre_summary_age_sp_22 <- pre_results_sp_ui_22$summary_cases_pre_age

prevacc_ui_go_22 <- pre_results_go_ui_22$summary_cases_pre
prevacc_ui_go_all_22 <- pre_results_go_ui_22$summary_cases_pre_all
pre_summary_age_go_22 <- pre_results_go_ui_22$summary_cases_pre_age

prevacc_ui_sg_22 <- pre_results_sg_ui_22$summary_cases_pre
prevacc_ui_sg_all_22 <- pre_results_sg_ui_22$summary_cases_pre_all
pre_summary_age_sg_22 <- pre_results_sg_ui_22$summary_cases_pre_age

prevacc_ui_ms_22 <- pre_results_ms_ui_22$summary_cases_pre
prevacc_ui_ms_all_22 <- pre_results_ms_ui_22$summary_cases_pre_all
pre_summary_age_ms_22 <- pre_results_ms_ui_22$summary_cases_pre_age

prevacc_ui_mh_22 <- pre_results_mh_ui_22$summary_cases_pre
prevacc_ui_mh_all_22 <- pre_results_mh_ui_22$summary_cases_pre_all
pre_summary_age_mh_22 <- pre_results_mh_ui_22$summary_cases_pre_age

# postsim 95% UI run ----------------------------------------------------------

postsim_sp_ui_22 <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_rg,
  observed           = observed_sp_22,
  N                 = N_sp,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "São Paulo",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_sp_22,
  posterior         = posterior_sp_22
)


postsim_go_ui_22 <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_ag,
  observed          = observed_go_22,
  N                 = N_go,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Goiás",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_go_22,
  posterior         = posterior_go_22
)

postsim_sg_ui_22 <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_sg_22,
  N                 = N_sg,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Sergipe",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_sg_22,
  posterior         = posterior_sg_22
)

postsim_ms_ui_22 <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_ms_22,
  N                 = N_ms,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Mato Grosso do Sul",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_ms_22,
  posterior         = posterior_ms_22
)

postsim_mh_ui_22 <- run_simulation_scenarios_ui(
  target_age_list   = target_age_list,
  #supply            = supply_mg,
  observed          = observed_mh_22,
  N                 = N_mh,
  bra_foi_state_summ = bra_foi_state_summ,
  age_groups        = age_groups,
  region_name       = "Maranhão",      # region name used for FOI-based immunity
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  age_gr_levels     = age_gr_levels,
  prevacc_ui        = prevacc_ui_mh_22,
  posterior         = posterior_mh_22
)

# summarise for post sim ------------------------------------------------------

postsim_all_sp_ui_22 <- postsim_all_ui(
  scenario_result      = postsim_sp_ui_22,
  observed             = observed_sp_22,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_sp_22,
  pre_summary_cases     = prevacc_ui_sp_22,
  pre_summary_cases_all     = prevacc_ui_sp_all_22,
  region                = "São Paulo"
)


postsim_all_go_ui_22 <- postsim_all_ui(
  scenario_result      = postsim_go_ui_22,
  observed             = observed_go_22,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_go_22,
  pre_summary_cases     = prevacc_ui_go_22,
  pre_summary_cases_all     = prevacc_ui_go_all_22,
  region                = "Goiás"
)

postsim_all_sg_ui_22 <- postsim_all_ui(
  scenario_result      = postsim_sg_ui_22,
  observed             = observed_sg_22,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_sg_22,
  pre_summary_cases     = prevacc_ui_sg_22,
  pre_summary_cases_all     = prevacc_ui_sg_all_22,
  region                = "Sergipe"
)

postsim_all_ms_ui_22 <- postsim_all_ui(
  scenario_result      = postsim_ms_ui_22,
  observed             = observed_ms_22,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_ms_22,
  pre_summary_cases     = prevacc_ui_ms_22,
  pre_summary_cases_all     = prevacc_ui_ms_all_22,
  region                = "Mato Grosso do Sul"
)

postsim_all_mh_ui_22 <- postsim_all_ui(
  scenario_result      = postsim_mh_ui_22,
  observed             = observed_mh_22,
  age_gr_levels        = age_gr_levels,
  pre_summary_cases_age = pre_summary_age_mh_22,
  pre_summary_cases     = prevacc_ui_mh_22,
  pre_summary_cases_all     = prevacc_ui_mh_all_22,
  region                = "Maranhão"
)

# pre-post graphs -------------------------------------------------------------

epi_graph(postsim_all_sp_ui_22,
                observed_sp_22)

epi_graph(postsim_all_go_ui_22,
                observed_go_22)

epi_graph(postsim_all_sg_ui_22,
                observed_sg_22)

epi_graph(postsim_all_ms_ui_22,
          observed_ms_22)

epi_graph(postsim_all_mh_ui_22,
          observed_mh_22)

combined_prepost_case_22_med <- bind_rows(
  postsim_all_sp_ui_22$summary_week_df,
  postsim_all_go_ui_22$summary_week_df,
  postsim_all_sg_ui_22$summary_week_df,
  postsim_all_ms_ui_22$summary_week_df,
  postsim_all_mh_ui_22$summary_week_df
)

global_impact_all <- combined_prepost_case_22_med %>%
  group_by(Scenario, region) %>%
  summarise(
    total_post_cases = sum(post_cases),
    total_pre_case   = sum(pre_cases),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_case - total_post_cases,
    impact = diff / total_pre_case * 100
  )

annotation_df <- global_impact_all %>%
  group_by(region) %>%
  summarise(ann = paste0(Scenario, ": ", round(impact, 1), "%", collapse = "\n"),
            .groups = "drop")

overall_max <- max(combined_prepost_case_22_med$hi95, na.rm = TRUE)

# Create an annotation data frame that has one row per region
vaccine_ann_df <- combined_prepost_case_22_med %>%
  distinct(region) %>%
  mutate(x = 4, 
         y = overall_max * 0.95, 
         label = "<----- Vaccine impact start (+ 2 wks after initiation)")

max_cases_22_med <- max(combined_prepost_case_22_med$hi95, na.rm = TRUE)
vacc_start_week_s1 <- postsim_all_sg_ui_22$vacc_weeks$scenario1$start
vacc_end_week_s1   <- postsim_all_sg_ui_22$vacc_weeks$scenario1$end
vacc_start_week_s2 <- postsim_all_sg_ui_22$vacc_weeks$scenario2$start
vacc_end_week_s2   <- postsim_all_sg_ui_22$vacc_weeks$scenario2$end
vacc_start_week_s3 <- postsim_all_sg_ui_22$vacc_weeks$scenario3$start
vacc_end_week_s3   <- postsim_all_sg_ui_22$vacc_weeks$scenario3$end

pre_post_graph_22_med <- 
  ggplot(combined_prepost_case_22_med) +
  # Scenario shading ribbons (existing)
  geom_ribbon(data = combined_prepost_case_22_med %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 2/3 * max_cases_22_med, ymax = max_cases_22_med, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_case_22_med %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
              aes(x = Week, ymin = 1/3 * max_cases_22_med, ymax = 2/3 * max_cases_22_med, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_case_22_med %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
              aes(x = Week, ymin = 0, ymax = 1/3 * max_cases_22_med, fill = Scenario),
              alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
              alpha = 0.2) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
              fill = "gray", alpha = 0.6) +
  # Post-vaccination case line (colored by scenario)
  geom_line(aes(x = Week, y = post_cases, color = Scenario, linetype = "Post Cases"), size = 0.5) +
  # Pre-vaccination case line (dashed, in black)
  geom_line(aes(x = Week, y = pre_cases, group = Scenario, linetype = "Pre Cases"), color = "black", size = 0.5) +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Scenario_1" = "<20 years only", 
                               "Scenario_2" = "20-59 years only", 
                               "Scenario_3" = ">60 years only")) +
  scale_color_brewer(palette = "Set1",
                     labels = c("Scenario_1" = "<20 years only", 
                                "Scenario_2" = "20-59 years only", 
                                "Scenario_3" = ">60 years only")) +
  scale_y_continuous(labels = scales::comma) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
       x = "Week", y = "Predicted symptomatic reported cases") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    plot.margin = margin(5, 30, 5, 5)
  ) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = Inf, label = ann), 
            inherit.aes = FALSE,
            hjust = 1.1, 
            vjust = 1.1, 
            size = 3) +
  #geom_text(data = vaccine_ann_df,
  #          aes(x = x, y = y, label = label),
  #          inherit.aes = FALSE,
  #          hjust = 0, vjust = 1, size = 2.5)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  geom_point(data = observed_all_22_med, aes(x = Week, y = Observed),size = 0.3) + 
  facet_wrap(~region, ncol= 3)

ggsave(filename = "02_Outputs/2_1_Figures/pre_post_region_med.jpg", pre_post_graph_22_med, width = 12, height = 7, dpi = 1200)


# vacc allocations -------------------------------------------------------------

vacc_alloc_sp_22 <- vacc_allocation(postsim_all_sp_ui_22, observed_sp_22, "São Paulo")

vacc_alloc_go_22 <- vacc_allocation(postsim_all_go_ui_22, observed_go_22, "Goiás")

vacc_alloc_sg_22 <- vacc_allocation(postsim_all_sg_ui_22, observed_sg_22, "Sergipe" )

vacc_alloc_ms_22 <- vacc_allocation(postsim_all_ms_ui_22, observed_ms_22, "Mato Grosso do Sul" )

vacc_alloc_mh_22 <- vacc_allocation(postsim_all_mh_ui_22, observed_mh_22, "Maranhão" )

# vacc allocations graph -------------------------------------------------------------
vacc_alloc_graph(vacc_alloc_sp_22)
vacc_alloc_graph(vacc_alloc_go_22)
vacc_alloc_graph(vacc_alloc_sg_22)
vacc_alloc_graph(vacc_alloc_ms_22)
vacc_alloc_graph(vacc_alloc_mh_22)

combined_vacc_alloc_22_med <- bind_rows(
  vacc_alloc_sp_22$weekly_allocation_long,
  vacc_alloc_go_22$weekly_allocation_long,
  vacc_alloc_sg_22$weekly_allocation_long,
  vacc_alloc_ms_22$weekly_allocation_long,
  vacc_alloc_mh_22$weekly_allocation_long
)

combined_vacc_alloc_summ_22_med <- combined_vacc_alloc_22_med %>%
  group_by(Scenario, AgeGroup, age_gr, Week) %>%
  summarise(
    Vaccinated = sum(Vaccinated, na.rm = TRUE),
    tot_sum = sum(tot_sum, na.rm = TRUE)
  ) %>%
  ungroup()

paired <- brewer.pal(12, "Paired")
extra <- brewer.pal(8, "Dark2")[1:6]
full_palette <- c(paired, extra)

p <- 
  ggplot(combined_vacc_alloc_summ_22_med, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=comma)+
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c("Scenario1" = "<20 years only", 
                                              "Scenario2" = "20-59 years only", 
                                              "Scenario3" = ">60 years only")))  +         # Facet by Scenario if you want separate plots for each
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/fig_vacc_alloc_tot.jpg", p, width = 12, height = 7, dpi = 1200)

# nnv  ------------------------------------------------------------------------

nnv_sp_22 <- nnv_list(vacc_alloc_sp_22, postsim_all_sp_ui_22, N_sp, "São Paulo", observed_sp_22)

nnv_go_22 <- nnv_list(vacc_alloc_go_22, postsim_all_go_ui_22, N_go, "Goiás", observed_go_22)

nnv_sg_22 <- nnv_list(vacc_alloc_sg_22, postsim_all_sg_ui_22, N_sg, "Sergipe", observed_sg_22)

nnv_ms_22 <- nnv_list(vacc_alloc_ms_22, postsim_all_ms_ui_22, N_ms, "Mato Grosso do Sul", observed_ms_22)

nnv_mh_22 <- nnv_list(vacc_alloc_mh_22, postsim_all_mh_ui_22, N_mh, "Maranhão", observed_mh_22)

# global total -----------------------------------------------------------------
## total pre-post infection graph (national level) ------------------------------
df_sp_22 <- postsim_all_sp_ui_22$summary_week_df
df_go_22 <- postsim_all_go_ui_22$summary_week_df
df_sg_22 <- postsim_all_sg_ui_22$summary_week_df
df_ms_22 <- postsim_all_ms_ui_22$summary_week_df
df_mh_22 <- postsim_all_mh_ui_22$summary_week_df

# Combine the two data frames
combined_df_22_med <- bind_rows(df_sp_22, df_go_22, 
                                df_sg_22, df_ms_22, df_mh_22) %>%
  group_by(Scenario, Week) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

# (2) Calculate global impact from your combined_df
global_impact_week <- combined_df_22_med %>%
  group_by(Scenario) %>%
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_case   = sum(pre_cases, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_case - total_post_cases,
    impact = diff / total_pre_case * 100
  )

# (3) Create annotation text from the global impact
annotation_text <- paste0(
  global_impact_week$Scenario, ": ",
  round(global_impact_week$impact, 1), "%",
  collapse = "\n"
)

postsim_all_global_ui_22_med <- list(
  postsim_all_sg_ui_22 = postsim_all_sg_ui_22,
  postsim_all_sg_ui_22 = postsim_all_go_ui_22,
  summary_week_df_22   = combined_df_22_med,
  global_impact_week = global_impact_week,
  annotation_text = annotation_text,
  vacc_start_week_s1 = postsim_all_sg_ui_22$vacc_weeks$scenario1,
  vacc_start_week_s2 = postsim_all_sg_ui_22$vacc_weeks$scenario2,
  vacc_start_week_s3 = postsim_all_sg_ui_22$vacc_weeks$scenario3
  
)

epi_graph_all_22_med <- 
  epi_graph_nat(postsim_all_global_ui_22_med, bra_all_sum_22_med)

## combined graph
pre_post_graph_comb_22_med <- pre_post_graph_22_med +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Vaccine impact at sub-national level")

# Remove legend from epi_graph_all (Scenario A)
epi_graph_all_22_med <- epi_graph_all_22_med + theme(legend.position = "none") + 
  ggtitle("Vaccine impact at national level")

# Arrange plots side-by-side with a common legend
combined_graph_22_med <- ggarrange(epi_graph_all_22_med, pre_post_graph_comb_22_med,
                               ncol = 2, nrow = 1,
                               labels = c("A", "B"),
                               common.legend = TRUE,
                               legend = "bottom",
                               align = "none")
combined_graph_22_med <- combined_graph_22_med +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )


ggsave(filename = "02_Outputs/2_1_Figures/combined_graph_22_med.jpg", combined_graph_22_med, width = 12, height = 6, dpi = 1200)


# table: pre-post cases -----------------------------------------------------
scenario_result_sg <- attach_target_column(postsim_all_sg_ui_22,
                                           observed = observed_sg_22,
                                           region   = "Sergipe")
scenario_result_sp <- attach_target_column(postsim_all_sp_ui_22,
                                           observed = observed_sp_22,
                                           region   = "São Paulo")
scenario_result_go <- attach_target_column(postsim_all_go_ui,
                                           observed = observed_go_22,
                                           region   = "Goiás")

combined_scenario <- bind_rows(scenario_result_sg, scenario_result_sp, scenario_result_gp)

combined_scenario_summ <- combined_scenario %>%
  group_by(target, Scenario) %>%
  summarise(
    pre_cases    = sum(pre_cases, na.rm = TRUE),
    post_cases   = sum(post_cases, na.rm = TRUE),
    post_low95   = sum(post_low95, na.rm = TRUE),
    post_hi95    = sum(post_hi95, na.rm = TRUE),
    pre_low95    = sum(pre_low95, na.rm = TRUE),
    pre_hi95     = sum(pre_hi95, na.rm = TRUE),
    post_fatal   = sum(post_fatal, na.rm = TRUE),
    pre_fatal    = sum(pre_fatal, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  mutate(
    diff_case    = pre_cases - post_cases,
    impact       = diff_case / pre_cases * 100,
    diff_fatal   = pre_fatal - post_fatal,
    imapact_fatal = diff_fatal / pre_fatal * 100
  )

## total nnv graph (national level) ------------------------------
n_scenarios = 3

df_nnv_sp_22 <- nnv_sp_22$final_summ_df
df_nnv_go_22 <- nnv_go_22$final_summ_df
df_nnv_sg_22 <- nnv_sg_22$final_summ_df
df_nnv_ms_22 <- nnv_ms_22$final_summ_df
df_nnv_mh_22 <- nnv_mh_22$final_summ_df


df_nnv_sp_22$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_sp_22$age_gr <- factor(df_nnv_sp_22$age_gr, levels = age_gr_levels)

df_nnv_go_22$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_go_22$age_gr <- factor(df_nnv_go_22$age_gr, levels = age_gr_levels)

df_nnv_sg_22$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_sg_22$age_gr <- factor(df_nnv_sg_22$age_gr, levels = age_gr_levels)

df_nnv_ms_22$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_ms_22$age_gr <- factor(df_nnv_ms_22$age_gr, levels = age_gr_levels)

df_nnv_mh_22$age_gr <- rep(age_gr[1:18],n_scenarios)
df_nnv_mh_22$age_gr <- factor(df_nnv_mh_22$age_gr, levels = age_gr_levels)

combined_nnv_df_region_22 <- bind_rows(df_nnv_sp_22, df_nnv_go_22, df_nnv_sg_22,
                                       df_nnv_ms_22, df_nnv_mh_22) 

combined_nnv_df_22 <- combined_nnv_df_region_22 %>%
  group_by(scenario, AgeGroup) %>%
  summarise(
    Median            = sum(Median, na.rm = TRUE),
    low95             = sum(low95, na.rm = TRUE),
    hi95              = sum(hi95, na.rm = TRUE),
    hospitalised      = sum(hospitalised, na.rm = TRUE),
    hospitalised_lo   = sum(hospitalised_lo, na.rm = TRUE),
    hospitalised_hi   = sum(hospitalised_hi, na.rm = TRUE),
    yld_acute         = sum(yld_acute, na.rm = TRUE),
    yld_acute_lo      = sum(yld_acute_lo, na.rm = TRUE),
    yld_acute_hi      = sum(yld_acute_hi, na.rm = TRUE),
    yld_subacute      = sum(yld_subacute, na.rm = TRUE),
    yld_subacute_lo   = sum(yld_subacute_lo, na.rm = TRUE),
    yld_subacute_hi   = sum(yld_subacute_hi, na.rm = TRUE),
    yld_chronic       = sum(yld_chronic, na.rm = TRUE),
    yld_chronic_lo    = sum(yld_chronic_lo, na.rm = TRUE),
    yld_chronic_hi    = sum(yld_chronic_hi, na.rm = TRUE),
    yld_total         = sum(yld_total, na.rm = TRUE),
    yld_total_lo      = sum(yld_total_lo, na.rm = TRUE),
    yld_total_hi      = sum(yld_total_hi, na.rm = TRUE),
    yll               = sum(yll, na.rm = TRUE),
    yll_lo            = sum(yll_lo, na.rm = TRUE),
    yll_hi            = sum(yll_hi, na.rm = TRUE),
    daly_tot          = sum(daly_tot, na.rm = TRUE),
    daly_tot_lo       = sum(daly_tot_lo, na.rm = TRUE),
    daly_tot_hi       = sum(daly_tot_hi, na.rm = TRUE),
    fatal             = sum(fatal, na.rm = TRUE),
    fatal_lo          = sum(fatal_lo, na.rm = TRUE),
    fatal_hi          = sum(fatal_hi, na.rm = TRUE),
    pre_vacc          = sum(pre_vacc, na.rm = TRUE),
    pre_vacc_low95    = sum(pre_vacc_low95, na.rm = TRUE),
    pre_vacc_hi95     = sum(pre_vacc_hi95, na.rm = TRUE),
    pre_yld_acute     = sum(pre_yld_acute, na.rm = TRUE),
    pre_yld_acute_low95 = sum(pre_yld_acute_low95, na.rm = TRUE),
    pre_yld_acute_hi  = sum(pre_yld_acute_hi, na.rm = TRUE),
    pre_yld_subacute  = sum(pre_yld_subacute, na.rm = TRUE),
    pre_yld_subacute_low95 = sum(pre_yld_subacute_low95, na.rm = TRUE),
    pre_yld_subacute_hi = sum(pre_yld_subacute_hi, na.rm = TRUE),
    pre_yld_chronic   = sum(pre_yld_chronic, na.rm = TRUE),
    pre_yld_chronic_low95 = sum(pre_yld_chronic_low95, na.rm = TRUE),
    pre_yld_chronic_hi = sum(pre_yld_chronic_hi, na.rm = TRUE),
    pre_yld_total     = sum(pre_yld_total, na.rm = TRUE),
    pre_yld_total_low95 = sum(pre_yld_total_low95, na.rm = TRUE),
    pre_yld_total_hi  = sum(pre_yld_total_hi, na.rm = TRUE),
    pre_yll           = sum(pre_yll, na.rm = TRUE),
    pre_yll_low95     = sum(pre_yll_low95, na.rm = TRUE),
    pre_yll_hi        = sum(pre_yll_hi, na.rm = TRUE),
    pre_daly          = sum(pre_daly, na.rm = TRUE),
    pre_daly_low95    = sum(pre_daly_low95, na.rm = TRUE),
    pre_daly_hi       = sum(pre_daly_hi, na.rm = TRUE),
    pre_fatal         = sum(pre_fatal, na.rm = TRUE),
    pre_fatal_low95   = sum(pre_fatal_low95, na.rm = TRUE),
    pre_fatal_hi      = sum(pre_fatal_hi, na.rm = TRUE),
    tot_vacc          = sum(tot_vacc, na.rm = TRUE)
  ) %>% 
  ungroup()

# Next, compute differences, impacts, and NNV values (with UI) for each metric.
final_nnv_df_22 <- combined_nnv_df_22 %>%
  mutate(
    # For Cases:
    diff           = pre_vacc - Median,
    diff_low       = pre_vacc_low95 - hi95,
    diff_hi        = pre_vacc_hi95 - low95,
    impact         = diff / pre_vacc * 100,
    impact_low     = diff_low / pre_vacc_hi95 * 100,
    impact_hi      = diff_hi / pre_vacc_low95 * 100,
    nnv            = tot_vacc / diff,
    nnv_lo         = tot_vacc / diff_hi,
    nnv_hi         = tot_vacc / diff_low,
    
    # For Fatal outcomes:
    diff_fatal     = pre_fatal - fatal,
    diff_fatal_low = pre_fatal_low95 - fatal_hi,
    diff_fatal_hi  = pre_fatal_hi - fatal_lo,
    impact_fatal   = diff_fatal / pre_fatal * 100,
    impact_fatal_low = diff_fatal_low / pre_fatal_hi * 100,
    impact_fatal_hi  = diff_fatal_hi / pre_fatal_low95 * 100,
    nnv_fatal      = tot_vacc / diff_fatal,
    nnv_fatal_lo   = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi   = tot_vacc / diff_fatal_low,
    
    # For YLD Acute:
    diff_yld_acute       = pre_yld_acute - yld_acute,
    diff_yld_acute_low   = pre_yld_acute_low95 - yld_acute_hi,
    diff_yld_acute_hi    = pre_yld_acute_hi - yld_acute_lo,
    impact_yld_acute     = diff_yld_acute / pre_yld_acute * 100,
    impact_yld_acute_low = diff_yld_acute_low / pre_yld_acute_hi * 100,
    impact_yld_acute_hi  = diff_yld_acute_hi / pre_yld_acute_low95 * 100,
    nnv_yld_acute        = tot_vacc / diff_yld_acute,
    nnv_yld_acute_lo     = tot_vacc / diff_yld_acute_hi,
    nnv_yld_acute_hi     = tot_vacc / diff_yld_acute_low,
    
    # For YLD Subacute:
    diff_yld_subacute       = pre_yld_subacute - yld_subacute,
    diff_yld_subacute_low   = pre_yld_subacute_low95 - yld_subacute_hi,
    diff_yld_subacute_hi    = pre_yld_subacute_hi - yld_subacute_lo,
    impact_yld_subacute     = diff_yld_subacute / pre_yld_subacute * 100,
    impact_yld_subacute_low = diff_yld_subacute_low / pre_yld_subacute_hi * 100,
    impact_yld_subacute_hi  = diff_yld_subacute_hi / pre_yld_subacute_low95 * 100,
    nnv_yld_subacute        = tot_vacc / diff_yld_subacute,
    nnv_yld_subacute_lo     = tot_vacc / diff_yld_subacute_hi,
    nnv_yld_subacute_hi     = tot_vacc / diff_yld_subacute_low,
    
    # For YLD Chronic:
    diff_yld_chronic       = pre_yld_chronic - yld_chronic,
    diff_yld_chronic_low   = pre_yld_chronic_low95 - yld_chronic_hi,
    diff_yld_chronic_hi    = pre_yld_chronic_hi - yld_chronic_lo,
    impact_yld_chronic     = diff_yld_chronic / pre_yld_chronic * 100,
    impact_yld_chronic_low = diff_yld_chronic_low / pre_yld_chronic_hi * 100,
    impact_yld_chronic_hi  = diff_yld_chronic_hi / pre_yld_chronic_low95 * 100,
    nnv_yld_chronic        = tot_vacc / diff_yld_chronic,
    nnv_yld_chronic_lo     = tot_vacc / diff_yld_chronic_hi,
    nnv_yld_chronic_hi     = tot_vacc / diff_yld_chronic_low,
    
    # For YLD Total:
    diff_yld_total       = pre_yld_total - yld_total,
    diff_yld_total_low   = pre_yld_total_low95 - yld_total_hi,
    diff_yld_total_hi    = pre_yld_total_hi - yld_total_lo,
    impact_yld_total     = diff_yld_total / pre_yld_total * 100,
    impact_yld_total_low = diff_yld_total_low / pre_yld_total_hi * 100,
    impact_yld_total_hi  = diff_yld_total_hi / pre_yld_total_low95 * 100,
    nnv_yld_total        = tot_vacc / diff_yld_total,
    nnv_yld_total_lo     = tot_vacc / diff_yld_total_hi,
    nnv_yld_total_hi     = tot_vacc / diff_yld_total_low,
    
    # For YLL:
    diff_yll       = pre_yll - yll,
    diff_yll_low   = pre_yll_low95 - yll_hi,
    diff_yll_hi    = pre_yll_hi - yll_lo,
    impact_yll     = diff_yll / pre_yll * 100,
    impact_yll_low = diff_yll_low / pre_yll_hi * 100,
    impact_yll_hi  = diff_yll_hi / pre_yll_low95 * 100,
    nnv_yll        = tot_vacc / diff_yll,
    nnv_yll_lo     = tot_vacc / diff_yll_hi,
    nnv_yll_hi     = tot_vacc / diff_yll_low,
    
    # For DALY:
    diff_daly       = pre_daly - daly_tot,
    diff_daly_low   = pre_daly_low95 - daly_tot_hi,
    diff_daly_hi    = pre_daly_hi - daly_tot_lo,
    impact_daly     = diff_daly / pre_daly * 100,
    impact_daly_low = diff_daly_low / pre_daly_hi * 100,
    impact_daly_hi  = diff_daly_hi / pre_daly_low95 * 100,
    nnv_daly        = tot_vacc / diff_daly,
    nnv_daly_lo     = tot_vacc / diff_daly_hi,
    nnv_daly_hi     = tot_vacc / diff_daly_low
  )


final_nnv_df_22$age_gr <- rep(age_gr[1:18],n_scenarios)
final_nnv_df_22$age_gr <- factor(final_nnv_df_22$age_gr, levels = age_gr_levels)


## nnv graphs
p1 <- 
  nnv_gg(
    final_nnv_df_22,
    y_var = "nnv",  # name of y variable as a string
    y_lab = "NNV to avert a single symptomatic case",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )

p2 <- 
  nnv_gg(
    final_nnv_df_22,
    y_var = "nnv_fatal",  # name of y variable as a string
    y_lab = "NNV to avert a single fatal case",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )

p3 <- 
  nnv_gg(
    final_nnv_df_22,
    y_var = "nnv_daly",  # name of y variable as a string
    y_lab = "NNV to avert a single DALY",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )


ggsave(filename = "02_Outputs/2_1_Figures/nnv_symp_nat_24.jpg", p1, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_fatal_nat_24.jpg", p2, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_daly_nat_24.jpg", p3, width = 10, height = 8, dpi = 1200)

## for faceting
nnv_nat_gg(combined_nnv_df_region_22, 
           y_var = "nnv_fatal",
           y_lab = "NNV to avert a single fatal case",
           x_lab = "Age group",
           title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

p1 <- p1 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single symptomatic case")

p2 <- p2 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single fatal case")

p3 <- p3 + theme(
  axis.text.y  = element_text(size = 10, color = "black"),  # restore y-axis text
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
) + 
  ggtitle("Number need to vaccinate to avert a single DALY")


combined_nnv <- ggarrange(p1, p2, p3,
                          ncol = 1, nrow = 3,
                          labels = c("A", "B", "C"),
                          common.legend = TRUE,
                          legend = "bottom",
                          align = "v")
combined_nnv <- annotate_figure(combined_nnv,
                                left = text_grob("Number needed to vaccinate to avert each outcome", rot = 90, size = 12),
                                bottom = text_grob("Age group", size = 12))
ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_24.jpg", combined_nnv, width = 9, height = 14, dpi = 1200)


# table: pre-post cases -----------------------------------------------------
scenario_result_go <- attach_target_column(df_nnv_go_22,
                                           observed = observed_go_22,
                                           region   = "Goiás")
scenario_result_sg <- attach_target_column(df_nnv_sg_22,
                                           observed = observed_sg_22,
                                           region   = "Sergipe")
scenario_result_sp <- attach_target_column(df_nnv_sp_22,
                                           observed = observed_sp_22,
                                           region   = "São Paulo")
scenario_result_ms <- attach_target_column(df_nnv_ms_22,
                                           observed = observed_ms_22,
                                           region   = "Mato Grosso do Sul")
scenario_result_mh <- attach_target_column(df_nnv_mh_22,
                                           observed = observed_mh_22,
                                           region   = "Maranhão")

combined_scenario_med <- bind_rows(scenario_result_go, scenario_result_sg, scenario_result_sp,
                               scenario_result_ms, scenario_result_mh)

combined_scenario_summ_med <- combined_scenario_med %>%
  group_by(target, scenario) %>%
  summarise(
    pre_cases    = sum(pre_cases, na.rm = TRUE),
    post_cases   = sum(post_cases, na.rm = TRUE),
    post_low95   = sum(post_low95, na.rm = TRUE),
    post_hi95    = sum(post_hi95, na.rm = TRUE),
    pre_low95    = sum(pre_low95, na.rm = TRUE),
    pre_hi95     = sum(pre_hi95, na.rm = TRUE),
    post_fatal   = sum(post_fatal, na.rm = TRUE),
    pre_fatal    = sum(pre_fatal, na.rm = TRUE),
    post_fatal_lo = sum(post_fatal_lo, na.rm = TRUE),
    post_fatal_hi = sum(post_fatal_hi, na.rm = TRUE),
    pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
    pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
    
    tot_vacc     = sum(tot_vacc),
    
    post_daly      = sum(post_daly, na.rm = TRUE),
    pre_daly       = sum(pre_daly, na.rm = TRUE),
    post_daly_lo   = sum(post_daly_lo, na.rm = TRUE),
    post_daly_hi   = sum(post_daly_hi, na.rm = TRUE),
    pre_daly_lo    = sum(pre_daly_lo, na.rm = TRUE),
    pre_daly_hi    = sum(pre_daly_hi, na.rm = TRUE),
    
    .groups      = "drop"
  ) %>%
  mutate(
    ## For Cases:
    diff_case      = pre_cases - post_cases,
    diff_case_low  = pre_low95 - post_low95,
    diff_case_hi   = pre_hi95 - post_hi95,
    impact         = diff_case / pre_cases * 100,
    impact_low     = diff_case_low / pre_low95 * 100,
    impact_hi      = diff_case_hi / pre_hi95 * 100,
    
    ## For Fatal outcomes:
    diff_fatal      = pre_fatal - post_fatal,
    diff_fatal_low  = pre_fatal_lo - post_fatal_lo,
    diff_fatal_hi   = pre_fatal_hi - post_fatal_hi,
    impact_fatal    = diff_fatal / pre_fatal * 100,
    impact_fatal_low = diff_fatal_low / pre_fatal_lo * 100,
    impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,
    
    ## For DALY:
    diff_daly       = pre_daly - post_daly,
    diff_daly_low   = pre_daly_lo - post_daly_lo,
    diff_daly_hi    = pre_daly_hi - post_daly_hi,
    impact_daly     = diff_daly / pre_daly * 100,
    impact_daly_low = diff_daly_low / pre_daly_lo * 100,
    impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100
  )

write_xlsx(combined_scenario_summ_med, path = "00_Data/0_2_Processed/combined_scenario_summ_med.xlsx")




# delay heatmap -----------------------------------------------------------
heatmap_ce <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_ce,
  posterior            = posterior_ce,
  N                    = N_ceara,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Ceará",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ce
)

heatmap_bh <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_bh,
  posterior            = posterior_bh,
  N                    = N_bahia,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Bahia",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_bh
)

heatmap_pn <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pn,
  posterior            = posterior_pn,
  N                    = N_pemam,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Pernambuco",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pn
)

heatmap_rg <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_rg,
  posterior            = posterior_rg,
  N                    = N_rg,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Rio Grande do Norte",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_rg
)

heatmap_pa <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pa,
  posterior            = posterior_pa,
  N                    = N_pa,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Paraíba",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pa
)

heatmap_pi <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_pi,
  posterior            = posterior_pi,
  N                    = N_pi,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Piauí",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_pi
)

heatmap_tc <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_tc,
  posterior            = posterior_tc,
  N                    = N_tc,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Tocantins",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_tc
)

heatmap_ag <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_ag,
  posterior            = posterior_ag,
  N                    = N_ag,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Alagoas",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_ag
)

heatmap_mg <- heatmap_delay(
  target_age_list      = target_age_list,
  observed             = observed_mg,
  posterior            = posterior_mg,
  N                    = N_mg,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  region_name          = "Minas Gerais",
  hosp                 = hosp,
  fatal                = fatal,
  nh_fatal             = nh_fatal,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  pre_summary_cases    = prevacc_ui_mg
)

combined_heatmap <- bind_rows(
  heatmap_ag,
  heatmap_bh,
  heatmap_ce,
  heatmap_mg,
  heatmap_pa,
  heatmap_pi,
  heatmap_pn,
  heatmap_rg,
  heatmap_tc
)

heatmap_summ <- combined_heatmap %>% group_by(scenario_id, Scenario, Supply, Delay) %>% summarise(
  tot_pre_cases  = sum(tot_pre_cases),
  tot_post_cases = sum(tot_post_cases),
  tot_pre_fatal  = sum(tot_pre_fatal),
  tot_post_fatal = sum(tot_post_fatal)
) %>% mutate(
  tot_diff = tot_pre_cases - tot_post_cases,
  impact = tot_diff / tot_pre_cases * 100,
  tot_pre_fatal = sum(tot_pre_fatal, na.rm = TRUE),
  tot_post_fatal = sum(tot_post_fatal, na.rm = TRUE),
  diff_fatal = tot_pre_fatal - tot_post_fatal,
  impact_fatal = diff_fatal / tot_pre_fatal * 100
)

p <- 
  ggplot(heatmap_summ, aes(x = Supply, y = Delay, fill = impact)) +
  geom_tile() +
  # Use a diverging palette with more pronounced color differences
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "% Cases averted") +  # Add color legend breaks
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Total vaccine coverage") +
  facet_wrap(~ Scenario, 
             labeller = as_labeller(c("Scenario_1" = "<20 years only", 
                                      "Scenario_2" = "20-59 years only", 
                                      "Scenario_3" = ">60 years only")))+
  # Adjust x-axis labels to show 10% steps
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format(accuracy = 1))

ggsave(filename = "02_Outputs/2_1_Figures/fig_tot_heatmap.jpg", p, width = 10, height = 4, dpi = 1200)


# table: pre-post cases -----------------------------------------------------
scenario_result_sg <- attach_target_column(df_nnv_sg_22,
                                           observed = observed_sg_22,
                                           region   = "Sergipe")
scenario_result_sp <- attach_target_column(df_nnv_sp_22,
                                           observed = observed_sp_22,
                                           region   = "São Paulo")
scenario_result_go <- attach_target_column(df_nnv_go_22,
                                           observed = observed_go_22,
                                           region   = "Goiás")

combined_scenario <- bind_rows(scenario_result_sg, scenario_result_sp, scenario_result_go)

combined_scenario_summ <- combined_scenario %>%
  group_by(target, scenario) %>%
  summarise(
    pre_cases    = sum(pre_cases, na.rm = TRUE),
    post_cases   = sum(post_cases, na.rm = TRUE),
    post_low95   = sum(post_low95, na.rm = TRUE),
    post_hi95    = sum(post_hi95, na.rm = TRUE),
    pre_low95    = sum(pre_low95, na.rm = TRUE),
    pre_hi95     = sum(pre_hi95, na.rm = TRUE),
    post_fatal   = sum(post_fatal, na.rm = TRUE),
    pre_fatal    = sum(pre_fatal, na.rm = TRUE),
    post_fatal_lo = sum(post_fatal_lo, na.rm = TRUE),
    post_fatal_hi = sum(post_fatal_hi, na.rm = TRUE),
    pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
    pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
    
    tot_vacc     = sum(tot_vacc),
    
    post_daly      = sum(post_daly, na.rm = TRUE),
    pre_daly       = sum(pre_daly, na.rm = TRUE),
    post_daly_lo   = sum(post_daly_lo, na.rm = TRUE),
    post_daly_hi   = sum(post_daly_hi, na.rm = TRUE),
    pre_daly_lo    = sum(pre_daly_lo, na.rm = TRUE),
    pre_daly_hi    = sum(pre_daly_hi, na.rm = TRUE),
    
    .groups      = "drop"
  ) %>%
  mutate(
    ## For Cases:
    diff_case      = pre_cases - post_cases,
    diff_case_low  = pre_low95 - post_low95,
    diff_case_hi   = pre_hi95 - post_hi95,
    impact         = diff_case / pre_cases * 100,
    impact_low     = diff_case_low / pre_low95 * 100,
    impact_hi      = diff_case_hi / pre_hi95 * 100,
    
    ## For Fatal outcomes:
    diff_fatal      = pre_fatal - post_fatal,
    diff_fatal_low  = pre_fatal_lo - post_fatal_lo,
    diff_fatal_hi   = pre_fatal_hi - post_fatal_hi,
    impact_fatal    = diff_fatal / pre_fatal * 100,
    impact_fatal_low = diff_fatal_low / pre_fatal_lo * 100,
    impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,
    
    ## For DALY:
    diff_daly       = pre_daly - post_daly,
    diff_daly_low   = pre_daly_lo - post_daly_lo,
    diff_daly_hi    = pre_daly_hi - post_daly_hi,
    impact_daly     = diff_daly / pre_daly * 100,
    impact_daly_low = diff_daly_low / pre_daly_lo * 100,
    impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100
  )

write_xlsx(combined_scenario_summ, path = "00_Data/0_2_Processed/combined_scenario_summ_v2.xlsx")

# ------------------------------------------------------------------------------

## fatal summary
cum_fatal_go <- fatal_summ_func(postsim_all_go_ui_22,
                                N_go,
                                observed_go_22,
                                "Goiás")
cum_fatal_sg <- fatal_summ_func(postsim_all_sg_ui_22,
                                N_sg,
                                observed_sg_22,
                                "Sergipe")
cum_fatal_sp <- fatal_summ_func(postsim_all_sp_ui_22,
                                N_sp,
                                observed_sp_22,
                                "São Paulo")

cum_fatal_ms <- fatal_summ_func(postsim_all_ms_ui_22,
                                N_ms,
                                observed_ms_22,
                                "Mato Grosso do Sul")

cum_fatal_mh <- fatal_summ_func(postsim_all_mh_ui_22,
                                N_mh,
                                observed_mh_22,
                                "Maranhão")

combined_prepost_fatal <- bind_rows(
  cum_fatal_go$fatal_week_df,
  cum_fatal_sg$fatal_week_df,
  cum_fatal_sp$fatal_week_df,
  cum_fatal_ms$fatal_week_df,
  cum_fatal_mh$fatal_week_df
)

global_impact_fatal <- combined_prepost_fatal %>%
  group_by(Scenario, region) %>%
  summarise(
    total_post_fatal = sum(fatal),
    total_pre_fatal  = sum(pre_fatal),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_fatal - total_post_fatal,
    impact = diff / total_pre_fatal * 100
  )

annotation_df <- global_impact_fatal %>%
  group_by(region) %>%
  summarise(ann = paste0(Scenario, ": ", round(impact, 1), "%", collapse = "\n"),
            .groups = "drop")

overall_max_med <- max(combined_prepost_fatal$pre_fatal_prop_hi, na.rm = TRUE)

# Create an annotation data frame that has one row per region
vaccine_ann_df <- combined_prepost_fatal %>%
  distinct(region) %>%
  mutate(x = 4, 
         y = overall_max * 0.95, 
         label = "<----- Vaccine impact start (+ 2 wks after initiation)")

pre_post_fatal_med <- 
  
  ggplot(combined_prepost_fatal) +
  # Scenario shading ribbons (existing)
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_1", Week >= vacc_start_week_s1 & Week <= vacc_end_week_s1),
              aes(x = Week, ymin = 2/3 * overall_max_med, ymax = overall_max_med, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_2", Week >= vacc_start_week_s2 & Week <= vacc_end_week_s2),
              aes(x = Week, ymin = 1/3 * overall_max_med, ymax = 2/3 * overall_max_med, fill = Scenario),
              alpha = 0.4) +
  geom_ribbon(data = combined_prepost_fatal %>% filter(Scenario == "Scenario_3", Week >= vacc_start_week_s3 & Week <= vacc_end_week_s3),
              aes(x = Week, ymin = 0, ymax = 1/3 * overall_max_med, fill = Scenario),
              alpha = 0.4) +
  # Post-vaccination UI ribbon (95% uncertainty interval)
  geom_ribbon(aes(x = Week, ymin = post_fatal_prop_lo, ymax = post_fatal_prop_hi, fill = Scenario),
              alpha = 0.2) +
  # Pre-vaccination UI ribbon (if available; omit if not)
  geom_ribbon(aes(x = Week, ymin = pre_fatal_prop_lo, ymax = pre_fatal_prop_hi),
              fill = "gray", alpha = 0.2) +
  # Post-vaccination case line (colored by scenario)
  geom_line(aes(x = Week, y = post_fatal_prop, color = Scenario, linetype = "Post Cases"), size = 0.5) +
  # Pre-vaccination case line (dashed, in black)
  geom_line(aes(x = Week, y = pre_fatal_prop, group = Scenario, linetype = "Pre Cases"), color = "black", size = 0.5) +
  scale_linetype_manual(values = c("Post Cases" = "solid", "Pre Cases" = "dashed")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("Scenario_1" = "<20 years only", 
                               "Scenario_2" = "20-59 years only", 
                               "Scenario_3" = ">60 years only")) +
  scale_color_brewer(palette = "Set1",
                     labels = c("Scenario_1" = "<20 years only", 
                                "Scenario_2" = "20-59 years only", 
                                "Scenario_3" = ">60 years only")) +
  #scale_y_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) + 
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 0.4) +
  labs(color = "Scenario", linetype = "Type", fill = "Scenario",
       title = "Coverage: 50%, Delivery Speed: 10%, Deployment: Week 2",
       x = "Week", y = "Predicted cumulative reported fatal cases (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(5, 30, 5, 5)
  ) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = Inf, label = ann), 
            inherit.aes = FALSE,
            hjust = 1.1, 
            vjust = 1.1, 
            size = 3) +
  #geom_text(data = vaccine_ann_df,
  #          aes(x = x, y = y, label = label),
  #          inherit.aes = FALSE,
  #          hjust = 0, vjust = 1, size = 2)+
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 30, 5, 5)) +
  #geom_point(data = observed_all, aes(x = Week, y = Observed), size = 0.5) + 
  facet_wrap(~region)

ggsave(filename = "02_Outputs/2_1_Figures/pre_post_fatal_med.jpg", pre_post_fatal_med, width = 12, height = 6, dpi = 1200)

combined_prepost_fatal_nat <- combined_prepost_fatal %>%
  group_by(Scenario, Week) %>%
  summarise(
    fatal    = sum(fatal),
    fatal_lo = sum(fatal_lo),
    fatal_hi = sum(fatal_hi),
    pre_fatal = sum(pre_fatal),
    pre_fatal_lo = sum(pre_fatal_lo),
    pre_fatal_hi = sum(pre_fatal_hi),
    tot_pop          = sum(tot_pop),
    .groups = "drop"
  ) %>% arrange(Scenario, Week) %>%
  group_by(Scenario) %>%
  mutate(
    # Now do cumsum within each Scenario in ascending Week order
    cum_fatal         = cumsum(fatal),
    cum_fatal_lo      = cumsum(fatal_lo),
    cum_fatal_hi      = cumsum(fatal_hi),
    cum_fatal_pre     = cumsum(pre_fatal),
    cum_fatal_pre_lo  = cumsum(pre_fatal_lo),
    cum_fatal_pre_hi  = cumsum(pre_fatal_hi)
  ) %>%
  ungroup() %>%
  mutate(
    pre_fatal_prop      = (cum_fatal_pre / tot_pop) * 100,
    pre_fatal_prop_lo   = (cum_fatal_pre_lo / tot_pop) * 100,
    pre_fatal_prop_hi   = (cum_fatal_pre_hi / tot_pop) * 100,
    post_fatal_prop     = (cum_fatal / tot_pop) * 100,
    post_fatal_prop_lo  = (cum_fatal_lo / tot_pop) * 100,
    post_fatal_prop_hi  = (cum_fatal_hi / tot_pop) * 100
  )

global_impact_fatal_nat <- combined_prepost_fatal_nat %>%
  group_by(Scenario) %>%
  summarise(
    total_post_fatal = sum(fatal),
    total_pre_fatal  = sum(pre_fatal),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_fatal - total_post_fatal,
    impact = diff / total_pre_fatal * 100
  )

global_impact_fatal_nat <- combined_prepost_fatal %>%
  group_by(Scenario) %>%
  summarise(
    total_post_fatal = sum(cum_fatal),  
    total_post_fatal_lo = sum(cum_fatal_lo),
    total_post_fatal_hi = sum(cum_fatal_hi),
    total_pre_fatal   = sum(cum_fatal_pre),
    total_pre_fatal_lo = sum(cum_fatal_pre_lo),
    total_pre_fatal_hi = sum(cum_fatal_pre_hi),
    .groups = "drop"
  ) %>%
  mutate(
    diff_fatal   = `total_pre_fatal` - `total_post_fatal`,
    impact_fatal = diff_fatal / `total_pre_fatal`  * 100,
    diff_fatal_lo = `total_pre_fatal_lo` - `total_post_fatal_lo`,
    impact_fatal_lo = diff_fatal_lo / `total_pre_fatal_lo`  * 100,
    diff_fatal_hi = `total_pre_fatal_hi` - `total_post_fatal_hi`,
    impact_fatal_hi = diff_fatal_hi / `total_pre_fatal_hi`  * 100
  )

annotation_text_fatal <- paste0(
  global_impact_fatal_nat$Scenario, ": ", 
  round(global_impact_fatal_nat$impact_fatal, 1), "%",
  collapse = "\n"
)

cum_fatal_nat <- list(
  fatal_week_df   = combined_prepost_fatal_nat,
  annotation_text_fatal = annotation_text_fatal
)


fatal_graph_nat_med <- 
  cum_fatal_plot(cum_fatal_nat, postsim_all_sg_ui_22)


# combining
pre_post_fatal_med <- pre_post_fatal_med +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Vaccine impact at sub-national level")

# Remove legend from epi_graph_all (Scenario A)
fatal_graph_nat_med <- fatal_graph_nat_med + theme(legend.position = "none") + 
  ggtitle("Vaccine impact at national level")

# Arrange plots side-by-side with a common legend
combined_fatal_graph <- ggarrange(fatal_graph_nat_med, pre_post_fatal_med,
                                  ncol = 2, nrow = 1,
                                  labels = c("C", "D"),
                                  common.legend = TRUE,
                                  legend = "bottom",
                                  align = "none")

combined_fatal_graph <- combined_fatal_graph +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(filename = "02_Outputs/2_1_Figures/combined_brazil_fatal.jpg", combined_fatal_graph, width = 12, height = 6, dpi = 1200)

