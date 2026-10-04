combined_df_2022 <- bind_rows(df_ce, df_bh, df_pa,
                              df_pn, df_rg, df_pi,
                              df_tc, df_ag, df_mg,
                              df_sp_22, df_go_22, 
                              df_sg_22, df_ms_22, df_mh_22) %>%
  group_by(Scenario, Week) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

global_impact_week <- combined_df_2022 %>%
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

postsim_all_global_ui_22 <- list(
  postsim_all_sg_ui_22 = postsim_all_sg_ui_22,
  postsim_all_sg_ui_22 = postsim_all_go_ui_22,
  summary_week_df_22   = combined_df_2022,
  global_impact_week = global_impact_week,
  annotation_text = annotation_text,
  vacc_start_week_s1 = postsim_all_ce_ui$vacc_weeks$scenario1,
  vacc_start_week_s2 = postsim_all_ce_ui$vacc_weeks$scenario2,
  vacc_start_week_s3 = postsim_all_ce_ui$vacc_weeks$scenario3
  
)

tot_sum <- bra_all_sum$Observed + bra_all_sum_22_med$Observed
bra_all_sum_22 <- bra_all_sum[,1:2]
bra_all_sum_22$Observed <- tot_sum

epi_graph_all_22 <- 
  epi_graph_nat(postsim_all_global_ui_22, bra_all_sum_22)


# combined fatal
combined_prepost_fatal_all <- bind_rows(
  cum_fatal_go$fatal_week_df,
  cum_fatal_sg$fatal_week_df,
  cum_fatal_sp$fatal_week_df,
  cum_fatal_ms$fatal_week_df,
  cum_fatal_mh$fatal_week_df,
  cum_fatal_ce$fatal_week_df,
  cum_fatal_bh$fatal_week_df,
  cum_fatal_ag$fatal_week_df,
  cum_fatal_mg$fatal_week_df,
  cum_fatal_pa$fatal_week_df,
  cum_fatal_pi$fatal_week_df,
  cum_fatal_pn$fatal_week_df,
  cum_fatal_rg$fatal_week_df,
  cum_fatal_tc$fatal_week_df
)


combined_prepost_fatal_2022 <- combined_prepost_fatal_all %>%
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

global_impact_fatal_2022 <- combined_prepost_fatal_2022 %>%
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

annotation_text_fatal_2022 <- paste0(
  global_impact_fatal_2022$Scenario, ": ", 
  round(global_impact_fatal_2022$impact_fatal, 1), "%",
  collapse = "\n"
)

cum_fatal_nat_2022 <- list(
  fatal_week_df   = combined_prepost_fatal_2022,
  annotation_text_fatal = annotation_text_fatal_2022
)

fatal_graph_nat_med <- 
  cum_fatal_plot(cum_fatal_nat_2022, postsim_all_sg_ui_22)

epi_graph_all_22 <- epi_graph_all_22 + ggtitle("A. Vaccine impact at national level (% symptomatic case reduction)")
fatal_graph_nat_med <- fatal_graph_nat_med + ggtitle("B. Vaccine impact at national level (% cumulative fatal reduction)")

comb_symp_fatal_2022 <- 
ggarrange(epi_graph_all_22, fatal_graph_nat_med,
          ncol = 2, nrow = 1,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "right",
          align = "none")

ggsave(filename = "02_Outputs/2_1_Figures/comb_symp_fatal_2022.jpg", comb_symp_fatal_2022, width = 12, height = 6, dpi = 1200)


## combined graph
pre_post_graph_comb_22_med <- pre_post_graph_22_med +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Vaccine impact at sub-national level")

pre_post_graph <- pre_post_graph +
  theme(
    axis.title.x = element_blank()
  )

pre_post_graph <- pre_post_graph + ggtitle("A. Vaccine impact in high-burden states (% symptomatic case reduction)")
pre_post_graph_comb_22_med <- pre_post_graph_comb_22_med + ggtitle("B. Vaccine impact in moderate-low burden states (% symptomatic case reduction)")

# Remove legend from epi_graph_all (Scenario A)
epi_graph_all_22_med <- epi_graph_all_22_med + theme(legend.position = "none") + 
  ggtitle("Vaccine impact at national level")

# Arrange plots side-by-side with a common legend
combined_graph_22_med <- ggarrange(pre_post_graph, pre_post_graph_comb_22_med,
                                   ncol = 2, nrow = 1,
                                   labels = c("C", "D"),
                                   common.legend = TRUE,
                                   legend = "right",
                                   align = "none")

combined_graph_22_med <- annotate_figure(
  combined_graph_22_med,
  bottom = text_grob("Week", size = 12)
)

ggsave(filename = "02_Outputs/2_1_Figures/combined_graph_22_med.jpg", combined_graph_22_med, width = 12, height = 6, dpi = 1200)


## combined fatal
pre_post_fatal_graph <- pre_post_fatal_graph + ggtitle("A. Vaccine impact in high-burden states (% cumulative fatal reduction)")
pre_post_fatal_med <- pre_post_fatal_med + ggtitle("B. Vaccine impact in moderate-low burden states (% cumulative fatal reduction)")

combined_fatal_22_med <- ggarrange(pre_post_fatal_graph, pre_post_fatal_med,
                                   ncol = 2, nrow = 1,
                                   labels = c("E", "F"),
                                   common.legend = TRUE,
                                   legend = "right",
                                   align = "none")
symp_fatal_comb_2022 <- 
ggarrange(pre_post_graph, pre_post_graph_comb_22_med,
          pre_post_fatal_graph, pre_post_fatal_med,
          ncol = 2, nrow = 2,
          labels = c("C", "D", "E", "F"),
          common.legend = TRUE,
          legend = "right",
          align = "none")

ggsave(filename = "02_Outputs/2_1_Figures/combined_fatal_22_med.jpg", combined_fatal_22_med, width = 12, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/symp_fatal_comb_2022.jpg", symp_fatal_comb_2022, width = 15, height = 10, dpi = 1200)



#### nnv 
## total nnv graph (national level) ------------------------------
n_scenarios = 3

combined_nnv_df_2022 <- bind_rows(df_nnv_ce, df_nnv_bh, df_nnv_pa,
                                    df_nnv_pn, df_nnv_rg, df_nnv_pi,
                                    df_nnv_tc, df_nnv_ag, df_nnv_mg,
                                    df_nnv_sp_22, df_nnv_go_22, df_nnv_sg_22,
                                    df_nnv_ms_22, df_nnv_mh_22) 


combined_nnv_region_2022 <- combined_nnv_df_2022 %>%
  group_by(scenario, region) %>%
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

write_xlsx(combined_nnv_region_all, path = "02_Outputs/2_2_Tables/combined_nnv_region_all.xlsx")


write_xlsx(combined_nnv_df_region, path = "02_Outputs/2_2_Tables/combined_nnv_df_region.xlsx")

combined_nnv_2022 <- combined_nnv_df_2022 %>%
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
final_nnv_df_2022 <- combined_nnv_2022 %>%
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


final_nnv_df_2022$age_gr <- rep(age_gr[1:18],n_scenarios)
final_nnv_df_2022$age_gr <- factor(final_nnv_df_2022$age_gr, levels = age_gr_levels)


save("final_nnv_df", file = "02_Outputs/2_2_Tables/final_nnv_df.RData")
write_xlsx(final_nnv_df, path = "02_Outputs/2_2_Tables/final_nnv_df.xlsx")


## nnv graphs
p1 <- 
  nnv_gg(
    final_nnv_df_2022,
    y_var = "nnv",  # name of y variable as a string
    y_lab = "NNV to avert a single symptomatic case",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )

p2 <- 
  nnv_gg(
    final_nnv_df_2022,
    y_var = "nnv_fatal",  # name of y variable as a string
    y_lab = "NNV to avert a single fatal case",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )

p3 <- 
  nnv_gg(
    final_nnv_df_2022,
    y_var = "nnv_daly",  # name of y variable as a string
    y_lab = "NNV to avert a single DALY",
    x_lab = "Age group",
    title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2"
  )

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


combined_nnv_2022 <- ggarrange(p1, p2, p3,
                          ncol = 1, nrow = 3,
                          labels = c("A", "B", "C"),
                          common.legend = TRUE,
                          legend = "bottom",
                          align = "v")
combined_nnv_2022 <- annotate_figure(combined_nnv_2022,
                                left = text_grob("Number needed to vaccinate to avert each outcome", rot = 90, size = 12),
                                bottom = text_grob("Age group", size = 12))


ggsave(filename = "02_Outputs/2_1_Figures/nnv_symp_nat.jpg", p1, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_fatal_nat.jpg", p2, width = 10, height = 8, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/nnv_daly_nat.jpg", p3, width = 10, height = 8, dpi = 1200)

## for faceting
g1 <- 
  nnv_nat_gg(combined_nnv_df_2022, 
             y_var = "nnv",
             y_lab = "NNV to avert a single symptomatic case",
             x_lab = "Age group",
             title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")
g2 <-
  nnv_nat_gg(combined_nnv_df_region, 
             y_var = "nnv_fatal",
             y_lab = "NNV to avert a single fatal case",
             x_lab = "Age group",
             title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

g3 <- 
  nnv_nat_gg(combined_nnv_df_region, 
             y_var = "nnv_daly",
             y_lab = "NNV to avert a single DALY",
             x_lab = "Age group",
             title = "Coverage:50%, Delivery speed: 10%, Deployment: Week 2")

g1 <- g1 +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single symptomatic case")

g2 <- g2 +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single fatal case")

g3 <- g3 +
  scale_x_discrete(labels = function(x) gsub(" years", "", x)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black")  # restore y-axis text
  ) + 
  ggtitle("Number needed to vaccinate to avert a single DALY")

# Arrange plots side-by-side with a common legend
combined_nnv_graph <- ggarrange(g1, g2, g3,
                                ncol = 1, nrow = 3,
                                labels = c("C", "D", "E"),
                                common.legend = TRUE,
                                legend = "bottom",
                                align = "v")

combined_nnv_graph <- annotate_figure(combined_nnv_graph,
                                      left = text_grob("Number needed to vaccinate to avert each outcome", rot = 90, size = 12),
                                      bottom = text_grob("Age group", size = 12))

ggsave(filename = "02_Outputs/2_1_Figures/combined_nnv_region.jpg", combined_nnv_graph, width = 11, height = 8, dpi = 1200)


## combined table
combined_scenario_2022 <- bind_rows(scenario_result_go, scenario_result_sg, scenario_result_sp,
                                   scenario_result_ms, scenario_result_mh, scenario_result_ce, scenario_result_bh, scenario_result_pa,
                                   scenario_result_pn, scenario_result_pi, scenario_result_rg,
                                   scenario_result_tc, scenario_result_ag, scenario_result_mg)

combined_scenario_summ_2022 <- combined_scenario_2022 %>%
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

write_xlsx(combined_scenario_summ_2022, path = "00_Data/0_2_Processed/combined_scenario_summ_2022.xlsx")


# total supply 
combined_vacc_alloc_2022 <- bind_rows(
  vacc_alloc_ce$weekly_allocation_long,
  vacc_alloc_bh$weekly_allocation_long,
  vacc_alloc_pa$weekly_allocation_long,
  vacc_alloc_pn$weekly_allocation_long,
  vacc_alloc_rg$weekly_allocation_long,
  vacc_alloc_pi$weekly_allocation_long,
  vacc_alloc_tc$weekly_allocation_long,
  vacc_alloc_ag$weekly_allocation_long,
  vacc_alloc_mg$weekly_allocation_long,
  vacc_alloc_sp_22$weekly_allocation_long,
  vacc_alloc_go_22$weekly_allocation_long,
  vacc_alloc_sg_22$weekly_allocation_long,
  vacc_alloc_ms_22$weekly_allocation_long,
  vacc_alloc_mh_22$weekly_allocation_long
)

combined_vacc_alloc_summ_2022 <- combined_vacc_alloc_2022 %>%
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
  ggplot(combined_vacc_alloc_summ_2022, aes(x = Week, y = Vaccinated, fill = factor(age_gr))) +
  geom_bar(stat = "identity") +   # Bars are automatically stacked
  labs(
    x = "Week",
    y = "Total doses of vaccines distributed per age group (millions)",
    fill = "Age Group"
  ) +
  scale_fill_manual(values = full_palette) + 
  scale_y_continuous(label=scales::label_number(scale_cut = scales::cut_short_scale()),
                     breaks = seq(0, 3000000, by = 500000))+
  facet_wrap(~ Scenario, 
             labeller = labeller(Scenario = c("Scenario1" = "<20 years only", 
                                              "Scenario2" = "20-59 years only", 
                                              "Scenario3" = ">60 years only")))  +         # Facet by Scenario if you want separate plots for each
  theme_minimal()

ggsave(filename = "02_Outputs/2_1_Figures/fig_combined_vacc_alloc_summ_2022.jpg", p, width = 12, height = 7, dpi = 1200)


