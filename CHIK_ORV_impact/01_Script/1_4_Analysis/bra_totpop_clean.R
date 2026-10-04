load("00_Data/0_2_Processed/ibge_agegr.RData")

tot_pop_age <- rowSums(ibge_agegr_wide[, 2:29])
ibge_agegr_wide$tot_pop <- tot_pop_age

pop_age_nat <- ibge_agegr_wide[, c(1, 30)]

N_bahia <- ibge_agegr_wide[,6]
N_ceara <- ibge_agegr_wide[,8]
N_mg    <- ibge_agegr_wide[,15]
N_pemam <- ibge_agegr_wide[,19]
N_pa    <- ibge_agegr_wide[,17]
N_rg    <- ibge_agegr_wide[,21]
N_pi    <- ibge_agegr_wide[,20]
N_ag    <- ibge_agegr_wide[,3]
N_tc    <- ibge_agegr_wide[,29]
N_se <- ibge_agegr_wide[,27]
N_go <- ibge_agegr_wide[,11]

age_group_levels <- c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
)

scenario <- c("s_0", 
              rep("s_1", 3),
              "s_2",
              rep("s_3",9),
              rep("s_4", 6)
              )

N_bra <- cbind(N_bahia, N_ceara, N_mg, N_pemam, N_pa, N_rg, N_pi, N_ag, N_tc, N_se, N_go)
N_bra <- N_bra %>% mutate(
  age_gr = age_gr_levels,
  scenario = scenario
)

region_totals <- N_bra %>%
  summarise(
    across(
      Bahia:Goiás,        
      ~ sum(.x, na.rm = TRUE)
    )
  )

scenario_region <- N_bra %>%
  group_by(scenario) %>%
  summarise(
    across(
      Bahia:Goiás,        # 같은 지역 컬럼을 지정
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  ungroup()

scenario_pct <- scenario_region %>%
  mutate(
    across(
      Bahia:Goiás,
      ~ .x / region_totals[[cur_column()]] * 100
    )
  )


row_s2 <- scenario_pct[scenario_pct$scenario == "s_2", ]
region_coverage <- as.numeric(row_s2[ , -1])

# 3) 이름 붙이기
names(region_coverage) <- colnames(scenario_pct)[-1]

save(region_coverage, file = "00_Data/0_2_Processed/region_coverage.RData")

save(N_bahia, N_ceara, N_mg, N_pemam, N_pa, N_rg, N_pi, N_ag, N_tc, N_se, N_go, file = "00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")

N_ac <- ibge_agegr_wide[,2]
N_ap <- ibge_agegr_wide[,4]
N_am <- ibge_agegr_wide[,5]
N_df <- ibge_agegr_wide[,9]
N_es <- ibge_agegr_wide[,10]
N_go <- ibge_agegr_wide[,11]
N_ma <- ibge_agegr_wide[,12]
N_mt <- ibge_agegr_wide[,13]
N_ms <- ibge_agegr_wide[,14]
N_pr <- ibge_agegr_wide[,16]
N_para <- ibge_agegr_wide[,18]
N_rs <- ibge_agegr_wide[,22]
N_rj <- ibge_agegr_wide[,23]
N_ro <- ibge_agegr_wide[,24]
N_rr <- ibge_agegr_wide[,25]
N_sc <- ibge_agegr_wide[,26]
N_se <- ibge_agegr_wide[,27]
N_sp <- ibge_agegr_wide[,28]
N_tot <- ibge_agegr_wide[,7]

save(N_ac, N_ap, N_am,
     N_df, N_es, N_go,
     N_ma, N_mt, N_ms,
     N_pr, N_para, N_rs,
     N_rj, N_ro, N_rr, 
     N_sc, N_se, N_sp,
     N_tot, 
     file = "00_Data/0_2_Processed/bra_pop_2022_cleaned_otherstates.RData")

tot_pop <- colSums(ibge_agegr_wide[, 2:29])
tot_pop_df <- data.frame(
  state_full = names(tot_pop),
  tot_pop = as.numeric(tot_pop)
)
tot_pop_df <- tot_pop_df[-6,]

save(tot_pop_df, file = "00_Data/0_2_Processed/tot_pop_df.RData")


tot_pop_df <- tot_pop_df %>%
  mutate(outbreak_status = if_else(
    state_full %in% unique(chikv_all$state_full),
    "outbreak",
    "no outbreak"
  ))

tot_pop_df_ob <- tot_pop_df %>% filter(outbreak_status == "outbreak")



