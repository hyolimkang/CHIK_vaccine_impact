# prevacc 95%UI simulation ----------------------------------------------------
preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups, 
                                N = N_ceara$Ceará, 
                                region = "Ceará",
                                observed = observed_ce,
                                lhs_sample = lhs_combined_finite$Ceará
)

preui_bh <- simulate_pre_ui_age(posterior = posterior_bh, bra_foi_state_summ, age_groups, 
                                N = N_bahia$Bahia, 
                                region = "Bahia",
                                observed = observed_bh,
                                lhs_sample = lhs_combined_finite$Bahia
)

preui_pa <- simulate_pre_ui_age(posterior = posterior_pa, bra_foi_state_summ, age_groups, 
                                N = N_pa$Paraíba, 
                                region = "Paraíba",
                                observed = observed_pa,
                                lhs_sample = lhs_combined_finite$Paraíba
)

preui_pn <- simulate_pre_ui_age(posterior = posterior_pn, bra_foi_state_summ, age_groups, 
                                N = N_pemam$Pernambuco, 
                                region = "Pernambuco",
                                observed = observed_pn,
                                lhs_sample = lhs_combined_finite$Pernambuco
)

preui_rg <- simulate_pre_ui_age(posterior = posterior_rg, bra_foi_state_summ, age_groups, 
                                N = N_rg$`Rio Grande do Norte`, 
                                region = "Rio Grande do Norte",
                                observed = observed_rg,
                                lhs_sample = lhs_combined_finite$`Rio Grande do Norte`
)

preui_pi <- simulate_pre_ui_age(posterior = posterior_pi, bra_foi_state_summ, age_groups, 
                                N = N_pi$Piauí, 
                                region = "Piauí",
                                observed = observed_pi,
                                lhs_sample = lhs_combined_finite$Piauí
)

preui_tc <- simulate_pre_ui_age(posterior = posterior_tc, bra_foi_state_summ, age_groups, 
                                N = N_tc$Tocantins, 
                                region = "Tocantins",
                                observed = observed_tc,
                                lhs_sample = lhs_combined_finite$Tocantins
)

preui_ag <- simulate_pre_ui_age(posterior = posterior_ag, bra_foi_state_summ, age_groups, 
                                N = N_ag$Alagoas, 
                                region = "Alagoas",
                                observed = observed_ag,
                                lhs_sample = lhs_combined_finite$Alagoas
)

preui_mg <- simulate_pre_ui_age(posterior = posterior_mg, bra_foi_state_summ, age_groups, 
                                N = N_mg$`Minas Gerais`, 
                                region = "Minas Gerais",
                                observed = observed_mg,
                                lhs_sample = lhs_combined_finite$`Minas Gerais`
)

preui_se <- simulate_pre_ui_age(posterior = posterior_se, bra_foi_state_summ, age_groups, 
                                N = N_se$Sergipe, 
                                region = "Sergipe",
                                observed = observed_se,
                                lhs_sample = lhs_combined_finite$Sergipe
)

preui_go <- simulate_pre_ui_age(posterior = posterior_go, bra_foi_state_summ, age_groups, 
                                N = N_go$Goiás, 
                                region = "Goiás",
                                observed = observed_go,
                                lhs_sample = lhs_combined_finite$Goiás
)

save(
  preui_ce, preui_bh, preui_pa, preui_pn, preui_rg, preui_pi,
  preui_ag, preui_tc, preui_mg, preui_se, preui_go,
  file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/preui_flat_all.RData"
)

save(
  preui_ce, preui_bh, preui_pa, preui_pn, preui_rg, preui_pi,
  preui_ag, preui_tc, preui_mg, preui_se, preui_go,
  file = "00_Data/0_2_Processed/preui_flat_all.RData"
)
