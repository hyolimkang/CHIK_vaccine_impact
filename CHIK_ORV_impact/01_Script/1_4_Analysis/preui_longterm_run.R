# prevacc 95%UI simulation ----------------------------------------------------
preui_ce <- simulate_pre_ui_age(posterior = posterior_ce, bra_foi_state_summ, age_groups,
                                N = N_ceara$Ceará,
                                region = "Ceará",
                                observed = observed_ce,
                                lhs_sample = lhs_combined$Ceará,
                                finite_history = FALSE
)

preui_bh <- simulate_pre_ui_age(posterior = posterior_bh, bra_foi_state_summ, age_groups,
                                N = N_bahia$Bahia,
                                region = "Bahia",
                                observed = observed_bh,
                                lhs_sample = lhs_combined$Bahia,
                                finite_history = FALSE
)

preui_pa <- simulate_pre_ui_age(posterior = posterior_pa, bra_foi_state_summ, age_groups,
                                N = N_pa$Paraíba,
                                region = "Paraíba",
                                observed = observed_pa,
                                lhs_sample = lhs_combined$Paraíba,
                                finite_history = FALSE
)

preui_pn <- simulate_pre_ui_age(posterior = posterior_pn, bra_foi_state_summ, age_groups,
                                N = N_pemam$Pernambuco,
                                region = "Pernambuco",
                                observed = observed_pn,
                                lhs_sample = lhs_combined$Pernambuco,
                                finite_history = FALSE
)

preui_rg <- simulate_pre_ui_age(posterior = posterior_rg, bra_foi_state_summ, age_groups,
                                N = N_rg$`Rio Grande do Norte`,
                                region = "Rio Grande do Norte",
                                observed = observed_rg,
                                lhs_sample = lhs_combined$`Rio Grande do Norte`,
                                finite_history = FALSE
)

preui_pi <- simulate_pre_ui_age(posterior = posterior_pi, bra_foi_state_summ, age_groups,
                                N = N_pi$Piauí,
                                region = "Piauí",
                                observed = observed_pi,
                                lhs_sample = lhs_combined$Piauí,
                                finite_history = FALSE
)

preui_tc <- simulate_pre_ui_age(posterior = posterior_tc, bra_foi_state_summ, age_groups,
                                N = N_tc$Tocantins,
                                region = "Tocantins",
                                observed = observed_tc,
                                lhs_sample = lhs_combined$Tocantins,
                                finite_history = FALSE
)

preui_ag <- simulate_pre_ui_age(posterior = posterior_ag, bra_foi_state_summ, age_groups,
                                N = N_ag$Alagoas,
                                region = "Alagoas",
                                observed = observed_ag,
                                lhs_sample = lhs_combined$Alagoas,
                                finite_history = FALSE
)

preui_mg <- simulate_pre_ui_age(posterior = posterior_mg, bra_foi_state_summ, age_groups,
                                N = N_mg$`Minas Gerais`,
                                region = "Minas Gerais",
                                observed = observed_mg,
                                lhs_sample = lhs_combined$`Minas Gerais`,
                                finite_history = FALSE
)

preui_se <- simulate_pre_ui_age(posterior = posterior_se, bra_foi_state_summ, age_groups,
                                N = N_se$Sergipe,
                                region = "Sergipe",
                                observed = observed_se,
                                lhs_sample = lhs_combined$Sergipe,
                                finite_history = FALSE
)

preui_go <- simulate_pre_ui_age(posterior = posterior_go, bra_foi_state_summ, age_groups,
                                N = N_go$Goiás,
                                region = "Goiás",
                                observed = observed_go,
                                lhs_sample = lhs_combined$Goiás,
                                finite_history = FALSE
)

save(
  preui_ce, preui_bh, preui_pa, preui_pn, preui_rg, preui_pi,
  preui_ag, preui_tc, preui_mg, preui_se, preui_go,
  file = "C:/Users/user/OneDrive/CHIK_benefit_risk/01_Data/preui_longterm_all.RData"
)
