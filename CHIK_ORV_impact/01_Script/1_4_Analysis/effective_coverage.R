### same dose
cf_dos1 <- sim_results_vc_dose_vimkun[["Ceará"]][["VE0"]][["cov5"]][[1]]$
  sim_out$coverage_frac  

total_eff_dos1 <- colSums(cf_dos1 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_dos1[52]

cf_dos2 <- sim_results_vc_dose_vimkun[["Ceará"]][["VE0"]][["cov5"]][[2]]$
  sim_out$coverage_frac  

total_eff_dos2 <- colSums(cf_dos2 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_dos2[52]


cf_dos3 <- sim_results_vc_dose_vimkun[["Ceará"]][["VE0"]][["cov5"]][[3]]$
  sim_out$coverage_frac  

total_eff_dos3 <- colSums(cf_dos3 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_dos3[52]


cf_dos4 <- sim_results_vc_dose_vimkun[["Ceará"]][["VE0"]][["cov5"]][[4]]$
  sim_out$coverage_frac  

total_eff_dos4 <- colSums(cf_dos4 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_dos4[52]

### coverage
cf_cov1 <- sim_results_vc_coverage_vimkun[["Ceará"]][["VE0"]][["cov50"]][[1]]$
  sim_out$coverage_frac  

total_eff_cov1 <- colSums(cf_cov1 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_cov1[52]

cf_cov2 <- sim_results_vc_coverage_vimkun[["Ceará"]][["VE0"]][["cov50"]][[2]]$
  sim_out$coverage_frac  

total_eff_cov2 <- colSums(cf_cov2 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_cov2[52]

cf_cov3 <- sim_results_vc_coverage_vimkun[["Ceará"]][["VE0"]][["cov50"]][[3]]$
  sim_out$coverage_frac  

total_eff_cov3 <- colSums(cf_cov3 * N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_cov3[52]

cf_cov4 <- sim_results_vc_coverage_vimkun[["Ceará"]][["VE0"]][["cov50"]][[4]]$
  sim_out$coverage_frac  

total_eff_cov4 <- colSums(cf_cov4* N_ceara$Ceará) / sum(N_ceara$Ceará)

total_eff_cov4[52]

coverage_strat <- c(total_eff_cov1[52], total_eff_cov2[52], total_eff_cov3[52], total_eff_cov4[52])
dose_strat <- c(total_eff_dos1[52], total_eff_dos2[52], total_eff_dos3[52], total_eff_dos4[52])


