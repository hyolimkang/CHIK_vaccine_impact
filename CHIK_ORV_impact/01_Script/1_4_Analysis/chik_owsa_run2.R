purrr::iwalk(sensitivity_ranges, function(x, nm) {
  message("Param: ", nm,
          " | length = ", length(x),
          " | names = ", paste(names(x), collapse = ", "))
  print(x)
})

sensitivity_ranges[["foi"]]
str(sensitivity_ranges[["foi"]])
label_lo_hi <- function(x) {
  x <- unname(x)
  if (length(x) != 2) stop("Expected length 2 for OWSA param, got ", length(x))
  ## 정렬해서 lo<hi 확실히 (원래 순서 쓰고 싶으면 주석 처리)
  ord <- order(x)
  x <- x[ord]
  setNames(x, c("lo","hi"))
}
sensitivity_ranges_labeled <- lapply(sensitivity_ranges, label_lo_hi)

test_res <- run_simulation_with_params(
  param_overrides      = setNames(list(sensitivity_ranges_labeled$foi["lo"]), "foi"),
  baseline_params      = baseline_params,
  posterior            = posterior_bh,
  bra_foi_state_summ   = bra_foi_state_summ,
  age_groups           = age_groups,
  N                    = N_bahia$Bahia,
  region               = region,
  observed             = observed_bh,
  age_gr_levels        = age_gr_levels,
  lhs_sample_young     = lhs_sample_young,
  lhs_old              = lhs_old,
  le_sample            = le_sample,
  target_age_list      = target_age_list
)
str(test_res, max.level = 1)


# lo/hi 이름 붙이기
label_lo_hi <- function(x) {
  x <- unname(x)
  stopifnot(length(x) == 2)
  setNames(x, c("lo","hi"))
}
sensitivity_ranges_labeled <- lapply(sensitivity_ranges, label_lo_hi)

sensitivity_results <- purrr::imap(sensitivity_ranges_labeled, function(param_values, param) {
  purrr::imap_dfr(param_values, function(val, scen) {
    res <- run_simulation_with_params(
      param_overrides      = setNames(list(val), param),
      baseline_params      = baseline_params,
      posterior            = posterior_bh,
      bra_foi_state_summ   = bra_foi_state_summ,
      age_groups           = age_groups,
      N                    = N_bahia$Bahia,
      region               = region,
      observed             = observed_bh,
      age_gr_levels        = age_gr_levels,
      lhs_sample_young     = lhs_sample_young,
      lhs_old              = lhs_old,
      le_sample            = le_sample,
      target_age_list      = target_age_list
    )
    tibble::tibble(
      Parameter                 = param,
      Scenario                  = scen,
      Value                     = val,
      pre_tot_cases_agegr1_owsa = res$pre_tot_cases_agegr1_owsa,
      pre_tot_cases_agegr2_owsa = res$pre_tot_cases_agegr2_owsa,
      pre_tot_cases_agegr3_owsa = res$pre_tot_cases_agegr3_owsa,
      pre_tot_cases_agegr4_owsa = res$pre_tot_cases_agegr4_owsa,
      pre_tot_cases_agegr1_mid  = res$pre_tot_cases_agegr1_mid,
      pre_tot_cases_agegr2_mid  = res$pre_tot_cases_agegr2_mid,
      pre_tot_cases_agegr3_mid  = res$pre_tot_cases_agegr3_mid,
      pre_tot_cases_agegr4_mid  = res$pre_tot_cases_agegr4_mid,
      post_vacc_age_by_week     = list(res$post_vacc_age_by_week),
      mid_post_vacc_age_by_week = list(res$mid_post_vacc_age_by_week)
    )
  })
})

tornado_bar_pct_bahia_age <- td_bahia %>%
  # 그래프 파이프라인에 맞추기 위해 기존 이름을 표준화
  mutate(
    Mid_Impact_pct = mid_impact,
    Lo_Impact_pct  = lo_impact,
    Hi_Impact_pct  = hi_impact
  ) %>%
  select(Region, Parameter, AgeGroup,
         Mid_Impact_pct, Lo_Impact_pct, Hi_Impact_pct) %>%
  pivot_longer(
    cols      = c(Lo_Impact_pct, Hi_Impact_pct),
    names_to  = "Type",
    values_to = "Value"
  ) %>%
  mutate(
    Direction = ifelse(Type == "Lo_Impact_pct", "Lower", "Upper"),
    Bar_start = Mid_Impact_pct,  # 기준선
    Bar_end   = Value           # Lo 또는 Hi
  ) %>%
  group_by(Region, Parameter, AgeGroup) %>%
  mutate(TotalImpact = sum(abs(Bar_end - Bar_start), na.rm = TRUE)) %>%
  ungroup() %>%
  # 파라미터별 정렬: 영향 큰 순
  mutate(
    Parameter = fct_reorder(Parameter, TotalImpact, .desc = TRUE)
  )

ggplot(tornado_bar_pct_bahia_age, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction), linewidth = 5, alpha = 0.7) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5) +
  facet_grid(Region ~ AgeGroup, scales = "free_y", switch = "y",
             labeller = labeller(
               Scenario = as_labeller(c(
                 "Scenario_1" = "<20 years", 
                 "Scenario_2" = "20–59 years", 
                 "Scenario_3" = ">60 years"
               ))
             )) +
  scale_color_manual(values = c("Lower" = "#E41A1C", "Upper" = "#377EB8")) +
  labs(
    x = "% Cases averted",
    y = "Parameter",
    color = "Bounds"  
  ) +
  theme_pubclean()+
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    strip.text.y = element_text(size = 6)
  )



###

out_dir <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_2_Processed"
dir.create(out_dir, showWarnings = FALSE)

tornado_by_region <- setNames(vector("list", length(regions )), regions )

for (reg in regions ) {
  message("Running ", reg, " ...")
  res <- tryCatch(
    run_owsa_for_region(
      region_name        = reg,
      bra_foi_state_summ = bra_foi_state_summ,
      age_groups         = age_groups,
      age_gr_levels      = age_gr_levels,
      lhs_sample_young   = lhs_sample_young,
      lhs_old            = lhs_old,
      le_sample          = le_sample,
      target_age_list    = target_age_list,
      hosp               = hosp,
      fatal              = fatal,
      nh_fatal           = nh_fatal
    ),
    error = function(e) {
      message("  ERROR in ", reg, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (!is.null(res)) {
    tornado_by_region[[reg]] <- res
    saveRDS(res$tornado, file = file.path(out_dir, paste0("tornado_", reg, ".rds")))
  }
}

td_bahia <- run_owsa_for_region(
  region_name       = "Bahia",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_ceara <- run_owsa_for_region(
  region_name       = "Ceará",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_pernam <- run_owsa_for_region(
  region_name       = "Pernambuco",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_mg <- run_owsa_for_region(
  region_name       = "Minas Gerais",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_pa <- run_owsa_for_region(
  region_name       = "Paraíba",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_rg <- run_owsa_for_region(
  region_name       = "Rio Grande do Norte",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_pi <- run_owsa_for_region(
  region_name       = "Piauí",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_ag <- run_owsa_for_region(
  region_name       = "Alagoas",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_tc <- run_owsa_for_region(
  region_name       = "Tocantins",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_se <- run_owsa_for_region(
  region_name       = "Sergipe",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

td_go <- run_owsa_for_region(
  region_name       = "Goiás",
  bra_foi_state_summ= bra_foi_state_summ,
  age_groups        = age_groups,
  age_gr_levels     = age_gr_levels,
  lhs_sample_young  = lhs_sample_young,
  lhs_old           = lhs_old,
  le_sample         = le_sample,
  target_age_list   = target_age_list,
  hosp              = hosp,
  fatal             = fatal,
  nh_fatal          = nh_fatal
)

save(td_bahia, td_ceara, td_mg, td_pa, td_pi, td_pernam, td_rg, 
     td_ag, td_tc, td_se, td_go, 
     file =  "00_Data/0_2_Processed/tornado_results.RData")


td_all <- bind_rows(td_bahia, td_ceara, td_mg, td_pa, td_pi, td_pernam, td_rg, 
                    td_ag, td_tc, td_se, td_go,)