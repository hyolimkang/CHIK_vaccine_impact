load("00_Data/0_2_Processed/postsim_vc_ixchiq_model.RData")
load("00_Data/0_2_Processed/postsim_vc_vimkun.RData")
load("00_Data/0_2_Processed/sim_results_vc_coverage_vimkun.RData")
load("00_Data/0_2_Processed/combined_nnv_two_vacc.RData")
load("00_Data/0_2_Processed/sim_results_vc_ixchiq_model.RData")
combined_prepost_case_all_ixchiq <- imap_dfr(                  # ─ region loop ─
  postsim_vc_ixchiq_model,
  #postsim_dose_ixchiq_model,
  function(region_list, region_nm) {
    
    imap_dfr(region_list,                               # ─ VE loop ─
             function(ve_list, ve_tag) {
               
               # ─ coverage loop ─
               map_dfr(ve_list, "summary_week_df", .id = "coverage") %>% 
                 mutate(
                   region   = region_nm,   
                   VE       = ve_tag       
                 )
             }
    )
  }
)

combined_prepost_case_all_vimkun <- imap_dfr(                  # ─ region loop ─
  postsim_vc_vimkun,
  #postsim_dose_vimkun,
  function(region_list, region_nm) {
    
    imap_dfr(region_list,                               # ─ VE loop ─
             function(ve_list, ve_tag) {
               
               # ─ coverage loop ─
               map_dfr(ve_list, "summary_week_df", .id = "coverage") %>% 
                 mutate(
                   region   = region_nm,   
                   VE       = ve_tag       
                 )
             }
    )
  }
)

combined_prepost_case_all_ixchiq$vaccine <- "Ixchiq"
combined_prepost_case_all_vimkun$vaccine <- "Vimkunya"

combined_prepost_case_all_twovacc <- bind_rows(combined_prepost_case_all_ixchiq,
                                               combined_prepost_case_all_vimkun)

# draw base impact
preui_all <- list(
  preui_ce,
  preui_ag,
  preui_bh,
  preui_go,
  preui_mg,
  preui_pa,
  preui_pi,
  preui_pn,
  preui_rg,
  preui_se,
  preui_tc
)

names(preui_all) <- c("Ceará","Alagoas", "Bahia", "Goiás",
                      "Minas Gerais", "Paraíba", "Piauí",
                      "Pernambuco", "Rio Grande do Norte", 
                      "Sergipe", "Tocantins")

calc_total_impact <- function(pre_list, post_arr) {
  # pre_list: list of length n_draws, each (age × week) matrix
  # post_arr: (age × week × draw) array
  
  # 1. pre 리스트 → 3D 배열로 변환
  pre_arr <- simplify2array(pre_list)
  
  # 차원 체크
  stopifnot(length(dim(pre_arr)) == 3, length(dim(post_arr)) == 3)
  
  # 2. draw별 총합
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  # 3. impact 계산
  impact <- (total_pre - total_post) / total_pre * 100
  
  tibble(
    impact_lo  = quantile(impact, 0.025, na.rm = TRUE),
    impact_mid = median(impact, na.rm = TRUE),
    impact_hi  = quantile(impact, 0.975, na.rm = TRUE)
  )
}
results <- purrr::imap_dfr(postsim_vc_vimkun, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        out <- calc_total_impact(pre_list, post_arr)
        
        mutate(out,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})


# national level impact
calc_total_impact_draws <- function(pre_list, post_arr) {
  pre_arr <- simplify2array(pre_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}
all_draws_vim <- purrr::imap_dfr(postsim_vc_vimkun, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        draw_totals <- calc_total_impact_draws(pre_list, post_arr)
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})
all_draws_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        draw_totals <- calc_total_impact_draws(pre_list, post_arr)
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})
all_draws_vim$vaccine <- "Vimkunya"
all_draws_ix$vaccine  <- "Ixchiq"
all_draws <- bind_rows(all_draws_ix, all_draws_vim)

## scaling rho
set.seed(1)

all_draws_ix_true <- all_draws_ix %>%
  group_by(Region) %>%
  mutate(
    rho = sample(rho_pool$rho[rho_pool$Region == first(Region)],
                 size = n(), replace = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    total_pre_true  = total_pre  / rho,
    total_post_true = total_post / rho
  )

global_impact_nat <- all_draws %>%
  group_by(vaccine, VE, Coverage, Scenario, draw_id) %>%
  summarise(
    total_pre_all  = sum(total_pre, na.rm = TRUE),
    total_post_all = sum(total_post, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(impact = (total_pre_all - total_post_all) / total_pre_all * 100) %>%
  group_by(vaccine, VE, Coverage, Scenario) %>%
  summarise(
    impact_lo  = quantile(impact, 0.025, na.rm = TRUE),
    impact_mid = median(impact, na.rm = TRUE),
    impact_hi  = quantile(impact, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

global_impact_subnat <- all_draws %>%
  group_by(vaccine, VE, Coverage, Scenario, draw_id, Region) %>%
  summarise(
    total_pre_all  = sum(total_pre, na.rm = TRUE),
    total_post_all = sum(total_post, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(impact = (total_pre_all - total_post_all) / total_pre_all * 100) %>%
  group_by(vaccine, VE, Coverage, Scenario, Region) %>%
  summarise(
    impact_lo  = quantile(impact, 0.025, na.rm = TRUE),
    impact_mid = median(impact, na.rm = TRUE),
    impact_hi  = quantile(impact, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(Coverage == "cov50")

### hospitalisation 95%UI



## fatal 95% UI 
make_fatal_hosp_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised       <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised   <- symp_matrix - hospitalised
    fatal              <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    list(
      hosp  = hospitalised,
      fatal = fatal
    )
  })
}
calc_total_fatal_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_fatal_list <- lapply(pre_conv, `[[`, "fatal")
  pre_arr <- simplify2array(pre_fatal_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_fatal_list <- lapply(post_conv, `[[`, "fatal")
  post_arr <- simplify2array(post_fatal_list)
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}
make_fatal_hosp_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised     <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised <- symp_matrix - hospitalised
    
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    list(
      hosp  = hospitalised,
      fatal = fatal
    )
  })
}

calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_hosp_list <- lapply(pre_conv, `[[`, "hosp")
  pre_arr <- simplify2array(pre_hosp_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_hosp_list <- lapply(post_conv, `[[`, "hosp")
  post_arr2 <- simplify2array(post_hosp_list)
  if (length(dim(post_arr2)) == 2) {
    post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  }
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr   <- pre_arr[ , , 1:n_draws, drop = FALSE]
  post_arr2 <- post_arr2[ , , 1:n_draws, drop = FALSE]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr2, 3, sum, na.rm = TRUE)
  
  tibble(
    draw_id    = seq_len(n_draws),
    total_pre  = total_pre,
    total_post = total_post
  )
}


all_draws_fatal_vim <- purrr::imap_dfr(postsim_vc_vimkun, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        draw_totals <- calc_total_fatal_draws(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal
        )
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})

### true
scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  # pre_list: length n_draws, each is (age x week) matrix
  # rho_vec: length n_draws
  Map(function(mat, r) mat / r, pre_list, rho_vec)
}

scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  # post_arr: (age x week x draw)
  # rho_vec: length draw
  for (i in seq_along(rho_vec)) post_arr[,,i] <- post_arr[,,i] / rho_vec[i]
  post_arr
}

calc_total_fatal_draws_rho <- function(pre_list, post_arr,
                                       hosp_rate, fatal_rate, nh_fatal_rate,
                                       rho_vec) {
  # 1) draw 길이 맞추기
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) 기존 함수 그대로 사용
  calc_total_fatal_draws(pre_list_true, post_arr_true,
                         hosp_rate, fatal_rate, nh_fatal_rate)
}


set.seed(1)

all_draws_fatal_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool (필수)
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        # 해당 조합에서 실제 사용 가능한 draw 수
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        # rho를 독립적으로 리샘플링 (500 rho pool -> 1000 draws도 OK)
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_fatal_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})


calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- lapply(pre_list, function(symp_matrix) {
    sweep(symp_matrix, 1, hosp_rate, `*`)   # hospitalised
  })
  pre_arr <- simplify2array(pre_conv)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- lapply(post_list, function(symp_matrix) {
    sweep(symp_matrix, 1, hosp_rate, `*`)
  })
  post_arr2 <- simplify2array(post_conv)
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr2 <- post_arr2[ , , 1:n_draws]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,   3, sum, na.rm = TRUE)
  total_post <- apply(post_arr2, 3, sum, na.rm = TRUE)
  
  tibble(
    draw_id = 1:n_draws,
    total_pre = total_pre,
    total_post = total_post
  )
}
calc_total_hosp_draws_rho <- function(pre_list, post_arr, hosp_rate, rho_vec) {
  
  # 1) draw 길이 맞추기
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling (true symptomatic)
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) hosp totals
  calc_total_hosp_draws(pre_list_true, post_arr_true, hosp_rate = hosp_rate)
}
set.seed(1)

all_draws_hosp_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        # 가능한 draw 수
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        # rho는 독립적으로 리샘플링
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_hosp_draws_rho(
          pre_list = pre_list,
          post_arr = post_arr,
          hosp_rate = hosp,
          rho_vec = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})



#####
all_draws_fatal_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        draw_totals <- calc_total_fatal_draws(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal
        )
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})

all_draws_fatal_vim$vaccine <- "Vimkunya"
all_draws_fatal_ix$vaccine  <- "Ixchiq"
all_draws_fatal <- bind_rows(all_draws_fatal_vim, all_draws_fatal_ix)

all_draws_hosp_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        draw_totals <- calc_total_hosp_draws(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal
        )
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})


global_impact_nat_fatal <- all_draws_fatal %>%
  group_by(vaccine, VE, Coverage, Scenario, draw_id) %>%
  summarise(
    total_pre_all  = sum(total_pre, na.rm = TRUE),
    total_post_all = sum(total_post, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(impact = (total_pre_all - total_post_all) / total_pre_all * 100) %>%
  group_by(vaccine, VE, Coverage, Scenario) %>%
  summarise(
    impact_lo  = quantile(impact, 0.025, na.rm = TRUE),
    impact_mid = median(impact, na.rm = TRUE),
    impact_hi  = quantile(impact, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

fatal_nat_draws <- all_draws_fatal_ix %>%
  group_by(draw_id, VE, Coverage, Scenario) %>%
  summarise(
    pre_fatal  = sum(total_pre,  na.rm = TRUE),
    post_fatal = sum(total_post, na.rm = TRUE),
    diff_fatal = pre_fatal - post_fatal,
    .groups = "drop"
  )

fatal_nat_draws_cov50 <- fatal_nat_draws %>%
  filter(Coverage == "cov50")

tot_vacc_nat <- combined_nnv_df_region_coverage_model %>%
  group_by(VE, VC, scenario) %>%
  summarise(
    tot_vacc_nat = sum(tot_vacc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(VC == "cov50")%>%
  rename(Scenario = scenario, Coverage = VC)


fatal_nat_draws_cov50 <- fatal_nat_draws_cov50 %>%
  mutate(Scenario = as.character(Scenario))

fatal_nat_draws_cov50 <- fatal_nat_draws_cov50 %>%
  mutate(
    Scenario = paste0("Scenario_", as.integer(Scenario))
  )

tot_vacc_nat <- tot_vacc_nat %>%
  mutate(Scenario = as.character(Scenario))

fatal_nat_draws_cov50 <- fatal_nat_draws_cov50 %>%
  left_join(
    tot_vacc_nat,
    by = c("VE", "Coverage", "Scenario")
  )


risk_draws <- lhs_sample %>%
  mutate(draw_id = row_number()) %>%
  transmute(
    draw_id,
    p_death_vacc_u65,
    p_death_vacc_65
  )

fatal_nat_draws_cov50 <- fatal_nat_draws_cov50 %>%
  left_join(risk_draws, by = "draw_id")

fatal_nat_draws_cov50 <- fatal_nat_draws_cov50 %>%
  mutate(
    p_death_vacc_all = 0.5 * p_death_vacc_u65 + 0.5 * p_death_vacc_65
  )
death_draw_df <- fatal_nat_draws_cov50 %>%
  filter(is.finite(tot_vacc_nat), tot_vacc_nat > 0) %>%
  mutate(
    y_10k   = (diff_fatal / tot_vacc_nat) * 1e4,   # averted deaths per 10k vaccinated
    x_10k   = 1e4 * p_death_vacc_all,              # vaccine AE deaths per 10k vaccinated
    net_10k = y_10k - x_10k
  )


## daly 95% UI 
make_daly_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate,
                            dw_hosp, dw_nonhosp, dw_chronic,
                            dur_acute, dur_subacute, dur_chronic,
                            dw_subacute, subac_prop, chr_prop,
                            le_left_vec) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised       <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised   <- symp_matrix - hospitalised
    fatal              <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    yld_acute    <- (hospitalised * dw_hosp     * dur_acute) +
      (non_hospitalised * dw_nonhosp * dur_acute)
    yld_subacute <- (hospitalised * subac_prop * dw_subacute * dur_subacute) +
      (non_hospitalised * chr_prop * dw_subacute * dur_subacute)
    yld_chronic  <- (hospitalised * chr_prop * dw_chronic * dur_chronic) +
      (non_hospitalised * chr_prop * dw_chronic * dur_chronic)
    
    yld_total <- yld_acute + yld_subacute + yld_chronic
    
    # precomputed 기대수명 벡터 사용
    yll <- sweep(fatal, 1, le_left_vec, `*`)
    
    daly_tot <- yld_total + yll
    
    list(
      hosp     = hospitalised,
      fatal    = fatal,
      yld_tot  = yld_total,
      yll      = yll,
      daly_tot = daly_tot
    )
  })
}
le_by_age <- function(age_numeric) {
  if (age_numeric <= 1) return(quantile(le_sample$le_1, 0.5))
  if (age_numeric < 20) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 30) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 40) return(quantile(le_sample$le_3, 0.5))
  if (age_numeric < 50) return(quantile(le_sample$le_4, 0.5))
  if (age_numeric < 60) return(quantile(le_sample$le_5, 0.5))
  if (age_numeric < 70) return(quantile(le_sample$le_6, 0.5))
  if (age_numeric < 80) return(quantile(le_sample$le_7, 0.5))
  return(quantile(le_sample$le_8, 0.5))
}

calc_total_daly_draws <- function(pre_list, post_arr,
                                  hosp_rate, fatal_rate, nh_fatal_rate,
                                  dw_hosp, dw_nonhosp, dw_chronic,
                                  dur_acute, dur_subacute, dur_chronic,
                                  dw_subacute, subac_prop, chr_prop,
                                  le_left_vec) {
  # Pre
  pre_conv <- make_daly_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate,
                              dw_hosp, dw_nonhosp, dw_chronic,
                              dur_acute, dur_subacute, dur_chronic,
                              dw_subacute, subac_prop, chr_prop,
                              le_left_vec)
  pre_daly_list <- lapply(pre_conv, `[[`, "daly_tot")
  pre_arr <- simplify2array(pre_daly_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # Post
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_daly_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate,
                               dw_hosp, dw_nonhosp, dw_chronic,
                               dur_acute, dur_subacute, dur_chronic,
                               dw_subacute, subac_prop, chr_prop,
                               le_left_vec)
  post_daly_list <- lapply(post_conv, `[[`, "daly_tot")
  post_arr <- simplify2array(post_daly_list)
  
  # Align
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

le_left_vec <- sapply(age_groups, le_by_age)
daly_draw_totals <- calc_total_daly_draws(
  pre_list = preui_all$Ceará$sim_results_list_rawsymp,
  post_arr = postsim_vc_vimkun$Ceará$VE97.8$cov50$scenario_result[[1]]$sim_result$age_array_raw_symp,
  hosp_rate     = hosp,
  fatal_rate    = fatal,
  nh_fatal_rate = nh_fatal,
  dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
  dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
  dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),
  dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
  dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
  dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),
  dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),
  subac_prop    = quantile(lhs_sample_young$subac, 0.5),
  chr_prop      = (quantile(lhs_sample_young$chr6m, 0.5) +
                     quantile(lhs_sample_young$chr12m, 0.5) +
                     quantile(lhs_sample_young$chr30m, 0.5)),
  le_left_vec   = le_left_vec   
)

all_draws_daly_vim <- purrr::imap_dfr(postsim_vc_vimkun, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        draw_totals <- calc_total_daly_draws(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
          dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
          dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),
          dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
          dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
          dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),
          dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),
          subac_prop    = quantile(lhs_sample_young$subac, 0.5),
          chr_prop      = (quantile(lhs_sample_young$chr6m, 0.5) +
                             quantile(lhs_sample_young$chr12m, 0.5) +
                             quantile(lhs_sample_young$chr30m, 0.5)),
          le_left_vec   = le_left_vec    
        )
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})
all_draws_daly_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        draw_totals <- calc_total_daly_draws(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
          dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
          dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),
          dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
          dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
          dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),
          dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),
          subac_prop    = quantile(lhs_sample_young$subac, 0.5),
          chr_prop      = (quantile(lhs_sample_young$chr6m, 0.5) +
                             quantile(lhs_sample_young$chr12m, 0.5) +
                             quantile(lhs_sample_young$chr30m, 0.5)),
          le_left_vec   = le_left_vec    
        )
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})

all_draws_daly_vim$vaccine <- "Vimkunya"
all_draws_daly_ix$vaccine  <- "Ixchiq"
all_draws_daly <- bind_rows(all_draws_daly_vim, all_draws_daly_ix)


## true
scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  # pre_list: list of (age x week) matrices
  Map(function(mat, r) mat / r, pre_list, rho_vec)
}

scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  # post_arr: (age x week x draw)
  for (i in seq_along(rho_vec)) post_arr[,,i] <- post_arr[,,i] / rho_vec[i]
  post_arr
}

calc_total_daly_draws_rho <- function(pre_list, post_arr,
                                      hosp_rate, fatal_rate, nh_fatal_rate,
                                      dw_hosp, dw_nonhosp, dw_chronic,
                                      dur_acute, dur_subacute, dur_chronic,
                                      dw_subacute, subac_prop, chr_prop,
                                      le_left_vec,
                                      rho_vec) {
  
  # 1) draw 길이 맞추기
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  # 3) 기존 DALY 함수 그대로
  calc_total_daly_draws(
    pre_list = pre_list_true,
    post_arr = post_arr_true,
    hosp_rate     = hosp_rate,
    fatal_rate    = fatal_rate,
    nh_fatal_rate = nh_fatal_rate,
    dw_hosp       = dw_hosp,
    dw_nonhosp    = dw_nonhosp,
    dw_chronic    = dw_chronic,
    dur_acute     = dur_acute,
    dur_subacute  = dur_subacute,
    dur_chronic   = dur_chronic,
    dw_subacute   = dw_subacute,
    subac_prop    = subac_prop,
    chr_prop      = chr_prop,
    le_left_vec   = le_left_vec
  )
}


set.seed(1)

all_draws_daly_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        # 해당 조합에서 가능한 draw 수
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        # rho는 독립적으로 리샘플링 (500 rho pool -> 1000 draw도 가능)
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_daly_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          
          dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
          dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
          dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),
          
          dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
          dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
          dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),
          
          dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),
          subac_prop    = quantile(lhs_sample_young$subac, 0.5),
          chr_prop      = (quantile(lhs_sample_young$chr6m, 0.5) +
                             quantile(lhs_sample_young$chr12m, 0.5) +
                             quantile(lhs_sample_young$chr30m, 0.5)),
          
          le_left_vec   = le_left_vec,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})

######
global_impact_nat_daly <- all_draws_daly %>%
  group_by(vaccine, VE, Coverage, Scenario, draw_id) %>%
  summarise(
    total_pre_all  = sum(total_pre, na.rm = TRUE),
    total_post_all = sum(total_post, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(impact = (total_pre_all - total_post_all) / total_pre_all * 100) %>%
  group_by(vaccine, VE, Coverage, Scenario) %>%
  summarise(
    impact_lo  = quantile(impact, 0.025, na.rm = TRUE),
    impact_mid = median(impact, na.rm = TRUE),
    impact_hi  = quantile(impact, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## only 50% 
combined_prepost_case_all_50 <- combined_prepost_case_all_twovacc %>% filter(coverage == "cov50")

global_impact_ve <- combined_prepost_case_all_50 %>%       
  group_by(VE, Scenario, region, vaccine) %>%               
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_post_cases_lo = sum(post_cases_lo, na.rm = TRUE),
    total_post_cases_hi = sum(post_cases_hi, na.rm = TRUE),
    
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    total_pre_cases_lo  = sum(pre_cases_lo,  na.rm = TRUE),
    total_pre_cases_hi  = sum(pre_cases_hi,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    diff_lo   = total_pre_cases_lo - total_post_cases_lo,
    diff_hi   = total_pre_cases_hi - total_post_cases_hi,
    
    impact    = diff / total_pre_cases * 100,
    impact_lo = diff_lo / total_pre_cases_lo * 100,
    impact_hi = diff_hi / total_pre_cases_hi * 100
    
  )

annotation_ve <- global_impact_ve %>%
  mutate(
    Strategy = dplyr:: recode(Scenario,
                              "Scenario_1" = "1-11 y",
                              "Scenario_2" = "12-17 y",
                              "Scenario_3" = "18-59 y",
                              "Scenario_4" = "60+ y",)
  ) %>%
  group_by(region, VE, vaccine) %>%                          # ★ VE별로 따로 주석
  summarise(
    ann = paste0(Strategy, ": ", round(impact, 1), "%", collapse = "\n"),
    .groups = "drop"
  )

# 2) 관측치도 VE별로 복제
observed_all_subnat <- tidyr::crossing(
  observed_all,                            # Week, Observed, region
  VE = unique(combined_prepost_case_all_50$VE)
)%>%
  mutate(
    setting = setting_key[region]
  )

observed_by_setting <- observed_all_subnat %>%
  group_by(setting, Week, VE) %>%       # grouping variables
  summarise(
    Observed = sum(Observed),        
    .groups = "drop"                     # ungroup
  )

observed_nat <- observed_all_subnat %>%
  group_by(Week, VE) %>%       # grouping variables
  summarise(
    Observed = sum(Observed),        
    .groups = "drop"                     # ungroup
  ) %>% mutate(
    setting = "National"
  )

combined_by_setting <- combined_prepost_case_all_50 %>% 
  mutate(
    setting = dplyr::recode(region, !!!setting_key)   # add new column
  )

setting_summ <- combined_by_setting %>% 
  group_by(setting, VE, Scenario, Week, vaccine) %>% 
  summarise(
    post_cases        = sum(post_cases),
    pre_cases         = sum(pre_cases),
    diff              = sum(diff),
    post_fatal        = sum(post_fatal),
    pre_fatal         = sum(pre_fatal),
    diff_fatal        = sum(diff_fatal),
    post_daly         = sum(post_daly),
    pre_daly          = sum(pre_daly),
    diff_daly         = sum(diff_daly),
    post_weekly_median = sum(post_weekly_median),
    post_weekly_low95  = sum(post_weekly_low95),
    post_weekly_hi95   = sum(post_weekly_hi95),
    lo95               = sum(lo95),
    hi95               = sum(hi95),
    .groups = "drop"
  )

combined_nat_setting <- combined_by_setting %>% 
  group_by(VE, Scenario, Week, vaccine) %>% 
  summarise(
    post_cases        = sum(post_cases),
    pre_cases         = sum(pre_cases),
    diff              = sum(diff),
    post_fatal        = sum(post_fatal),
    pre_fatal         = sum(pre_fatal),
    diff_fatal        = sum(diff_fatal),
    post_daly         = sum(post_daly),
    pre_daly          = sum(pre_daly),
    diff_daly         = sum(diff_daly),
    post_weekly_median = sum(post_weekly_median),
    post_weekly_low95  = sum(post_weekly_low95),
    post_weekly_hi95   = sum(post_weekly_hi95),
    lo95               = sum(lo95),
    hi95               = sum(hi95),
    .groups = "drop"
  )

combined_nat_setting$setting <- "National"

combined_all <- bind_rows(combined_nat_setting, setting_summ) %>%
  mutate(
    setting = factor(setting,
                     levels = c("Low", "Moderate", "High", "National"))  # bottom → top
  )


observed_all_full <- bind_rows(observed_nat, observed_by_setting)
observed_all_full <- observed_all_full %>%
  tidyr::crossing(
    vaccine = unique(combined_all$vaccine)
  )
observed_all_full <- observed_all_full %>%
  filter(
    !(vaccine == "Ixchiq"   & VE == "VE97.8"),
    !(vaccine == "Vimkunya" & VE == "VE98.9")
  )


vacc_start_week_s1 = 2
vacc_end_week_s1   = 11 

ve_subnat<-
  ggplot(combined_prepost_case_all_50) +
  # ❶ 백신 투입 기간 회색 음영 (facet 전체 공통)
  geom_ribbon(data = combined_prepost_case_all_50 %>%
                filter(between(Week, vacc_start_week_s1, vacc_end_week_s1)),
              aes(x = Week, ymin = 0, ymax = Inf),
              fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  
  # ❷ 포스트-백신 UI 리본 (시나리오별 색)
  geom_ribbon(aes(x = Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = 0.25) +
  
  # ❸ 프리-백신 UI 리본 (연한 회색)
  geom_ribbon(aes(x = Week, ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = 0.5) +
  
  # ❹ 포스트-백신 선 (Scenario별 색)
  geom_line(aes(x = Week, y = post_cases, color = Scenario), size = 0.25) +
  
  # ❺ 프리-백신 선 (검정 점선)
  geom_line(aes(x = Week, y = pre_cases, group = Scenario),
            color = "black", linetype = "dashed", size = 0.25) +
  
  # ❻ 관측치 점
  geom_point(data = observed_all_subnat,
             aes(x = Week, y = Observed),
             size = 0.2, inherit.aes = FALSE) +
  
  # ❼ 지역-별 주석
  geom_text(data = annotation_ve,
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE) +
  
  scale_fill_brewer(name   = "Vaccination strategy",  
                    palette = "Set1",
                    labels = c("Scenario_1" = "1–19 y",
                               "Scenario_2" = "20–59 y",
                               "Scenario_3" = "≥60 y")) +
  scale_color_brewer(name   = "Vaccination strategy",  
                     palette = "Set1",
                     labels = c("Scenario_1" = "1–19 y",
                                "Scenario_2" = "20–59 y",
                                "Scenario_3" = "≥60 y")) +
  scale_y_continuous(labels = comma) +
  labs(x = "Week", y = "Predicted symptomatic cases",
       fill = "Scenario", color = "Scenario") +
  
  # ❽ facet: 행 = 지역, 열 = VE
  facet_grid(region ~ VE, scales = "free_y") +
  #facet_wrap(vars(region, VE), ncol = 6, scales = "free_y")+
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.margin  = margin(5, 20, 5, 5),
    legend.position = "bottom"
  )

make_prepost_plot <- function(df_plot, ann_df, vacc_title = NULL) {
  
  ggplot(df_plot) +
    # ① 백신 투입 기간 음영
    geom_ribbon(data = df_plot %>% filter(between(Week, vacc_start_week_s1, vacc_end_week_s1)),
                aes(Week, ymin = 0, ymax = Inf),
                fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
    
    # ② 포스트-백신 UI 리본
    geom_ribbon(aes(Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                    fill = Scenario), alpha = 0.25) +
    
    # ③ 프리-백신 UI 리본
    geom_ribbon(aes(Week, ymin = lo95, ymax = hi95),
                fill = "lightgrey", alpha = 0.5) +
    
    # ④ 포스트-백신 선
    geom_line(aes(Week, post_cases, colour = Scenario), size = 0.25) +
    
    # ⑤ 프리-백신 선
    geom_line(aes(Week, pre_cases, group = Scenario),
              colour = "black", linetype = "dashed", size = 0.25) +
    
    # ⑥ 관측치 점
    #geom_point(data = observed_all_subnat,
    #           aes(Week, Observed), size = 0.2, inherit.aes = FALSE) +
    
    # ⑦ 지역별 주석
    geom_text(data = ann_df,
              aes(x = Inf, y = Inf, label = ann),
              hjust = 1.05, vjust = 1.1, size = 2.4,
              inherit.aes = FALSE) +
    
    scale_fill_brewer(name = "Vaccination strategy",
                      palette = "Set1",
                      labels = c("Scenario_1" = "1–11 y",
                                 "Scenario_2" = "12-17 y",
                                 "Scenario_3" = "18-59 y",
                                 "Scenario_4" = "60+ y")) +
    scale_colour_brewer(name = "Vaccination strategy",
                        palette = "Set1",
                        labels = c("Scenario_1" = "1–11 y",
                                   "Scenario_2" = "12-17 y",
                                   "Scenario_3" = "18-59 y",
                                   "Scenario_4" = "60+ y")) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = vacc_title,
         x = "Week",
         y = "Predicted symptomatic cases",
         fill = "Scenario",
         colour = "Scenario") +
    
    ## ── (행) region × (열) VE ──
    facet_wrap(~ region, ncol = 4, scales = "free_y") +
    theme_pubclean(base_size = 9) +
    theme(
      strip.text.y    = element_text(angle = 0),
      plot.margin     = margin(5, 20, 5, 5),
      legend.position = "bottom"
    )
}

make_prepost_plot_select <- function(df_plot, ann_df,
                                     vacc_title = NULL,
                                     ncol_pair  = 4) {  # 한 행에 들어갈 패널 수(짝수 권장)
  
  library(dplyr)
  library(ggplot2)
  library(tidyr)      # uncount()
  
  ## ────────────────────────────────────────────────────────────────
  ## 1) setting 우선순위 고정 + 기본 strip("지역 (Setting)") 생성
  ## ────────────────────────────────────────────────────────────────
  df_plot <- df_plot %>%
    mutate(
      setting   = factor(setting, levels = c("High", "Moderate", "Low")),
      facet_lab = sprintf("%s (%s)", region, setting)
    )
  
  ann_df  <- ann_df %>%
    mutate(
      setting   = factor(setting, levels = c("High", "Moderate", "Low")),
      facet_lab = sprintf("%s (%s)", region, setting)
    )
  
  ## ────────────────────────────────────────────────────────────────
  ## 2) VE → 설명 텍스트 매핑
  ## ────────────────────────────────────────────────────────────────
  ve_map <- function(x) ifelse(x == "VE0",
                               "Disease blocking only",
                               "Disease & infection blocking")
  
  df_plot <- df_plot %>% mutate(VE_desc = ve_map(VE))
  ann_df  <- ann_df  %>% mutate(VE_desc = ve_map(VE))
  
  ## ────────────────────────────────────────────────────────────────
  ## 3) panel 수준(levels) 생성 :  High → Moderate → Low  정렬 후
  ##    각 지역마다 VE 두 칸(블로킹 / 블로킹+감염) 붙이기
  ## ────────────────────────────────────────────────────────────────
  panel_levels <- df_plot %>%
    distinct(region, setting) %>%           # 한 지역당 한 줄
    arrange(setting, region) %>%            # High → Mod → Low
    mutate(base = sprintf("%s (%s)", region, setting)) %>%
    uncount(2, .id = "ve_ix") %>%           # VE 두 칸 확보
    mutate(ve_lab = ifelse(ve_ix == 1,
                           "Disease blocking only",
                           "Disease & infection blocking"),
           level  = sprintf("%s\n%s", base, ve_lab)) %>%
    pull(level)
  
  ## panel_lab 컬럼 최종 생성 (factor 수준 강제)
  df_plot <- df_plot %>%
    mutate(panel_lab = sprintf("%s\n%s", facet_lab, VE_desc),
           panel_lab = factor(panel_lab, levels = panel_levels))
  
  ann_df  <- ann_df %>%
    mutate(panel_lab = sprintf("%s\n%s", facet_lab, VE_desc),
           panel_lab = factor(panel_lab, levels = panel_levels))
  
  ## ────────────────────────────────────────────────────────────────
  ## 4) 공통 ggplot 레이어
  ## ────────────────────────────────────────────────────────────────
  p_base <- ggplot(df_plot) +
    ## (1) 백신 투입 기간 음영
    geom_ribbon(
      data = df_plot %>%
        filter(between(Week, vacc_start_week_s1, vacc_end_week_s1)),
      aes(Week, ymin = 0, ymax = Inf),
      fill = "grey70", alpha = 0.2, inherit.aes = FALSE
    ) +
    ## (2) 포스트-백신 UI 리본
    geom_ribbon(
      aes(Week, ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = Scenario),
      alpha = 0.25
    ) +
    ## (3) 프리-백신 UI 리본
    geom_ribbon(
      aes(Week, ymin = lo95, ymax = hi95),
      fill = "lightgrey", alpha = 0.5
    ) +
    ## (4) 포스트-백신 선
    geom_line(aes(Week, post_cases, colour = Scenario), size = 0.25) +
    ## (5) 프리-백신 선
    geom_line(aes(Week, pre_cases, group = Scenario),
              colour = "black", linetype = "dashed", size = 0.25) +
    ## (6) 패널별 주석
    geom_text(
      data = ann_df,
      aes(x = Inf, y = Inf, label = ann),
      hjust = 1.05, vjust = 1.1, size = 2.4, inherit.aes = FALSE
    ) +
    ## 스케일·레이블
    scale_fill_brewer(
      palette = "Set1",
      labels = c(`Scenario_1` = "1–11 y", `Scenario_2` = "12–17 y",
                 `Scenario_3` = "18–59 y", `Scenario_4` = "60+ y")
    ) +
    scale_colour_brewer(
      palette = "Set1",
      labels = c(`Scenario_1` = "1–11 y", `Scenario_2` = "12–17 y",
                 `Scenario_3` = "18–59 y", `Scenario_4` = "60+ y")
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(title  = vacc_title,
         x      = "Week",
         y      = "Predicted symptomatic cases",
         fill   = "Scenario",
         colour = "Scenario") +
    theme_pubclean(base_size = 9) +
    theme(
      strip.text      = element_text(size = 7),
      plot.margin     = margin(5, 20, 5, 5),
      legend.position = "bottom"
    )
  
  ## ────────────────────────────────────────────────────────────────
  ## 5) facet_wrap :  짝수(ncol_pair) 유지 → 지역×2(VE) 한 행에 배치
  ## ────────────────────────────────────────────────────────────────
  p_final <- p_base +
    facet_wrap(~ panel_lab,
               ncol   = ncol_pair,   # 예: 4 → (VE쌍) 두 지역씩 한 행
               scales = "free_y")
  
  return(p_final)
}

make_prepost_plot_pair <- function(df_plot, ann_df,
                                   vacc_title = NULL,
                                   layout = c("grid", "wrap")) {
  
  layout <- match.arg(layout)          # "grid" (기본) 또는 "wrap"
  
  library(dplyr); library(ggplot2); library(tidyr)
  
  ## ─ 1. 지역·setting 라벨 --------------------------------------------------
  levels_setting <- c("High", "Moderate", "Low")
  
  df_plot <- df_plot %>% 
    mutate(
      setting    = factor(setting, levels = levels_setting),
      region_lab = paste0(region, " (", setting, ")")
    )
  
  region_levels <- df_plot %>% 
    distinct(region_lab, setting, region) %>% 
    arrange(setting, region) %>% 
    pull(region_lab)
  
  df_plot <- df_plot %>% 
    mutate(region_f = factor(region_lab, levels = region_levels))
  
  ann_df <- ann_df %>% 
    mutate(
      setting    = factor(setting, levels = levels_setting),
      region_lab = paste0(region, " (", setting, ")"),
      region_f   = factor(region_lab, levels = region_levels)
    )
  
  ## ─ 2. VE 라벨 -----------------------------------------------------------
  ve_map <- function(x) ifelse(x == "VE0",
                               "Disease blocking only",
                               "Disease & infection blocking")
  df_plot <- df_plot %>% mutate(VE_desc = ve_map(VE))
  ann_df  <- ann_df  %>% mutate(VE_desc = ve_map(VE))
  
  ## ─ 3. 주석 y 위치 -------------------------------------------------------
  ann_df <- ann_df %>% 
    group_by(region_f) %>% 
    arrange(VE_desc) %>% 
    mutate(
      y_pos = Inf - row_number() * 0.07 *
        diff(range(c(df_plot$post_cases,
                     df_plot$hi95,
                     df_plot$post_weekly_hi95),
                   na.rm = TRUE))
    ) %>% 
    ungroup()
  
  ## ─ 4. 기본 ggplot --------------------------------------------------------
  p <- ggplot(df_plot) +
    ## ❶ 백신 투여 기간 배경
    geom_ribbon(data = df_plot %>%
                  filter(between(Week, vacc_start_week_s1, vacc_end_week_s1)),
                aes(Week, ymin = 0, ymax = Inf),
                fill = "grey70", alpha = .20, inherit.aes = FALSE) +
    ## ❷ Post-vaccine 95% CI
    geom_ribbon(aes(Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                    fill = Scenario), alpha = .25) +
    ## ❸ Pre-vaccine 95% CI
    geom_ribbon(aes(Week, ymin = lo95, ymax = hi95),
                fill = "lightgrey", alpha = .50) +
    
    ## ❹ 곡선: with / without vaccination
    geom_line(aes(Week, post_cases,
                  colour   = Scenario,
                  linetype = "With vaccination"),
              size = .25) +
    geom_line(aes(Week, pre_cases,
                  linetype = "Without vaccination",
                  group    = Scenario),
              colour = "black", size = .25) +
    
    ## ❺ 주석
    geom_text(data = ann_df,
              aes(x = Inf, y = y_pos, label = ann),
              hjust = 1.05, vjust = 1, size = 2.4, inherit.aes = FALSE) +
    
    ## ❻ 범례 설정 -------------------------------------------------------
  scale_fill_brewer(
    palette = "Set1",
    labels  = c(Scenario_1 = "1–11 y", Scenario_2 = "12–17 y",
                Scenario_3 = "18–59 y", Scenario_4 = "60+ y"),
    name    = "Vaccination strategy") +
    scale_colour_brewer(
      palette = "Set1",
      labels  = c(Scenario_1 = "1–11 y", Scenario_2 = "12–17 y",
                  Scenario_3 = "18–59 y", Scenario_4 = "60+ y"),
      name    = "Vaccination strategy") +
    scale_linetype_manual(
      name   = NULL,
      values = c("With vaccination"    = "solid",
                 "Without vaccination" = "dashed")) +
    
    scale_y_continuous(labels = scales::comma) +
    labs(title = vacc_title,
         x     = "Week",
         y     = "Predicted symptomatic cases") +
    theme_pubclean(base_size = 9) +
    theme(
      axis.text.x       = element_text(angle = 0, hjust = .5),  # x축 숫자 0°
      legend.position   = "bottom",
      plot.margin       = margin(5, 20, 5, 5)
    )
  
  ## ─ 5. facet 배치 ---------------------------------------------------------
  if (layout == "grid") {
    # (행) VE ─ (열) 지역; strip.x 기본값 = 위쪽
    p <- p +
      facet_grid(
        rows   = vars(VE_desc),     # 2행
        cols   = vars(region_f),    # 3열
        scales = "free_y"
      ) +
      theme(
        strip.text.x      = element_text(angle = 0, hjust = .5), # 지역 strip 0°
        strip.text.y.left = element_text(angle = 0)
      )
    
  } else {  # layout == "wrap": 3행 × 2열 고정
    df_plot$facet_wrap_lab <- interaction(df_plot$region_f, df_plot$VE_desc,
                                          sep = " • ", lex.order = TRUE)
    ann_df$facet_wrap_lab  <- interaction(ann_df$region_f, ann_df$VE_desc,
                                          sep = " • ", lex.order = TRUE)
    
    p <- p %+% df_plot +
      geom_text(data = ann_df,
                aes(x = Inf, y = y_pos, label = ann),
                hjust = 1.05, vjust = 1, size = 2.4, inherit.aes = FALSE) +
      facet_wrap(~ facet_wrap_lab, ncol = 2, scales = "free_y") +
      theme(
        strip.text       = element_text(size = 8, margin = margin(b = 2)),
        strip.background = element_rect(fill = "#F0F0F0", colour = NA)
      )
  }
  
  return(p)
}

## -----------------------------------------
setting_key <- c(
  "Ceará"             = "High",
  "Bahia"             = "Low",
  "Paraíba"           = "High",
  "Pernambuco"        = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"             = "High",
  "Tocantins"         = "Moderate",
  "Alagoas"           = "High",
  "Minas Gerais"      = "Low",
  "Sergipe"           = "Low",
  "Goiás"             = "Low"
)

setting_key_df <- tibble(
  region  = names(setting_key),
  setting = unname(setting_key)
)

ix_df  <- combined_prepost_case_all_50 %>% 
  filter(vaccine == "Ixchiq", VE %in% c("VE0", "VE98.9")) %>%
  mutate(VE = factor(VE, levels = c("VE0", "VE98.9")))

ix_ann <- annotation_ve %>% 
  filter(vaccine == "Ixchiq", VE %in% levels(ix_df$VE))
ix_ann <- ix_ann %>%
  left_join(setting_key_df, by = c("region" = "region"))

ix_df_0    <- ix_df  %>% filter(VE == "VE0")
ix_df_98.9 <- ix_df  %>% filter(VE == "VE98.9")

ix_ann_0    <- ix_ann %>% filter(VE == "VE0")
ix_ann_98.9 <- ix_ann %>% filter(VE == "VE98.9")

p_ix_0    <- 
  make_prepost_plot(ix_df_0,    ix_ann_0, vacc_title = "Ixchiq") + theme(legend.position = "none")
p_ix_98.9 <- 
  make_prepost_plot(ix_df_98.9, ix_ann_98.9, vacc_title = "Ixchiq") + theme(axis.title.y = element_blank())

pw <- 
(p_ix_0 + p_ix_98.9) +
  plot_layout(guides = "collect") 
ix <- 
(pw & theme(
  legend.position = "bottom",
  legend.justification = "center")) +
  plot_annotation(tag_levels = list(c("B","C")))

ggsave(filename = "02_Outputs/2_1_Figures/fig1_bc.jpg", ix, width = 14, height = 6, dpi = 1200)

## -----------------------------------------
vim_df  <- combined_prepost_case_all_50 %>% 
  filter(vaccine == "Vimkunya", VE %in% c("VE0", "VE97.8")) %>%
  mutate(VE = factor(VE, levels = c("VE0", "VE97.8")))

vim_ann <- annotation_ve %>% 
  filter(vaccine == "Vimkunya", VE %in% levels(vim_df$VE))
vim_ann <- vim_ann %>%
  left_join(setting_key_df, by = c("region" = "region"))

vim_df_0    <- vim_df  %>% filter(VE == "VE0")
vim_df_97.8 <- vim_df  %>% filter(VE == "VE97.8")

vim_ann_0    <- vim_ann %>% filter(VE == "VE0")
vim_ann_97.8 <- vim_ann %>% filter(VE == "VE97.8")

p_vim_0    <- make_prepost_plot(vim_df_0,    vim_ann_0, vacc_title = "Vimkunya") + theme(legend.position = "none")
p_vim_97.8 <- make_prepost_plot(vim_df_97.8, vim_ann_97.8, vacc_title = "Vimkunya") + theme(axis.title.y = element_blank())

ve_max <- 
  (p_ix_98.9 + p_vim_97.8) +
  plot_layout(guides = "collect") 

ve_max <- 
(ve_max & theme(
  legend.position = "bottom",
  legend.justification = "center")) 
 # + plot_annotation(tag_levels = list(c("D","E")))

ggsave(filename = "02_Outputs/2_1_Figures/figs8_vemax.jpg", ve_max, width = 14, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs8_vemax.pdf", ve_max, width = 14, height = 6, dpi = 1200)

# ixchiq and vimkunya max VE
ve_min <- 
  (p_ix_0 + p_vim_0) +
  plot_layout(guides = "collect") 

vim <- 
  (pw & theme(
    legend.position = "bottom",
    legend.justification = "center")) +
  plot_annotation(tag_levels = list(c("C","D")))

ve_min <- 
  (p_ix_0 + p_vim_0) +
  plot_layout(guides = "collect") 

ve_min <- 
  (ve_min & theme(
    legend.position = "bottom",
    legend.justification = "center")) 
 #+plot_annotation(tag_levels = list(c("D","E")))

ggsave(filename = "02_Outputs/2_1_Figures/figs8_vemin.jpg", ve_min, width = 14, height = 6, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figs8_vemin.pdf", ve_min, width = 14, height = 6, dpi = 1200)



## selected----------------------------------------------------------------------

# 2) ix_df 에 setting 컬럼 추가
ix_df2 <- ix_df %>%
  left_join(setting_key_df, by = c("region" = "region"))
set.seed(2025)  # 재현 가능하도록 시드 고정
selected_regions <- ix_df2 %>%
  distinct(region, setting) %>%
  group_by(setting) %>%
  slice_sample(n = 2) %>%
  pull(region)


vim_df2 <- vim_df %>%
  left_join(setting_key_df, by = c("region" = "region"))
set.seed(2025)  # 재현 가능하도록 시드 고정

low_fixed <- c("Sergipe", "Goiás")
set.seed(2025)                                 
other_selected <- ix_df2 %>%                  
  distinct(region, setting) %>%                
  filter(setting != "Low") %>%                 
  group_by(setting) %>% 
  slice_sample(n = 2) %>%                     
  pull(region)
selected_regions <- c(low_fixed, other_selected)

ix_df_sel <- ix_df2 %>%
  filter(region %in% selected_regions)

# 5) ix_ann 도 동일하게 처리
ix_ann_sel <- ix_ann %>%
  left_join(setting_key_df, by = c("region" = "region")) %>%
  filter(region %in% selected_regions)

vim_df_sel <- vim_df2 %>%
  filter(region %in% selected_regions)

vim_ann_sel <- vim_ann %>%
  left_join(setting_key_df, by = c("region" = "region")) %>%
  filter(region %in% selected_regions)

ix  <- make_prepost_plot_pair(ix_df2, ix_ann) + plot_annotation(tag_levels = list(c("B")))
vim <- make_prepost_plot_pair(vim_df2, vim_ann, vacc_title = "Vimkunya") + theme(axis.title.y = element_blank())

pw <- 
(ix + vim) + 
  plot_annotation(tag_levels = list(c("A","B")))+
  plot_layout(guides = "collect") 

comb <- 
  (pw & theme(
    legend.position = "bottom",
    legend.justification = "center")) 

ggsave(filename = "02_Outputs/2_1_Figures/comb_fig1.jpg", ix, width = 14, height = 7, dpi = 1200)
ggsave(filename = "02_Outputs/2_1_Figures/figS_vim_epicurves.jpg", vim, width = 14, height = 7, dpi = 1200)

## -----------------------------------------------------------------------------
vaccine_sel  <- "Ixchiq"   
coverage_sel <- "cov50"     
ve_sel <- "VE0"

global_impact_ve <- combined_prepost_case_all_50 %>%
  group_by(VE, Scenario, region, vaccine) %>%
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    impact = diff / total_pre_cases * 100    # % 감소
  )

pre_tbl <- global_impact_ve %>%
  group_by(region, vaccine, VE) %>%
  summarise(
    impact = 0,    
    .groups = "drop"
  ) %>%
  mutate(
    Scenario = "Pre"
  )

plot_long <- bind_rows(global_impact_ve, pre_tbl) %>%
  mutate(
    Scenario = dplyr::recode(Scenario,
                      "Scenario_1" = "1-11y",
                      "Scenario_2" = "12-17y",
                      "Scenario_3" = "18-59y",
                      "Scenario_4" = "60+y")
  ) %>%
  mutate(
    Scenario = factor(Scenario, levels = c("Pre", "1-11y", "12-17y", "18-59y", "60+y")),
    lbl = ifelse(Scenario == "Pre", NA, paste0("-", sprintf("%.1f", impact), "%"))
  )

ix_ve0 <- 
ggplot(plot_long %>% filter(vaccine == vaccine_sel, VE == ve_sel),
       aes(x = Scenario, y = 100 - impact, fill = Scenario)) +
  geom_col(alpha = 0.5, width = 0.75) +
  geom_text(aes(label = lbl), vjust = -0.25, na.rm = TRUE, size = 2.6) +
  facet_wrap(~ region, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, NA)
  ) +
  scale_fill_manual(
    values = c(
      "Pre"     = "grey50",
      "1-11y"   = "#E41A1C",
      "12-17y"  = "#377EB8",
      "18-59y"  = "#4DAF4A",
      "60+y"    = "#984EA3"
    ),
    name = "Vaccination strategy"
  ) +
  labs(
    x = "Vaccination strategy",
    y = "Cumulative predicted symptomatic cases (% of pre-vaccination)"
    #y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_rect(fill = "grey95"),
    axis.text.x       = element_text(angle = 45, hjust = 1)
  )

ve0_all <- 
  (ix_ve0 + vim_ve0) +
  plot_layout(guides = "collect") 

ve0_all <- 
  (ve0_all & theme(
    legend.position = "bottom",
    legend.justification = "center")) +
  plot_annotation(tag_levels = list(c("B","C")))

ggsave(filename = "02_Outputs/2_1_Figures/fig1_ve0.jpg", ve0_all, width = 15, height = 6, dpi = 1200)


### summarised daly pre-post----------------------------------------------------
pre_tbl_daly <- global_impact_nat_daly %>%
  filter(Coverage == "cov50") %>%
  group_by(vaccine, VE) %>%
  summarise(
    impact_lo  = 0,
    impact_mid = 0,
    impact_hi  = 0,
    .groups = "drop"
  ) %>%
  mutate(Scenario = "No vaccination")

global_impact_nat_daly <- global_impact_nat_daly %>% mutate(
  Scenario = dplyr::recode(Scenario,
                           "1" = "1-11y",
                           "2" = "12-17y",
                           "3" = "18-59y",
                           "4" = "60+y")
) %>% filter(Coverage == "cov50")

# 3. Combined long table
plot_long_daly <- bind_rows(global_impact_nat_daly, pre_tbl_daly) %>%
  mutate(
    Scenario = factor(Scenario, levels = c("No vaccination", "1-11y", "12-17y", "18-59y", "60+y")),
    lbl = ifelse(
      Scenario == "Pre", 
      NA, 
      paste0("-", sprintf("%.1f", impact_mid), "% ",
             "(95% UI ", sprintf("%.1f", impact_lo), "–", sprintf("%.1f", impact_hi), "%)")
    )
  )%>%
  mutate(
    facet_label = case_when(
      VE == "VE0" ~ paste(vaccine, "(Disease blocking only)"),
      VE %in% c("VE97.8", "VE98.9") ~ paste(vaccine, "(Disease & infection blocking)"),
      TRUE ~ paste(vaccine, VE)  # fallback
    )
  ) 

v <- sort(unique(plot_long_daly$vaccine))
facet_levels <- as.vector(rbind(
  paste(v, "(Disease blocking only)"),
  paste(v, "(Disease & infection blocking)")
))
plot_long_daly$facet_label <- factor(plot_long_daly$facet_label, levels = facet_levels)


daly_impact <- 
ggplot(plot_long_daly,
       aes(x = Scenario, y = 100 - impact_mid, fill = Scenario)) +
  geom_col(alpha = 0.6, width = 0.75) +
  geom_errorbar(aes(
    ymin = 100 - impact_hi,   
    ymax = 100 - impact_lo   
  ),
  width = 0.12,               
  colour = "black"            
  )+
  #geom_text(aes(label = lbl), vjust = -2, na.rm = TRUE, size = 3) +
  facet_wrap(~ facet_label) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, NA)
  ) +
  scale_fill_manual(
    values = c(
      "No vaccination"     = "grey60",
      "1-11y"   = "#E41A1C",
      "12-17y"  = "#377EB8",
      "18-59y"  = "#4DAF4A",
      "60+y"    = "#984EA3"
    ),
    name = "Age specific vaccination strategy"
  ) +
  labs(
    x = NULL,
    y = "Cumulative DALYs (% of no vaccination)"
  ) +
  theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    #axis.text.y  = element_text(size = 8),
    axis.title.x = element_text(size = 12),                          # 축 제목
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

## pre post fatal ----------------------------------------------------------------------------
pre_tbl_fatal <- global_impact_nat_fatal %>%
  filter(Coverage == "cov50") %>%
  group_by(vaccine, VE) %>%
  summarise(
    impact_lo  = 0,
    impact_mid = 0,
    impact_hi  = 0,
    .groups = "drop"
  ) %>%
  mutate(Scenario = "No vaccination")

global_impact_nat_fatal <- global_impact_nat_fatal %>% mutate(
  Scenario = dplyr::recode(Scenario,
                           "1" = "1-11y",
                           "2" = "12-17y",
                           "3" = "18-59y",
                           "4" = "60+y")
)%>% filter(Coverage == "cov50")

plot_long_fatal <- bind_rows(global_impact_nat_fatal, pre_tbl_fatal) %>%
  mutate(
    Scenario = factor(Scenario, levels = c("No vaccination", "1-11y", "12-17y", "18-59y", "60+y")),
    lbl = ifelse(
      Scenario == "No vaccination", 
      NA, 
      paste0("-", sprintf("%.1f", impact_mid), "% ",
             "(95% UI ", sprintf("%.1f", impact_lo), "–", sprintf("%.1f", impact_hi), "%)")
    )
  ) %>%
  mutate(
    facet_label = case_when(
      VE == "VE0" ~ paste(vaccine, "(Disease blocking only)"),
      VE %in% c("VE97.8", "VE98.9") ~ paste(vaccine, "(Disease & infection blocking)"),
      TRUE ~ paste(vaccine, VE)
    )
  )

v <- sort(unique(plot_long_fatal$vaccine))
facet_levels <- as.vector(rbind(
  paste(v, "(Disease blocking only)"),
  paste(v, "(Disease & infection blocking)")
))
plot_long_fatal$facet_label <- factor(plot_long_fatal$facet_label, levels = facet_levels)


death_impact <- 
ggplot(plot_long_fatal,
       aes(x = Scenario, y = 100 - impact_mid, fill = Scenario)) +
  geom_col(alpha = 0.6, width = 0.75) +
  geom_errorbar(aes(
    ymin = 100 - impact_hi,   
    ymax = 100 - impact_lo   
  ),
  width = 0.12,               
  colour = "black"            
  ) + 
  #geom_text(aes(label = lbl), vjust = -1.5, na.rm = TRUE, size = 3) +
  facet_wrap(~ facet_label) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, NA)
  ) +
  scale_fill_manual(
    values = c(
      "No vaccination"     = "grey60",
      "1-11y"   = "#E41A1C",
      "12-17y"  = "#377EB8",
      "18-59y"  = "#4DAF4A",
      "60+y"    = "#984EA3"
    ),
    name = "Age specific vaccination strategy"
  ) +
  labs(
    x = "Age specific vaccination strategy",
    y = "Cumulative deaths (% of no vaccination)"
  ) +
  theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),                        
    axis.title.y = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

daly_death <- 
  (daly_impact + death_impact) +
  plot_layout(guides = "collect", ncol = 1) 

daly_death_all <- 
  (daly_death & theme(
    legend.position = "right",
    legend.justification = "center")) +
  plot_annotation(tag_levels = list(c("B","C")))

ggsave(filename = "02_Outputs/2_1_Figures/fig1_BC.jpg", daly_death_all, width = 9, height = 9, dpi = 1200)
ggsave(
  filename = "02_Outputs/2_1_Figures/fig1_BC.pdf",  
  plot     = daly_death_all, 
  width    = 10.5, 
  height   = 9.5, 
  dpi      = 300   
)
## -----------------------------------------------------------------------------

observed_all_full <- bind_rows(observed_nat, observed_by_setting)
observed_all_full <- observed_all_full %>%
  tidyr::crossing(
    vaccine = unique(combined_all$vaccine)
  )
observed_all_full <- observed_all_full %>%
  filter(
    !(vaccine == "Ixchiq"   & VE == "VE97.8"),
    !(vaccine == "Vimkunya" & VE == "VE98.9")
  )


## 2 ── one-row data frame for the grey vaccination window
vacc_window_df <- data.frame(
  xmin = vacc_start_week_s1,
  xmax = vacc_end_week_s1
)

global_impact_ve <- combined_all %>%       # ← VE 컬럼 포함된 데이터
  group_by(VE, Scenario, setting, vaccine) %>%                 # ★ VE도 그룹에 추가
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    impact = diff / total_pre_cases * 100
  )

annotation_setting_ve <- global_impact_ve %>%      # already Low / Moderate / …
  mutate(
    Strategy = dplyr::recode(Scenario,
                             "Scenario_1" = "Strategy 1",
                             "Scenario_2" = "Strategy 2",
                             "Scenario_3" = "Strategy 3",
                             "Scenario_4" = "Strategy 4"),
    setting  = factor(setting,                     # make sure order matches plot
                      levels = c("Low","Moderate","High","National"))
  ) %>% 
  group_by(setting, VE, vaccine) %>%                        # exactly the facet grid
  summarise(
    ann = paste0(Strategy, ": ", round(impact, 1), "%", collapse = "\n"),
    .groups = "drop"
  )

## 3 ── plot
library(forcats)
ve_epi_graph_all <- 
  ggplot(combined_all) +
  
  ## ❶ Vaccination window (all facets)
  geom_rect(data = vacc_window_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = .20, inherit.aes = FALSE) +
  
  ## ❷ Post-vaccine UI ribbon
  geom_ribbon(aes(Week,
                  ymin = post_weekly_low95,
                  ymax = post_weekly_hi95,
                  fill = Scenario),
              alpha = .25) +
  
  ## ❸ Pre-vaccine UI ribbon
  geom_ribbon(aes(Week,
                  ymin = lo95,
                  ymax = hi95),
              fill = "lightgrey", alpha = .50) +
  
  ## ❹ Post-vaccine median
  geom_line(aes(Week, post_cases, colour = Scenario), size = .25) +
  
  ## ❺ Pre-vaccine median
  geom_line(aes(Week, pre_cases, group = Scenario),
            colour = "black", linetype = "dashed", size = .25) +
  
  ## ❻ Facet-specific annotation (setting × VE × vaccine)
  geom_text(data = annotation_setting_ve,
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE) +
  geom_point(
    data        = observed_all_full,
    aes(x = Week, y = Observed),        
    size        = 0.2,
    inherit.aes = FALSE
  )+
  
  ## ❼ 색상 / 채우기
  scale_fill_brewer(palette = "Set1",
                    name   = "Vaccination strategy",
                    labels = c("Scenario_1" = "1–11 y",
                               "Scenario_2" = "12–17 y",
                               "Scenario_3" = "18–59 y",
                               "Scenario_4" = "60+ y")) +
  scale_colour_brewer(palette = "Set1",
                      name   = "Vaccination strategy",
                      labels = c("Scenario_1" = "1–11 y",
                                 "Scenario_2" = "12–17 y",
                                 "Scenario_3" = "18–59 y",
                                 "Scenario_4" = "60+ y")) +
  
  scale_y_continuous(labels = comma) +
  labs(x = "Week", y = "Predicted symptomatic cases") +
  
  ## ❽ Facet grid  (행 = setting, 열 = vaccine → VE)
  facet_grid(
    rows = vars(fct_rev(setting)),   # National이 맨 위
    cols = vars(vaccine, VE),        # ← vaccine 먼저, VE는 그 안쪽
    scales   = "free_y",
    drop   = TRUE, 
    labeller = labeller(
      .rows   = label_value,         # 행 strip: 값만 (National, High…)
      vaccine = label_value,         # 열 strip: 백신 이름만
      VE      = c(
        "VE0"   = "Disease blocking only",
        "VE98.9" = "Disease & infection blocking",
        "VE97.8" = "Disease & infection blocking"
      )
    )
  ) +
  
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y    = element_text(angle = 0),
    plot.margin     = margin(5, 20, 5, 5),
    legend.position = "bottom"
  )

# two vacc
two_vacc_epi <- 
ggplot(combined_all) +
  
  ## ❶ vaccination window (all facets)
  geom_rect(data = vacc_window_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = .20, inherit.aes = FALSE) +
  
  ## ❷ post-vaccine UI ribbon
  geom_ribbon(aes(Week, ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = .25) +
  
  ## ❸ pre-vaccine UI ribbon
  geom_ribbon(aes(Week, ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = .50) +
  
  ## ❹ post-vaccine median
  geom_line(aes(Week, post_cases, colour = Scenario), size = .25) +
  
  ## ❺ pre-vaccine median
  geom_line(aes(Week, pre_cases, group = Scenario),
            colour = "black", linetype = "dashed", size = .25) +
  geom_point(
    data         = observed_all_full,
    aes(x = Week, y = Observed),
    shape        = 16,      # 동그라미
    size         = 1.5,
    colour       = "black", # 원한다면 그룹별로 색을 바꿀 수도 있습니다
    inherit.aes  = FALSE
  ) +
  
  ## ❻ facet-specific annotation (setting × VE × vaccine)
  geom_text(data = annotation_setting_ve,
            aes(x = Inf, y = Inf, label = ann),
            hjust = 1.05, vjust = 1.1, size = 3,
            inherit.aes = FALSE) +
  
  ## ❼ 색상 / 채우기
  scale_fill_brewer(palette = "Set1",
                    name   = "Vaccination strategy",
                    labels = c("Scenario_1" = "1–11 y",
                               "Scenario_2" = "12–17 y",
                               "Scenario_3" = "18–59 y",
                               "Scenario_4" = "60+ y")) +
  scale_colour_brewer(palette = "Set1",
                      name   = "Vaccination strategy",
                      labels = c("Scenario_1" = "1–11 y",
                                 "Scenario_2" = "12–17 y",
                                 "Scenario_3" = "18–59 y",
                                 "Scenario_4" = "60+ y")) +
  
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Week", y = "Predicted symptomatic cases") +
  
  ## ❽ facet grid  (행 = setting, 열 = VE → vaccine)
  facet_grid(
    rows = vars(fct_rev(setting)),      # National이 맨 위
    cols = vars(VE, vaccine),           # VE → vaccine 순
    scales   = "free_y",
    labeller = labeller(
      .rows   = label_value,            # 행 strip: 값만 출력 (National, High…)
      vaccine = label_value,
      VE      = c(
        "VE0"   = "Disease blocking only",
        "VE98.9" = "Disease & infection blocking"
      )            # 열 strip: “97.5%” 처럼 값만
                   # 열 strip: 백신 이름만
    )
  ) +
  
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y    = element_text(angle = 0),
    plot.margin     = margin(5, 20, 5, 5),
    legend.position = "bottom"
  )

ggsave(filename = "02_Outputs/2_1_Figures/ve_epi_graph_all.jpg", ve_epi_graph_all, width = 12, height = 7, dpi = 1200)

combined_nat <- combined_prepost_case_all_50 %>%          # ← 이미 VE 포함
  group_by(VE, Scenario, Week, vaccine) %>%                 # ★ VE 추가
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),  # post_cases, pre_cases, hi95 …
    .groups = "drop"
  )

global_impact_nat <- combined_nat %>%          # ← VE, Scenario, Week, post/pre_cases
  group_by(VE, Scenario, vaccine) %>%                  # ★ VE도 포함
  summarise(
    total_post_cases = sum(post_cases, na.rm = TRUE),
    total_post_cases_lo = sum(post_cases_lo, na.rm = TRUE),
    total_post_cases_hi = sum(post_cases_hi, na.rm = TRUE),
    total_pre_cases  = sum(pre_cases,  na.rm = TRUE),
    total_pre_cases_lo = sum(pre_cases_lo, na.rm = TRUE),
    total_pre_cases_hi = sum(pre_cases_hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff   = total_pre_cases - total_post_cases,
    diff_lo = total_pre_cases_lo - total_post_cases_lo,
    diff_hi = total_pre_cases_hi - total_post_cases_hi,
    impact = diff / total_pre_cases * 100,
    impact_lo = diff_lo / total_pre_cases_lo * 100,
    impact_hi = diff_hi / total_pre_cases_hi * 100
  )

global_impact_nat <- combined_nat %>% 
  group_by(VE, Scenario, vaccine) %>% 
  summarise(
    total_post_cases     = sum(post_cases,     na.rm = TRUE),
    total_post_cases_lo  = sum(post_cases_lo,  na.rm = TRUE),
    total_post_cases_hi  = sum(post_cases_hi,  na.rm = TRUE),
    total_pre_cases      = sum(pre_cases,      na.rm = TRUE),
    total_pre_cases_lo   = sum(pre_cases_lo,   na.rm = TRUE),
    total_pre_cases_hi   = sum(pre_cases_hi,   na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    impact_mid = (1 - total_post_cases     / total_pre_cases)    * 100,
    impact_lo  = (1 - total_post_cases_hi  / total_pre_cases_hi) * 100,  # 최소
    impact_hi  = (1 - total_post_cases_lo  / total_pre_cases_lo) * 100   # 최대
  )

## 2) 주석 텍스트 만들기 ----------------------------------------------------
annotation_nat <- global_impact_nat %>%
  mutate(
    Strategy = dplyr::recode(Scenario,
                             "1" = "1-11y",
                             "2" = "12-17y",
                             "3" = "18-59y",
                             "4" = "60+y")
  )  %>%
  group_by(VE, vaccine, Strategy, Scenario, Coverage) %>%
  summarise(
    ann = paste0(Strategy, ": ",
                 round(impact_mid, 1), "% ",
                 "(95% UI ", round(impact_lo, 1), "–", round(impact_hi, 1), ")",
                 collapse = "\n"),
    .groups = "drop"
  ) 

vacc_tbl <- imap_dfr(                                   # ─ region loop ─
  postsim_vc_ixchiq_model,
  function(region_list, region_nm) {
    
    imap_dfr(region_list,                               # ─ VE loop ─
             function(ve_list, ve_tag) {
               
               imap_dfr(ve_list,                               # ─ coverage loop ─
                        function(res, cov_tag) {
                          tibble(
                            region   = region_nm,                     # 예: "Bahia"
                            VE       = ve_tag,                        # 예: "VE98.9"
                            coverage = cov_tag,                       # 예: "cov5"
                            start    = res$vacc_weeks$scenario1$start,
                            end      = res$vacc_weeks$scenario1$end
                          )
                        }
               )
             }
    )
  }
)

vacc_tbl <- vacc_tbl %>% filter(coverage == "cov50")

observed_nat <- observed_all %>% 
  group_by(Week) %>%                         # 주차별
  summarise(Observed = sum(Observed, na.rm = TRUE),
            .groups = "drop")

ve_levels <- unique(combined_nat$VE)         # "VE0" "VE50" "VE98.9"

observed_nat_full <- crossing(observed_nat, VE = ve_levels)

observed_nat_full <- observed_nat_full %>%
  tidyr::crossing(
    vaccine = unique(combined_all$vaccine)
  )
observed_nat_full <- observed_nat_full %>%
  filter(
    !(vaccine == "Ixchiq"   & VE == "VE97.8"),
    !(vaccine == "Vimkunya" & VE == "VE98.9")
  )

ve_lab <- c(
  "VE0"    = "Disease-blocking only",
  "VE97.8" = "Disease & infection-blocking",
  "VE98.9" = "Disease & infection-blocking"
)


ve_nat<- collapse_ve <- function(df) {
  df %>%
    mutate(
      VE_group = ifelse(VE == "VE0",
                        "Disease-blocking only",
                        "Disease & infection-blocking")
    )
}

combined_nat   <- collapse_ve(combined_nat)
vacc_tbl       <- collapse_ve(vacc_tbl)
annotation_nat <- collapse_ve(annotation_nat)

ve_nat<- 

ggplot(combined_nat, aes(x = Week)) +
  
  ## (a) 백신 투입 회색 음영
  geom_rect(data = vacc_tbl,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "lightgrey", alpha = .05,
            show.legend = FALSE) +
  
  ## (b) 포스트-백신 UI 리본
  geom_ribbon(aes(ymin = post_weekly_low95, ymax = post_weekly_hi95,
                  fill = Scenario), alpha = .1) +
  
  ## (c) 프리-백신 UI 리본
  geom_ribbon(aes(ymin = lo95, ymax = hi95),
              fill = "lightgrey", alpha = .5) +
  
  ## (d) 포스트-백신 선
  geom_line(aes(y = post_cases, colour = Scenario, linetype = "With"), size = .35) +
  
  ## (e) 프리-백신 선
  geom_line(aes(y = pre_cases, colour = Scenario, linetype = "Without"),
            color = "black", size = .35) +
  
  ## (f) 관측치 점
  geom_point(data = observed_nat_full,
             aes(y = Observed, shape = "Reported symptomatic cases"), size = .3) +
  
  ## (g) fill 범례 (리본 + 선을 함께 표현)
  scale_fill_manual(
    name = "Vaccination strategy",
    values = c(
      "Scenario_1" = "#E41A1C",
      "Scenario_2" = "#377EB8",
      "Scenario_3" = "#4DAF4A",
      "Scenario_4" = "#984EA3"
    ),
    labels = c(
      "Scenario_1" = "1–11 y",
      "Scenario_2" = "12–17 y",
      "Scenario_3" = "18–59 y",
      "Scenario_4" = "60+ y"
    ),
    guide = guide_legend(
      override.aes = list(
        color    = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        fill     = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        alpha    = 0.3,
        linetype = "solid",
        size     = 0.4
      )
    )
  ) +
  
  ## (h) 선 색 범례 숨김
  scale_colour_manual(
    values = c(
      "Scenario_1" = "#E41A1C",
      "Scenario_2" = "#377EB8",
      "Scenario_3" = "#4DAF4A",
      "Scenario_4" = "#984EA3"
    ),
    guide = "none"
  ) +
  
  ## (i) 선 종류 범례
  scale_linetype_manual(
    values = c(With = "solid", Without = "dashed"),
    labels = c(With = "With vaccination", Without = "Without vaccination"),
    name = ""
  ) +
  
  ## (j) 관측치 점 범례
  scale_shape_manual(
    name = "",
    values = c("Reported symptomatic cases" = 16)
  ) +
  
  ## (k) Y축 포맷
  scale_y_continuous(labels = comma) +
  
  ## (l) 축 라벨
  labs(x = "Week", y = "Predicted symptomatic cases") +
  
  ## (m) Facet
  facet_grid(rows = vars(vaccine),
             cols = vars(VE_group), scales = "free_y") +
  
  ## (n1) 숫자 박스 주석 (e.g. "32.7%")
  geom_label(data = annotation_nat,
             aes(x = Inf, y = Inf, label = ann),
             hjust = 1.05, vjust = 1.6,
             size = 3, fill = "white", label.size = 0.2, alpha = 0.95,
             inherit.aes = FALSE) +
  
  ## (n2) 상단 제목 주석
  geom_text(data = annotation_nat,
            aes(x = Inf, y = Inf, label = "% Reduction in predicted symptomatic cases"),
            hjust = 1.05, vjust = 2.8,
            size = 2.9, fontface = "plain",
            inherit.aes = FALSE) +
  
  ## (o) 테마
  theme_pubclean(base_size = 9) +
  theme(
    strip.text.y   = element_text(angle = 0),
    legend.key     = element_rect(fill = "white", colour = NA),
    legend.box     = "vertical",
    legend.position = "right",
    plot.margin    = margin(5, 25, 5, 5)
  )

ve_nat <- 
ve_nat + plot_annotation(tag_levels = "A")

ggsave(filename = "02_Outputs/2_1_Figures/fig1_a.jpg", ve_nat, width = 10, height = 6, dpi = 1200)


fig1_a <- 
  ggarrange(ve_nat ,
            ncol = 1, nrow = 1,
            labels = c("A"),
            common.legend = TRUE,
            legend = "bottom",
            align = "none")

fig1_b <- 
  ggarrange(ve_subnat ,
            ncol = 1, nrow = 1,
            labels = c("B"),
            common.legend = TRUE,
            legend = "bottom",
            align = "none")


### final VE nat ---------------------------------------------------------------

# ─────────────────────────────────────────────────────────────
# 1. Re-define VE_group (VE0 = disease-blocking, others = infection-blocking)
# ─────────────────────────────────────────────────────────────
combined_nat <- combined_nat %>% 
  mutate(
    VE_group = if_else(VE == "VE0",
                       "Disease-blocking only",
                       "Disease & infection-blocking")
  )

# ─────────────────────────────────────────────────────────────
# 2. Viridis palette WITHOUT yellow
#    option = "C": take indices 1,2,3,4   (5=lime, 6=yellow -> excluded)
# ─────────────────────────────────────────────────────────────
pal_vals <- setNames(
  viridis(6, option = "C")[c(1, 2, 3, 4)],
  c("1–11 y", "12–17 y", "18–59 y", "60+ y")
)

ve_pal <- setNames(
  viridis(2, option = "C"),    
  c("VE0", "VE98.9")
)

label_k <- function(x) {
  paste0(formatC(x / 1000, format = "f", digits = 0), " K")
}

scenario_labels <- c(
  "Scenario_1" = "1–11 y",
  "Scenario_2" = "12–17 y",
  "Scenario_3" = "18–59 y",
  "Scenario_4" = "60+ y"
)

annotation_nat_revised <- global_impact_nat %>% 
  mutate(
    Strategy  = dplyr::recode(
      Scenario,
      "1" = "1–11 y",
      "2" = "12–17 y",
      "3" = "18–59 y",
      "4" = "60+ y"
    ),
    VE_group = if_else(VE == "VE0",
                       "Disease-blocking only",
                       "Disease & infection-blocking")
  ) %>%                       # 백신·VE_group·Strategy마다 1행
  group_by(vaccine, VE_group, VE, Strategy, Coverage) %>% 
  summarise(
    ann    = paste0(round(impact_mid, 1), " %"),
    ann_lo = paste0(round(impact_lo, 1), " %"), 
    ann_hi = paste0(round(impact_hi, 1), " %"), 
    .groups = "drop"
  )%>% filter(Coverage == "cov50")

# ─────────────────────────────────────────────────────────────
# 3. Plot-building function (4 Scenario panels -> 2×2 layout)
# ─────────────────────────────────────────────────────────────
make_ve_plot <- function(vacc, ve_group_label, ve_code) {
  
  # subset for the given vaccine & VE level
  df <- combined_nat %>% 
    filter(vaccine  == vacc,
           VE_group == ve_group_label)
  
  
  # observed points – filter by vaccine only
  obs_df <- observed_nat_full %>% 
    filter(vaccine == vacc)
  
  ann_df <- annotation_nat_revised %>% 
    filter(vaccine  == vacc,
           VE_group == ve_group_label,
           VE       == ve_code) %>% 
    dplyr::rename(Scenario = Strategy)%>% 
    mutate(                                  
      label = paste0(
        ann,                       
        "(95% UI: ",             
        ann_lo, " – ", ann_hi,    
        ")"                        
      )
    )
  
  fill_vals <- c("Without vaccination" = "grey85", pal_vals)
  
  ggplot(df, aes(x = Week)) +
    
    # baseline (no vaccination) ribbon
    geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = "Without vaccination"),
                #fill = "grey85", 
                colour = NA) +
    
    # strategy ribbon
    geom_ribbon(aes(ymin = post_weekly_low95, ymax = post_weekly_hi95,
                    fill = Scenario),
                alpha = .35, colour = NA) +
    
    # strategy line
    geom_line(aes(y = post_cases, colour = Scenario), size = .35) +
    
    # observed symptomatic cases (points)
    geom_point(data = obs_df,
               aes(y = Observed,
                   shape = "Reported symptomatic cases"),
               colour = "black",
               size = 1) +
    
    # viridis colours
    scale_fill_manual(values = fill_vals, name = "Vaccination strategy") +
    scale_colour_manual(values = pal_vals, guide = "none") +
    scale_shape_manual(values = c("Reported symptomatic cases" = 16),
                       name   = "") +
    scale_y_continuous(
      labels = label_k,
      expand = expansion(mult = c(0, .02))
    ) +
    
    labs(
      x = "Week",
      y = "Predicted symptomatic cases"
    ) +
    
    # 2×2 layout: 4 scenarios
    facet_wrap(~ Scenario, ncol = 2, scales = "free_y") +
    
    # clean theme
    theme_pubclean(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 8)
    ) + 
    geom_label(
      data = ann_df,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.04, vjust = 1.55,
      size = 3, 
      label.size = 0, 
      alpha = .95,
      fill = "white", 
      inherit.aes = FALSE
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = "% Reduction in predicted symptomatic cases",
      hjust = 1.05, vjust = 1.1,
      size = 2.9, fontface = "plain"
    )
}

# ─────────────────────────────────────────────────────────────
# 4. Build two rows (VE0 top, VE98.9 bottom) for Ixchiq
# ─────────────────────────────────────────────────────────────
make_ve_plot("Ixchiq",
             "Disease-blocking only",   "VE0")

make_ve_plot("Ixchiq",
             "Disease & infection-blocking", "VE98.9")  # or "VE97.8"

make_ve_plot("Vimkunya",
             "Disease & infection-blocking", "VE97.8") 

#### 1 vacicne with 2 ves 
### 2 VEs within strategy 
ve_pal_ix <- c(
  VE0    = "#00468B",  # disease-blocking only
  VE98.9 = "#ED0000"   # disease + infection-blocking
)

ve_pal_vim <- c(
  VE0    = "#00468B",  # disease-blocking only
  VE97.8 = "#ED0000"   # disease + infection-blocking
)

vacc_tbl_nat <- vacc_tbl %>% filter(region == "Ceará")%>%
  crossing(Scenario = unique(combined_nat$Scenario)) 

make_ve_plot_contrast_ix <- function(vacc) {
  
  df <- combined_nat %>%
    filter(vaccine == vacc)
  
  obs_df <- observed_nat_full %>%
    filter(vaccine == vacc)
  
  ann_df <- annotation_nat_revised %>%
    filter(vaccine == vacc) %>%
    dplyr::rename(Scenario = Strategy) %>% 
    mutate(
      VE_lab = case_when(
        VE == "VE0"     ~ "Disease-blocking only",
        VE == "VE98.9"  ~ "Disease & infection-blocking",
        VE == "VE97.8"  ~ "Disease & infection-blocking", 
        TRUE            ~ VE
      ),
      label  = paste0(VE_lab, ": ", ann,
                      " (95% UI ", ann_lo, "–", ann_hi, ")")
    )
  
  ann_lab <- ann_df %>% 
    arrange(Scenario, VE) %>%                       
    group_by(Scenario) %>% 
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  
  scenarios <- unique(df$Scenario)
  
  vacc_tbl_nat <- vacc_tbl %>%
    filter(region == "Ceará") %>%
    tidyr::crossing(Scenario = scenarios)
  
  
  ggplot(df, aes(x = Week)) +
    
    geom_rect(
      data = vacc_tbl_nat,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill  = "lightgrey",
      alpha = .2,
      show.legend = FALSE
    ) + 
    
    geom_ribbon(aes(ymin = lo95, ymax = hi95,
                    fill = "Without vaccination"),
                colour = NA, alpha = 0.7) +   # alpha 복원
    
    geom_line(aes(y = pre_cases,
                  linetype = "Without vaccination",
                  colour   = "Without vaccination",
                  group = Scenario),
              size = 0.7) +
    
    geom_ribbon(aes(ymin = post_weekly_low95, ymax = post_weekly_hi95,
                    fill = VE),
                alpha = 0.3, colour = NA) +
    
    geom_line(aes(y = post_cases,
                  linetype = VE, colour = VE),
              size = 0.9) +
    
    geom_point(data = obs_df,
               aes(y = Observed,
                   shape = "Reported symptomatic cases"),
               colour = "black", size = 0.8) +
    
    scale_fill_manual(
      values = c("Without vaccination" = "grey", ve_pal_ix),
      breaks = c("Without vaccination", names(ve_pal_ix)),
      labels = c("Without vaccination",
                 "Disease-blocking only",
                 "Disease & infection-blocking"),
      name   = ""
    ) +
    scale_colour_manual(
      values = c("Without vaccination" = "black", ve_pal_ix),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("Without vaccination" = "solid",
                 "VE0"               = "dotted",
                 "VE98.9"            = "dashed"),
      breaks = c("Without vaccination", names(ve_pal_ix)),
      labels = c("Without vaccination",
                 "Disease-blocking only",
                 "Disease & infection-blocking"),
      name   = "Vaccination strategy"
    ) +
    scale_shape_manual(values = c("Reported symptomatic cases" = 4),
                       name = "") +
    scale_y_continuous(
      labels = label_k,
      expand = expansion(mult = c(0, .02))
    ) +  
    
    ylab("Predicted symptomatic cases") +
    facet_wrap(~Scenario, ncol = 2, scales = "free_y") +
    
    theme_pubclean(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8)) +
    
    geom_label(data = ann_lab,
               aes(x = Inf, y = Inf, label = label),
               hjust = 1.04, vjust = 1.55,
               size  = 2.7, fill = NA, label.size = 0,
               inherit.aes = FALSE) +
    annotate("text", x = Inf, y = Inf, size = 3,
             label = "% Reduction in predicted symptomatic cases",
             hjust = 1.05, vjust = 1.1)
}

make_ve_plot_contrast_ix("Ixchiq")

make_ve_plot_contrast_vim <- function(vacc) {
  
  df <- combined_nat %>%
    filter(vaccine == vacc)
  
  obs_df <- observed_nat_full %>%
    filter(vaccine == vacc)
  
  ann_df <- annotation_nat_revised %>%
    filter(vaccine == vacc) %>%
    dplyr::rename(Scenario = Strategy) %>% 
    mutate(
      VE_lab = case_when(
        VE == "VE0"     ~ "Disease-blocking only",
        VE == "VE98.9"  ~ "Disease & infection-blocking",
        VE == "VE97.8"  ~ "Disease & infection-blocking", 
        TRUE            ~ VE
      ),
      label  = paste0(VE_lab, ": ", ann,
                      " (95% UI ", ann_lo, "–", ann_hi, ")")
    )
  
  ann_lab <- ann_df %>% 
    arrange(Scenario, VE) %>%                       
    group_by(Scenario) %>% 
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  
  ggplot(df, aes(x = Week)) +
    
    # 백신 공급 시기 음영
    geom_rect(
      data = vacc_tbl_nat,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill  = "lightgrey",
      alpha = .2,
      show.legend = FALSE
    ) + 
    
    geom_ribbon(aes(ymin = lo95, ymax = hi95,
                    fill = "Without vaccination"),
                colour = NA, alpha = 0.7) +   # alpha 복원
    
    geom_line(aes(y = pre_cases,
                  linetype = "Without vaccination",
                  colour   = "Without vaccination",
                  group = Scenario),
              size = 0.7) +
    
    geom_ribbon(aes(ymin = post_weekly_low95, ymax = post_weekly_hi95,
                    fill = VE),
                alpha = 0.3, colour = NA) +
    
    geom_line(aes(y = post_cases,
                  linetype = VE, colour = VE),
              size = 0.9) +
    
    geom_point(data = obs_df,
               aes(y = Observed,
                   shape = "Reported symptomatic cases"),
               colour = "black", size = 0.8) +
    
    scale_fill_manual(
      values = c("Without vaccination" = "grey85", ve_pal_vim),
      breaks = c("Without vaccination", names(ve_pal_vim)),
      labels = c("Without vaccination",
                 "Disease-blocking only",
                 "Disease & infection-blocking"),
      name   = ""
    ) +
    scale_colour_manual(
      values = c("Without vaccination" = "black", ve_pal_vim),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c("Without vaccination" = "solid",
                 "VE0"               = "dotted",
                 "VE97.8"            = "dashed"),
      breaks = c("Without vaccination", names(ve_pal_vim)),
      labels = c("Without vaccination",
                 "Disease-blocking only",
                 "Disease & infection-blocking"),
      name   = "Vaccination strategy"
    ) +
    scale_shape_manual(values = c("Reported symptomatic cases" = 4),
                       name = "") +
    scale_y_continuous(
      labels = label_k,
      expand = expansion(mult = c(0, .02))
    ) +
    
    ylab("Predicted symptomatic cases") +
    facet_wrap(~Scenario, ncol = 2, scales = "free_y") +
    
    theme_pubclean(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8)) +
    
    geom_label(data = ann_lab,
               aes(x = Inf, y = Inf, label = label),
               hjust = 1.04, vjust = 1.55,
               size  = 2.7, fill = NA, label.size = 0,
               inherit.aes = FALSE) +
    annotate("text", x = Inf, y = Inf, size = 3,
             label = "% Reduction in predicted symptomatic cases",
             hjust = 1.05, vjust = 1.1)
}

make_ve_plot_contrast_vim("Vimkunya")

## make both
draw_key_rectline <- function(data, params, size) {
  grid::grobTree(
    ggplot2::draw_key_rect(data, params, size), 
    grid::segmentsGrob(                         
      x0 = unit(0.15, "npc"), x1 = unit(0.85, "npc"),
      y0 = unit(0.50, "npc"), y1 = unit(0.50, "npc"),
      gp = grid::gpar(
        col   = data$colour  %||% "black",
        lwd   = (data$size   %||% 0.6) * .pt,
        lty   = data$linetype%||% 1,
        lineend = "butt"
      )
    )
  )
}

# original
make_ve_plot_contrast_both <- function(vaccines = c("Ixchiq", "Vimkunya")) {
  library(ggnewscale)
  
  ## 1) 데이터 전처리
  df <- combined_nat %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    dplyr::mutate(
      Scenario = dplyr::case_when(
        Scenario == "Scenario_1" ~ "1–11 y",
        Scenario == "Scenario_2" ~ "12–17 y",
        Scenario == "Scenario_3" ~ "18–59 y",
        Scenario == "Scenario_4" ~ "60+ y",
        TRUE ~ Scenario
      ),
      VE_lab = dplyr::if_else(VE == "VE0", "Disease blocking only",
                              "Disease & infection blocking"),
      VE_lab = factor(VE_lab,
                      levels = c("Disease blocking only",
                                 "Disease & infection blocking"))
    )
  
  scenarios <- unique(df$Scenario)
  
  pre_base <- df %>%
    dplyr::filter(VE == "VE0", Scenario == "1–11 y") %>%
    dplyr::distinct(vaccine, Week, .keep_all = TRUE) %>%
    dplyr::select(-Scenario)
  pre_df <- tidyr::crossing(pre_base, Scenario = scenarios)
  
  obs_df <- observed_nat_full %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    tidyr::crossing(Scenario = scenarios)
  
  vacc_tbl_nat <- vacc_tbl %>%
    dplyr::filter(region == "Ceará") %>%
    tidyr::crossing(Scenario = scenarios)
  
  ann_df <- annotation_nat_revised %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    dplyr::rename(Scenario = Strategy) %>%
    dplyr::mutate(
      VE_lab = dplyr::if_else(VE == "VE0", "Disease blocking only",
                              "Disease & infection blocking"),
      VE_lab = factor(VE_lab,
                      levels = c("Disease blocking only",
                                 "Disease & infection blocking")),
      label = paste0(VE_lab, ": ", ann, " (95% UI ", ann_lo, "–", ann_hi, ")")
    )
  ann_lab <- ann_df %>%
    dplyr::arrange(vaccine, Scenario, VE_lab) %>%
    dplyr::group_by(vaccine, Scenario) %>%
    dplyr::summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  
  ## 팔레트
  fill_pal <- c(
    "No vaccination"               = "grey50",
    "Disease blocking only"        = "#B5DBFF",
    "Disease & infection blocking" = "#F9C2A8"
  )
  col_pal <- c(
    "No vaccination"               = "black",
    "Disease blocking only"        = "#0072B2",
    "Disease & infection blocking" = "#E53935"
  )
  lty_pal <- c(
    "No vaccination"               = "solid",
    "Disease blocking only"        = "dotted",
    "Disease & infection blocking" = "23"
  )
  
  ## 첫 줄 병합용 공통 변수
  pre_df_legend <- pre_df %>% dplyr::mutate(legend_first = "No vaccination")
  obs_df_legend <- obs_df %>% dplyr::mutate(legend_first = "Reported symptomatic cases")
  
  ## 2) 플롯
  ggplot(df, aes(x = Week)) +
    # 공급기간 음영
    geom_rect(
      data = vacc_tbl_nat,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#C8E6C9", alpha = .1
    ) +
    
    ## ── [첫째 줄] "No vaccination"(제목) → 회색 박스 → × Reported symptomatic cases
    geom_ribbon(
      data = pre_df_legend,
      aes(x = Week, ymin = lo95, ymax = hi95, fill = legend_first),
      colour = NA, alpha = 0.30, inherit.aes = FALSE
    ) +
    geom_point(
      data = obs_df_legend,
      aes(x = Week, y = Observed, shape = legend_first, fill = legend_first),
      colour = "black", size = .9, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name   = "No vaccination",   # ← 제목을 텍스트로 사용 (라벨을 비워 상자만 보이게)
      breaks = c("No vaccination", "Reported symptomatic cases"),
      values = c("No vaccination" = fill_pal[["No vaccination"]],
                 "Reported symptomatic cases" = NA),
      labels = c("", "Reported symptomatic cases"),  # ← 첫 키 라벨 비움
      guide  = guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        key_glyph = draw_key_rectline,
        override.aes = list(
          fill     = c(fill_pal[["No vaccination"]], NA),
          colour   = c(col_pal[["No vaccination"]], "black"),
          linetype = c(lty_pal[["No vaccination"]], "blank"),
          shape    = c(NA, 4),
          alpha    = c(.4, 1),
          size     = c(.6, 3)  # 두 번째 값이 '×' 크기
        )
      )
    ) +
    # 같은 제목으로 두되, 가이드는 숨김(키는 위 fill 가이드에서 그림)
    scale_shape_manual(
      name   = "No vaccination",
      breaks = c("No vaccination", "Reported symptomatic cases"),
      values = c("No vaccination" = NA, "Reported symptomatic cases" = 4),
      guide  = "none"
    ) +
    
    ## 둘째 줄(보호효과) 만들기 위해 fill 스케일 리셋
    ggnewscale::new_scale("fill") +
    
    ## ── [둘째 줄] Vaccine protection (파랑/주황 박스)
    geom_ribbon(
      aes(ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = VE_lab),
      alpha = 1, colour = NA
    ) +
    geom_line(
      aes(y = post_cases, linetype = VE_lab, colour = VE_lab),
      size = 1.3
    ) +
    scale_fill_manual(
      name   = "Vaccine protection",
      breaks = c("Disease blocking only", "Disease & infection blocking"),
      values = c(
        "Disease blocking only"        = fill_pal[["Disease blocking only"]],
        "Disease & infection blocking" = fill_pal[["Disease & infection blocking"]]
      ),
      guide  = guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        key_glyph = draw_key_rectline,
        title.position = "left",
        title.theme    = element_text(margin = margin(r = 10)), # 살짝 오른쪽으로
        override.aes = list(
          fill     = unname(fill_pal[c("Disease blocking only",
                                       "Disease & infection blocking")]),
          colour   = unname(col_pal[c("Disease blocking only",
                                      "Disease & infection blocking")]),
          linetype = unname(lty_pal[c("Disease blocking only",
                                      "Disease & infection blocking")]),
          size     = .6, alpha = 1
        )
      )
    ) +
    
    # 라인 범례 숨김(박스로 충분)
    scale_colour_manual(values = col_pal, guide = "none") +
    scale_linetype_manual(values = lty_pal, guide = "none") +
    
    # pre-vacc 선
    geom_line(
      data = pre_df,
      aes(x = Week, y = pre_cases,
          linetype = "No vaccination", colour = "No vaccination"),
      size = .7, inherit.aes = FALSE
    ) +
    
    ## 축/패싯/테마
    scale_y_continuous(labels = label_k, expand = expansion(mult = c(0, .02))) +
    ylab("Predicted symptomatic cases") +
    facet_grid(Scenario ~ vaccine, scales = "free_y") +
    theme_pubclean(base_size = 10) +
    theme(
      legend.position      = "bottom",
      legend.box           = "vertical",     # 두 줄
      legend.justification = "left",
      legend.box.just      = "left",
      legend.direction     = "horizontal",
      legend.title         = element_text(size = 12),
      legend.text          = element_text(size = 12),
      legend.spacing.x     = unit(6, "pt"),
      strip.text           = element_text(size = 11, face = "bold"),
      axis.text.y          = element_text(size = 13),
      axis.title.y         = element_text(size = 13),
      axis.title.x         = element_text(size = 13)
    ) +
    
    ## 주석
    geom_label(
      data = ann_lab,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.04, vjust = 1.55,
      size = 4, fill = NA, label.size = 0
    ) +
    annotate(
      "text", x = Inf, y = Inf, size = 4,
      label = "% Reduction in predicted symptomatic cases",
      hjust = 1.05, vjust = 1.1, fontface = "bold"
    )
}

# updatd
make_ve_plot_contrast_both <- function(vaccines = c("Ixchiq", "Vimkunya")) {
  library(ggnewscale)
  
  ## 1) 데이터 전처리
  df <- combined_nat %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    dplyr::mutate(
      Scenario = dplyr::case_when(
        Scenario == "Scenario_1" ~ "1–11 y",
        Scenario == "Scenario_2" ~ "12–17 y",
        Scenario == "Scenario_3" ~ "18–59 y",
        Scenario == "Scenario_4" ~ "60+ y",
        TRUE ~ Scenario
      ),
      VE_lab = dplyr::if_else(VE == "VE0", "Disease blocking only",
                              "Disease & infection blocking"),
      VE_lab = factor(VE_lab,
                      levels = c("Disease blocking only",
                                 "Disease & infection blocking"))
    )
  
  scenarios <- unique(df$Scenario)
  
  pre_base <- df %>%
    dplyr::filter(VE == "VE0", Scenario == "1–11 y") %>%
    dplyr::distinct(vaccine, Week, .keep_all = TRUE) %>%
    dplyr::select(-Scenario)
  pre_df <- tidyr::crossing(pre_base, Scenario = scenarios)
  
  obs_df <- observed_nat_full %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    tidyr::crossing(Scenario = scenarios)
  
  vacc_tbl_nat <- vacc_tbl %>%
    dplyr::filter(region == "Ceará") %>%
    tidyr::crossing(Scenario = scenarios)
  
  ann_df <- annotation_nat_revised %>%
    dplyr::filter(vaccine %in% vaccines) %>%
    dplyr::rename(Scenario = Strategy) %>%
    dplyr::mutate(
      VE_lab = dplyr::if_else(VE == "VE0", "Disease blocking only",
                              "Disease & infection blocking"),
      VE_lab = factor(VE_lab,
                      levels = c("Disease blocking only",
                                 "Disease & infection blocking")),
      label = paste0(VE_lab, ": ", ann, " (95% UI ", ann_lo, "–", ann_hi, ")")
    )
  ann_lab <- ann_df %>%
    dplyr::arrange(vaccine, Scenario, VE_lab) %>%
    dplyr::group_by(vaccine, Scenario) %>%
    dplyr::summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  
  ## 팔레트
  fill_pal <- c(
    "No vaccination"               = "grey50",
    "Disease blocking only"        = "#B5DBFF",
    "Disease & infection blocking" = "#F9C2A8"
  )
  col_pal <- c(
    "No vaccination"               = "black",
    "Disease blocking only"        = "#0072B2",
    "Disease & infection blocking" = "#E53935"
  )
  lty_pal <- c(
    "No vaccination"               = "solid",
    "Disease blocking only"        = "dotted",
    "Disease & infection blocking" = "23"
  )
  
  ## 첫 줄 병합용 공통 변수
  pre_df_legend <- pre_df %>% dplyr::mutate(legend_first = "No vaccination")
  obs_df_legend <- obs_df %>% dplyr::mutate(legend_first = "Reported symptomatic cases")
  
  ## 2) 플롯
  ggplot(df, aes(x = Week)) +
    # 공급기간 음영
    geom_rect(
      data = vacc_tbl_nat,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#C8E6C9", alpha = .1
    ) +
    
    ## ── [첫째 줄] "No vaccination"(제목) → 회색 박스 → × Reported symptomatic cases
    geom_ribbon(
      data = pre_df_legend,
      aes(x = Week, ymin = lo95, ymax = hi95, fill = legend_first),
      colour = NA, alpha = 0.30, inherit.aes = FALSE
    ) +
    geom_point(
      data = obs_df_legend,
      aes(x = Week, y = Observed, shape = legend_first, fill = legend_first),
      colour = "black", size = .9, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name   = "No vaccination",
      breaks = c("No vaccination", "Reported symptomatic cases"),
      values = c("No vaccination" = fill_pal[["No vaccination"]],
                 "Reported symptomatic cases" = NA),
      labels = c("", "Reported symptomatic cases"),
      guide  = guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        key_glyph = draw_key_rectline,
        override.aes = list(
          fill     = c(fill_pal[["No vaccination"]], NA),
          colour   = c(col_pal[["No vaccination"]], "black"),
          linetype = c(lty_pal[["No vaccination"]], "blank"),
          shape    = c(NA, 4),
          alpha    = c(.4, 1),
          size     = c(.6, 3)
        )
      )
    ) +
    scale_shape_manual(
      name   = "No vaccination",
      breaks = c("No vaccination", "Reported symptomatic cases"),
      values = c("No vaccination" = NA, "Reported symptomatic cases" = 4),
      guide  = "none"
    ) +
    
    ## 둘째 줄(보호효과) 만들기 위해 fill 스케일 리셋
    ggnewscale::new_scale("fill") +
    
    ## ── [둘째 줄] Vaccine protection (파랑/주황 박스)
    geom_ribbon(
      aes(ymin = post_weekly_low95, ymax = post_weekly_hi95, fill = VE_lab),
      alpha = 1, colour = NA
    ) +
    geom_line(
      aes(y = post_cases, linetype = VE_lab, colour = VE_lab),
      size = 1.3
    ) +
    scale_fill_manual(
      name   = "Vaccine protection",
      breaks = c("Disease blocking only", "Disease & infection blocking"),
      values = c(
        "Disease blocking only"        = fill_pal[["Disease blocking only"]],
        "Disease & infection blocking" = fill_pal[["Disease & infection blocking"]]
      ),
      guide  = guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        key_glyph = draw_key_rectline,
        title.position = "left",
        title.theme    = element_text(margin = margin(r = 10)),
        override.aes = list(
          fill     = unname(fill_pal[c("Disease blocking only",
                                       "Disease & infection blocking")]),
          colour   = unname(col_pal[c("Disease blocking only",
                                      "Disease & infection blocking")]),
          linetype = unname(lty_pal[c("Disease blocking only",
                                      "Disease & infection blocking")]),
          size     = .6, alpha = 1
        )
      )
    ) +
    
    # 라인 범례 숨김(박스로 충분)
    scale_colour_manual(values = col_pal, guide = "none") +
    scale_linetype_manual(values = lty_pal, guide = "none") +
    
    # pre-vacc 선
    geom_line(
      data = pre_df,
      aes(x = Week, y = pre_cases,
          linetype = "No vaccination", colour = "No vaccination"),
      size = .7, inherit.aes = FALSE
    ) +
    
    ## 축/패싯/테마
    scale_y_continuous(labels = label_k, expand = expansion(mult = c(0, .02))) +
    ylab("Predicted symptomatic cases") +
    facet_grid(Scenario ~ vaccine, scales = "free_y") +
    theme_pubclean(base_size = 10) +
    theme(
      legend.position      = "bottom",
      legend.box           = "vertical",
      legend.justification = "left",
      legend.box.just      = "left",
      legend.direction     = "horizontal",
      legend.title         = element_text(size = 12),
      legend.text          = element_text(size = 12),
      legend.spacing.x     = unit(6, "pt"),
      strip.text           = element_text(size = 11, face = "bold"),
      axis.text.y          = element_text(size = 13),
      axis.title.y         = element_text(size = 13),
      axis.title.x         = element_text(size = 13)
    ) +
    
    ## 주석
    geom_label(
      data = ann_lab,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.04, vjust = 1.55,
      size = 4, fill = NA, label.size = 0
    ) +
    annotate(
      "text", x = Inf, y = Inf, size = 4,
      label = "% Reduction in predicted symptomatic cases",
      hjust = 1.05, vjust = 1.1, fontface = "bold"
    ) +
    
    ## ⬇️ 맨 위로 한 번 더 점을 그려 가림 방지 (범례에는 미포함)
    geom_point(
      data = obs_df,
      aes(x = Week, y = Observed),
      shape = 4, colour = "black", size = .9, stroke = .6,
      inherit.aes = FALSE, show.legend = FALSE
    )
}


fig1_A <- 
make_ve_plot_contrast_both() +
  plot_annotation(tag_levels = list(c("A")))

ggsave(filename = "02_Outputs/2_1_Figures/fig1_A.jpg", fig1_A, width = 10, height = 8, dpi = 1200)
ggsave(
  filename = "02_Outputs/2_1_Figures/fig1_A.pdf",
  plot     = fig1_A,
  width    = 12,
  height   = 10
)
#-------------------------------------------------------------------------------

# ❶ region–VE–coverage 결합 (그대로)
combined_prepost_case_all_vevc <- imap_dfr(
  postsim_vc_ixchiq_model,                       # ─ region loop ─
  function(region_list, region_nm) {
    imap_dfr(region_list,                          # ─ VE loop ─
             function(ve_list, ve_tag) {
               map_dfr(ve_list, "summary_week_df", .id = "coverage") %>% 
                 mutate(region = region_nm, VE = ve_tag)
             })
  })

# ❷ setting·national 합치기 (coverage 포함 그대로)
combined_by_setting <- combined_prepost_case_all_vevc %>% 
  mutate(setting = dplyr::recode(region, !!!setting_key))

setting_summ <- combined_by_setting %>% 
  group_by(setting, VE, coverage, Scenario, Week) %>%     
  summarise(across(c(post_cases:hi95), sum), .groups = "drop")

combined_nat_setting <- combined_by_setting %>% 
  group_by(VE, coverage, Scenario, Week) %>%               
  summarise(across(c(post_cases:hi95), sum), .groups = "drop") %>% 
  mutate(setting = "National")

combined_all <- bind_rows(combined_nat_setting, setting_summ) %>% 
  mutate(setting = factor(setting, levels = c("Low","Moderate","High","National")),
         coverage = factor(coverage,                     
                           levels = c("cov10","cov50","cov90"),
                           labels = c("10 %","50 %","90 %")))

global_impact <- combined_all %>% 
  group_by(VE, coverage, Scenario, setting) %>%           
  summarise(total_post = sum(post_cases),
            total_pre  = sum(pre_cases), .groups = "drop") %>% 
  mutate(diff = total_pre - total_post,
         impact = diff / total_pre * 100)

annotation_setting <- global_impact %>% 
  mutate(Strategy = dplyr::recode(Scenario,
                           Scenario_1 = "Strategy 1",
                           Scenario_2 = "Strategy 2",
                           Scenario_3 = "Strategy 3",
                           Scenario_4 = "Strategy 4")) %>% 
  group_by(setting, VE, coverage) %>%                     
  summarise(ann = paste0(Strategy, ": ", round(impact,1), "%", collapse = "\n"),
            .groups = "drop")

