# cleaning pop
N_ce <- N_ceara
N_bh <- N_bahia
N_ce <- N_ce$Ceará
N_bh <- N_bh$Bahia
N_mg <- N_mg$`Minas Gerais`
N_pn <- N_pemam
N_pn <- N_pn$Pernambuco
N_pa <- N_pa$Paraíba
N_rg <- N_rg$`Rio Grande do Norte`
N_pi <- N_pi$Piauí
N_ag <- N_ag$Alagoas
N_tc <- N_tc$Tocantins
N_se <- N_se$Sergipe
N_go <- N_go$Goiás

region_abbr <- c(
  "Ceará" = "ce",
  "Pernambuco" = "pn",
  "Minas Gerais" = "mg",
  "Bahia" = "bh",
  "Paraíba" = "pa",
  "Rio Grande do Norte" = "rg",
  "Piauí" = "pi",
  "Alagoas" = "ag",
  "Tocantins" = "tc",
  "Sergipe"   = "se",
  "Goiás"     = "go"
)

region_abbr <- c("Ceará" = "ce")

label_lo_hi <- function(x) {
  x <- unname(x)
  stopifnot(length(x) == 2)
  ord <- order(x)
  setNames(x[ord], c("lo","hi"))
}

safe_pop <- function(N_obj, region_name) {
  if (is.list(N_obj)) {
    if (!is.null(N_obj[[region_name]])) return(N_obj[[region_name]])
    # 첫 원소 사용 (경고)
    warning("Population object is list; region not found; using first.")
    return(N_obj[[1]])
  }
  N_obj
}
run_owsa_for_region <- function(region_name,
                                bra_foi_state_summ,
                                age_groups,
                                age_gr_levels,
                                lhs_sample_young,
                                lhs_old,
                                le_sample,
                                target_age_list,
                                hosp,
                                fatal,
                                nh_fatal) {
  ## ---- region lookup ----
  region_abbr <- c(
    "Ceará" = "ce", "Pernambuco" = "pn", "Minas Gerais" = "mg", "Bahia" = "bh",
    "Paraíba" = "pa", "Rio Grande do Norte" = "rg", "Piauí" = "pi", "Alagoas" = "ag",
    "Tocantins" = "tc", "Sergipe" = "se", "Goiás" = "go"
  )
  abbr <- region_abbr[[region_name]]
  if (is.null(abbr)) stop("Unknown region_name: ", region_name)
  
  # dynamic objects by naming convention
  observed  <- get(paste0("observed_", abbr), envir = .GlobalEnv)
  N_region  <- get(paste0("N_", abbr),         envir = .GlobalEnv)
  posterior <- get(paste0("posterior_", abbr), envir = .GlobalEnv)
  
  # ensure population scalar
  N_use <- safe_pop(N_region, region_name)
  
  ## ---- assemble region params ----
  params_list <- make_params(
    region             = region_name,
    bra_foi_state_summ = bra_foi_state_summ,
    posterior          = posterior,
    lhs_sample_young   = lhs_sample_young,
    hosp               = hosp,
    fatal              = fatal,
    nh_fatal           = nh_fatal
  )
  baseline_params    <- params_list$baseline_params
  sensitivity_ranges <- params_list$sensitivity_ranges
  
  ## lo/hi labeling across all params
  sensitivity_ranges_labeled <- lapply(sensitivity_ranges, function(x) {
    # skip non-length-2 params
    if (length(x) != 2) {
      warning("Param ", deparse(substitute(x)), " not length 2; skipping.")
      return(x)
    }
    label_lo_hi(x)
  })
  
  ## ---- OWSA loops ----
  sensitivity_results <- purrr::imap(sensitivity_ranges_labeled, function(param_values, param) {
    purrr::imap_dfr(param_values, function(val, scen) {
      res <- run_simulation_with_params(
        param_overrides      = setNames(list(val), param),
        baseline_params      = baseline_params,
        posterior            = posterior,
        bra_foi_state_summ   = bra_foi_state_summ,
        age_groups           = age_groups,
        N                    = N_use,
        region               = region_name,
        observed             = observed,
        age_gr_levels        = age_gr_levels,
        lhs_sample_young     = lhs_sample_young,
        lhs_old              = lhs_old,
        le_sample            = le_sample,
        target_age_list      = target_age_list
      )
      tibble::tibble(
        Parameter                 = param,
        Scenario                  = scen,  # "lo"/"hi"
        Value                     = val,
        pre_tot_cases_agegr1_owsa = res$pre_tot_cases_agegr1_owsa,
        pre_tot_cases_agegr2_owsa = res$pre_tot_cases_agegr2_owsa,
        pre_tot_cases_agegr3_owsa = res$pre_tot_cases_agegr3_owsa,
        pre_tot_cases_agegr4_owsa = res$pre_tot_cases_agegr4_owsa,
        pre_tot_cases_agegr1_mid  = res$pre_tot_cases_agegr1_mid,
        pre_tot_cases_agegr2_mid  = res$pre_tot_cases_agegr2_mid,
        pre_tot_cases_agegr3_mid  = res$pre_tot_cases_agegr3_mid,
        pre_tot_cases_agegr4_mid  = res$pre_tot_cases_agegr4_mid,
        post_vacc_age_totals      = list(res$post_vacc_age_by_week),     # rename for clarity
        mid_post_vacc_age_totals  = list(res$mid_post_vacc_age_by_week)
      )
    })
  })
  
  ## ---- extract OWSA (lo/hi split) ----
  # Returns nested list param -> (lo list, hi list)
  extract_owsa <- purrr::imap(sensitivity_results, function(df, param) {
    stopifnot(nrow(df) == 2)  # expecting lo, hi
    # i=1 => lo, i=2=>hi (because label_lo_hi sorted)
    out <- vector("list", 2)
    names(out) <- df$Scenario  # should be c("lo","hi")
    
    for (i in seq_len(2)) {
      pre_owsa_vec <- as.numeric(df[i, paste0("pre_tot_cases_agegr", 1:4, "_owsa")])
      pre_mid_vec  <- as.numeric(df[i, paste0("pre_tot_cases_agegr", 1:4, "_mid")])
      
      # post_vacc_age_totals is list of scenario tibbles
      scen_list_owsa <- df$post_vacc_age_totals[[i]]
      scen_list_mid  <- df$mid_post_vacc_age_totals[[i]]
      
      # assume 1:4 map to agegrp1~4 vaccination scenarios
      post_vals <- vapply(1:4, function(j) {
        scen_list_owsa[[j]][[paste0("tot_post_cases_agegr", j)]]
      }, numeric(1))
      
      mid_vals <- vapply(1:4, function(j) {
        scen_list_mid[[j]][[paste0("tot_post_cases_agegr", j)]]
      }, numeric(1))
      
      out[[i]] <- list(
        parameter = df$Parameter[i],
        value     = df$Value[i],
        pre_owsa  = setNames(pre_owsa_vec, paste0("agegrp", 1:4)),
        pre_mid   = setNames(pre_mid_vec,  paste0("agegrp", 1:4)),
        post_vacc = setNames(post_vals,    paste0("agegrp", 1:4)),
        mid_vacc  = setNames(mid_vals,     paste0("agegrp", 1:4))
      )
    }
    out
  })
  
  ## ---- tidy tornado-style summary ----
  tornado_data <- purrr::map_dfr(names(extract_owsa), function(param) {
    ex <- extract_owsa[[param]]
    lo <- ex$lo
    hi <- ex$hi
    
    agegrps <- names(lo$pre_owsa)  # "agegrp1"..."agegrp4"
    
    tibble::tibble(
      Region        = region_name,
      Parameter     = param,
      AgeGroup      = agegrps,
      Value_lo      = lo$value,
      Value_hi      = hi$value,
      pre_owsa_lo   = as.numeric(lo$pre_owsa),
      pre_owsa_hi   = as.numeric(hi$pre_owsa),
      pre_mid       = as.numeric(lo$pre_mid),  # same for lo/hi; mid = baseline
      post_owsa_lo  = as.numeric(lo$post_vacc),
      post_owsa_hi  = as.numeric(hi$post_vacc),
      mid_vacc      = as.numeric(lo$mid_vacc)  # baseline vacc
    ) %>%
      dplyr::mutate(
        mid_diff   = pre_mid     - mid_vacc,
        lo_diff    = pre_owsa_lo - post_owsa_lo,
        hi_diff    = pre_owsa_hi - post_owsa_hi,
        mid_impact = ifelse(pre_mid     > 0, mid_diff / pre_mid     * 100, NA_real_),
        lo_impact  = ifelse(pre_owsa_lo > 0, lo_diff  / pre_owsa_lo * 100, NA_real_),
        hi_impact  = ifelse(pre_owsa_hi > 0, hi_diff  / pre_owsa_hi * 100, NA_real_)
      )
  })
  
  ## Optionally return raw as attr or list
  attr(tornado_data, "sensitivity_results") <- sensitivity_results
  tornado_data
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


regions <- c("Ceará", "Pernambuco",  "Minas Gerais","Bahia","Paraíba",            
             "Rio Grande do Norte", "Piauí", "Alagoas", "Tocantins", "Sergipe", "Goiás")


regions <- names(region_abbr)

tornado_results_all <- bind_rows(lapply(regions, run_owsa_for_region))

save("tornado_results_all", file = "00_Data/0_2_Processed/tornado_results_all.RData")

load("00_Data/0_2_Processed/tornado_results_all.RData")

tornado_bar_pct <- td_all %>%
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
    Bar_start = Mid_Impact_pct,  
    Bar_end   = Value           
  ) %>%
  group_by(Region, Parameter, AgeGroup) %>%
  mutate(TotalImpact = sum(abs(Bar_end - Bar_start), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Parameter = fct_reorder(Parameter, TotalImpact, .desc = TRUE)
  )

tornado_bar_pct <- tornado_bar_pct %>%
  mutate(Parameter = dplyr::recode(Parameter,
                            "foi"     = "FOI",
                            "VE_inf"  = "VE (disease & infection)",
                            "VE_block" = "VE (disease)",
                            "weekly_delivery_speed"  = "Weekly delivery rate",
                            "total_coverage"     = "Coverage",
                            "delay"   = "Delay",
                            "time_until_immunity" = "Time to immunity"
                           
  ))

tornado_bar_pct$Parameter <- factor(tornado_bar_pct$Parameter,
                                        levels = c("FOI", "Coverage", "Weekly delivery rate",  "Delay",  "Time to immunity",
                                                   "VE (disease & infection)", "VE (disease)")
)

p <- 
ggplot(tornado_bar_pct, aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction), linewidth = 5, alpha = 0.7) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5) +
  facet_grid(Region ~ AgeGroup, scales = "free_y", switch = "y",
             labeller = labeller(
               AgeGroup = as_labeller(c(
                 "agegrp1" = "Vaccination strategy 1 (1-11y)",
                 "agegrp2" = "Vaccination strategy 2 (12-17y)",
                 "agegrp3" = "Vaccination strategy 3 (18-59y)",
                 "agegrp4" = "Vaccination strategy 4 (60+y)"
               ))
             )) +
  scale_color_manual(values = c("Lower" = "#E41A1C", "Upper" = "#377EB8")) +
  labs(
    x = "% Symptomatic cases averted",
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

ggplot(tornado_bar_pct %>% filter(AgeGroup == "agegrp4"),
       aes(y = fct_rev(Parameter))) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction), linewidth = 5, alpha = 0.7) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5)  + facet_grid(Region ~ .)



ggsave(filename = "02_Outputs/2_1_Figures/figs9.jpg", p, width = 10, height = 12, dpi = 1200)

# national level---------------------------------------------------------------
td_nat <- td_all %>%
  group_by(AgeGroup, Parameter) %>%
  summarise(
    pre_mid      = sum(pre_mid,      na.rm = TRUE),
    mid_vacc     = sum(mid_vacc,     na.rm = TRUE),
    pre_owsa_lo  = sum(pre_owsa_lo,  na.rm = TRUE),
    post_owsa_lo = sum(post_owsa_lo, na.rm = TRUE),
    pre_owsa_hi  = sum(pre_owsa_hi,  na.rm = TRUE),
    post_owsa_hi = sum(post_owsa_hi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_mid = ifelse(pre_mid     > 0, 100 * (pre_mid     - mid_vacc)     / pre_mid,     NA_real_),
    pct_lo  = ifelse(pre_owsa_lo > 0, 100 * (pre_owsa_lo - post_owsa_lo) / pre_owsa_lo, NA_real_),
    pct_hi  = ifelse(pre_owsa_hi > 0, 100 * (pre_owsa_hi - post_owsa_hi) / pre_owsa_hi, NA_real_)
  )

tornado_bar_pct_nat <- td_nat %>%
  select(AgeGroup, Parameter, pct_mid, pct_lo, pct_hi) %>%
  pivot_longer(
    cols      = c(pct_lo, pct_hi),
    names_to  = "Direction",
    values_to = "Bar_end"
  ) %>%
  mutate(
    Direction = dplyr::recode(Direction, pct_lo = "Lower", pct_hi = "Upper"),
    Bar_start = pct_mid
  )


param_order <- td_nat %>%
  mutate(
    range_width = pmax(abs(pct_hi - pct_mid), abs(pct_lo - pct_mid))
  ) %>%
  group_by(AgeGroup) %>%
  arrange(desc(range_width), .by_group = TRUE) %>%
  mutate(ord = row_number()) %>%
  ungroup() %>%
  select(AgeGroup, Parameter, ord)

param_labels <- c(
  "total_coverage" = "Vaccine coverage",
  "foi"   = "FoI",
  "delay" = "Delay in deployment",
  "weekly_delivery_speed" = "Weekly delivery speed",
  "VE_inf" = "VE (infection blocking)",
  "VE_block" = "VE (disease blocking)",
  "time_until_immunity" = "Time until immunity acquisition"
)
  
tornado_bar_pct_nat <- tornado_bar_pct_nat %>%
  left_join(param_order, by = c("AgeGroup", "Parameter")) %>%
  group_by(AgeGroup) %>%
  mutate(Parameter = fct_reorder(Parameter, ord, .desc = TRUE)) %>%
  #mutate(Parameter = recode(Parameter, !!!param_labels)) %>% 
  ungroup()


strategy_lab <- c(
  "agegrp1" = "Strategy 1 (1–11y)",
  "agegrp2" = "Strategy 2 (12–17y)",
  "agegrp3" = "Strategy 3 (18–59y)",
  "agegrp4" = "Strategy 4 (60+y)"
)

lab_fn <- function(x) {
  out <- as.character(x)
  idx <- out %in% names(param_labels)
  out[idx] <- param_labels[out[idx]]
  out
}

ggplot(tornado_bar_pct_nat,
       aes(y = Parameter)) +
  geom_segment(aes(x = Bar_start, xend = Bar_end, yend = fct_rev(Parameter),
                   color = Direction),
               linewidth = 5, alpha = 0.7) +
  geom_point(aes(x = Bar_start), shape = 21, fill = "white",
             color = "black", size = 2.5) +
  facet_grid(. ~ AgeGroup, switch = "y",
             labeller = labeller(AgeGroup = as_labeller(strategy_lab))) +
  scale_color_manual(values = c("Lower" = "#E41A1C", "Upper" = "#377EB8")) +
  scale_y_discrete(labels = lab_fn) +
  labs(
    x = "% symptomatic cases averted",
    y = "Parameter",
    color = "Bounds"
  ) +
  theme_bw() +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.grid.major.y= element_blank(),
    legend.position   = "bottom",
    strip.text.y      = element_text(size = 6)
  )



ggsave(filename = "02_Outputs/2_1_Figures/owsa_all_v2.jpg", owsa_all, width = 10, height = 9, dpi = 1200)
write_xlsx(tornado_bar_pct_nat, path = "02_Outputs/2_2_Tables/tornado_bar_pct_nat_v2.xlsx")

