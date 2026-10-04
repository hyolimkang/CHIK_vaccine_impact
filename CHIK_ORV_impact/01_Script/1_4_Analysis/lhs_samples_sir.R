library(MASS)
library(truncnorm)
## 1. LHS sampling process for each parameter --------------------------------------

# check histograms visually
hist(posterior_ag$I0)  # left skewed, always should be > 0 
hist(posterior_ag$base_beta) # left skewed, always should be > 0 
hist(posterior_ag$rho) # normally distributed, always should be > 0  
hist(posterior_ag$gamma) # normally distributed, always should be > 0  
hist(posterior_ag$sigma) # normally distributed, always should be > 0  

# posterior list 
posterior_list <- list(
  "Ceará" = posterior_ce,
  "Bahia"= posterior_bh,
  "Paraíba" = posterior_pa,
  "Pernambuco" = posterior_pn,
  "Rio Grande do Norte" = posterior_rg,
  "Piauí" = posterior_pi,
  "Alagoas" = posterior_ag,
  "Tocantins" = posterior_tc,
  "Minas Gerais" = posterior_mg,
  "Sergipe" = posterior_se,
  "Goiás" = posterior_go
)

set.seed(123)
runs <- 1000
D <- randomLHS(runs, 12)

lhs_all_regions <- lapply(names(posterior_list), function(reg) {
  
  posterior <- posterior_list[[reg]]
  
  # fit distributions
  fit_I0    <- fitdistr(posterior$I0, "lognormal")
  fit_beta  <- fitdistr(posterior$base_beta, "lognormal")
  fit_rho   <- fitdistr(posterior$rho, "lognormal")
  fit_gamma <- fitdistr(posterior$gamma, "lognormal")
  fit_sigma <- fitdistr(posterior$sigma, "lognormal")
  
  # LHS sampling
  data.frame(
    I0        = qlnorm(D[,1], meanlog = fit_I0$estimate["meanlog"], 
                       sdlog   = fit_I0$estimate["sdlog"]),
    base_beta = qlnorm(D[,2], meanlog = fit_beta$estimate["meanlog"], 
                       sdlog   = fit_beta$estimate["sdlog"]),
    rho       = qlnorm(D[,3], meanlog = fit_rho$estimate["meanlog"], 
                       sdlog   = fit_rho$estimate["sdlog"]),
    gamma     = qlnorm(D[,4], meanlog = fit_gamma$estimate["meanlog"], 
                       sdlog   = fit_gamma$estimate["sdlog"]),
    sigma     = qlnorm(D[,5], meanlog = fit_sigma$estimate["meanlog"], 
                       sdlog   = fit_sigma$estimate["sdlog"]),
    
    foi       = foi_draws_list[[reg]],   
    
    ve        = qtruncnorm(D[,6], a = 0.967, b = 0.998,
                           mean = 0.989,
                           sd   = (0.998 - 0.967) / (2 * qnorm(0.975))),
    
    vc = qtruncnorm(D[,7], a = 0.4, b = 0.6,
                    mean = 0.5,
                    sd   = (0.6 - 0.4) / (2 * qnorm(0.975))),
    
    vc10      = qtruncnorm(D[,7], a = 0.05, b = 0.15, mean = 0.10,
                           sd = (0.15 - 0.05) / (2 * qnorm(0.975))),
    
    vc50      = qtruncnorm(D[,7], a = 0.40, b = 0.60, mean = 0.50,
                           sd = (0.60 - 0.40) / (2 * qnorm(0.975))),
    
    vc90      = qtruncnorm(D[,7], a = 0.85, b = 0.95, mean = 0.90,
                           sd = (0.95 - 0.85) / (2 * qnorm(0.975))),
    
    wd = qtruncnorm(D[,8], a = 0.1*0.9, b = 0.1*1.1,
                    mean = 0.1,
                    sd   = (0.1*1.1 - 0.1*0.9) / (2 * qnorm(0.975)))
  )
})

names(lhs_all_regions) <- names(posterior_list)

set.seed(123)
n_draws <- 5000
lhs_idx_ixchiq <- lapply(lhs_all_regions, function(x) sample(seq_len(nrow(x)), size = n_draws, replace = FALSE))
names(lhs_idx_ixchiq) <- names(lhs_all_regions)


### vimkunya
lhs_vimkunya <- lapply(names(posterior_list), function(reg) {
  
  posterior <- posterior_list[[reg]]  # posterior 추출
  
  # fit distributions
  fit_I0    <- fitdistr(posterior$I0, "lognormal")
  fit_beta  <- fitdistr(posterior$base_beta, "lognormal")
  fit_rho   <- fitdistr(posterior$rho, "lognormal")
  fit_gamma <- fitdistr(posterior$gamma, "lognormal")
  fit_sigma <- fitdistr(posterior$sigma, "lognormal")
  
  # LHS 
  data.frame(
    I0        = qlnorm(D[,1], meanlog = fit_I0$estimate["meanlog"], 
                       sdlog   = fit_I0$estimate["sdlog"]),
    base_beta = qlnorm(D[,2], meanlog = fit_beta$estimate["meanlog"], 
                       sdlog   = fit_beta$estimate["sdlog"]),
    rho       = qlnorm(D[,3], meanlog = fit_rho$estimate["meanlog"], 
                       sdlog   = fit_rho$estimate["sdlog"]),
    gamma     = qlnorm(D[,4], meanlog = fit_gamma$estimate["meanlog"], 
                       sdlog   = fit_gamma$estimate["sdlog"]),
    sigma     = qlnorm(D[,5], meanlog = fit_sigma$estimate["meanlog"], 
                       sdlog   = fit_sigma$estimate["sdlog"]),
    
    foi       = foi_draws_list[[reg]],   
    
    ve        = qtruncnorm(D[,6], a = 0.972, b = 0.983,
                           mean = 0.978,
                           sd   = (0.983 - 0.972) / (2 * qnorm(0.975))),
    
    vc = qtruncnorm(D[,7], a = 0.4, b = 0.6,
                    mean = 0.5,
                    sd   = (0.6 - 0.4) / (2 * qnorm(0.975))),
    
    vc10      = qtruncnorm(D[,7], a = 0.05, b = 0.15, mean = 0.10,
                           sd = (0.15 - 0.05) / (2 * qnorm(0.975))),
    
    vc50      = qtruncnorm(D[,7], a = 0.40, b = 0.60, mean = 0.50,
                           sd = (0.60 - 0.40) / (2 * qnorm(0.975))),
    
    vc90      = qtruncnorm(D[,7], a = 0.85, b = 0.95, mean = 0.90,
                           sd = (0.95 - 0.85) / (2 * qnorm(0.975))),
    
    wd = qtruncnorm(D[,8], a = 0.1*0.9, b = 0.1*1.1,
                    mean = 0.1,
                    sd   = (0.1*1.1 - 0.1*0.9) / (2 * qnorm(0.975)))
  )
})
names(lhs_vimkunya) <- names(posterior_list)


## combined lhs #----------------------------------------------------------------
# 1. combined LHS sample list 
set.seed(123)
n_draws <- 1000
lhs_combined <- lapply(names(posterior_list), function(reg) {
  
  posterior <- posterior_list[[reg]]
  
  # fit distributions
  fit_I0    <- fitdistr(posterior$I0, "lognormal")
  fit_beta  <- fitdistr(posterior$base_beta, "lognormal")
  fit_rho   <- fitdistr(posterior$rho, "lognormal")
  fit_gamma <- fitdistr(posterior$gamma, "lognormal")
  fit_sigma <- fitdistr(posterior$sigma, "lognormal")
  
  # LHS sampling
  data.frame(
    I0        = qlnorm(D[,1], meanlog = fit_I0$estimate["meanlog"], 
                       sdlog   = fit_I0$estimate["sdlog"]),
    base_beta = qlnorm(D[,2], meanlog = fit_beta$estimate["meanlog"], 
                       sdlog   = fit_beta$estimate["sdlog"]),
    rho       = qlnorm(D[,3], meanlog = fit_rho$estimate["meanlog"], 
                       sdlog   = fit_rho$estimate["sdlog"]),
    gamma     = qlnorm(D[,4], meanlog = fit_gamma$estimate["meanlog"], 
                       sdlog   = fit_gamma$estimate["sdlog"]),
    sigma     = qlnorm(D[,5], meanlog = fit_sigma$estimate["meanlog"], 
                       sdlog   = fit_sigma$estimate["sdlog"]),
    
    foi       = foi_draws_list[[reg]],
    
    # Ixchiq VE
    ve_ix     = qtruncnorm(D[,6], a = 0.967, b = 0.998,
                           mean = 0.989,
                           sd   = (0.998 - 0.967) / (2 * qnorm(0.975))),
    
    # Vimkunya VE
    ve_vimkun = qtruncnorm(D[,6], a = 0.972, b = 0.983,
                           mean = 0.978,
                           sd   = (0.983 - 0.972) / (2 * qnorm(0.975))),
    
    vc   = qtruncnorm(D[,7], a = 0.4, b = 0.6,
                      mean = 0.5,
                      sd   = (0.6 - 0.4) / (2 * qnorm(0.975))),
    
    vc10 = qtruncnorm(D[,7], a = 0.05, b = 0.15, mean = 0.10,
                      sd = (0.15 - 0.05) / (2 * qnorm(0.975))),
    
    vc50 = qtruncnorm(D[,7], a = 0.40, b = 0.60, mean = 0.50,
                      sd = (0.60 - 0.40) / (2 * qnorm(0.975))),
    
    vc90 = qtruncnorm(D[,7], a = 0.85, b = 0.95, mean = 0.90,
                      sd = (0.95 - 0.85) / (2 * qnorm(0.975))),
    
    wd   = qtruncnorm(D[,8], a = 0.1*0.9, b = 0.1*1.1,
                      mean = 0.1,
                      sd   = (0.1*1.1 - 0.1*0.9) / (2 * qnorm(0.975)))
  )
})

names(lhs_combined) <- names(posterior_list)

# 2. for reproducibility, create LHS sample ids 
set.seed(123)  # reproducibility
n_draws <- 1000

lhs_idx_list <- lapply(lhs_combined, function(df) {
  sample(seq_len(nrow(df)), size = n_draws, replace = FALSE)
})

names(lhs_idx_list) <- names(lhs_combined)

# 3. for reproducibility, create posterior sapmle ids 
set.seed(123)
posterior_idx_list <- lapply(posterior_list, function(posterior) {
  sample(seq_len(nrow(posterior$base_beta)), size = 1000, replace = TRUE)
})

names(posterior_idx_list) <- names(posterior_list)


save(
  posterior_idx_list, lhs_idx_list, lhs_combined, 
  file = "00_Data/0_2_Processed/lhs_orv.RData"
)



## debugging -----------------------------------------------------------------------------------------------------------------
# why negative values in pre-vacc calculation for lo values? 

# very upper function (prevacc-sim)
has_negative <- sapply(preui_ce, function(x) {
  if (is.data.frame(x)) {
    any(sapply(x, function(col) is.numeric(col) && any(col < 0, na.rm = TRUE)))
  } else if (is.numeric(x)) {
    any(x < 0, na.rm = TRUE)
  } else {
    FALSE
  }
}) # ok so far so good -- there is no negative values here in the raw simulation 

# prevacc summ 
neg_cols <- names(prevacc_ui_ce)[sapply(prevacc_ui_ce, function(x) is.numeric(x) && any(x < 0, na.rm = TRUE))]
# ok now fixed so far so good 

# next check if there is any neg value in post-vacc sim
has_negative <- sapply(postsim_ce_ui , function(x) {
  if (is.data.frame(x)) {
    any(sapply(x, function(col) is.numeric(col) && any(col < 0, na.rm = TRUE)))
  } else if (is.numeric(x)) {
    any(x < 0, na.rm = TRUE)
  } else {
    FALSE
  }
}) # ok so far so good 

# now check if postsim summarisation is okay
neg_cols <- names(postsim_all_ce_ui)[sapply(postsim_all_ce_ui, function(x) is.numeric(x) && any(x < 0, na.rm = TRUE))]
# ok non negative values 


# now check if df_nnv
has_negative <- sapply(nnv_ce , function(x) {
  if (is.data.frame(x)) {
    any(sapply(x, function(col) is.numeric(col) && any(col < 0, na.rm = TRUE)))
  } else if (is.numeric(x)) {
    any(x < 0, na.rm = TRUE)
  } else {
    FALSE
  }
}) 

neg_cols <- names(postsim_all_bh_ui$summary_list_df)[sapply(postsim_all_bh_ui$summary_list_df , function(x) is.numeric(x) && any(x < 0, na.rm = TRUE))]
neg_cols
# ok so far so good    

neg_cols <- names(df_nnv_bh)[sapply(df_nnv_bh , function(x) is.numeric(x) && any(x < 0, na.rm = TRUE))]
neg_cols
# now everything cleared 


neg_cols <- names(postsim_vc_ixchiq_model$Ceará$VE0$cov90$final_summ[[1]])[sapply(postsim_vc_ixchiq_model$Ceará$VE0$cov50$final_summ[[1]], function(x) is.numeric(x) && any(x < 0, na.rm = TRUE))]
neg_cols


## check again
# test for convergence using nnv

vim_nnv_check <- df_nnv_bh %>%
  group_by(scenario) %>%            # VE·VC·scenario 기준으로 그룹화
  summarise(
    across(c(tot_vacc, pre_inf:diff_daly_hi), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ## NNV
    nnv_inf      = tot_vacc / diff_inf,
    nnv_inf_lo   = tot_vacc / diff_inf_hi,
    nnv_inf_hi   = tot_vacc / diff_inf_lo,
    nnv          = tot_vacc / diff,
    nnv_lo       = tot_vacc / diff_hi,
    nnv_hi       = tot_vacc / diff_low,
    nnv_fatal    = tot_vacc / diff_fatal,
    nnv_fatal_lo = tot_vacc / diff_fatal_hi,
    nnv_fatal_hi = tot_vacc / diff_fatal_low,
    nnv_daly     = tot_vacc / diff_daly,
    nnv_daly_lo  = tot_vacc / diff_daly_hi,
    nnv_daly_hi  = tot_vacc / diff_daly_low,
    
    ## per 100k 예방 효과
    per1M_inf     = diff_inf       / tot_vacc * 1e5,
    per1M_inf_lo  = diff_inf_lo    / tot_vacc * 1e5,
    per1M_inf_hi  = diff_inf_hi    / tot_vacc * 1e5,
    
    per1M_diff     = diff       / tot_vacc * 1e5,
    per1M_diff_lo  = diff_low    / tot_vacc * 1e5,
    per1M_diff_hi  = diff_hi    / tot_vacc * 1e5,
    
    per1M_fatal     = diff_fatal       / tot_vacc * 1e5,
    per1M_fatal_lo  = diff_fatal_low    / tot_vacc * 1e5,
    per1M_fatal_hi  = diff_fatal_hi    / tot_vacc * 1e5,
    
    per1M_daly     = diff_daly       / tot_vacc * 1e5,
    per1M_daly_lo  = diff_daly_low    / tot_vacc * 1e5,
    per1M_daly_hi  = diff_daly_hi    / tot_vacc * 1e5
  )

ggplot(vim_nnv_check) +
  geom_bar(aes(x = scenario, y = nnv), stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(x = scenario, ymin = nnv_lo, ymax = nnv_hi),
                width = 0.2)+
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))+
  scale_y_log10()


ggplot(vim_nnv_check) +
  geom_bar(aes(x = scenario, y = nnv_fatal), stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(x = scenario, ymin = nnv_fatal_lo, ymax = nnv_fatal_hi),
                width = 0.2)+
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))+
  scale_y_log10()


ggplot(vim_nnv_check) +
  geom_bar(aes(x = scenario, y = nnv_daly), stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(x = scenario, ymin = nnv_daly_lo, ymax = nnv_daly_hi),
                width = 0.2)+
  scale_y_log10()




ve_waning_curve <- function(T = 52, ve_week3, ve_week26) {
  slope <- (ve_week26 - ve_week3) / (26 - 3)
  
  VE_time <- sapply(1:T, function(t) {
    if (t <= 3) {
      ve_week3
    } else if (t <= 26) {
      ve_week3 + slope * (t - 3)
    } else {
      ve_week26 + slope * (t - 26)
    }
  })
  
  return(VE_time)
}

ve_waning_curve <- function(T = 520 * 1.5, ve_week3, ve_week26) {
  
  lambda <- log(ve_week3 / ve_week26) / (26 - 3)
  
  VE_time <- numeric(T)
  
  for (t in 1:T) {
    if (t <= 3) {
      VE_time[t] <- ve_week3
    } else {
      VE_time[t] <- ve_week3 * exp(-lambda * (t - 3))
    }
  }
  
  return(VE_time)
}
VE_time <- ve_waning_curve(T = 520 * 1.5,
                           ve_week3  = 0.978,
                           ve_week26 = 0.841)

VE_time <- ve_waning_curve(T = 52,
                           ve_week3  = 0.978,
                           ve_week26 = 0.948)

# debug ve x vc function with using one region
regions <- list(
  "Ceará" = list(observed = observed_ce, N = N_ceara$Ceará,
                 prevacc_ui = prevacc_ui_ce, posterior = posterior_ce,
                 lhs_sample = lhs_all_regions$ce)
)

make_prevacc <- function(posterior, observed, N, region, lhs_sample, sero_vec = NULL) {

  pre_ui <- simulate_pre_ui_age(
    posterior          = posterior,
    bra_foi_state_summ = bra_foi_state_summ,
    age_groups         = age_groups,
    N                  = N,
    region             = region,
    observed           = observed,
    lhs_sample         = lhs_sample,
    sero_vec           = sero_vec
  )
  
  presum <- summarise_presim_ui(
    sim_result       = pre_ui,
    observed         = observed,
    age_gr_levels    = age_gr_levels,
    lhs_sample_young = lhs_sample_young,
    lhs_old          = lhs_old,
    le_sample        = le_sample,
    hosp             = hosp,
    fatal            = fatal,
    nh_fatal         = nh_fatal,
    region           = region
  )
  
  list(
    prevacc_ui        = presum$summary_cases_pre,
    prevacc_ui_all    = presum$summary_cases_pre_all,
    pre_summary_age   = presum$summary_cases_pre_age
  )
}

prevacc_sets <- imap(regions, function(reg_args, region_name) {
  make_prevacc(
    posterior = reg_args$posterior,
    observed  = reg_args$observed,
    N         = reg_args$N,
    region    = region_name,
    lhs_sample = reg_args$lhs_sample
  )
})

ve_list <- postsim_ce_ui[[4]]$VE_time_list
weeks <- 1:52

summary_df <- imap_dfr(ve_list, function(df, idx) {
  tibble(
    Week   = weeks,
    Median = apply(df, 2, median, na.rm = TRUE),
    Low95  = apply(df, 2, quantile, probs = 0.025, na.rm = TRUE),
    High95 = apply(df, 2, quantile, probs = 0.975, na.rm = TRUE),
    Scenario = paste0("VE", idx)
  )
})

# ggplot
ggplot(summary_df, aes(x = Week, y = Median, color = Scenario, fill = Scenario)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Low95, ymax = High95), alpha = 0.2, color = NA) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")