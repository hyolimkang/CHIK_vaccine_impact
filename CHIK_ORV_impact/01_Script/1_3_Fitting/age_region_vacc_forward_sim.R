posterior_prevacc <- rstan::extract(fit_prevacc)

sirv_sim <- function(
    # model dimensions
  T, A,
  N,                # dimension: [A]
  r,                # length A
  
  # params from the 'no-vacc' posterior
  base_beta,        # length T
  I0_draw,          # vector [A]
  rho,              # detection probability
  gamma,            # recovery rate
  
  # vaccine scenario
  vaccination_rate,
  delay,
  VE_block,
  vaccine_target_age,   # length A (0/1)
  vaccine_duration
){
  # 1) Initialize compartments: S, I, R, V, IV, RV
  S <- matrix(0, nrow = A, ncol = T)
  I <- matrix(0, nrow = A, ncol = T)
  R_ <- matrix(0, nrow = A, ncol = T)
  V <- matrix(0, nrow = A, ncol = T)
  IV <- matrix(0, nrow = A, ncol = T)
  RV <- matrix(0, nrow = A, ncol = T)
  age_stratified_cases <- matrix(0, nrow = A, ncol = T)
  
  # initial conditions at t=1
  for(a in seq_len(A)) {
    S[a, 1] <- N[a] - I0_draw[a]
    I[a, 1] <- I0_draw[a]
    R_[a, 1] <- 0
    V[a, 1]  <- 0
    IV[a, 1] <- 0
    RV[a, 1] <- 0
    age_stratified_cases[a, 1] <- rho * (I[a, 1] + IV[a, 1])
  }
  
  vaccine_end <- delay + vaccine_duration
  
  # 2) Forward simulation from t=2..T
  for(t_ in 2:T) {
    # Compute total infected in (t_-1)
    infected_prev <- sum(I[, t_-1] + IV[, t_-1])
    region_pop <- sum(N)
    
    # Force of infection
    phi <- base_beta[t_-1] * (infected_prev / region_pop)
    
    # Update compartments by age
    for(a in seq_len(A)) {
      # Decide vaccination_effect for age 'a'
      vacc_eff <- ifelse((t_ >= delay) && (t_ <= vaccine_end) && (vaccine_target_age[a] == 1),
                         vaccination_rate, 0)
      
      if(a == 1) {
        # First age group (no inflow from a-1)
        S[a, t_] <- S[a, t_-1] - phi * S[a, t_-1] -
          vacc_eff * S[a, t_-1] -
          r[a] * S[a, t_-1]
        I[a, t_] <- I[a, t_-1] + phi * S[a, t_-1] -
          gamma * I[a, t_-1] -
          r[a] * I[a, t_-1]
        R_[a, t_] <- R_[a, t_-1] + gamma * I[a, t_-1] -
          r[a] * R_[a, t_-1]
        V[a, t_] <- V[a, t_-1] + vacc_eff * S[a, t_-1] -
          (1 - VE_block) * phi * V[a, t_-1] -
          r[a] * V[a, t_-1]
        IV[a, t_] <- IV[a, t_-1] + (1 - VE_block) * phi * V[a, t_-1] -
          gamma * IV[a, t_-1] -
          r[a] * IV[a, t_-1]
        RV[a, t_] <- RV[a, t_-1] + gamma * IV[a, t_-1] -
          r[a] * RV[a, t_-1]
        
      } else {
        # Age group a >= 2 (inflow from a-1)
        S[a, t_] <- S[a, t_-1] - phi * S[a, t_-1] -
          vacc_eff * S[a, t_-1] -
          r[a] * S[a, t_-1] +
          r[a-1] * S[a-1, t_-1]
        I[a, t_] <- I[a, t_-1] + phi * S[a, t_-1] -
          gamma * I[a, t_-1] -
          r[a] * I[a, t_-1] +
          r[a-1] * I[a-1, t_-1]
        R_[a, t_] <- R_[a, t_-1] + gamma * I[a, t_-1] -
          r[a] * R_[a, t_-1] +
          r[a-1] * R_[a-1, t_-1]
        V[a, t_] <- V[a, t_-1] + vacc_eff * S[a, t_-1] -
          (1 - VE_block) * phi * V[a, t_-1] -
          r[a] * V[a, t_-1] +
          r[a-1] * V[a-1, t_-1]
        IV[a, t_] <- IV[a, t_-1] + (1 - VE_block) * phi * V[a, t_-1] -
          gamma * IV[a, t_-1] -
          r[a] * IV[a, t_-1] +
          r[a-1] * IV[a-1, t_-1]
        RV[a, t_] <- RV[a, t_-1] + gamma * IV[a, t_-1] -
          r[a] * RV[a, t_-1] +
          r[a-1] * RV[a-1, t_-1]
      }
      
      # Optionally clamp compartments at 0
      S[a, t_]  <- pmax(0, S[a, t_])
      I[a, t_]  <- pmax(0, I[a, t_])
      R_[a, t_] <- pmax(0, R_[a, t_])
      V[a, t_]  <- pmax(0, V[a, t_])
      IV[a, t_] <- pmax(0, IV[a, t_])
      RV[a, t_] <- pmax(0, RV[a, t_])
    } # end a
    for(a in seq_len(A)) {
      age_stratified_cases[a, t_] <- rho * (I[a, t_] + IV[a, t_])
    }
  }     # end t_
  
  # Return compartments
  list(
    S=S, I=I, R=R_, 
    V=V, IV=IV, RV=RV,
    age_stratified_cases = age_stratified_cases
  )
}

base_beta       <- apply(posterior_prevacc$base_beta, 2, median)       # length T
I0              <- apply(posterior_prevacc$I0, 2, median)           
gamma           <- median(posterior_prevacc$gamma)
rho             <- median(posterior_prevacc$rho)

sim_result <- sirv_sim(
  T = 52, 
  A = 18,
  N = N[,2],
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  rho = rho,
  gamma = gamma,
  
  # Vaccination scenario
  vaccination_rate = 1 - (0.6)^(1/6),  # Weekly vaccination rate
  delay = 3,                   # Vaccination starts from Week 1
  VE_block = 0.75,             # Vaccine efficacy
  vaccine_target_age = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0),  # Targeting specific age groups
  vaccine_duration = 6         # Vaccination duration in weeks
)

sim_result_novacc <- sirv_sim(
  T = 52, A = 18,
  N = N[,2],
  r = rep(0, 18),
  
  base_beta = base_beta,
  I0 = I0,
  rho = rho,
  gamma = gamma,
  
  # Vaccination scenario
  vaccination_rate = 0,
  delay = 1,
  VE_block = 0.75,
  vaccine_target_age = rep(0, 18),
  vaccine_duration = 0
)

summary_cases <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Median")

# Rename columns for clarity
colnames(summary_cases) <- c("AgeGroup", "Week", "Median")

# Convert columns to numeric for correct ordering
summary_cases <- summary_cases %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )
df_post <- summary_cases

summary_cases_pre <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Median")
summary_R_pre <- as.data.frame.table(sim_result_novacc$R, responseName = "R")
summary_RV_pre <- as.data.frame.table(sim_result_novacc$RV, responseName = "RV")

# Rename columns for clarity
colnames(summary_cases_pre) <- c("AgeGroup", "Week","Median")

# Convert columns to numeric for correct ordering
summary_cases_pre <- summary_cases_pre %>%
  mutate(
    AgeGroup = as.numeric(AgeGroup),
    Week = as.numeric(Week)
  )

summary_cases_pre$Scenario <- "Pre-Vaccination"
summary_cases$Scenario <- "Post-Vaccination"

# Combine both datasets
all_cases <- rbind(summary_cases_pre, summary_cases)

ggplot()+
  geom_line(data = summary_cases, aes(x = Week, y = Median, color = "blue"))+
  geom_line(data = summary_cases_pre, aes(x = Week, y = Median, color = "red"))+
  facet_wrap(~AgeGroup)
#-------------------------------------------------------------------------------
num_iters <- nrow(posterior_prevacc$base_beta)  # Number of iterations
all_results <- list() 
all_results_pre <- list()

# post-vacc
for (i_iter in seq_len(num_iters)) {
  
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Extract parameters for the current iteration
  base_beta <- posterior_prevacc$base_beta[i_iter, ]
  I0        <- posterior_prevacc$I0[i_iter, ]
  gamma     <- posterior_prevacc$gamma[i_iter]
  rho       <- posterior_prevacc$rho[i_iter]
  
  # Simulate for Post-Vaccination
  sim_result <- sirv_sim(
    T = 52, 
    A = 18,
    N = N[,1],
    r = rep(0, 18),
    base_beta = base_beta,
    I0_draw = I0,
    rho = rho,
    gamma = gamma,
    vaccination_rate = 1 - (0.6)^(1/6),
    delay = 3,
    VE_block = 0.75,
    vaccine_target_age = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
    vaccine_duration = 6
  )
  
  # Store results in a combined data frame
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
  sim_df <- sim_df %>%
    mutate(
      Iteration = i_iter,
      Scenario = "Post-Vaccination",
      AgeGroup = as.numeric(Var1),
      Week = as.numeric(Var2)
    )
  
  all_results[[length(all_results) + 1]] <- sim_df
}

# pre-vacc
for (i_iter in seq_len(num_iters)) {
  
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Extract parameters for the current iteration
  base_beta <- posterior_prevacc$base_beta[i_iter, ]
  I0        <- posterior_prevacc$I0[i_iter, ]
  gamma     <- posterior_prevacc$gamma[i_iter]
  rho       <- posterior_prevacc$rho[i_iter]
  
  # Simulate for Post-Vaccination
  sim_result_novacc <- sirv_sim(
    T = 52, A = 18,
    N = N[,1],
    r = rep(0, 18),
    
    base_beta = base_beta,
    I0 = I0,
    rho = rho,
    gamma = gamma,
    
    # Vaccination scenario
    vaccination_rate = 0,
    delay = 52,
    VE_block = 0.75,
    vaccine_target_age = rep(0, 18),
    vaccine_duration = 0
  )
  
  # Store results in a combined data frame
  sim_df <- as.data.frame.table(sim_result_novacc$age_stratified_cases, responseName = "Cases")
  sim_df <- sim_df %>%
    mutate(
      Iteration = i_iter,
      Scenario = "Pre-Vaccination",
      AgeGroup = as.numeric(Var1),
      Week = as.numeric(Var2)
    )
  
  all_results_pre[[length(all_results_pre) + 1]] <- sim_df
}

combined_results <- do.call(rbind, all_results)
combined_results_pre <- do.call(rbind, all_results_pre)

# Summarize results: median and 95% credible intervals
summary_results <- combined_results %>%
  group_by(AgeGroup, Week, Scenario) %>%
  summarise(
    Median = median(Cases),
    Lower = quantile(Cases, probs = 0.025),
    Upper = quantile(Cases, probs = 0.975),
    .groups = "drop"
  )

summary_results_pre <- combined_results_pre %>%
  group_by(AgeGroup, Week, Scenario) %>%
  summarise(
    Median = median(Cases),
    Lower = quantile(Cases, probs = 0.025),
    Upper = quantile(Cases, probs = 0.975),
    .groups = "drop"
  )

combined_sim <- rbind(summary_results, summary_results_pre)

ggplot(summary_results, aes(x = Week, y = Median, color = Scenario)) +
  geom_line() +
  facet_wrap(~ AgeGroup)

ggplot(summary_results_pre, aes(x = Week, y = Median, color = Scenario)) +
  geom_line() +
  facet_wrap(~ AgeGroup)

ggplot()+
  geom_line(data = summary_results, aes(x = Week, y = Median, color = Scenario))+
  facet_wrap(~ AgeGroup, ncol = 3)+
  geom_line(data = summary_results_pre, aes(x = Week, y = Median, color = Scenario))+
  theme_bw()


impact <- summary_results_pre[, 4:6] - summary_results[, 4:6]
impact_df <- cbind(summary_results_pre[,1:3], impact)

impact_summ <- impact_df %>% group_by(AgeGroup) %>% summarise(impact_all = sum(Median))
pre_vacc <- summary_results_pre %>% group_by(AgeGroup) %>% summarise(prevacc_med = sum(Median))
pre_vacc <- pre_vacc[,2]
impact_perc <- impact_summ$impact_all/pre_vacc * 100


#-------------------------------------------------------------------------------

# Create a sequence of vaccination rates
x_values <- seq(1, 0.1, by = -0.1)

# Calculate vaccination_rate_seq using the formula 1 - (x)^(1/6)
vaccination_rate_seq <- 1 - (x_values)^(1/6)

# Create an empty list to store results for all combinations
all_results_by_delay_rate <- list()

# Loop through delays from 1 to 52
for (delay in seq_len(52)) {
  
  # Loop through vaccination rates
  for (vaccination_rate in vaccination_rate_seq) {
    
    # Print current delay and vaccination rate for tracking
    print(sprintf("Running simulations for delay: %d and vaccination_rate: %.4f", delay, vaccination_rate))
    
    # Create a temporary list to store results for the current combination
    combination_results <- list()
    
    for (i_iter in seq_len(num_iters)) {
      
      print(sprintf("  Running iteration: %d of %d for delay: %d and vaccination_rate: %.4f", 
                    i_iter, num_iters, delay, vaccination_rate))
      
      # Extract parameters for the current iteration
      base_beta <- posterior_cases_prevacc$base_beta[i_iter, ]
      I0        <- posterior_cases_prevacc$I0[i_iter, ]
      gamma     <- posterior_cases_prevacc$gamma[i_iter]
      rho       <- posterior_cases_prevacc$rho[i_iter]
      
      # Simulate for Post-Vaccination with current delay and vaccination rate
      sim_result <- sirv_sim(
        T = 52, 
        A = 18,
        N = N[,1],
        r = rep(0, 18),
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        vaccination_rate = vaccination_rate,  # Varying vaccination rate
        delay = delay,  # Varying delay
        VE_block = 0.75,
        vaccine_target_age = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
        vaccine_duration = 6
      )
      
      # Store results in a combined data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = "Post-Vaccination",
          Delay = delay,  # Add delay as a variable
          VaccinationRate = vaccination_rate,  # Add vaccination rate as a variable
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      combination_results[[length(combination_results) + 1]] <- sim_df
    }
    
    # Combine all results for the current combination and add to the main list
    all_results_by_delay_rate[[paste(delay, vaccination_rate, sep = "_")]] <- do.call(rbind, combination_results)
  }
}

combinations <- expand.grid(
  delay = seq_len(52),
  vaccination_rate = vaccination_rate_seq
)

for (i in seq_len(nrow(combinations))) {
  delay <- combinations$delay[i]
  vaccination_rate <- combinations$vaccination_rate[i]
  
  # Print current combination for tracking
  print(sprintf("Running simulations for delay: %d and vaccination_rate: %.4f", delay, vaccination_rate))
  
  # Simulate for Post-Vaccination with current delay and vaccination rate
  sim_result <- sirv_sim(
    T = 52, 
    A = 18,
    N = N[,1],
    r = rep(0, 18),
    base_beta = base_beta,
    I0_draw = I0,
    rho = rho,
    gamma = gamma,
    vaccination_rate = vaccination_rate,  # Varying vaccination rate
    delay = delay,  # Varying delay
    VE_block = 0.75,
    vaccine_target_age = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
    vaccine_duration = 6
  )
  
  # Store results in a combined data frame
  sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
  sim_df <- sim_df %>%
    mutate(
      Scenario = "Post-Vaccination",
      Delay = delay,  # Add delay as a variable
      VaccinationRate = vaccination_rate,  # Add vaccination rate as a variable
      AgeGroup = as.numeric(Var1),
      Week = as.numeric(Var2)
    )
  
  # Combine results for the current combination
  all_results_by_delay_rate[[paste(delay, vaccination_rate, sep = "_")]] <- sim_df
}

# Initialize a data frame to store total cases and % cases averted
# Step 1: Calculate and save baseline cases (VaccinationRate = 0)
baseline_cases <- NULL
for (combination_key in names(all_results_by_delay_rate)) {
  sim_results <- all_results_by_delay_rate[[combination_key]]
  
  # Extract baseline cases (VaccinationRate == 0)
  baseline_cases_combination <- sim_results %>%
    filter(VaccinationRate == 0) %>%
    summarise(BaselineCases = sum(Cases), .groups = "drop") %>%
    pull(BaselineCases)
  
  # Assign the baseline cases (they should be consistent across all combinations)
  if (is.null(baseline_cases)) {
    baseline_cases <- baseline_cases_combination
    break  # Exit loop as baseline cases are now calculated
  }
}

# Step 2: Loop through combinations to calculate Cases Averted Percent
summary_results <- data.frame()

for (combination_key in names(all_results_by_delay_rate)) {
  # Extract results for the combination
  sim_results <- all_results_by_delay_rate[[combination_key]]
  
  # Calculate total cases for each combination of delay and vaccination rate
  total_cases <- sim_results %>%
    group_by(Delay, VaccinationRate) %>%
    summarise(TotalCases = sum(Cases), .groups = "drop")
  
  # Add BaselineCases and calculate CasesAvertedPercent
  total_cases <- total_cases %>%
    mutate(
      BaselineCases = baseline_cases,
      CasesAvertedPercent = (BaselineCases - TotalCases) / BaselineCases * 100
    )
  
  # Store results for the current combination
  summary_results <- bind_rows(summary_results, total_cases)
}

vaccination_rate_labels <- seq_along(vaccination_rate_seq)  
summary_results <- summary_results %>%
  mutate(
    CasesAvertedPercent = round(CasesAvertedPercent, 1),  # Round to 1 decimal point
    VaccinationRateStep = factor(
      match(VaccinationRate, vaccination_rate_seq),  # Map actual rates to equal steps
      labels = paste0(1:10, " (", round(vaccination_rate_seq * 100, 1), "%)")  # Combine with vaccination rate in %
    )
  )

ggplot(summary_results, aes(x = VaccinationRateStep, y = Delay, fill = CasesAvertedPercent)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "% cases averted") +
  #scale_fill_gradient(low = "grey", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )+
  labs(x = "weekly vaccination rate")


#-------------------------------------------------------------------------------
# vaccination supply: 30% of total population (but we can vary from 0 to 50%)
# total vacc supply = total_pop x 0.3
# scenario 1: prioritised age < 20 years old
# roll out speed (r): vary from 0.2% to 1% 
# n_vax_total = N x r 
# n_vax_agegr20 = N_20 / N x r (total duration)

# population data 
N_region1 <- N[,1]

N_20 <- sum(N_region1[1:4])
N_20_49 <- sum(N_region1[5:12])
N_20plus <- sum(N_region1[5:18])
N_60plus <- sum(N_region1[15:18])
N_all <- sum(N_region1)
total_coverage <- 0.3
daily_rate <- 0.01
total_supply <- N_all * total_coverage
tot_vacc_dur <- total_supply / (N_all * daily_rate) # 30 days 
tot_vacc_dur_week <- tot_vacc_dur / 7

n_vax_tot <- tot_supply / tot_vacc_dur_week # weekly 

n_vax_20 <- N_20 / n_vax_tot /7
n_vax_2049 <- N_20_49 / n_vax_tot /7
n_vax_20plus <- N_20plus / n_vax_tot/7
N_vax_60plus <- N_60plus / n_vax_tot/7

sim_result <- sirv_sim(
  T = 52, 
  A = 18,
  N = N[,1],
  r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
  
  base_beta = base_beta,
  I0_draw = I0,
  rho = rho,
  gamma = gamma,
  
  # Vaccination scenario
  vaccination_rate = 0.2*7,  # Weekly vaccination rate
  delay = 0,                   # Vaccination starts from Week 1
  VE_block = 0.75,             # Vaccine efficacy
  vaccine_target_age = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),  # Targeting specific age groups
  vaccine_duration = round(N_20 / n_vax_tot / 7)         # Vaccination duration in weeks
)

target_age_list <- list(c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                        c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1))

vax_duration <- c(round(N_20 / n_vax_tot),
                  round(N_20_49 / n_vax_tot),
                  round(N_20plus / n_vax_tot),
                  round(N_60plus / n_vax_tot))
                        
vax_result <- list()

for(i in seq_along(target_age_list)){
   
    target <- target_age_list[[i]]
    duration <- vax_duration[i]
    
    print(paste("Scenario", i, ": Target Age =", toString(target), 
                ", Duration =", duration))

    sim_result <- sirv_sim(
      T = 52, 
      A = 18,
      N = N[,2],
      r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
      
      base_beta = base_beta,
      I0_draw = I0,
      rho = rho,
      gamma = gamma,
      
      # Vaccination scenario
      vaccination_rate = 0.2*7,  # Weekly vaccination rate
      delay = 0,                   # Vaccination starts from Week 1
      VE_block = 0.75,             # Vaccine efficacy
      vaccine_target_age = target,  # Targeting specific age groups
      vaccine_duration = duration        # Vaccination duration in weeks
    )
    
   vax_result[[i]] <- sim_result
}

all_results <- list()

for (i_iter in seq_len(num_iters)) {
  
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  vax_result <- list()
  
  # Inner loop for vaccination scenarios
  for (i in seq_along(target_age_list)) {
    
    target <- target_age_list[[i]]
    duration <- vax_duration[i]
    
    print(paste("Iteration", i_iter, "Scenario", i, ": Target Age =", toString(target), 
                ", Duration =", duration))
    
    base_beta <- posterior_prevacc$base_beta[i_iter, ]
    I0        <- posterior_prevacc$I0[i_iter, ]
    gamma     <- posterior_prevacc$gamma[i_iter]
    rho       <- posterior_prevacc$rho[i_iter]
    
    # Simulate the SIRV model
    sim_result <- sirv_sim(
      T = 52, 
      A = 18,
      N = N[,1],
      r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
      
      base_beta = base_beta,
      I0_draw = I0,
      rho = rho,
      gamma = gamma,
      
      # Vaccination scenario
      vaccination_rate = 1 * 7 /100,  # Weekly vaccination rate
      delay = 3,                   # Vaccination starts from Week 1
      VE_block = 0.75,             # Vaccine efficacy
      vaccine_target_age = target, # Targeting specific age groups
      vaccine_duration = duration  # Vaccination duration in weeks
    )
    
    # Store results of this scenario
    vax_result[[i]] <- sim_result
  }
  
  # Store all results from this iteration
  all_results[[i_iter]] <- vax_result
}

scenario1 <- lapply(all_results, function(x) x[[1]])

sim_df_list <- lapply(seq_along(scenario1), function(i_iter){
  as.data.frame.table(scenario1[[i_iter]]$age_stratified_cases, responseName = "Cases")%>% 
    mutate(
      Iteration = i_iter,
      Scenario = "Post-Vaccination",
      AgeGroup = as.numeric(Var1),
      Week = as.numeric(Var2)
    )
})

comb_sim <- do.call(rbind, sim_df_list)
summary_cases <- comb_sim %>%
  group_by(AgeGroup, Week, Scenario) %>%
  summarise(
    Median = median(Cases),
    Lower = quantile(Cases, probs = 0.025),
    Upper = quantile(Cases, probs = 0.975),
    .groups = "drop"
  )

summary_cases_med <- summary_cases[,1:4]
summary_cases_pre_med <- summary_results_pre[,1:4]
all_cases <- rbind(summary_cases_pre_med, summary_cases_med)

ggplot(all_cases)+
  geom_line(aes(x = Week, y = Median, colour = Scenario))+
  facet_wrap(~AgeGroup)+
  theme_bw()

all_cases$Color <- ifelse(all_cases$AgeGroup %in% 1:4, "Highlight", "Grey")
all_cases_postvacc <- all_cases %>% filter(Scenario == "Post-Vaccination")
all_cases_postvacc$tot_pop <- sum(N_region1)
all_cases_postvacc <- all_cases_postvacc %>% mutate(prop_infec = Median/tot_pop * 100)
all_cases_postvacc$prevacc <- summary_cases_pre_med$Median
all_cases_postvacc <- all_cases_postvacc %>% mutate(diff = prevacc - Median)
total_impact <- sum(all_cases_postvacc$diff) / sum(all_cases_postvacc$prevacc) * 100
  
custom_colors <- c("Highlight" = "magenta", "Grey" = "grey")


ggplot(all_cases_postvacc, aes(x = Week, y = prop_infec, group = AgeGroup, color = Color)) +
  geom_line(size = 0.8) +  # Line plot for each age group
  annotate(
    "rect", xmin = 1, xmax = 21, ymin = 0, ymax = Inf,  # Shaded area for vaccination duration
    alpha = 0.2, fill = "grey"
  ) +
  scale_color_manual(values = custom_colors) +  # Custom colors for highlighted groups
  labs(
    x = "Week",
    y = "infected (%)",
    title = "Age-Stratified Infections",
    color = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

summary_cases_pre_all <- summary_results_pre %>% group_by(Week) %>%
                                summarise(
                                  Median   = sum(Median)
                                ) %>% mutate(
                                  Scenario = "Pre-vaccination"
                                )
summary_cases_post_all <- summary_cases_med %>% group_by(Week) %>%
  summarise(
    Median   = sum(Median)
  ) %>% mutate(
    Scenario = "Post-vaccination"
  )

all_cases_summ <- rbind(summary_cases_pre_all, summary_cases_post_all)

# Define the number of scenarios
num_scenarios <- 4

# Create an empty list to store the results for all scenarios
all_scenario_summaries <- vector("list", length = num_scenarios)

# Loop over all scenarios
for (scen_idx in seq_len(num_scenarios)) {
  
  # Extract the current scenario data from all iterations
  scenario_data <- lapply(all_results, function(x) x[[scen_idx]])
  
  # Process each iteration for the current scenario
  sim_df_list <- lapply(seq_along(scenario_data), function(i_iter) {
    as.data.frame.table(
      scenario_data[[i_iter]]$age_stratified_cases,
      responseName = "Cases"
    ) %>%
      mutate(
        Iteration = i_iter,
        Scenario = paste0("Scenario_", scen_idx),  # Label by scenario index
        AgeGroup = as.numeric(Var1),
        Week     = as.numeric(Var2)
      )
  })
  
  # Combine all iterations for the current scenario
  comb_sim <- do.call(rbind, sim_df_list)
  
  # Summarize the results (Median, Lower, Upper)
  summary_cases <- comb_sim %>%
    group_by(AgeGroup, Week, Scenario) %>%
    summarise(
      Median = median(Cases),
      Lower  = quantile(Cases, probs = 0.025),
      Upper  = quantile(Cases, probs = 0.975),
      .groups = "drop"
    )
  
  # Summarize across all age groups for the current scenario
  summary_cases_all <- summary_cases %>%
    group_by(Week, Scenario) %>%
    summarise(
      Median = sum(Median),
      Lower  = sum(Lower),
      Upper  = sum(Upper),
      .groups = "drop"
    )
  
  # Store the result for the current scenario
  all_scenario_summaries[[scen_idx]] <- summary_cases_all
}

# Combine all scenarios into one data frame
final_summary_all_scenarios <- bind_rows(all_scenario_summaries, summary_cases_pre_all)

ggplot(final_summary_all_scenarios) +
  geom_line(aes(x = Week, y = Median, colour = Scenario, linetype = Scenario)) +
  scale_linetype_manual(
    values = c("Pre-vaccination" = "dashed", "Scenario_1" = "solid", 
               "Scenario_2" = "solid", "Scenario_3" = "solid", "Scenario_4" = "solid")
  ) +
  theme_bw()


#------------------------------------------------------------------------------------------
# vary delay with fixed rolling speed across different scenarios under "30% total coverage"
# fixed roll out speed as (0.2% per day) 1.4% per week but varied duration by total coverage goal 
# for each scenario of a) vaccinate under 20, b) 20-40, c) 20+, d) 60+, e) all 
# if coverage goal is 30%, duration = 
# for a current scenario (under 20), we are varying delay (0-52) 
# let's test within a single scenario, vary duration and delay
#------------------------------------------------------------------------------------------
total_coverage <- seq(0.1, 1, by = 0.1)
daily_rate <- 0.01
total_supply <- N_all * total_coverage
tot_vacc_dur <- total_supply / (N_all * daily_rate) # 30 days 
tot_vacc_dur_week <- tot_vacc_dur / 7

n_vax_tot <- N_all * daily_rate * 7 # total weekly supply
n_vax_tot * tot_vacc_dur_week

all_results_comb <- list()  # Initialize storage for all iterations
summary_dfs <- list()  # Initialize storage for summarized data frames

for (i_iter in seq_len(num_iters)) {
  
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  # Inner loop for vaccination scenarios
  for (i in seq_along(target_age_list)) {
    
    target <- target_age_list[[i]]
    duration <- vax_duration[i]
    
    print(paste("Iteration", i_iter, "Scenario", i, ": Target Age =", toString(target), 
                ", Duration =", duration))
    
    # Initialize storage for results within this scenario
    scenario_result <- list()
    
    # Loop through delay values (e.g., 0 to 52 weeks)
    for (delay in seq(0, 52, by = 1)) {
      print(paste("Iteration", i_iter, "Scenario", i, "Delay:", delay))
      
      base_beta <- posterior_prevacc$base_beta[i_iter, ]
      I0        <- posterior_prevacc$I0[i_iter, ]
      gamma     <- posterior_prevacc$gamma[i_iter]
      rho       <- posterior_prevacc$rho[i_iter]
      
      # Simulate the SIRV model
      sim_result <- sirv_sim(
        T = 52, 
        A = 18,
        N = N[,1],
        r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        # Vaccination scenario
        vaccination_rate = 1 * 7 / 100,    # Weekly vaccination rate
        delay = delay,                     # Varying delay
        VE_block = 0.75,                   # Vaccine efficacy
        vaccine_target_age = target,       # Targeting specific age groups
        vaccine_duration = duration        # Vaccination duration in weeks
      )
      
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = paste0("Scenario_", i, "_Delay_", delay),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs[[length(summary_dfs) + 1]] <- sim_df
      
      # Store results for this delay
      scenario_result[[as.character(delay)]] <- sim_result
    }
    
    # Store results for this scenario
    iteration_result[[i]] <- scenario_result
  }
  
  # Store all results from this iteration
  all_results[[i_iter]] <- iteration_result
}

# Define parameters
coverage_rate_combinations <- expand.grid(
  total_coverage = seq(0.1, 1, by = 0.1),  # 10 values
  daily_rate = seq(0.01, 0.05, by = 0.01)  # 5 values
)

# Initialize storage for results
all_results_comb <- list()
summary_dfs <- list()

# Simulation loop
# sample run
target <- target_age_list[[1]]
for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  iteration_result <- list()  # Initialize storage for this iteration
  
  for (comb in seq_len(nrow(coverage_rate_combinations))) {
    # Extract combination values
    total_coverage <- coverage_rate_combinations$total_coverage[comb]
    daily_rate <- coverage_rate_combinations$daily_rate[comb]
    
    # Calculate total supply and vaccine duration
    tot_supply <- N_all * total_coverage
    daily_capacity <- N_all * daily_rate
    vaccine_duration <- tot_supply / daily_capacity
    vaccine_duration_weeks <- vaccine_duration / 7
    
    print(sprintf("Iteration: %d, Coverage: %.1f%%, Daily Rate: %.2f%%, Vaccine Duration: %.2f weeks", 
                  i_iter, total_coverage * 100, daily_rate * 100, vaccine_duration_weeks))
    
    scenario_result <- list()  # Initialize results for delays
    
    for (delay in seq(0, 52, by = 1)) {
      print(paste("Iteration", i_iter, "Delay:", delay))
      
      base_beta <- posterior_prevacc$base_beta[i_iter, ]
      I0        <- posterior_prevacc$I0[i_iter, ]
      gamma     <- posterior_prevacc$gamma[i_iter]
      rho       <- posterior_prevacc$rho[i_iter]
      
      # Simulate the SIRV model
      sim_result <- sirv_sim(
        T = 52, 
        A = 18,
        N = N[,1],
        r = rep(0, 18),
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        # Adjusted vaccination scenario
        vaccination_rate = daily_rate * 7,
        delay = delay,
        VE_block = 0.75,
        vaccine_target_age = target,
        vaccine_duration = vaccine_duration_weeks
      )
      
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = paste0("Scenario_1_Delay_", delay, 
                            "_Coverage_", total_coverage, 
                            "_Rate_", daily_rate),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs[[length(summary_dfs) + 1]] <- sim_df
      
      # Store results for this delay
      scenario_result[[as.character(delay)]] <- sim_result
    }
    
    # Store results for this combination
    iteration_result[[paste0("Scenario_1_Coverage_", total_coverage, "_Rate_", daily_rate)]] <- scenario_result
  }
  
  # Store all results for this iteration
  all_results_comb[[i_iter]] <- iteration_result
}

# loop through all scenarios
for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %d of %d", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  for (comb in seq_len(nrow(coverage_rate_combinations))) {
    # Extract combination values
    total_coverage <- coverage_rate_combinations$total_coverage[comb]
    daily_rate <- coverage_rate_combinations$daily_rate[comb]
    
    # Calculate total supply and vaccine duration
    tot_supply <- N_all * total_coverage
    daily_capacity <- N_all * daily_rate
    vaccine_duration <- tot_supply / daily_capacity  # Duration in days
    vaccine_duration_weeks <- vaccine_duration / 7  # Convert to weeks
    
    print(sprintf("Iteration: %d, Coverage: %.1f%%, Daily Rate: %.2f%%, Vaccine Duration: %.2f weeks", 
                  i_iter, total_coverage * 100, daily_rate * 100, vaccine_duration_weeks))
    
    # Inner loop for vaccination scenarios
    for (i in seq_along(target_age_list)) {
      target <- target_age_list[[i]]
      
      print(paste("Iteration", i_iter, "Scenario", i, ": Target Age =", toString(target)))
      
      # Initialize storage for results within this scenario
      scenario_result <- list()
      
      # Loop through delay values (e.g., 0 to 52 weeks)
      for (delay in seq(0, 52, by = 1)) {
        print(paste("Iteration", i_iter, "Scenario", i, "Delay:", delay))
        
        base_beta <- posterior_prevacc$base_beta[i_iter, ]
        I0        <- posterior_prevacc$I0[i_iter, ]
        gamma     <- posterior_prevacc$gamma[i_iter]
        rho       <- posterior_prevacc$rho[i_iter]
        
        # Simulate the SIRV model
        sim_result <- sirv_sim(
          T = 52, 
          A = 18,
          N = N[,1],
          r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
          
          base_beta = base_beta,
          I0_draw = I0,
          rho = rho,
          gamma = gamma,
          
          # Adjusted vaccination scenario
          vaccination_rate = daily_rate * 7,         # Weekly vaccination rate
          delay = delay,                             # Varying delay
          VE_block = 0.75,                           # Vaccine efficacy
          vaccine_target_age = target,              # Targeting specific age groups
          vaccine_duration = vaccine_duration_weeks # Adjusted duration
        )
        
        # Summarize the simulation result into a data frame
        sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
        sim_df <- sim_df %>%
          mutate(
            Iteration = i_iter,
            Scenario = paste0("Scenario_", i, "_Delay_", delay, 
                              "_Coverage_", total_coverage, 
                              "_Rate_", daily_rate),
            AgeGroup = as.numeric(Var1),
            Week = as.numeric(Var2)
          )
        
        # Store the summarized data frame
        summary_dfs[[length(summary_dfs) + 1]] <- sim_df
        
        # Store results for this delay
        scenario_result[[as.character(delay)]] <- sim_result
      }
      
      # Store results for this scenario
      iteration_result[[paste0("Scenario_", i, "_Coverage_", total_coverage, "_Rate_", daily_rate)]] <- scenario_result
    }
  }
  
  # Store all results from this iteration
  all_results_comb[[i_iter]] <- iteration_result
}

# test the code
all_results_delay <- list()  # Initialize storage for all iterations
summary_dfs <- list()  # Initialize storage for summarized data frames

# Fixed to the first scenario
target <- target_age_list[[1]]
duration <- tot_vacc_dur_week

# Loop through iterations
for (i_iter in seq_len(num_iters)) {
  print(sprintf("Running iteration: %d of %d for the first scenario", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  # Loop through vaccination durations
  for (duration in tot_vacc_dur_week) {  # Loop over varied durations
    print(sprintf("Testing Duration: %.2f weeks", duration))
    
    # Loop through delay values (e.g., 0 to 52 weeks)
    for (delay in seq(0, 52, by = 1)) {
      print(paste("Iteration", i_iter, "Delay:", delay, "Duration:", duration))
      
      # Extract posterior parameters for the current iteration
      base_beta <- posterior_prevacc$base_beta[i_iter, ]
      I0        <- posterior_prevacc$I0[i_iter, ]
      gamma     <- posterior_prevacc$gamma[i_iter]
      rho       <- posterior_prevacc$rho[i_iter]
      
      # Simulate the SIRV model
      sim_result <- sirv_sim(
        T = 52, 
        A = 18,
        N = N[,1],
        r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        # Vaccination scenario
        vaccination_rate = 7 / 100,       # Weekly vaccination rate
        delay = delay,                    # Varying delay
        VE_block = 0.75,                  # Vaccine efficacy
        vaccine_target_age = target,      # Targeting specific age groups
        vaccine_duration = duration       # Varying duration
      )
      
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Scenario = paste0("Scenario_1_Delay_", delay, "_Duration_", round(duration, 2)),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs[[length(summary_dfs) + 1]] <- sim_df
      
      # Store results for this delay and duration
      scenario_key <- paste0("Delay_", delay, "_Duration_", round(duration, 2))
      iteration_result[[scenario_key]] <- sim_result
    }
  }
  
  # Store all results from this iteration
  all_results[[i_iter]] <- iteration_result
}


# Extract unique scenarios from the nested structure
unique_scenarios <- unique(unlist(lapply(all_results, function(iteration) names(iteration))))

# Function to extract data for a given scenario
extract_scenario_data <- function(scenario, all_results) {
  scenario_data <- lapply(all_results, function(iteration) {
    if (scenario %in% names(iteration)) {
      iteration[[scenario]]  # Extract the data frame
    } else {
      NULL  # Handle missing scenario gracefully
    }
  })
  
  # Remove NULL values
  scenario_data <- scenario_data[!sapply(scenario_data, is.null)]
  
  # Combine all iterations for this scenario into one data frame
  combined_scenario_data <- do.call(rbind, lapply(scenario_data, function(sim_result) {
    as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases") %>%
      mutate(
        Scenario = scenario,
        AgeGroup = as.numeric(Var1),
        Week = as.numeric(Var2)
      )
  }))
  
  return(combined_scenario_data)
}

# Test for one scenario
test_scenario <- unique_scenarios[1]
test_data <- extract_scenario_data(test_scenario, all_results)
print(head(test_data))  # Check the extracted data

compute_median <- function(scenario_data) {
  scenario_data %>%
    group_by(Week, AgeGroup) %>%  # Group by Week and AgeGroup
    summarise(
      MedianCases = median(Cases, na.rm = TRUE),  # Compute median
      .groups = "drop"
    ) %>%
    mutate(Scenario = unique(scenario_data$Scenario))  # Add Scenario column
}

# Test median computation for one scenario
median_data <- compute_median(test_data)

# Initialize a list to store results for all scenarios
scenario_summary_list <- list()

for (scenario in unique_scenarios) {
  print(paste("Processing:", scenario))
  
  # Extract data for the current scenario
  scenario_data <- extract_scenario_data(scenario, all_results)
  
  # Compute median across iterations
  median_summary <- compute_median(scenario_data)
  
  # Store the summarized data
  scenario_summary_list[[scenario]] <- median_summary
}

final_summ <- lapply(scenario_summary_list, function(df) {
  df %>%
    group_by(Week, Scenario) %>%  # Include Scenario in grouping
    summarise(
      Median = sum(MedianCases, na.rm = TRUE),  # Summarize MedianCases
      .groups = "drop"  # Drop grouping for clarity in the output
    ) %>% mutate(
      pre_vacc = summary_cases_pre_all$Median,
      diff     = pre_vacc - Median,
      impact   = diff / pre_vacc * 100
    )
})
total_impact <- lapply(final_summ, function(df) {
  df %>% group_by(Scenario) %>%
    summarise(tot_diff = sum(diff),
              pre_vacc = sum(pre_vacc))%>% mutate(
      pre_vacc = sum(pre_vacc),
      tot_impact = tot_diff / pre_vacc * 100
    )
})

final_total_impact_df <- do.call(rbind, total_impact)

summary_table <- final_total_impact_df %>%
  # Extract Delay and Duration using regex
  mutate(
    Delay = as.numeric(sub("Delay_", "", str_extract(Scenario, "Delay_\\d+"))),
    Duration = as.numeric(sub("Duration_", "", str_extract(Scenario, "Duration_\\d+(\\.\\d+)?")))
  )
summary_table$tot_impact <- as.numeric(summary_table$tot_impact)
summary_table$Duration <- as.factor(summary_table$Duration)
unique_durations <- unique(summary_table$Duration)
coverage_mapping <- data.frame(
  Duration = unique_durations,
  vaccine_coverage = total_coverage
)
summary_table <- summary_table %>%
  left_join(coverage_mapping, by = "Duration")
summary_table$vaccine_coverage <- as.factor(summary_table$vaccine_coverage)

ggplot(summary_table, aes(x = vaccine_coverage, y = Delay, fill = tot_impact)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "% cases averted") +
  #scale_fill_gradient(low = "grey", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )+
  labs(x = "Total vaccine coverage")


#-------------------------------------------------------------------------------
# varying duration by coverage 0-100% (with fixed weekly rate 1.4%) and delay 1-52
all_results <- list()  # Initialize storage for all iterations
summary_dfs <- list()  # Initialize storage for summarized data frames

# Fixed to the first scenario
target <- target_age_list[[1]]

# Define total coverage levels (from 0% to 100% by 10%)
coverage_levels <- seq(0, 1, by = 0.1)

# Define weekly vaccination speed (percentage of the population vaccinated per week)
weekly_speed <- 0.014  # 1.4%

for (i_iter in seq_len(num_iters)) {
  
  print(sprintf("Running iteration: %d of %d for the first scenario", i_iter, num_iters))
  
  # Initialize storage for results within this iteration
  iteration_result <- list()
  
  # Loop through delay values (e.g., 0 to 52 weeks)
  for (delay in seq(0, 52, by = 1)) {
    print(paste("Iteration", i_iter, "Delay:", delay))
    
    delay_results <- list()  # Storage for results at this delay
    
    # Loop through coverage levels
    for (coverage in coverage_levels) {
      print(paste("Coverage:", coverage * 100, "%"))
      
      # Calculate total vaccine supply and weekly vaccination capacity
      total_population <- sum(N[, 1])  # Total population across all age groups
      total_vaccine_supply <- coverage * total_population
      n_vax_tot <- weekly_speed * total_population  # Weekly vaccination capacity
      vaccination_duration <- (total_vaccine_supply / n_vax_tot) / 7  # Duration in weeks
      
      base_beta <- posterior_prevacc$base_beta[i_iter, ]
      I0        <- posterior_prevacc$I0[i_iter, ]
      gamma     <- posterior_prevacc$gamma[i_iter]
      rho       <- posterior_prevacc$rho[i_iter]
      
      # Simulate the SIRV model
      sim_result <- sirv_sim(
        T = 52, 
        A = 18,
        N = N[, 1],
        r = rep(0, 18),  # Aging rate, here set to 0 for all age groups
        
        base_beta = base_beta,
        I0_draw = I0,
        rho = rho,
        gamma = gamma,
        
        # Vaccination scenario
        vaccination_rate = weekly_speed,  # Normalized weekly rate
        delay = delay,                     # Varying delay
        VE_block = 0.75,                   # Vaccine efficacy
        vaccine_target_age = target,       # Targeting specific age groups
        vaccine_duration = vaccination_duration  # Calculated vaccination duration
      )
      
      # Summarize the simulation result into a data frame
      sim_df <- as.data.frame.table(sim_result$age_stratified_cases, responseName = "Cases")
      sim_df <- sim_df %>%
        mutate(
          Iteration = i_iter,
          Delay = delay,
          Coverage = coverage * 100,  # Coverage as a percentage
          WeeklySpeed = weekly_speed * 100,  # Weekly vaccination speed as a percentage
          Duration = vaccination_duration,
          Scenario = paste0("Scenario_1_Delay_", delay, "_Coverage_", coverage * 100),
          AgeGroup = as.numeric(Var1),
          Week = as.numeric(Var2)
        )
      
      # Store the summarized data frame
      summary_dfs[[length(summary_dfs) + 1]] <- sim_df
      
      # Store results for this coverage level
      delay_results[[as.character(coverage)]] <- sim_result
    }
    
    # Store results for this delay
    iteration_result[[as.character(delay)]] <- delay_results
  }
  
  # Store all results from this iteration
  all_results[[i_iter]] <- iteration_result
}


