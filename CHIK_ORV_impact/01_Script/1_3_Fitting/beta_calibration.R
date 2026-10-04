library(deSolve)

bra_foi <- allfoi %>% filter(country == "Brazil")
lon_target <- -43.7908
lat_target <- -17.9302
tol <- 0.05

mean(bra_foi$foi_mid, na.rm = TRUE)

# Filter the dataframe using subset
bra_mg <- subset(bra_foi, abs(x - lon_target) <= tol & abs(y - lat_target) <= tol)

foi_bra_mg <- mean(bra_mg$foi_mid)

# state
bra_foi_sf <- st_as_sf(bra_foi, coords = c("x", "y"), crs = 4326)
br_states <- st_read("00_Data/0_1_Raw/country_shape/gadm41_BRA_shp/gadm41_BRA_1.shp")
br_states <- st_transform(br_states, crs = st_crs(bra_foi_sf))
bra_foi_states <- st_join(bra_foi_sf, br_states, join = st_intersects)

bra_foi_state_summ <- bra_foi_states %>%
  group_by(NAME_1) %>%
  summarise(avg_foi = mean(foi_mid, na.rm = TRUE))%>%
  filter(!is.na(NAME_1))

# 1. Define the SIR model function (using standard incidence and proportions)
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # With normalized populations, the standard SIR equations become:
    dS <- - beta * S * I       # since S + I + R = 1
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# 2. Set up parameters and initial conditions in PROPORTION terms
# Here, we normalize the total population to 1.
N <- 1       # normalized total population
gamma_val <- 1/7            # recovery rate (per day); infectious period = 7 days
FOI_target <- foi_bra_mg   # target annual FOI (1% per year)

## foi target list (loop through all possible `fois`)
foi_list <- quantile(bra_foi$foi_mid[bra_foi$foi_mid != 0],
                     probs = c(0.25, 0.5, 0.975, 1),
                     na.rm = TRUE)

foi_state_list <- bra_foi_state_summ$avg_foi

age_groups <- c(mean(0:4),
                mean(5:9),
                mean(10:14),
                mean(15:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)
w_age <- N_mg / sum(N_mg)
S_age <- exp(-FOI_target * age_groups)
I_age <- (FOI_target / 365) * exp(-FOI_target * age_groups) * (1 / gamma_val)
R_age <- 1 - exp(-FOI_target * age_groups)

S_avg <- sum(S_age * w_age)
I_avg <- sum(I_age * w_age)
R_avg <- sum(R_age * w_age)

init_state <- c(
  S = S_avg,
  I = I_avg,
  R = R_avg
)

## loop
init_state_list <- vector("list", length(foi_state_list))

for(i in seq_along(foi_state_list)){
  FOI_target <- foi_state_list[[i]]
  # Age-specific S, I, R
  S_age <- exp(-FOI_target * age_groups)
  I_age <- (FOI_target / 365) * exp(-FOI_target * age_groups) * (1 / gamma_val)
  R_age <- 1 - exp(-FOI_target * age_groups)
  
  # Weighted averages
  S_avg <- sum(S_age * w_age)
  I_avg <- sum(I_age * w_age)
  R_avg <- sum(R_age * w_age)
  
  # Create your init_state vector (or list)
  init_state <- c(S = S_avg, I = I_avg, R = R_avg)
  
  # Store in results_list
  init_state_list[[i]] <- list(
    init_state = init_state
  )
}



# 3. Time vector: simulate for 365 days (daily steps)
times <- seq(0, 365, by = 1)

# 4. Define the objective function that computes (and prints) daily FOI values.
#    In a normalized model, the daily FOI is computed as:
#         FOI_daily = S * beta_candidate * I * (1/gamma_val)
#    (Note: since N = 1, S/N becomes S)
objective_function <- function(beta_candidate, init_state, FOI_target) {
  # 1. Set up parameters for the SIR model
  params <- c(beta = beta_candidate, gamma = gamma_val)  # N=1 if normalized
  
  # 2. Run the model for 365 days using the provided init_state
  sim <- ode(y = init_state, times = times, func = sir_model, parms = params)
  sim_df <- as.data.frame(sim)
  
  # 3. Get S(0) and S(365)
  S0_sim  <- sim_df$S[1]
  S365_sim <- tail(sim_df$S, 1)
  
  # 4. Compute the annual FOI as the fraction of S lost over the year
  FOI_simulated <- (S0_sim - S365_sim) / S0_sim
  
  cat(sprintf(
    "Beta candidate = %.4f | FOI_target = %.6f | Simulated Annual FOI = %.6f\n",
    beta_candidate, FOI_target, FOI_simulated
  ))
  
  # 5. Compute error
  error <- (FOI_simulated - FOI_target)^2
  
  return(list(
    beta_candidate = beta_candidate,
    error = error,
    FOI_simulated = FOI_simulated,
    FOI_target = FOI_target
  ))
}


# 5. Grid search for beta candidates (increments of 0.05)
beta_values <- seq(0, 1, by = 0.001)
results_list <- vector("list", length(beta_values))

for (i in seq_along(beta_values)) {
  res <- objective_function(beta_values[i], init_state, FOI_target)
  results_list[[i]] <- res
}

# Convert the results list to a data frame for easier inspection
results_df <- do.call(rbind, lapply(results_list, as.data.frame))
print(results_df)

all_results <- vector("list", length(init_state_list))

for (j in seq_along(init_state_list)) {
  # Set the current initial state
  current_init_state <- as.numeric(unlist(init_state_list[[j]]))
  names(current_init_state) <- c("S", "I", "R")
  
  # Create an empty list for the beta grid search results for this init state
  results_list <- vector("list", length(beta_values))
  
  # Loop over beta candidates
  for (i in seq_along(beta_values)) {
    results_list[[i]] <- objective_function(beta_values[i], current_init_state, FOI_target)
  }
  
  # Convert the results list to a data frame
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  
  # Identify the best beta for this initial state (minimizing error)
  best_idx <- which.min(results_df$error)
  best_beta <- results_df$beta_candidate[best_idx]
  
  # Store results for the jth init state
  all_results[[j]] <- list(
    init_state = current_init_state,
    results_df = results_df,
    best_beta = best_beta
  )
}

# Identify the beta with the smallest error
best_idx <- which.min(results_df$error)
best_beta <- results_df$beta_candidate[best_idx]
cat("\nBest beta from grid search:", best_beta, "\n")

best_beta_list <- list()
for(i in seq_along(all_results)){
  result <- all_results[[i]]$best_beta
  best_beta_list[[i]] <- result
}

params_final <- list(beta = best_beta, gamma = gamma_val)

# Run the simulation with the calibrated beta
sim_result <- ode(y = init_state, times = times, func = sir_model, parms = params_final)
sim_result <- as.data.frame(sim_result)

# --- Plot the epidemic (infection) curve using ggplot2 --- #

ggplot(sim_result, aes(x = time, y = I)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Epidemic Curve (Proportion Infectious Over Time)",
       x = "Time (days)",
       y = "Proportion Infectious (I)") +
  theme_minimal()

# post-process
load("00_Data/0_2_Processed/bra_combined_indexP.RData")
bra_indexp_mg <- bra_combined_output[[27]]
indexP_raw <- bra_indexp_mg %>% filter(grepl("^2023", date))
indexP_mid <- indexP_raw$mid
indexP_normalized <- indexP_mid / mean(indexP_mid)

beta_time <- best_beta * indexP_normalized

indexP_normalized_list <- lapply(bra_combined_output, function(df) {
  # Filter rows where the date starts with "2023"
  indexP_raw <- df %>% filter(grepl("^2023", date))
  
  # Extract the 'mid' column
  indexP_mid <- indexP_raw$mid
  
  # Normalize: divide by the mean (ignoring NAs if needed)
  indexP_normalized <- indexP_mid / mean(indexP_mid, na.rm = TRUE)
  
  return(indexP_normalized)
})

beta_time_list <- list()

for(i in seq_along(best_beta_list)){
  result <- best_beta_list[[i]] * indexP_normalized_list[[i]]
  beta_time_list[[i]] <- result
}

df_beta <- data.frame(day = 1:365, beta_time = beta_time)

weekly_beta <- df_beta %>%
  mutate(week = ceiling(day / 7)) %>%
  group_by(week) %>%
  summarise(weekly_beta = mean(beta_time))

beta_week <- weekly_beta$weekly_beta[1:52]

