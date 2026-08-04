/// underying immunity incorporated with single rho

data {
  int<lower=1> T;                           // Number of weeks
  int<lower=1> A;                           // Number of age groups
  int<lower=0> observed_cases_by_age[A, T];   // Age-specific observed cases
  
  real<lower=0> N[A];                       // Population per age group
  real<lower=0> r[A];                       // Aging rate per age group
  
  real<lower=0> indexP[T];                  // Capacity index for each week
  real<lower=0, upper=1> eta;               // Desired vaccine coverage (unused directly)
  real<lower=0> vaccination_rate;           // Weekly vaccination rate
  int<lower=1> delay;                       // First week of vaccination
  real<lower=0, upper=1> VE_block;           // Vaccine efficacy
  
  // Binary indicator: which ages get vaccinated (0 or 1)
  int<lower=0, upper=1> vaccine_target_age[A];
  
  // Informative priors for initial infections (e.g., based on external data)
  // prior_I0[a] is the prior mean for I0 in age group a
  real<lower=0> prior_I0[A];
  // A common standard deviation for the initial infection priors
  real<lower=0> prior_sd_I0;
}

parameters {
  // Age-specific initial infected count, with bounds and informative prior
  real<lower=10, upper=500> I0[A];
  
  // Overall time-varying baseline transmission rate (applied to all ages)
  real<lower=0, upper=6> base_beta[T];
  
  // Overall detection probability (uniform across ages)
  real<lower=0, upper=1> rho;
  
  // Overdispersion parameter for the negative binomial likelihood
  real<lower=15, upper=20> shape;
  
  // Recovery rate (assumed the same for all ages)
  real<lower=0, upper=1> gamma;
}

transformed parameters {
  // Normalize the indexP values to [0, 1]
  real<lower=0, upper=1> indexP_normalised[T];
  // Final transmission rate for each week
  real<lower=0> beta[T];
  
  {
    real max_index = indexP[1];
    real min_index = indexP[1];
    for (t in 2:T) {
      if (indexP[t] > max_index) max_index = indexP[t];
      if (indexP[t] < min_index) min_index = indexP[t];
    }
    for (t in 1:T) {
      if (max_index > min_index)
        indexP_normalised[t] = (indexP[t] - min_index) / (max_index - min_index);
      else
        indexP_normalised[t] = 0;
    }
  }
  
  for (t in 1:T) {
    beta[t] = base_beta[t] * (1 + 0.5 * indexP_normalised[t]);
  }
}

model {
  // Declare state variables for each age and week:
  real S[A, T];
  real I[A, T];
  real R[A, T];
  real V[A, T];   // Vaccinated and protected
  real IV[A, T];  // Vaccinated but infected (vaccine failure)
  real RV[A, T];  // Recovered from vaccinated infections
  
  // --- Priors ---
  // Prior for base_beta[t]
  for (t in 1:T) {
    base_beta[t] ~ lognormal(-1, 0.5);
  }
  // Prior for detection probability
  rho ~ beta(30, 70);
  // Prior for recovery rate
  gamma ~ normal(0.67, 0.02);
  // Prior for overdispersion parameter
  shape ~ uniform(15, 20);
  
  // Age-specific informative priors for I0:
  for (a in 1:A) {
    I0[a] ~ normal(prior_I0[a], prior_sd_I0) T[10, 1000];
  }
  
  // --- 1) Initial conditions at t=1 ---
  for (a in 1:A) {
    S[a, 1] = N[a] - I0[a];
    I[a, 1] = I0[a];
    R[a, 1] = 0;
    V[a, 1] = 0;
    IV[a, 1] = 0;
    RV[a, 1] = 0;
  }
  
  // --- 2) SIRV Dynamics ---
  for (t in 2:T) {
    // Compute total infected from previous week (across ages)
    real total_infected = 0;
    for (a in 1:A) {
      total_infected += I[a, t-1] + IV[a, t-1];
    }
    
    // Compute total population (across ages)
    real pop_total = 0;
    for (a in 1:A) {
      pop_total += N[a];
    }
    
    // Calculate force of infection for time t-1
    real phi = beta[t-1] * total_infected / pop_total;
    
    for (a in 1:A) {
      // Vaccination is applied only if:
      // (t >= delay) and (t <= delay+6) and if the age is targeted.
      real vacc_effect = ((t >= delay) && (t <= delay + 6) && (vaccine_target_age[a] == 1))
                          ? vaccination_rate : 0;
      
      if (a == 1) {
        S[a, t] = fmax(0, S[a, t-1] - phi * S[a, t-1]
                            - vacc_effect * S[a, t-1]
                            - r[a] * S[a, t-1]);
        I[a, t] = fmax(0, I[a, t-1] + phi * S[a, t-1]
                            - gamma * I[a, t-1]
                            - r[a] * I[a, t-1]);
        R[a, t] = fmax(0, R[a, t-1] + gamma * I[a, t-1]
                            - r[a] * R[a, t-1]);
        V[a, t] = fmax(0, V[a, t-1] + vacc_effect * S[a, t-1]
                            - (1 - VE_block) * phi * V[a, t-1]
                            - r[a] * V[a, t-1]);
        IV[a, t] = fmax(0, IV[a, t-1] + (1 - VE_block) * phi * V[a, t-1]
                            - gamma * IV[a, t-1]
                            - r[a] * IV[a, t-1]);
        RV[a, t] = fmax(0, RV[a, t-1] + gamma * IV[a, t-1]
                            - r[a] * RV[a, t-1]);
      } else {
        S[a, t] = fmax(0, S[a, t-1] - phi * S[a, t-1]
                           - vacc_effect * S[a, t-1]
                           - r[a] * S[a, t-1]
                           + r[a-1] * S[a-1, t-1]);
        I[a, t] = fmax(0, I[a, t-1] + phi * S[a, t-1]
                           - gamma * I[a, t-1]
                           - r[a] * I[a, t-1]
                           + r[a-1] * I[a-1, t-1]);
        R[a, t] = fmax(0, R[a, t-1] + gamma * I[a, t-1]
                           - r[a] * R[a, t-1]
                           + r[a-1] * R[a-1, t-1]);
        V[a, t] = fmax(0, V[a, t-1] + vacc_effect * S[a, t-1]
                           - (1 - VE_block) * phi * V[a, t-1]
                           - r[a] * V[a, t-1]
                           + r[a-1] * V[a-1, t-1]);
        IV[a, t] = fmax(0, IV[a, t-1] + (1 - VE_block) * phi * V[a, t-1]
                           - gamma * IV[a, t-1]
                           - r[a] * IV[a, t-1]
                           + r[a-1] * IV[a-1, t-1]);
        RV[a, t] = fmax(0, RV[a, t-1] + gamma * IV[a, t-1]
                           - r[a] * RV[a, t-1]
                           + r[a-1] * RV[a-1, t-1]);
      }
    }
  }
  
  // --- 3) Observation Model ---
  // Use the age-specific observed data to calibrate the model.
  for (t in 1:T) {
    for (a in 1:A) {
      real mu = rho * (I[a, t] + IV[a, t]);
      observed_cases_by_age[a, t] ~ neg_binomial_2(mu, shape);
    }
  }
}

generated quantities {
  // Posterior predictive of compartments
  real S_pred[A, T];
  real I_pred[A, T];
  real R_pred[A, T];
  real V_pred[A, T];
  real IV_pred[A, T];
  real RV_pred[A, T];
  
  // Predicted weekly cases (aggregated across ages)
  real pred_cases[T];
  // Aggregated infections (sum of I and IV across ages)
  real aggregated_infections[T];
  // Age-specific predicted cases
  real age_stratified_cases[A, T];
  real R_eff[T];
  real phi_pred[T];
  real N_pred[T];
  
  // Calculate total population ONCE at the start
  real pop_total = 0;
  for (a in 1:A)
    pop_total += N[a];
  
  // Copy over initial states
  for (a in 1:A) {
    S_pred[a, 1] = N[a] - I0[a];
    I_pred[a, 1] = I0[a];
    R_pred[a, 1] = 0;
    V_pred[a, 1] = 0;
    IV_pred[a, 1] = 0;
    RV_pred[a, 1] = 0;
    // Age-specific cases at t=1
    age_stratified_cases[a, 1] = rho * (I_pred[a, 1] + IV_pred[a, 1]);
  }
  
  // Calculate initial N_pred[1]
  {
    real current_pop = 0;
    for (a in 1:A) {
      current_pop += S_pred[a, 1] + I_pred[a, 1] + R_pred[a, 1] + V_pred[a, 1] + IV_pred[a, 1] + RV_pred[a, 1];
    }
    N_pred[1] = current_pop;
  }
  
  // Compute R_eff at t=1
  {
    real total_susceptible_unvacc = 0;
    real total_susceptible_vacc = 0;
    for (a in 1:A) {
      total_susceptible_unvacc += S_pred[a, 1];
      total_susceptible_vacc += V_pred[a, 1];
    }
    R_eff[1] = (beta[1] / gamma) * ((total_susceptible_unvacc + (1 - VE_block) * total_susceptible_vacc) / pop_total);
  }
  
  // Initial phi calculation for t=1
  {
    real total_infected_1 = 0;
    for (a in 1:A) {
      total_infected_1 += I_pred[a, 1] + IV_pred[a, 1];
    }
    phi_pred[1] = beta[1] * total_infected_1 / pop_total;
  }
  
  // Recompute dynamics for predictions
  for (t in 2:T) {
    // First, calculate total infected from previous timestep
    real total_infected_prev = 0;
    for (a in 1:A) {
      total_infected_prev += I_pred[a, t-1] + IV_pred[a, t-1];
    }
    
    // Calculate force of infection using total population
    phi_pred[t] = beta[t-1] * total_infected_prev / pop_total;
    real phi = phi_pred[t];
    // Vaccination window: assume vaccination from delay to delay+6
    int vaccine_end = delay + 6;
    for (a in 1:A) {
      real vaccination_effect = ((t >= delay) && (t <= vaccine_end) && (vaccine_target_age[a] == 1))
                                  ? vaccination_rate : 0;
      if (a == 1) {
        S_pred[a, t] = fmax(0, S_pred[a, t-1] - phi * S_pred[a, t-1]
                                 - vaccination_effect * S_pred[a, t-1]
                                 - r[a] * S_pred[a, t-1]);
        I_pred[a, t] = fmax(0, I_pred[a, t-1] + phi * S_pred[a, t-1]
                                 - gamma * I_pred[a, t-1]
                                 - r[a] * I_pred[a, t-1]);
        R_pred[a, t] = fmax(0, R_pred[a, t-1] + gamma * I_pred[a, t-1]
                                 - r[a] * R_pred[a, t-1]);
        V_pred[a, t] = fmax(0, V_pred[a, t-1] + vaccination_effect * S_pred[a, t-1]
                                 - (1 - VE_block) * phi * V_pred[a, t-1]
                                 - r[a] * V_pred[a, t-1]);
        IV_pred[a, t] = fmax(0, IV_pred[a, t-1] + (1 - VE_block) * phi * V_pred[a, t-1]
                                  - gamma * IV_pred[a, t-1]
                                  - r[a] * IV_pred[a, t-1]);
        RV_pred[a, t] = fmax(0, RV_pred[a, t-1] + gamma * IV_pred[a, t-1]
                                  - r[a] * RV_pred[a, t-1]);
      } else {
        S_pred[a, t] = fmax(0, S_pred[a, t-1] - phi * S_pred[a, t-1]
                                 - vaccination_effect * S_pred[a, t-1]
                                 - r[a] * S_pred[a, t-1]
                                 + r[a-1] * S_pred[a-1, t-1]);
        I_pred[a, t] = fmax(0, I_pred[a, t-1] + phi * S_pred[a, t-1]
                                 - gamma * I_pred[a, t-1]
                                 - r[a] * I_pred[a, t-1]
                                 + r[a-1] * I_pred[a-1, t-1]);
        R_pred[a, t] = fmax(0, R_pred[a, t-1] + gamma * I_pred[a, t-1]
                                 - r[a] * R_pred[a, t-1]
                                 + r[a-1] * R_pred[a-1, t-1]);
        V_pred[a, t] = fmax(0, V_pred[a, t-1] + vaccination_effect * S_pred[a, t-1]
                                 - (1 - VE_block) * phi * V_pred[a, t-1]
                                 - r[a] * V_pred[a, t-1]
                                 + r[a-1] * V_pred[a-1, t-1]);
        IV_pred[a, t] = fmax(0, IV_pred[a, t-1] + (1 - VE_block) * phi * V_pred[a, t-1]
                                  - gamma * IV_pred[a, t-1]
                                  - r[a] * IV_pred[a, t-1]
                                  + r[a-1] * IV_pred[a-1, t-1]);
        RV_pred[a, t] = fmax(0, RV_pred[a, t-1] + gamma * IV_pred[a, t-1]
                                  - r[a] * RV_pred[a, t-1]
                                  + r[a-1] * RV_pred[a-1, t-1]);
      }
      // Age-specific cases for time t
      age_stratified_cases[a, t] = rho * (I_pred[a, t] + IV_pred[a, t]);
    }
    
    // Calculate total population at time t
    {
      real current_pop = 0;
      for (a in 1:A) {
        current_pop += S_pred[a, t] + I_pred[a, t] + R_pred[a, t] +
                       V_pred[a, t] + IV_pred[a, t] + RV_pred[a, t];
      }
      N_pred[t] = current_pop;
    }
    
    // Calculate R_eff using population at t-1
    {
      real susceptible_unvacc_t = 0;
      real susceptible_vacc_t = 0;
      for (a in 1:A) {
        susceptible_unvacc_t += S_pred[a, t-1];
        susceptible_vacc_t += V_pred[a, t-1];
      }
      R_eff[t] = (beta[t] / gamma) * ((susceptible_unvacc_t + (1 - VE_block)*susceptible_vacc_t) / pop_total);
    }
  }
  
  // Summaries for each week
  for (t in 1:T) {
    real total_infections_t = 0;
    for (a in 1:A) {
      total_infections_t += (I_pred[a, t] + IV_pred[a, t]);
    }
    pred_cases[t] = rho * total_infections_t;
    aggregated_infections[t] = total_infections_t;
  }
}
