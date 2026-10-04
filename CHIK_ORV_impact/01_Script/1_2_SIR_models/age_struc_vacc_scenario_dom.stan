data {
  int<lower=1> T;                    // Number of weeks
  int<lower=1> A;                    // Number of age groups
  
  int<lower=0> observed_cases[T];    // Weekly observed cases
  
  real<lower=0> N[A];                // Population per age group
  real<lower=0> r[A];                // Aging rate per age group
  
  real<lower=0> indexP[T];           // Vector capacity index for each week
  
  real<lower=0, upper=1> eta;        // Desired vaccine coverage (may be unused directly)
  real<lower=0> vaccination_rate;    // Weekly vaccination rate
  int<lower=1> delay;                // First week of vaccination
  real<lower=0, upper=1> VE_block;   // Vaccine efficacy
  int<lower=0, upper=1> vaccine_target_age[A];  // Which ages get vaccinated (0 or 1)
}

parameters {
  // Initial infected (absolute count) in each age group, forced between 1 and 100
  // Adjust if you want to allow zero or larger upper bound
  real<lower=1, upper=1000> I0[A];
  
  // Time-varying baseline transmission rate
  real<lower=0, upper=5> base_beta[T];
  
  // Other parameters
  real<lower=0,   upper=1>   rho;    // Detection probability
  real<lower=30,  upper=31>  shape;  // Overdispersion
  real<lower=0,   upper=1>   gamma;  // Recovery rate
}

transformed parameters {
  // Normalized indexP in [0,1]
  real<lower=0, upper=1> indexP_normalised[T];
  
  // Final transmission rate each week
  real<lower=0> beta[T];

  // 1) Normalize indexP
  {
    real max_index = indexP[1];
    real min_index = indexP[1];
    for (t in 2:T) {
      if (indexP[t] > max_index) max_index = indexP[t];
      if (indexP[t] < min_index) min_index = indexP[t];
    }
    for (t in 1:T) {
      if (max_index > min_index) {
        indexP_normalised[t] = (indexP[t] - min_index) / (max_index - min_index);
      } else {
        indexP_normalised[t] = 0;  // If all indexP are the same
      }
    }
  }
  
  // 2) Compute final beta[t] using normalized indexP
  for (t in 1:T) {
    beta[t] = base_beta[t] * (1 + 0.5 * indexP_normalised[t]);
  }
}

model {
  int vaccine_end = delay + 6;  // Example: vaccinate from 'delay' up to 'delay+24'
  
  // State variables (age x time)
  real S[A, T];
  real I[A, T];
  real R[A, T];
  real V[A, T];
  real IV[A, T];
  real RV[A, T];
  
  // Expected weekly cases
  real expected_cases[T];
  
  //------------------
  // 1) Priors
  //------------------
  // base_beta[t]
  for (t in 1:T) {
    base_beta[t] ~ lognormal(0.5, 0.5);  // Example prior
  }
  
  // detection
  rho   ~ beta(5, 10);
  
  // recovery
  gamma ~ normal(0.1, 0.05);
  
  // overdispersion
  shape ~ uniform(30, 31);  // originally (17,18)
  
  //------------------
  // 2) Initial conditions at t=1
  //------------------
  for (a in 1:A) {
    S[a, 1] = N[a] - I0[a];
    I[a, 1] = I0[a];
    R[a, 1] = 0;
    V[a, 1] = 0;
    IV[a, 1] = 0;
    RV[a, 1] = 0;
  }
  
  //------------------
  // 3) SIRV dynamics
  //------------------
  for (t in 2:T) {
    // Sum total infected across ages (previous week)
    real total_infected = 0;
    for (a in 1:A) {
      total_infected += I[a, t-1] + IV[a, t-1];
    }
    
    // Force of infection = beta[t-1] * (total_infected / sum(N)) 
    // (One region => sum(N[a]) is total population)
    real pop_total = 0;
    for (a in 1:A) pop_total += N[a];
    
    real phi = beta[t-1] * total_infected / pop_total;
    
    for (a in 1:A) {
      // Only vaccinate targeted ages, after 'delay' week and before 'vaccine_end'
      real vaccination_effect = ((t >= delay) && (t <= vaccine_end) && (vaccine_target_age[a] == 1))
                                ? vaccination_rate 
                                : 0;
      
      if (a == 1) {
        // Age group 1 (no in-flow from previous age)
        S[a, t]  = fmax(0, S[a, t-1] - phi * S[a, t-1] 
                                  - vaccination_effect * S[a, t-1] 
                                  - r[a] * S[a, t-1]);
        I[a, t]  = fmax(0, I[a, t-1] + phi * S[a, t-1] 
                                  - gamma * I[a, t-1] 
                                  - r[a] * I[a, t-1]);
        R[a, t]  = fmax(0, R[a, t-1] + gamma * I[a, t-1] 
                                  - r[a] * R[a, t-1]);
        V[a, t]  = fmax(0, V[a, t-1] + vaccination_effect * S[a, t-1]
                                  - (1 - VE_block) * phi * V[a, t-1]
                                  - r[a] * V[a, t-1]);
        IV[a, t] = fmax(0, IV[a, t-1] + (1 - VE_block) * phi * V[a, t-1]
                                  - gamma * IV[a, t-1]
                                  - r[a] * IV[a, t-1]);
        RV[a, t] = fmax(0, RV[a, t-1] + gamma * IV[a, t-1]
                                  - r[a] * RV[a, t-1]);
        
      } else {
        // Age group 2..A (in-flow from previous age group due to aging)
        S[a, t]  = fmax(0, S[a, t-1] - phi * S[a, t-1] 
                                  - vaccination_effect * S[a, t-1] 
                                  - r[a] * S[a, t-1]
                                  + r[a-1] * S[a-1, t-1]);
        I[a, t]  = fmax(0, I[a, t-1] + phi * S[a, t-1] 
                                  - gamma * I[a, t-1]
                                  - r[a] * I[a, t-1]
                                  + r[a-1] * I[a-1, t-1]);
        R[a, t]  = fmax(0, R[a, t-1] + gamma * I[a, t-1]
                                  - r[a] * R[a, t-1]
                                  + r[a-1] * R[a-1, t-1]);
        V[a, t]  = fmax(0, V[a, t-1] + vaccination_effect * S[a, t-1]
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
  
  //------------------
  // 4) Observation model
  //------------------
  for (t in 1:T) {
    // total infected across all ages in week t
    real total_infections_t = 0;
    for (a in 1:A) {
      total_infections_t += I[a, t] + IV[a, t];
    }
    expected_cases[t] = rho * total_infections_t;
    
    // Negative binomial observation
    observed_cases[t] ~ neg_binomial_2(expected_cases[t], shape);
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
  
  real pred_cases[T];             // Predicted weekly cases
  real aggregated_infections[T];  // Sum of I and IV across ages
  real age_stratified_cases[A, T]; // Age-specific cases
  real R_eff[T];
  real phi_pred[T];
  real N_pred[T];
  
  // Calculate total population ONCE at the start
  real pop_total = 0;
  for (a in 1:A) pop_total += N[a];
  
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
      current_pop += S_pred[a, 1] + I_pred[a, 1] + R_pred[a, 1] + 
                    V_pred[a, 1] + IV_pred[a, 1] + RV_pred[a, 1];
    }
    N_pred[1] = current_pop;
  }
  
  // Compute R_eff at t=1
  {
    real total_susceptible_unvacc = 0;
    real total_susceptible_vacc   = 0;
    for (a in 1:A) {
      total_susceptible_unvacc += S_pred[a, 1];
      total_susceptible_vacc   += V_pred[a, 1];
    }
    R_eff[1] = (beta[1] / gamma) 
                * ( (total_susceptible_unvacc 
                     + (1 - VE_block)*total_susceptible_vacc)
                    / pop_total );
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
    // First calculate total infected from previous timestep
    real total_infected_prev = 0;
    for (a in 1:A) {
      total_infected_prev += I_pred[a, t-1] + IV_pred[a, t-1];
    }
    
    // Calculate force of infection using correct total population
    phi_pred[t] = beta[t-1] * total_infected_prev / pop_total;
    real phi = phi_pred[t];
    // Vaccination window
    int vaccine_end = delay + 6; // same logic as in model
    for (a in 1:A) {
      real vaccination_effect = ((t >= delay) && (t <= vaccine_end) && (vaccine_target_age[a] == 1))
                                ? vaccination_rate
                                : 0;
      
      if (a == 1) {
        // Age group 1
        S_pred[a, t]  = fmax(0, S_pred[a, t-1] - phi * S_pred[a, t-1]
                                       - vaccination_effect * S_pred[a, t-1]
                                       - r[a] * S_pred[a, t-1]);
        I_pred[a, t]  = fmax(0, I_pred[a, t-1] + phi * S_pred[a, t-1]
                                       - gamma * I_pred[a, t-1]
                                       - r[a] * I_pred[a, t-1]);
        R_pred[a, t]  = fmax(0, R_pred[a, t-1] + gamma * I_pred[a, t-1]
                                       - r[a] * R_pred[a, t-1]);
        V_pred[a, t]  = fmax(0, V_pred[a, t-1] + vaccination_effect * S_pred[a, t-1]
                                       - (1 - VE_block)*phi * V_pred[a, t-1]
                                       - r[a] * V_pred[a, t-1]);
        IV_pred[a, t] = fmax(0, IV_pred[a, t-1] + (1 - VE_block)*phi * V_pred[a, t-1]
                                       - gamma * IV_pred[a, t-1]
                                       - r[a] * IV_pred[a, t-1]);
        RV_pred[a, t] = fmax(0, RV_pred[a, t-1] + gamma * IV_pred[a, t-1]
                                       - r[a] * RV_pred[a, t-1]);
      } else {
        // Age group 2..A
        S_pred[a, t]  = fmax(0, S_pred[a, t-1] - phi * S_pred[a, t-1]
                                       - vaccination_effect * S_pred[a, t-1]
                                       - r[a] * S_pred[a, t-1]
                                       + r[a-1] * S_pred[a-1, t-1]);
        I_pred[a, t]  = fmax(0, I_pred[a, t-1] + phi * S_pred[a, t-1]
                                       - gamma * I_pred[a, t-1]
                                       - r[a] * I_pred[a, t-1]
                                       + r[a-1] * I_pred[a-1, t-1]);
        R_pred[a, t]  = fmax(0, R_pred[a, t-1] + gamma * I_pred[a, t-1]
                                       - r[a] * R_pred[a, t-1]
                                       + r[a-1] * R_pred[a-1, t-1]);
        V_pred[a, t]  = fmax(0, V_pred[a, t-1] + vaccination_effect * S_pred[a, t-1]
                                       - (1 - VE_block)*phi * V_pred[a, t-1]
                                       - r[a] * V_pred[a, t-1]
                                       + r[a-1] * V_pred[a-1, t-1]);
        IV_pred[a, t] = fmax(0, IV_pred[a, t-1] + (1 - VE_block)*phi * V_pred[a, t-1]
                                       - gamma * IV_pred[a, t-1]
                                       - r[a] * IV_pred[a, t-1]
                                       + r[a-1] * IV_pred[a-1, t-1]);
        RV_pred[a, t] = fmax(0, RV_pred[a, t-1] + gamma * IV_pred[a, t-1]
                                       - r[a] * RV_pred[a, t-1]
                                       + r[a-1] * RV_pred[a-1, t-1]);
      }
      // Age-specific cases
      age_stratified_cases[a, t] = rho * (I_pred[a, t] + IV_pred[a, t]);
    }
    
    // Calculate N_pred for current timestep
    {
      real current_pop = 0;
      for (a in 1:A) {
        current_pop += S_pred[a, t] + I_pred[a, t] + R_pred[a, t] + 
                      V_pred[a, t] + IV_pred[a, t] + RV_pred[a, t];
      }
      N_pred[t] = current_pop;
    }
    
    // Calculate R_eff using the same pop_total
    {
      real susceptible_unvacc_t = 0;
      real susceptible_vacc_t = 0;
      for (a in 1:A) {
        susceptible_unvacc_t += S_pred[a, t-1];
        susceptible_vacc_t += V_pred[a, t-1];
      }
      R_eff[t] = (beta[t] / gamma) 
                * ((susceptible_unvacc_t + (1 - VE_block)*susceptible_vacc_t)
                   / pop_total);
    }
  }
  
  // Summaries for each week
  for (t in 1:T) {
    real total_infections_t = 0;
    for (a in 1:A) {
      total_infections_t += (I_pred[a, t] + IV_pred[a, t]);
    }
    pred_cases[t]           = rho * total_infections_t; 
    aggregated_infections[t] = total_infections_t;
  }
}
