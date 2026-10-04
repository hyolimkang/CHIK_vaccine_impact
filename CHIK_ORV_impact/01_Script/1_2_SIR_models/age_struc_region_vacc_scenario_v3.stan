data {
  int<lower=1> T;                     // Number of weeks
  int<lower=1> A;                     // Number of age groups
  int<lower=1> nRegions;              // Number of regions
  int<lower=0> observed_cases[T];     // Weekly observed case data (aggregated)
  real<lower=0> N[A, nRegions];       // Population per age group and region
  real<lower=0> r[A];                 // Aging rate per age group
  real<lower=0> indexP[T, nRegions];  // Vector capacity index per region and week
  real<lower=0, upper=1> eta;         // Total desired vaccine coverage 
  real<lower=0> vaccination_rate;     // Weekly vaccination rate
  int<lower=1> delay;                 // Vaccine start time
  real<lower=0, upper=1> VE_block;    // Vaccine efficacy 
  int<lower=0, upper=1> vaccine_target_age[A]; // Which age groups get vaccine
}

parameters {
  // (1) We treat I0 as an *absolute count* of initially infected
  //     for each age group and region, but we constrain it between 1 and 100.
  //     - You may also do <lower=0, upper=100> if you allow zero initial infection
  real<lower=1, upper=100> I0[A, nRegions];

  // (2) Single time-varying baseline for beta
  real<lower=0, upper=6> base_beta[T];

  // (3) Region-specific multiplier to capture differences
  real<lower=0, upper=2> region_beta_mult[nRegions];

  // (4) Other parameters
  real<lower=0, upper=0.1> rho;   // Detection probability
  real<lower=17, upper=18> shape; // Overdispersion param
  real<lower=0, upper=1> gamma;   // Recovery rate
}

transformed parameters {
  real<lower=0> beta[T, nRegions];           
  real<lower=0, upper=1> indexP_normalised[T, nRegions];

  // -- (A) Normalize indexP
  for (reg in 1:nRegions) {
    real max_indexP = indexP[1, reg];
    real min_indexP = indexP[1, reg];
    for (t in 2:T) {
      if (indexP[t, reg] > max_indexP) max_indexP = indexP[t, reg];
      if (indexP[t, reg] < min_indexP) min_indexP = indexP[t, reg];
    }
    for (t in 1:T) {
      if (max_indexP > min_indexP) {
        indexP_normalised[t, reg] = 
          (indexP[t, reg] - min_indexP) / (max_indexP - min_indexP);
      } else {
        indexP_normalised[t, reg] = 0;
      }
    }
  }

  // -- (B) Compute region-specific beta[t, reg]
  //        = base_beta[t] * region_beta_mult[reg] * (0.5 + 0.5 * indexP_normalised[t, reg])
  for (t in 1:T) {
    for (reg in 1:nRegions) {
      beta[t, reg] = base_beta[t]
                     * region_beta_mult[reg]
                     * (0.5 + 0.5 * indexP_normalised[t, reg]);
    }
  }
}

model {
  int vaccine_end = delay + 6;
  real S[A, T, nRegions];             // Susceptible individuals
  real I[A, T, nRegions];             // Infectious individuals
  real R[A, T, nRegions];             // Recovered individuals
  real V[A, T, nRegions];             // Vaccinated individuals
  real IV[A, T, nRegions];            // Vaccinated but infected individuals (leaky vaccine)
  real RV[A, T, nRegions];            // Vaccinated and recovered
  real expected_cases[T];             // Expected weekly cases (aggregated)

  // Initial conditions
  for (a in 1:A) {
  for (reg in 1:nRegions) {
    S[a, 1, reg] = N[a, reg] - I0[a, reg];
    V[a, 1, reg] = 0;  // Always start with zero vaccinated individuals
    I[a, 1, reg] = I0[a, reg];
    R[a, 1, reg] = 0;
    IV[a, 1, reg] = 0;
    RV[a, 1, reg] = 0;
  }
}

  // (i) Priors for base_beta[t]
  for (t in 1:T) {
    // e.g. a lognormal prior centered near ~ e^-1 = 0.367
    base_beta[t] ~ lognormal(-1, 0.5);  
  }
  
  // (ii) Prior for region_beta_mult
  //      We expect region multipliers near 1.0, with some variation
  for (reg in 1:nRegions) {
    region_beta_mult[reg] ~ normal(1, 0.3); 
  }

  // (iii) Other priors
  rho   ~ beta(2, 18);       // or as you had
  gamma ~ normal(0.1, 0.05);
  shape ~ uniform(17, 18);

  // SIR dynamics with ageing
  for (t in 2:T) {
    for (reg in 1:nRegions) {
      real total_infected = 0;
      
      for (a in 1:A) {
        total_infected += I[a, t-1, reg] + IV[a, t-1, reg];
      }
      
      real phi = beta[t-1, reg] * total_infected / sum(N[, reg]);

      for (a in 1:A) {
        // Vaccination effect only after the delay and for targeted age groups
        real vaccination_effect = ((t >= delay) && (t <= vaccine_end) && (vaccine_target_age[a] == 1)) ? vaccination_rate : 0;

        if (a == 1) {
          // First age group (no ageing into it)
          S[a, t, reg] = fmax(0, S[a, t-1, reg] - phi * S[a, t-1, reg] - vaccination_effect * S[a, t-1, reg] - r[a] * S[a, t-1, reg]);
          I[a, t, reg] = fmax(0, I[a, t-1, reg] + phi * S[a, t-1, reg] - gamma * I[a, t-1, reg] - r[a] * I[a, t-1, reg]);
          R[a, t, reg] = fmax(0, R[a, t-1, reg] + gamma * I[a, t-1, reg] - r[a] * R[a, t-1, reg]);
          V[a, t, reg] = fmax(0, V[a, t-1, reg] + vaccination_effect * S[a, t-1, reg] - (1 - VE_block) * phi * V[a, t-1, reg] - r[a] * V[a, t-1, reg]);
          IV[a, t, reg] = fmax(0, IV[a, t-1, reg] + (1 - VE_block) * phi * V[a, t-1, reg] - gamma * IV[a, t-1, reg] - r[a] * IV[a, t-1, reg]);
          RV[a, t, reg] = fmax(0, RV[a, t-1, reg] + gamma * IV[a, t-1, reg] - r[a] * RV[a, t-1, reg]);
        } else {
          // Remaining age groups
          S[a, t, reg] = fmax(0, S[a, t-1, reg] - phi * S[a, t-1, reg] - vaccination_effect * S[a, t-1, reg] - r[a] * S[a, t-1, reg] + r[a-1] * S[a-1, t-1, reg]);
          I[a, t, reg] = fmax(0, I[a, t-1, reg] + phi * S[a, t-1, reg] - gamma * I[a, t-1, reg] - r[a] * I[a, t-1, reg] + r[a-1] * I[a-1, t-1, reg]);
          R[a, t, reg] = fmax(0, R[a, t-1, reg] + gamma * I[a, t-1, reg] - r[a] * R[a, t-1, reg] + r[a-1] * R[a-1, t-1, reg]);
          V[a, t, reg] = fmax(0, V[a, t-1, reg] + vaccination_effect * S[a, t-1, reg] - (1 - VE_block) * phi * V[a, t-1, reg] - r[a] * V[a, t-1, reg] + r[a-1] * V[a-1, t-1, reg]);
          IV[a, t, reg] = fmax(0, IV[a, t-1, reg] + (1 - VE_block) * phi * V[a, t-1, reg] - gamma * IV[a, t-1, reg] - r[a] * IV[a, t-1, reg] + r[a-1] * IV[a-1, t-1, reg]);
          RV[a, t, reg] = fmax(0, RV[a, t-1, reg] + gamma * IV[a, t-1, reg] - r[a] * RV[a, t-1, reg] + r[a-1] * RV[a-1, t-1, reg]);
        }
      }
    }
  }


  // Observation model: aggregate across regions and age groups
  for (t in 1:T) {
    real total_infections = 0;
    for (reg in 1:nRegions) {
      for (a in 1:A) {
        total_infections += I[a, t, reg] + IV[a, t, reg];
      }
    }
    expected_cases[t] = fmax(0, rho * total_infections);
    observed_cases[t] ~ neg_binomial_2(expected_cases[t], shape);
  }
}

generated quantities {
  real pred_cases[T];                         // Predicted weekly cases (national)
  real S_pred[A, T, nRegions];                // Predicted susceptible individuals
  real I_pred[A, T, nRegions];                // Predicted infectious individuals
  real R_pred[A, T, nRegions];                // Predicted recovered individuals
  real V_pred[A, T, nRegions];                // Predicted vaccinated individuals
  real IV_pred[A, T, nRegions];               // Predicted vaccinated and infected individuals
  real RV_pred[A, T, nRegions];               // Predicted vaccinated, infected, and recovered individuals
  real age_stratified_cases[A, T, nRegions];  // Age-specific predicted cases by region
  real region_specific_cases[T, nRegions];    // Region-specific total cases
  real aggregated_infections[T];              // National total infections
  real R_eff[T, nRegions];                    // Effective reproduction number
  int vaccine_end = delay + 6;
  
  // Initialize
  for (a in 1:A) {
  for (reg in 1:nRegions) {
    S_pred[a, 1, reg] = N[a, reg] - I0[a, reg];
    V_pred[a, 1, reg] = 0;  // No initial vaccinated individuals
    I_pred[a, 1, reg] = I0[a, reg];
    R_pred[a, 1, reg] = 0;
    IV_pred[a, 1, reg] = 0;
    RV_pred[a, 1, reg] = 0;
    age_stratified_cases[a, 1, reg] = rho * (I_pred[a, 1, reg] + IV_pred[a, 1, reg]);
  }
}
  // Initialize R_eff for t=1
  for (reg in 1:nRegions) {
    real total_susceptible_unvaccinated = 0;
    real total_susceptible_vaccinated = 0;
    for (a in 1:A) {
      total_susceptible_unvaccinated += S_pred[a, 1, reg];
      total_susceptible_vaccinated   += V_pred[a, 1, reg];
    }
    R_eff[1, reg] = (beta[1, reg] / gamma) 
                   * ((total_susceptible_unvaccinated + (1 - VE_block) * total_susceptible_vaccinated) 
                      / sum(N[, reg]));
  }
  
  // Recompute dynamics for predictions
  for (t in 2:T) {
    for (reg in 1:nRegions) {
      // Compute total infections across all age groups
      real total_infected = 0;
      real total_susceptible_unvaccinated = 0;
      real total_susceptible_vaccinated = 0;
      
      for (a in 1:A) {
        total_infected += I_pred[a, t-1, reg] + IV_pred[a, t-1, reg];
        total_susceptible_unvaccinated += S_pred[a, t-1, reg];
        total_susceptible_vaccinated += V_pred[a, t-1, reg];
      }
      
      real phi = beta[t-1, reg] * total_infected / sum(N[, reg]);

      for (a in 1:A) {
        // Vaccination effect only after the delay and for targeted age groups
        real vaccination_effect = ((t >= delay) && (t <= vaccine_end) &&  (vaccine_target_age[a] == 1)) ? vaccination_rate : 0;

        if (a == 1) {
          // First age group
          S_pred[a, t, reg] = fmax(0, S_pred[a, t-1, reg] - phi * S_pred[a, t-1, reg] - vaccination_effect * S_pred[a, t-1, reg] - r[a] * S_pred[a, t-1, reg]);
          I_pred[a, t, reg] = fmax(0, I_pred[a, t-1, reg] + phi * S_pred[a, t-1, reg] - gamma * I_pred[a, t-1, reg] - r[a] * I_pred[a, t-1, reg]);
          R_pred[a, t, reg] = fmax(0, R_pred[a, t-1, reg] + gamma * I_pred[a, t-1, reg] - r[a] * R_pred[a, t-1, reg]);
          V_pred[a, t, reg] = fmax(0, V_pred[a, t-1, reg] + vaccination_effect * S_pred[a, t-1, reg] - (1 - VE_block) * phi * V_pred[a, t-1, reg] - r[a] * V_pred[a, t-1, reg]);
          IV_pred[a, t, reg] = fmax(0, IV_pred[a, t-1, reg] + (1 - VE_block) * phi * V_pred[a, t-1, reg] - gamma * IV_pred[a, t-1, reg] - r[a] * IV_pred[a, t-1, reg]);
          RV_pred[a, t, reg] = fmax(0, RV_pred[a, t-1, reg] + gamma * IV_pred[a, t-1, reg] - r[a] * RV_pred[a, t-1, reg]);
        } else {
          // Other age groups
          S_pred[a, t, reg] = fmax(0, S_pred[a, t-1, reg] - phi * S_pred[a, t-1, reg] - vaccination_effect * S_pred[a, t-1, reg] - r[a] * S_pred[a, t-1, reg] + r[a-1] * S_pred[a-1, t-1, reg]);
          I_pred[a, t, reg] = fmax(0, I_pred[a, t-1, reg] + phi * S_pred[a, t-1, reg] - gamma * I_pred[a, t-1, reg] - r[a] * I_pred[a, t-1, reg] + r[a-1] * I_pred[a-1, t-1, reg]);
          R_pred[a, t, reg] = fmax(0, R_pred[a, t-1, reg] + gamma * I_pred[a, t-1, reg] - r[a] * R_pred[a, t-1, reg] + r[a-1] * R_pred[a-1, t-1, reg]);
          V_pred[a, t, reg] = fmax(0, V_pred[a, t-1, reg] + vaccination_effect * S_pred[a, t-1, reg] - (1 - VE_block) * phi * V_pred[a, t-1, reg] - r[a] * V_pred[a, t-1, reg] + r[a-1] * V_pred[a-1, t-1, reg]);
          IV_pred[a, t, reg] = fmax(0, IV_pred[a, t-1, reg] + (1 - VE_block) * phi * V_pred[a, t-1, reg] - gamma * IV_pred[a, t-1, reg] + r[a-1] * IV_pred[a-1, t-1, reg]- r[a] * IV_pred[a, t-1, reg]);
          RV_pred[a, t, reg] = fmax(0, RV_pred[a, t-1, reg] + gamma * IV_pred[a, t-1, reg] + r[a-1] * RV_pred[a-1, t-1, reg]- r[a] * RV_pred[a, t-1, reg]);
        }

        age_stratified_cases[a, t, reg] = rho * (I_pred[a, t, reg] + IV_pred[a, t, reg]);
      }

      R_eff[t, reg] = (beta[t, reg] / gamma) * ((total_susceptible_unvaccinated + (1 - VE_block) * total_susceptible_vaccinated) / sum(N[, reg]));
    }
  }

  // Aggregate predictions
  for (t in 1:T) {
    real total_infections = 0;

    // Regional totals
    for (reg in 1:nRegions) {
      real region_total = 0;
      for (a in 1:A) {
        region_total += (I_pred[a, t, reg] + IV_pred[a, t, reg]);
      }
      region_specific_cases[t, reg] = rho * region_total;
      total_infections += region_total;
    }

    // National totals
    pred_cases[t] = rho * total_infections;
    aggregated_infections[t] = total_infections;
  }
}

