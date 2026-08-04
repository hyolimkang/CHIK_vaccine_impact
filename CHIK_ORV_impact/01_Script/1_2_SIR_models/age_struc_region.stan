data {
  int<lower=1> T;                     // Number of weeks (52)
  int<lower=1> A;                     // Number of age groups
  int<lower=1> nRegions;              // Number of regions (27)
  int<lower=0> observed_cases[T];     // Weekly observed case data (aggregated)
  real<lower=0> N[A, nRegions];       // Population per age group and region
  real<lower=0> r[A];                 // Ageing rate per age group
  real<lower=0> indexP[T, nRegions];  // Vector capacity index per region and week
}

parameters {
  real<lower=1, upper=100> I0[A, nRegions];  // Initial infectious individuals per age group and region
  real<lower=0, upper=6> base_beta[T];       // Baseline transmission parameter
  real<lower=0, upper=0.1> rho;              // Detection probability
  real<lower=17, upper=18> shape;            // Overdispersion parameter
  real<lower=0, upper=1> gamma;              // Recovery rate
}

transformed parameters {
  real<lower=0> beta[T, nRegions];            // Region-specific transmission parameters
  real<lower=0, upper=1> indexP_normalised[T, nRegions];  // Normalized vector capacity index

  // Normalize indexP to [0, 1] for each region
  for (reg in 1:nRegions) {
    real max_indexP = indexP[1, reg];  // Initialize maximum
    real min_indexP = indexP[1, reg];  // Initialize minimum

    // Find the maximum and minimum indexP for the region
    for (t in 2:T) {
      max_indexP = fmax(max_indexP, indexP[t, reg]);
      min_indexP = fmin(min_indexP, indexP[t, reg]);
    }

    // Normalize indexP for this region to the range [0, 1]
    for (t in 1:T) {
      if (max_indexP > min_indexP) { // Avoid division by zero
        indexP_normalised[t, reg] = (indexP[t, reg] - min_indexP) / (max_indexP - min_indexP);
      } else {
        indexP_normalised[t, reg] = 0; // In case all values are the same
      }
    }
  }

  // Scale beta using normalized indexP
  for (t in 1:T) {
    for (reg in 1:nRegions) {
      beta[t, reg] = base_beta[t] * (0.5 + 0.5 * indexP_normalised[t, reg]);
    }
  }
}

model {
  real S[A, T, nRegions];              // Susceptible individuals
  real I[A, T, nRegions];              // Infectious individuals
  real R[A, T, nRegions];              // Recovered individuals
  real expected_cases[T];              // Expected weekly cases (aggregated)

  // Initial conditions
  for (a in 1:A) {
    for (reg in 1:nRegions) {
      S[a, 1, reg] = N[a, reg] - I0[a, reg];
      I[a, 1, reg] = I0[a, reg];
      R[a, 1, reg] = 0;
    }
  }

  // Priors
  for (t in 1:T) {
    base_beta[t] ~ lognormal(-1, 0.5);
  }
  rho ~ beta(2, 18);
  gamma ~ normal(0.1, 0.05);
  shape ~ uniform(17, 18);

  // SIR dynamics with ageing
  for (t in 2:T) {
    for (reg in 1:nRegions) {
      real phi = beta[t-1, reg] * sum(I[, t-1, reg]) / sum(N[, reg]);  // Region-specific FOI

      for (a in 1:A) {
        if (a == 1) {
          // First age group (no ageing into it)
          S[a, t, reg] = fmax(0, S[a, t-1, reg] - phi * S[a, t-1, reg] - r[a] * S[a, t-1, reg]);
          I[a, t, reg] = fmax(0, I[a, t-1, reg] + phi * S[a, t-1, reg] - gamma * I[a, t-1, reg] - r[a] * I[a, t-1, reg]);
          R[a, t, reg] = fmax(0, R[a, t-1, reg] + gamma * I[a, t-1, reg] - r[a] * R[a, t-1, reg]);
        } else {
          // Remaining age groups
          S[a, t, reg] = fmax(0, S[a, t-1, reg] - phi * S[a, t-1, reg] - r[a] * S[a, t-1, reg] + r[a-1] * S[a-1, t-1, reg]);
          I[a, t, reg] = fmax(0, I[a, t-1, reg] + phi * S[a, t-1, reg] - gamma * I[a, t-1, reg] - r[a] * I[a, t-1, reg] + r[a-1] * I[a-1, t-1, reg]);
          R[a, t, reg] = fmax(0, R[a, t-1, reg] + gamma * I[a, t-1, reg] - r[a] * R[a, t-1, reg] + r[a-1] * R[a-1, t-1, reg]);
        }
      }
    }
  }

  // Observation model: aggregate across regions and age groups
  for (t in 1:T) {
    real total_infections = 0;
    for (reg in 1:nRegions) {
      for (a in 1:A) {
        total_infections += I[a, t, reg];
      }
    }
    expected_cases[t] = fmax(0, rho * total_infections);
    observed_cases[t] ~ neg_binomial_2(expected_cases[t], shape);
  }
}

generated quantities {
  real pred_cases[T];                         // Predicted weekly cases (national)
  real S_pred[A, T, nRegions];               // Predicted susceptible individuals
  real I_pred[A, T, nRegions];               // Predicted infectious individuals
  real R_pred[A, T, nRegions];               // Predicted recovered individuals
  real age_stratified_cases[A, T, nRegions]; // Age-specific predicted cases by region
  real region_specific_cases[T, nRegions];   // Region-specific total cases
  real aggregated_infections[T];             // National total infections
  
  // Initial conditions
  for (a in 1:A) {
    for (reg in 1:nRegions) {
      S_pred[a, 1, reg] = N[a, reg] - I0[a, reg];
      I_pred[a, 1, reg] = I0[a, reg];
      R_pred[a, 1, reg] = 0;
      age_stratified_cases[a, 1, reg] = rho * I_pred[a, 1, reg];
    }
  }

  // Recompute SIR dynamics for predictions
  for (t in 2:T) {
    for (reg in 1:nRegions) {
      real phi = beta[t-1, reg] * sum(I_pred[, t-1, reg]) / sum(N[, reg]);

      for (a in 1:A) {
        if (a == 1) {
          S_pred[a, t, reg] = fmax(0, S_pred[a, t-1, reg] - phi * S_pred[a, t-1, reg] - r[a] * S_pred[a, t-1, reg]);
          I_pred[a, t, reg] = fmax(0, I_pred[a, t-1, reg] + phi * S_pred[a, t-1, reg] - gamma * I_pred[a, t-1, reg] - r[a] * I_pred[a, t-1, reg]);
          R_pred[a, t, reg] = fmax(0, R_pred[a, t-1, reg] + gamma * I_pred[a, t-1, reg] - r[a] * R_pred[a, t-1, reg]);
        } else {
          S_pred[a, t, reg] = fmax(0, S_pred[a, t-1, reg] - phi * S_pred[a, t-1, reg] - r[a] * S_pred[a, t-1, reg] + r[a-1] * S_pred[a-1, t-1, reg]);
          I_pred[a, t, reg] = fmax(0, I_pred[a, t-1, reg] + phi * S_pred[a, t-1, reg] - gamma * I_pred[a, t-1, reg] - r[a] * I_pred[a, t-1, reg] + r[a-1] * I_pred[a-1, t-1, reg]);
          R_pred[a, t, reg] = fmax(0, R_pred[a, t-1, reg] + gamma * I_pred[a, t-1, reg] - r[a] * R_pred[a, t-1, reg] + r[a-1] * R_pred[a-1, t-1, reg]);
        }
        age_stratified_cases[a, t, reg] = rho * I_pred[a, t, reg];
      }
    }
  }

  // Generate predictions at different levels of aggregation
  for (t in 1:T) {
    real total_infections = 0;
    
    // Regional totals
    for (reg in 1:nRegions) {
      real region_total = 0;
      for (a in 1:A) {
        region_total += I_pred[a, t, reg];
      }
      region_specific_cases[t, reg] = rho * region_total;
      total_infections += region_total;
    }
    
    // National totals
    pred_cases[t] = rho * total_infections;
    aggregated_infections[t] = total_infections;
  }
}
