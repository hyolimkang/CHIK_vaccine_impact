//--------------------------------------------------------------------
//  AGE-STRUCTURED SEIR (NO VACCINATION, WITH E COMPARTMENT)
//--------------------------------------------------------------------
data {
  int<lower=1> T;                           // Number of weeks
  int<lower=1> A;                           // Number of age groups
  int<lower=0> observed_cases_by_age[A, T]; // Age-specific observed cases
  real<lower=0> N[A];                       // Population per age group
  real<lower=0> r[A];                       // Aging rate per age group
  real<lower=0> indexP[T];                  // Capacity index for each week
  real<lower=0> prior_I0[A];
  real<lower=0> prior_sd_I0;
  real<lower=0, upper=1> sero[A];
}

parameters {
  real<lower=10, upper=500> I0[A];
  real<lower=0, upper=6> base_beta[T];
  real<lower=0, upper=1> rho;
  real<lower=15, upper=20> shape;
  real<lower=0, upper=1> gamma;
  real<lower=0, upper=1> sigma;  // Rate from E to I (1/incubation period)
}

transformed parameters {
  real<lower=0, upper=1> indexP_normalised[T];
  real<lower=0> beta[T];

  {
    real max_index = indexP[1];
    real min_index = indexP[1];
    for (t in 2:T) {
      if (indexP[t] > max_index) max_index = indexP[t];
      if (indexP[t] < min_index) min_index = indexP[t];
    }
    for (t in 1:T) {
      indexP_normalised[t] = (max_index > min_index) ?
        (indexP[t] - min_index) / (max_index - min_index) : 0;
    }
  }

  for (t in 1:T) {
    beta[t] = base_beta[t] * (1 + 0.5 * indexP_normalised[t]);
  }
}

model {
  real S[A, T];
  real E[A, T];
  real I[A, T];
  real R[A, T];

  for (t in 1:T)
    base_beta[t] ~ lognormal(-1, 0.5);
  rho ~ beta(20, 60);
  gamma ~ normal(0.67, 0.02);
  sigma ~ normal(1.0, 0.2);  // Mean latent = 1 week (i.e. 2-6 days)
  shape ~ uniform(15, 20);
  for (a in 1:A)
    I0[a] ~ normal(prior_I0[a], prior_sd_I0) T[10, 500];

  for (a in 1:A) {
    R[a,1] = sero[a] * N[a];
    I[a,1] = I0[a];
    E[a,1] = 0;
    S[a,1] = N[a] - I0[a] - R[a,1];
  }

  for (t in 2:T) {
    real total_infected = 0;
    for (a in 1:A)
      total_infected += I[a,t-1];
    real phi = beta[t-1] * total_infected / sum(N);

    for (a in 1:A) {
      if (a == 1) {
        S[a,t] = fmax(0, S[a,t-1] - phi*S[a,t-1] - r[a]*S[a,t-1]);
        E[a,t] = fmax(0, E[a,t-1] + phi*S[a,t-1] - sigma*E[a,t-1] - r[a]*E[a,t-1]);
        I[a,t] = fmax(0, I[a,t-1] + sigma*E[a,t-1] - gamma*I[a,t-1] - r[a]*I[a,t-1]);
        R[a,t] = fmax(0, R[a,t-1] + gamma*I[a,t-1] - r[a]*R[a,t-1]);
      } else {
        S[a,t] = fmax(0, S[a,t-1] - phi*S[a,t-1] - r[a]*S[a,t-1] + r[a-1]*S[a-1,t-1]);
        E[a,t] = fmax(0, E[a,t-1] + phi*S[a,t-1] - sigma*E[a,t-1] - r[a]*E[a,t-1] + r[a-1]*E[a-1,t-1]);
        I[a,t] = fmax(0, I[a,t-1] + sigma*E[a,t-1] - gamma*I[a,t-1] - r[a]*I[a,t-1] + r[a-1]*I[a-1,t-1]);
        R[a,t] = fmax(0, R[a,t-1] + gamma*I[a,t-1] - r[a]*R[a,t-1] + r[a-1]*R[a-1,t-1]);
      }
    }
  }

  for (t in 1:T) {
    for (a in 1:A) {
      real mu = rho * I[a,t];
      observed_cases_by_age[a,t] ~ neg_binomial_2(mu, shape);
    }
  }
}

generated quantities {
  real S_pred[A, T];
  real E_pred[A, T];
  real I_pred[A, T];
  real R_pred[A, T];

  real pred_cases[T];
  real aggregated_infections[T];
  real age_stratified_cases[A, T];
  real R_eff[T];
  real phi_pred[T];
  real N_pred[T];

  real pop_total = 0;
  for (a in 1:A) pop_total += N[a];

  for (a in 1:A) {
    S_pred[a,1] = N[a] - I0[a] - sero[a] * N[a];
    E_pred[a,1] = 0;
    I_pred[a,1] = I0[a];
    R_pred[a,1] = sero[a] * N[a];
    age_stratified_cases[a,1] = rho * I_pred[a,1];
  }

  {
    real current_pop = 0;
    for (a in 1:A)
      current_pop += S_pred[a,1] + E_pred[a,1] + I_pred[a,1] + R_pred[a,1];
    N_pred[1] = current_pop;
  }

  {
    real susceptible_total = 0;
    for (a in 1:A) susceptible_total += S_pred[a,1];
    R_eff[1] = (beta[1] / gamma) * (susceptible_total / pop_total);
  }

  {
    real total_infected_1 = 0;
    for (a in 1:A)
      total_infected_1 += I_pred[a,1];
    phi_pred[1] = beta[1] * total_infected_1 / pop_total;
  }

  for (t in 2:T) {
    real total_infected_prev = 0;
    for (a in 1:A)
      total_infected_prev += I_pred[a,t-1];

    phi_pred[t] = beta[t-1] * total_infected_prev / pop_total;
    real phi = phi_pred[t];

    for (a in 1:A) {
      if (a == 1) {
        S_pred[a,t] = fmax(0, S_pred[a,t-1] - phi * S_pred[a,t-1] - r[a] * S_pred[a,t-1]);
        E_pred[a,t] = fmax(0, E_pred[a,t-1] + phi * S_pred[a,t-1] - sigma * E_pred[a,t-1] - r[a] * E_pred[a,t-1]);
        I_pred[a,t] = fmax(0, I_pred[a,t-1] + sigma * E_pred[a,t-1] - gamma * I_pred[a,t-1] - r[a] * I_pred[a,t-1]);
        R_pred[a,t] = fmax(0, R_pred[a,t-1] + gamma * I_pred[a,t-1] - r[a] * R_pred[a,t-1]);
      } else {
        S_pred[a,t] = fmax(0, S_pred[a,t-1] - phi * S_pred[a,t-1] - r[a] * S_pred[a,t-1] + r[a-1] * S_pred[a-1,t-1]);
        E_pred[a,t] = fmax(0, E_pred[a,t-1] + phi * S_pred[a,t-1] - sigma * E_pred[a,t-1] - r[a] * E_pred[a,t-1] + r[a-1] * E_pred[a-1,t-1]);
        I_pred[a,t] = fmax(0, I_pred[a,t-1] + sigma * E_pred[a,t-1] - gamma * I_pred[a,t-1] - r[a] * I_pred[a,t-1] + r[a-1] * I_pred[a-1,t-1]);
        R_pred[a,t] = fmax(0, R_pred[a,t-1] + gamma * I_pred[a,t-1] - r[a] * R_pred[a,t-1] + r[a-1] * R_pred[a-1,t-1]);
      }
      age_stratified_cases[a,t] = rho * I_pred[a,t];
    }

    {
      real current_pop = 0;
      for (a in 1:A)
        current_pop += S_pred[a,t] + E_pred[a,t] + I_pred[a,t] + R_pred[a,t];
      N_pred[t] = current_pop;
    }

    {
      real susceptible_t = 0;
      for (a in 1:A)
        susceptible_t += S_pred[a,t-1];
      R_eff[t] = (beta[t] / gamma) * (susceptible_t / pop_total);
    }
  }

  for (t in 1:T) {
    real total_infections_t = 0;
    for (a in 1:A)
      total_infections_t += I_pred[a,t];
    pred_cases[t] = rho * total_infections_t;
    aggregated_infections[t] = total_infections_t;
  }
}

