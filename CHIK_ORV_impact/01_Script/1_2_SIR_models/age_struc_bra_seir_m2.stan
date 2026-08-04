//--------------------------------------------------------------------
// AGE-STRUCTURED DISCRETE-TIME SEIR MODEL
//--------------------------------------------------------------------
// - No vaccination
// - Weekly time step
// - Unobserved initialization period before the first observed week
// - Incidence-based observation model
// - First-order random walk for the weekly transmission rate
// - Baseline seroprevalence supplied externally
// - No external transmission-capacity index
//--------------------------------------------------------------------


//--------------------------------------------------------------------
// DATA BLOCK
//--------------------------------------------------------------------

data {

  // Number of observed epidemiological weeks.
  int<lower=1> T;

  // Number of unobserved epidemiological weeks simulated
  // before the first observed epidemiological week.
  int<lower=1> B;

  // Number of age groups.
  int<lower=1> A;

  // Weekly reported cases by age group.
  // The first dimension represents age groups.
  // The second dimension represents observed epidemiological weeks.
  int<lower=0> observed_cases_by_age[A, T];

  // Population size in each age group.
  vector<lower=0>[A] N;

  // Weekly probability of aging out of each age group.
  // These values must already be expressed on a weekly scale.
  // The final age group should normally have an aging rate of zero.
  vector<lower=0, upper=1>[A] r;

  // Prior mean for the initial infectious population
  // in each age group.
  vector<lower=0>[A] prior_I0;

  // Common prior standard deviation for the initial
  // infectious population.
  real<lower=1e-9> prior_sd_I0;

  // Prespecified baseline seroprevalence in each age group.
  // This is supplied before model fitting.
  // It is not estimated from the epidemic case data.
  vector<lower=0, upper=1>[A] sero;
}


//--------------------------------------------------------------------
// TRANSFORMED DATA BLOCK
//--------------------------------------------------------------------

transformed data {

  // Total number of simulated epidemiological weeks.
  // This includes both unobserved and observed weeks.
  int<lower=2> K;

  // Calculate the full simulation duration.
  K = T + B;
}


//--------------------------------------------------------------------
// PARAMETERS BLOCK
//--------------------------------------------------------------------

parameters {

  // Initial infectious population in each age group.
  // No separate initial exposed population is estimated.
  vector<lower=0>[A] I0;

  // Overall average transmission rate on the logarithmic scale.
  real alpha_log_beta;

  // Standard deviation of week-to-week random-walk innovations.
  // Smaller values imply a smoother transmission trajectory.
  real<lower=0> sigma_beta_rw;

  // Standard-normal innovations used to construct
  // the first-order random walk.
  vector[K - 1] z_beta;

  // Probability that a new infectious onset is observed
  // as a reported case.
  //
  // Unless symptomatic probability is modelled separately,
  // this parameter combines symptomatic probability and
  // case-reporting probability.
  real<lower=1e-6, upper=1> rho;

  // Negative-binomial overdispersion parameter.
  // Larger values imply less extra-Poisson variation.
  real<lower=1e-6> shape;

  // Weekly probability of leaving the infectious compartment.
  // This is a discrete-time transition probability.
  real<lower=1e-6, upper=1> gamma;

  // Weekly probability of leaving the exposed compartment
  // and becoming infectious.
  // This is also a discrete-time transition probability.
  real<lower=1e-6, upper=1> sigma;
}


//--------------------------------------------------------------------
// TRANSFORMED PARAMETERS BLOCK
//--------------------------------------------------------------------

transformed parameters {

  // Total population across all age groups.
  real<lower=0> pop_total;

  // Uncentred random-walk trajectory on the log-transmission scale.
  vector[K] beta_rw_raw;

  // Mean-centred random-walk trajectory.
  vector[K] beta_rw_centered;

  // Weekly transmission-rate trajectory.
  vector<lower=0>[K] beta;

  // Weekly force of infection.
  vector<lower=0>[K] lambda;

  // Weekly probability of infection among susceptible individuals.
  vector<lower=0, upper=1>[K] p_infection;

  // Susceptible population at the start of each model week.
  //
  // There are K transmission intervals and therefore
  // K + 1 compartment states.
  matrix[A, K + 1] S;

  // Exposed population at the start of each model week.
  matrix[A, K + 1] E;

  // Infectious population at the start of each model week.
  matrix[A, K + 1] I;

  // Recovered or previously immune population
  // at the start of each model week.
  matrix[A, K + 1] R;

  // New infections entering the exposed compartment
  // during each model week.
  matrix[A, K] new_exposed;

  // New infectious onsets during each model week.
  matrix[A, K] new_infectious;

  // New recoveries during each model week.
  matrix[A, K] new_recovered;

  // Expected reported cases during each observed week.
  matrix[A, T] expected_reported_cases;


  //------------------------------------------------------------------
  // CALCULATE THE TOTAL POPULATION
  //------------------------------------------------------------------

  // Sum population sizes across all age groups.
  pop_total = sum(N);


  //------------------------------------------------------------------
  // CONSTRUCT THE RANDOM-WALK TRANSMISSION TRAJECTORY
  //------------------------------------------------------------------

  // Anchor the first random-walk value at zero.
  beta_rw_raw[1] = 0;

  // Construct all subsequent random-walk values.
  for (k in 2:K) {

    // Add a scaled standard-normal innovation
    // to the previous week's value.
    beta_rw_raw[k] =
      beta_rw_raw[k - 1] +
      sigma_beta_rw * z_beta[k - 1];
  }

  // Remove the temporal mean from the random-walk component.
  //
  // This separates the overall transmission level,
  // alpha_log_beta, from temporal variation.
  beta_rw_centered =
    beta_rw_raw -
    rep_vector(mean(beta_rw_raw), K);

  // Transform the log-scale trajectory
  // to the natural transmission-rate scale.
  for (k in 1:K) {

    // Calculate the weekly transmission rate.
    beta[k] =
      exp(
        alpha_log_beta +
        beta_rw_centered[k]
      );
  }


  //------------------------------------------------------------------
  // SET THE INITIAL COMPARTMENT SIZES
  //------------------------------------------------------------------

  for (a in 1:A) {

    // Set the initially immune population using
    // externally supplied baseline seroprevalence.
    R[a, 1] =
      sero[a] * N[a];

    // Set the initial exposed population to zero.
    //
    // The unobserved initialization period allows
    // infectious individuals to generate exposed individuals
    // before the observation period begins.
    E[a, 1] = 0;

    // Set the initial infectious population.
    I[a, 1] =
      I0[a];

    // Assign the remaining population
    // to the susceptible compartment.
    S[a, 1] =
      N[a] -
      R[a, 1] -
      I[a, 1];
  }


  //------------------------------------------------------------------
  // SIMULATE THE FULL SEIR TRAJECTORY
  //------------------------------------------------------------------

  // Simulate both the unobserved initialization period
  // and the observed epidemic period.
  for (k in 1:K) {

    // Total infectious population at the start of the week.
    real total_infectious;

    // Compartment values after epidemiological transitions
    // but before aging transitions.
    vector[A] S_after_epi;
    vector[A] E_after_epi;
    vector[A] I_after_epi;
    vector[A] R_after_epi;

    // Numbers aging out of each compartment.
    vector[A] aging_out_S;
    vector[A] aging_out_E;
    vector[A] aging_out_I;
    vector[A] aging_out_R;

    // Initialise the total infectious population.
    total_infectious = 0;

    // Sum infectious individuals across all age groups.
    for (a in 1:A) {

      // Add age-specific infectious prevalence.
      total_infectious +=
        I[a, k];
    }

    // Calculate the force of infection during the current week.
    lambda[k] =
      beta[k] *
      total_infectious /
      pop_total;

    // Convert the force of infection
    // to a weekly infection probability.
    //
    // This prevents the number of new infections
    // from exceeding the susceptible population.
    p_infection[k] =
      1 -
      exp(-lambda[k]);


    //----------------------------------------------------------------
    // CALCULATE EPIDEMIOLOGICAL TRANSITIONS
    //----------------------------------------------------------------

    for (a in 1:A) {

      // Calculate new infections entering
      // the exposed compartment.
      new_exposed[a, k] =
        S[a, k] *
        p_infection[k];

      // Calculate exposed individuals
      // becoming infectious.
      new_infectious[a, k] =
        sigma *
        E[a, k];

      // Calculate infectious individuals recovering.
      new_recovered[a, k] =
        gamma *
        I[a, k];

      // Update susceptible individuals before aging.
      S_after_epi[a] =
        S[a, k] -
        new_exposed[a, k];

      // Update exposed individuals before aging.
      E_after_epi[a] =
        E[a, k] +
        new_exposed[a, k] -
        new_infectious[a, k];

      // Update infectious individuals before aging.
      I_after_epi[a] =
        I[a, k] +
        new_infectious[a, k] -
        new_recovered[a, k];

      // Update recovered individuals before aging.
      R_after_epi[a] =
        R[a, k] +
        new_recovered[a, k];

      // Calculate susceptible individuals aging out.
      aging_out_S[a] =
        r[a] *
        S_after_epi[a];

      // Calculate exposed individuals aging out.
      aging_out_E[a] =
        r[a] *
        E_after_epi[a];

      // Calculate infectious individuals aging out.
      aging_out_I[a] =
        r[a] *
        I_after_epi[a];

      // Calculate recovered individuals aging out.
      aging_out_R[a] =
        r[a] *
        R_after_epi[a];
    }


    //----------------------------------------------------------------
    // APPLY AGING TRANSITIONS
    //----------------------------------------------------------------

    for (a in 1:A) {

      // Handle the youngest age group.
      // This group receives no inflow from a younger age group.
      if (a == 1) {

        // Update susceptible individuals.
        S[a, k + 1] =
          S_after_epi[a] -
          aging_out_S[a];

        // Update exposed individuals.
        E[a, k + 1] =
          E_after_epi[a] -
          aging_out_E[a];

        // Update infectious individuals.
        I[a, k + 1] =
          I_after_epi[a] -
          aging_out_I[a];

        // Update recovered individuals.
        R[a, k + 1] =
          R_after_epi[a] -
          aging_out_R[a];
      }

      // Handle all remaining age groups.
      else {

        // Update susceptible individuals.
        S[a, k + 1] =
          S_after_epi[a] -
          aging_out_S[a] +
          aging_out_S[a - 1];

        // Update exposed individuals.
        E[a, k + 1] =
          E_after_epi[a] -
          aging_out_E[a] +
          aging_out_E[a - 1];

        // Update infectious individuals.
        I[a, k + 1] =
          I_after_epi[a] -
          aging_out_I[a] +
          aging_out_I[a - 1];

        // Update recovered individuals.
        R[a, k + 1] =
          R_after_epi[a] -
          aging_out_R[a] +
          aging_out_R[a - 1];
      }
    }
  }


  //------------------------------------------------------------------
  // MAP MODEL WEEKS TO OBSERVED WEEKS
  //------------------------------------------------------------------

  for (t in 1:T) {

    // Model-week index corresponding to observed week t.
    int model_week;

    // Observed week 1 corresponds to model week B + 1.
    model_week =
      B + t;

    for (a in 1:A) {

      // Link reported cases to new infectious onsets.
      expected_reported_cases[a, t] =
        rho *
        new_infectious[a, model_week] +
        1e-9;
    }
  }
}


//--------------------------------------------------------------------
// MODEL BLOCK
//--------------------------------------------------------------------

model {

  //------------------------------------------------------------------
  // PRIOR FOR THE INITIAL INFECTIOUS POPULATION
  //------------------------------------------------------------------

  for (a in 1:A) {

    // Apply the supplied prior to the initial infectious population.
    I0[a] ~ normal(
      prior_I0[a],
      prior_sd_I0
    );
  }


  //------------------------------------------------------------------
  // PRIORS FOR THE TRANSMISSION TRAJECTORY
  //------------------------------------------------------------------

  // Prior for the overall transmission level
  // on the logarithmic scale.
  alpha_log_beta ~ normal(
    -1,
    0.5
  );

  // Half-normal prior for the standard deviation
  // of weekly random-walk innovations.
  sigma_beta_rw ~ normal(
    0,
    0.15
  );

  // Standard-normal priors
  // for the random-walk innovations.
  z_beta ~ std_normal();


  //------------------------------------------------------------------
  // PRIORS FOR THE OBSERVATION MODEL
  //------------------------------------------------------------------

  // Prior mean reporting probability is 0.25.
  rho ~ beta(
    20,
    60
  );

  // Prior for negative-binomial overdispersion.
  shape ~ lognormal(
    log(17.5),
    0.5
  );


  //------------------------------------------------------------------
  // PRIORS FOR DISEASE-PROGRESSION PARAMETERS
  //------------------------------------------------------------------

  // Prior for the weekly probability
  // of leaving the infectious compartment.
  gamma ~ normal(
    0.67,
    0.02
  );

  // Prior for the weekly probability
  // of becoming infectious after exposure.
  sigma ~ normal(
    0.90,
    0.10
  );


  //------------------------------------------------------------------
  // OBSERVATION LIKELIHOOD
  //------------------------------------------------------------------

  for (t in 1:T) {

    for (a in 1:A) {

      // Fit age-specific weekly reported cases
      // using a negative-binomial observation model.
      observed_cases_by_age[a, t] ~
        neg_binomial_2(
          expected_reported_cases[a, t],
          shape
        );
    }
  }
}


//--------------------------------------------------------------------
// GENERATED QUANTITIES BLOCK
//--------------------------------------------------------------------

generated quantities {

  // Full susceptible trajectory.
  matrix[A, K + 1] S_pred;

  // Full exposed trajectory.
  matrix[A, K + 1] E_pred;

  // Full infectious trajectory.
  matrix[A, K + 1] I_pred;

  // Full recovered trajectory.
  matrix[A, K + 1] R_pred;

  // Expected total reported cases during each observed week.
  vector[T] pred_cases;

  // Total incident infections during each observed week.
  vector[T] incident_infections;

  // Total new infectious onsets during each observed week.
  vector[T] incident_infectious_onsets;

  // Backward-compatible alias for incident infections.
  vector[T] aggregated_infections;

  // Infectious prevalence at the start of each observed week.
  vector[T] infectious_prevalence;

  // Effective reproduction number during each observed week.
  vector[T] R_eff;

  // Force of infection during each observed week.
  vector[T] phi_pred;

  // Weekly transmission rate during each observed week.
  vector[T] beta_observed;

  // Total model population at the start of each observed week
  // and at the end of the final observed week.
  vector[T + 1] N_pred;

  // Expected reported cases by age group and observed week.
  matrix[A, T] age_stratified_cases;

  // Posterior-predictive replicated case counts.
  int replicated_cases_by_age[A, T];

  // Pointwise log-likelihood values.
  matrix[A, T] log_lik;

  // Susceptible fraction at the start of each observed week.
  vector[T] susceptible_fraction;


  //------------------------------------------------------------------
  // COPY THE FULL COMPARTMENT TRAJECTORIES
  //------------------------------------------------------------------

  // Copy the susceptible trajectory.
  S_pred = S;

  // Copy the exposed trajectory.
  E_pred = E;

  // Copy the infectious trajectory.
  I_pred = I;

  // Copy the recovered trajectory.
  R_pred = R;


  //------------------------------------------------------------------
  // CALCULATE OBSERVED-PERIOD OUTPUTS
  //------------------------------------------------------------------

  for (t in 1:T) {

    // Model-week index corresponding to observed week t.
    int model_week;

    // Total susceptible population at the start of the week.
    real susceptible_total;

    // Map the observed week to the full model week.
    model_week =
      B + t;

    // Initialise weekly expected reported cases.
    pred_cases[t] = 0;

    // Initialise weekly incident infections.
    incident_infections[t] = 0;

    // Initialise weekly infectious onsets.
    incident_infectious_onsets[t] = 0;

    // Initialise weekly infectious prevalence.
    infectious_prevalence[t] = 0;

    // Initialise weekly susceptible population.
    susceptible_total = 0;

    // Store the force of infection.
    phi_pred[t] =
      lambda[model_week];

    // Store the weekly transmission rate.
    beta_observed[t] =
      beta[model_week];

    // Aggregate age-specific quantities.
    for (a in 1:A) {

      // Store expected age-specific reported cases.
      age_stratified_cases[a, t] =
        expected_reported_cases[a, t];

      // Sum expected reported cases.
      pred_cases[t] +=
        expected_reported_cases[a, t];

      // Sum incident infections.
      incident_infections[t] +=
        new_exposed[a, model_week];

      // Sum new infectious onsets.
      incident_infectious_onsets[t] +=
        new_infectious[a, model_week];

      // Sum infectious prevalence.
      infectious_prevalence[t] +=
        I[a, model_week];

      // Sum susceptible individuals.
      susceptible_total +=
        S[a, model_week];

      // Generate posterior-predictive reported cases.
      replicated_cases_by_age[a, t] =
        neg_binomial_2_rng(
          expected_reported_cases[a, t],
          shape
        );

      // Calculate the pointwise log likelihood.
      log_lik[a, t] =
        neg_binomial_2_lpmf(
          observed_cases_by_age[a, t] |
          expected_reported_cases[a, t],
          shape
        );
    }

    // Retain the original output name for compatibility.
    aggregated_infections[t] =
      incident_infections[t];

    // Calculate the susceptible fraction.
    susceptible_fraction[t] =
      susceptible_total /
      pop_total;

    // Calculate the approximate effective reproduction number.
    R_eff[t] =
      (
        beta[model_week] /
        gamma
      ) *
      susceptible_fraction[t];
  }


  //------------------------------------------------------------------
  // CHECK POPULATION CONSERVATION
  //------------------------------------------------------------------

  for (t in 1:(T + 1)) {

    // Full model-state index corresponding
    // to the observed-period state.
    int state_week;

    // Map the observed state to the full trajectory.
    state_week =
      B + t;

    // Initialise the population total.
    N_pred[t] = 0;

    // Sum all compartments across age groups.
    for (a in 1:A) {

      // Add the age-specific compartment sizes.
      N_pred[t] +=
        S[a, state_week] +
        E[a, state_week] +
        I[a, state_week] +
        R[a, state_week];
    }
  }
}
