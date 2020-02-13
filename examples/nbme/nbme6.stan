functions {
  real nbme(int start, int end, int[] y_slice, matrix X_subset, vector b,
	    real sigma_state, vector a_state, int[] state,
	    real sigma_specialty, vector a_specialty, int[] specialty,
	    real sigma_state_specialty, matrix a_state_specialty) {
    int Ns = size(y_slice);
    vector[Ns] mu = X_subset[start:end] * b;
    
    for (n in start:end) {
      mu[n - start + 1] = mu[n - start + 1] +
	sigma_state * a_state[state[n]] +
	sigma_specialty * a_specialty[specialty[n]] +
	sigma_state_specialty * a_state_specialty[state[n], specialty[n]];
    }

    return bernoulli_logit_lpmf(y_slice | mu);
  }
}

data{
  int<lower=1> N;
  int<lower=1> N_state;
  int<lower=1> N_specialty;
  int<lower=1,upper=N_state> state[N];
  int<lower=1,upper=N_specialty> specialty[N];
  int<lower=0,upper=1> y[N];
  int female[N];
  real S2[N];
  real year[N];
  int<lower=1> K;
  matrix[N,K] X;
  int<lower=1> N_subset;
}

transformed data {
  int y_subset[N_subset] = y[1:N_subset];
  matrix[N_subset, K] X_subset = X[1:N_subset];
}

parameters {
  real<lower=0> sigma_state;
  real<lower=0> sigma_specialty;
  real<lower=0> sigma_state_specialty;
  vector[K] b;
  vector[N_state] a_state;
  vector[N_specialty] a_specialty;
  matrix[N_state,N_specialty] a_state_specialty;
}

model {
  sigma_state ~ exponential(1);
  sigma_specialty ~ exponential(1);
  sigma_state_specialty ~ exponential(1);
  a_state ~ normal(0,1);
  a_specialty ~ normal(0,1);

  to_vector(a_state_specialty) ~ normal(0,1);

  target += reduce_sum(nbme, y_subset, 125,
		       X_subset, b,
		       sigma_state, a_state, state,
		       sigma_specialty, a_specialty, specialty,
		       sigma_state_specialty, a_state_specialty);
}
