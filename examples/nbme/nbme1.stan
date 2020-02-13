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
  {
    vector[N_subset] y_pred_subset = X[1:N_subset] * b;
    for (n in 1:N_subset){
      y_pred_subset[n] = sigma_state*a_state[state[n]] + sigma_specialty*a_specialty[specialty[n]] + sigma_state_specialty*a_state_specialty[state[n],specialty[n]];
    }
    head(y,N_subset) ~ bernoulli_logit(y_pred_subset);
  }
}
