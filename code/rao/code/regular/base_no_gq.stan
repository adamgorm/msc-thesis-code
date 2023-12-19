#include ../functions.stan

data {
  int n_obs;
  int n_group;
  array[n_obs] real t_obs;
  matrix[n_obs, n_group] y_obs;
  int n_pred;
  array[n_pred] real t_pred;
}

transformed data {
  int n_obs_total = n_obs * n_group;
  vector[n_obs_total] zero_n_obs_total = rep_vector(0, n_obs_total);
  vector[n_obs_total] y_obs_vec = to_vector(y_obs);
}

parameters {
  real<lower = 0> magnitude_mu;
  real<lower = 0> length_scale_mu;
  real<lower = 0> magnitude_eta;
  real<lower = 0> length_scale_eta;
  real<lower = 0> sigma;
}

model {
  length_scale_mu ~ inv_gamma(2, 2);
  magnitude_mu ~ normal(0, 2);
  length_scale_eta ~ inv_gamma(2, 1);
  magnitude_eta ~ normal(0, 2);
  sigma ~ std_normal();
  
  target += log_lik_reg_base(y_obs_vec, t_obs,
                             n_obs, n_group, n_obs_total,
                             magnitude_mu, length_scale_mu,
                             magnitude_eta, length_scale_eta,
                             sigma, zero_n_obs_total);
}

