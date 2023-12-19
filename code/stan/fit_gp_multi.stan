#include ../rao/code/functions.stan

data {
  int n_obs;
  int n_group;
  array[n_obs] real t_obs;
  matrix[n_obs, n_group] y_obs;
  int n_pred;
  array[n_pred] real t_pred;
  real<lower = 0> magnitude_mu;
  real<lower = 0> length_scale_mu;
  real<lower = 0> magnitude_eta;
  real<lower = 0> length_scale_eta;
  real<lower = 0> sigma;  
}

transformed data {
  int n_obs_total = n_obs * n_group;
  matrix[n_group, n_group] one_mat_n_group = rep_matrix(1, n_group, n_group);
  vector[n_obs_total] y_obs_vec = to_vector(y_obs);
}

generated quantities {
  vector[n_pred] mu_pred;
  matrix[n_pred, n_group] eta_pred;
  matrix[n_pred, n_group] f_pred;
  {
    matrix[n_pred, n_group+1] mu_eta_pred =
      mu_eta_reg_smart_iter_rng(y_obs_vec, y_obs, t_obs, t_pred,
                                n_obs, n_group, n_obs_total, n_pred,
                                magnitude_mu, length_scale_mu,
                                magnitude_eta, length_scale_eta,
                                sigma, one_mat_n_group);
    mu_pred = mu_eta_pred[, 1];
    eta_pred = mu_eta_pred[, 2:];
    f_pred = f_draws(mu_eta_pred);
  }
}
