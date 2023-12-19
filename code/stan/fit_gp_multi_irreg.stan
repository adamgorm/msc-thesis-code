#include ../rao/code/functions.stan

data {
  int n;
  array[n] int J;
  array[sum(J)] real t;
  int J_pred;
  array[J_pred] real t_pred;
  vector[sum(J)] y_vec;
  real<lower = 0> magnitude_mu;
  real<lower = 0> length_scale_mu;
  real<lower = 0> magnitude_eta;
  real<lower = 0> length_scale_eta;
  real<lower = 0> sigma;  
}

transformed data {
  array[n+1] int t_is = append_array({0}, cumulative_sum(J));
}

generated quantities {
  vector[J_pred] mu_pred;
  matrix[J_pred, n] eta_pred;
  matrix[J_pred, n] f_pred;
  {
    matrix[J_pred, n+1] mu_eta_pred =
      mu_eta_completely_irreg_rng(y_vec,
                                  t, t_pred,
                                  n, J_pred,
                                  J, t_is,
                                  magnitude_mu, length_scale_mu,
                                  magnitude_eta, length_scale_eta,
                                  sigma);
    mu_pred = mu_eta_pred[, 1];
    eta_pred = mu_eta_pred[, 2:];
    f_pred = f_draws(mu_eta_pred);
  }
}
