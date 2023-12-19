functions {
  real log_lik_smart_no_chol(array[] real t_obs, array[] real t_pred,
                     int n_obs, int n_group, int n_obs_total,
                     matrix y_obs, vector y_obs_vec,
                     real magnitude_mu, real length_scale_mu,
                     real magnitude_eta, real length_scale_eta,
                     real sigma)
  {
    matrix[n_obs, n_obs] K_mu =
      gp_exp_quad_cov(t_obs, magnitude_mu, length_scale_mu);
    matrix[n_obs, n_obs] K_eta =
      gp_exp_quad_cov(t_obs, magnitude_eta, length_scale_eta);
    matrix[n_obs, n_obs] Sigma0 = add_diag(n_group / (n_group - 1.0) * K_eta, square(sigma));
    matrix[n_obs, n_obs] Sigma1 = add_diag(n_group * K_mu, square(sigma));
    real log_det0 = log_determinant(Sigma0);
    real log_det1 = log_determinant(Sigma1);
    real log_det = (n_group - 1.0) * log_det0 + log_det1;
    matrix[n_obs, n_group] prod0 = Sigma0 \ y_obs;
    matrix[n_obs, n_group] prod1 = Sigma1 \ y_obs;
    real quad_form = y_obs_vec' *
      to_vector(prod0 + 1.0 / n_group * (prod1 - prod0) * rep_matrix(1.0, n_group, n_group));

    // Note: Here I include the normalization constant, for comparability with the lpdf function where we have tu include the normalization constant.
    return -0.5 * (n_obs_total * log(2*pi()) + log_det + quad_form);
  }
}
