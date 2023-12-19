data {
  int n_obs;
  int n_pred;
  array[n_obs] real t_obs;
  vector[n_obs] y_obs;
  array[n_pred] real t_pred;
  real magnitude;
  real length_scale;
  real sigma;
}

transformed data {
  vector[n_pred] mean_post;
  matrix[n_pred, n_pred] L_post;
  {
  int N = n_obs + n_pred;  
  matrix[N, N] K = add_diag(
    gp_exp_quad_cov(append_array(t_obs, t_pred), magnitude, length_scale),
    1e-8
  );
  matrix[n_obs, n_obs] K_obs = K[1:n_obs, 1:n_obs];
  matrix[n_pred, n_obs] K_pred_obs = K[(n_obs + 1):N, 1:n_obs];
  matrix[n_pred, n_pred] K_pred = K[(n_obs + 1):N, (n_obs + 1):N];
  matrix[n_pred, n_obs] inverse_prod = K_pred_obs / add_diag(K_obs, square(sigma));
  mean_post = inverse_prod * y_obs;
  matrix[n_pred, n_pred] K_post = K_pred - inverse_prod * K_pred_obs';
  L_post = cholesky_decompose(K_post);
  }
}

generated quantities {
  vector[n_pred] f_pred;
  f_pred = multi_normal_cholesky_rng(mean_post, L_post);
}

