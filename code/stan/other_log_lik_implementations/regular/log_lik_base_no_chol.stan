functions {
  /* Very plain calculation of the log-lik without any fancy tricks,
     for testing purposes */
  real log_lik_base_no_chol(vector y_obs_vec, array[] real t_obs,
               int n_obs, int n_group, int n_obs_total,
               real magnitude_mu, real length_scale_mu,
               real magnitude_eta, real length_scale_eta,
               real sigma,
               vector zero_n_obs_total)
  {
    matrix[n_obs, n_obs] K_mu =
      gp_exp_quad_cov(t_obs, magnitude_mu, length_scale_mu);
    matrix[n_obs, n_obs] K_eta =
      gp_exp_quad_cov(t_obs, magnitude_eta, length_scale_eta);
    matrix[n_obs_total, n_obs_total] cov_mat;
    matrix[n_obs, n_obs] diagonal_block = add_diag(K_mu + K_eta, max([square(sigma), 1e-8]));
    matrix[n_obs, n_obs] off_diagonal_block = K_mu - K_eta / (n_group - 1);
    int row_start, row_end, col_start, col_end;
    for (i_col in 1:n_group) {
      col_start = (i_col-1) * n_obs + 1;
      col_end = col_start + n_obs - 1;
      for (i_row in i_col:n_group) {
        row_start = (i_row-1) * n_obs + 1;
        row_end = row_start + n_obs - 1;
        if (i_row == i_col) {
          cov_mat[row_start:row_end, col_start:col_end] = diagonal_block;
        } else {
          cov_mat[row_start:row_end, col_start:col_end] = off_diagonal_block;
          cov_mat[col_start:col_end, row_start:row_end] = off_diagonal_block;
        }
      }
    }
    // Note: Here we have to use the lpdf version (instead of the lupdf), since we can only use lupdf if we call our own function _lpdf, but then we can't expose the function to Stan.
    return multi_normal_lpdf(y_obs_vec | zero_n_obs_total, cov_mat);
  }
}
