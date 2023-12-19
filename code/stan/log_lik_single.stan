functions {
  real
  log_lik_single(vector y,
                 int J,
                 array[] real t,
                 real magnitude, real length_scale,
                 real sigma)
  {
    matrix[J, J] K = add_diag(gp_exp_quad_cov(t, magnitude, length_scale),
                              max([square(sigma), 1e-8]));
    return multi_normal_lpdf(y | rep_vector(0, J), K);
  }
}
