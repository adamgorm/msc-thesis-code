data {
  int<lower = 1> Ny;
  int<lower = 1> Nf;
  array[Ny] real ty;
  array[Nf] real tf;
  real<lower = 0> magnitude;
  real<lower = 0> length_scale;
  real<lower = 0> sigma;
}

transformed data {
  int N = Ny + Nf;
  array[N] real t = append_array(tf, ty);
  matrix[N, N] K;
  matrix[N, N] L;
  vector[N] zero_N;
  zero_N = rep_vector(0, N);
  // We add some jitter to the diagonal for numerical stability.
  K = add_diag(
    gp_exp_quad_cov(t, magnitude, length_scale),
    1e-8
  );
  L = cholesky_decompose(K);
}

generated quantities {
  vector[Nf] f;
  array[Ny] real y;
  {
  vector[N] sim;
  sim = multi_normal_cholesky_rng(zero_N, L);
  f = sim[1:Nf];
  y = normal_rng(sim[(Nf+1):N], sigma);
  }
}
