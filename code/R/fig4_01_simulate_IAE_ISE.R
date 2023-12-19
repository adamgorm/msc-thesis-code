## This script must be run interactively - otherwise one gets a segfault (likely
## due to bug in R interface to Stan tuples).

## OBS: Note that the table will not be completely reproducible, since we only
## set R seeds and not Stan seeds (since this would require that we set a Stan
## seed in each simulation run, which would be cumbersome).

library(tidyverse)
theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
source("code/rao/code/irregular/sim_irregular.R")
source("code//R/sim_list_to_tib.R")
sim_mod <- cmdstanr::cmdstan_model("code//rao/code//sim.stan")
fit_mod <- cmdstanr::cmdstan_model("code/stan//fit_gp_multi.stan")
fit_mod_irreg <- cmdstanr::cmdstan_model("code/stan//fit_gp_multi_irreg.stan")

stan_funs <- cmdstanr::cmdstan_model("code/rao//code/functions.stan",
                                     force_recompile = TRUE,
                                     compile_standalone = TRUE)$functions

n <- n_a <- 10
n_b <- 1
J <- 50 # note: final number of observations will be J + 1
scale_factor_f <- 4 # factor of extra observations for the latent functions to
                                        # make them look smooth
J_a <- J * scale_factor_f + 1 # add 1 to be able to include end point
J_b <- rep(J_a, n_b)
t_a <- seq(-10, 10, length.out = J_a)
Jpred <- J_a
tpred <- t_a
## t_a <- c(seq(-5, 1, length.out = J_a/2),
##          seq(2, 5, length.out = J_a/2))
## t_a <- runif(J_a, -5, 5)
t_b <- list()
for (i in 1:n_b) {
    t_b <- c(t_b, list(t_a))
}
magnitude_mu <- 1
length_scale_mu <- 1
magnitude_eta <- 0.5
length_scale_eta <- 0.5
sigma <- 0.2

## Number of different data sets
n_sim <- 100
## Number of draws from the posterior for each data set
n_draw <- 100

get_ise_mu <- function(mu_mat, mu_true)
{
    mu_ise_vec <- apply((mu_mat - mu_true)^2, 2, \(y) pracma::trapz(tpred, y))
    as.vector(mu_ise_vec)
}

get_ise_f <- function(f_array, f_true)
{
    f_ise_mat <- matrix(NA, nrow = n_draw, ncol = n)
    for (i in 1:n_draw) {
        f_ise_mat[i, ] <- apply((f_array[, , i] - f_true)^2, 2, \(y) pracma::trapz(tpred, y))
    }
    as.vector(f_ise_mat)    
}

get_iae_mu <- function(mu_mat, mu_true)
{
    mu_iae_vec <- apply(abs(mu_mat - mu_true), 2, \(y) pracma::trapz(tpred, y))
    as.vector(mu_iae_vec)
}

get_iae_f <- function(f_array, f_true)
{
    f_iae_mat <- matrix(NA, nrow = n_draw, ncol = n)
    for (i in 1:n_draw) {
        f_iae_mat[i, ] <- apply(abs(f_array[, , i] - f_true), 2, \(y) pracma::trapz(tpred, y))
    }
    as.vector(f_iae_mat)
}

## Function to either load data from a file or initialize to NULL
load_or_initialize <- function(file_path) {
  if (file.exists(file_path)) {
    return(as.vector(read.csv(file_path, row.names = NULL)[[1]]))
  } else {
    return(NULL)
  }
}

file_paths <- list(
  mu_ise = "code/R/IAE_ISE_results/mu_ise.csv",
  mu_iae = "code/R/IAE_ISE_results/mu_iae.csv",
  f_ise = "code/R/IAE_ISE_results/f_ise.csv",
  f_iae = "code/R/IAE_ISE_results/f_iae.csv",
  mu_ise_oracle = "code/R/IAE_ISE_results/mu_ise_oracle.csv",
  mu_iae_oracle = "code/R/IAE_ISE_results/mu_iae_oracle.csv",
  f_ise_oracle = "code/R/IAE_ISE_results/f_ise_oracle.csv",
  f_iae_oracle = "code/R/IAE_ISE_results/f_iae_oracle.csv",
  mu_ise_exclude = "code/R/IAE_ISE_results/mu_ise_exclude.csv",
  mu_iae_exclude = "code/R/IAE_ISE_results/mu_iae_exclude.csv",
  f_ise_exclude = "code/R/IAE_ISE_results/f_ise_exclude.csv",
  f_iae_exclude = "code/R/IAE_ISE_results/f_iae_exclude.csv",
  mu_ise_include = "code/R/IAE_ISE_results/mu_ise_include.csv",
  mu_iae_include = "code/R/IAE_ISE_results/mu_iae_include.csv",
  f_ise_include = "code/R/IAE_ISE_results/f_ise_include.csv",
  f_iae_include = "code/R/IAE_ISE_results/f_iae_include.csv"
)

mu_ise <- load_or_initialize(file_paths$mu_ise)
mu_iae <- load_or_initialize(file_paths$mu_iae)
f_ise <- load_or_initialize(file_paths$f_ise)
f_iae <- load_or_initialize(file_paths$f_iae)
mu_ise_oracle <- load_or_initialize(file_paths$mu_ise_oracle)
mu_iae_oracle <- load_or_initialize(file_paths$mu_iae_oracle)
f_ise_oracle <- load_or_initialize(file_paths$f_ise_oracle)
f_iae_oracle <- load_or_initialize(file_paths$f_iae_oracle)
mu_ise_exclude <- load_or_initialize(file_paths$mu_ise_exclude)
mu_iae_exclude <- load_or_initialize(file_paths$mu_iae_exclude)
f_ise_exclude <- load_or_initialize(file_paths$f_ise_exclude)
f_iae_exclude <- load_or_initialize(file_paths$f_iae_exclude)
mu_ise_include <- load_or_initialize(file_paths$mu_ise_include)
mu_iae_include <- load_or_initialize(file_paths$mu_iae_include)
f_ise_include <- load_or_initialize(file_paths$f_ise_include)
f_iae_include <- load_or_initialize(file_paths$f_iae_include)

for (i in 1:n_sim) {
    print("Starting next simulation")
    print(i)
    sim_res <- sim_irregular(n_a, J_a, t_a,
                             n_b, J_b, t_b,
                             magnitude_mu, length_scale_mu,
                             magnitude_eta, length_scale_eta,
                             sigma, sim_mod)
    sim_tib_full <- sim_list_to_tib(sim_res[[2]])
    sim_tib_full <- filter(sim_tib_full, type == "a") |>
        select(-type, -group_id, -observed)

    sim_res[[2]]$f[1, , 1:n] |>
        apply(1, mean) ->
        mu_true

    eta_true <- matrix(NA,
                       nrow = J_a,
                       ncol = n)

    for (i in 1:n) {
        eta_true[, i] <- sim_res[[2]]$f[1, , i] - mu_true
    }

    f_true <- sim_res[[2]]$f[1, , 1:n]

    f1_true <- sim_res[[2]]$f[1, , 1]

    t_obs <- t_a[(0:J) * scale_factor_f + 1]

    sim_tib <- filter(sim_tib_full, t %in% t_obs)

    censoring_val <- 0.75
    sim_tib |>
        mutate(censored = y > censoring_val,
               y_obs = ifelse(censored, censoring_val, y),
               group = as.factor(group)) ->
        sim_tib
    
    sim_tib

    sim_tib |>
        filter(!censored) |>
        count(t) |>
        rename(t_count = n) |>
        filter(t_count == n) |>
        (\(x) x$t)() ->
                 to1
    Jo1 <- length(to1)

    to2 <- c()
    Jo2 <- c()
    tc <- c()
    Jc <- c()
    to <- c()
    Jo <- c()
    for (i in 1:n) {
        sim_tib |>
            filter(group == i, !censored, ! (t %in% to1)) |>
            (\(x) x$t)() ->
                     to2_i
        to2 <- c(to2, to2_i)
        Jo2 <- c(Jo2, length(to2_i))
        sim_tib |>
            filter(group == i, censored, ! (t %in% to1)) |>
            (\(x) x$t)() ->
                     tc_i
        tc <- c(tc, tc_i)
        Jc <- c(Jc, length(tc_i))
        sim_tib |>
            filter(group == i, !censored) |>
            (\(x) x$t)() ->
                     to_i
        to <- c(to, to_i)
        Jo <- c(Jo, length(to_i))
    }
    to2_is <- c(0, cumsum(Jo2))
    tc_is <- c(0, cumsum(Jc))
    to_is <- c(0, cumsum(Jo))


    yo1 <- filter(sim_tib, !censored, t %in% to1)$y_obs
    yo2 <- filter(sim_tib, !censored, ! (t %in% to1))$y_obs
    yo1o2 <- c(yo1, yo2)

    yc <- filter(sim_tib, censored)$y_obs

    yo <- filter(sim_tib, !censored)$y_obs

    marg_log_lik <- function(params)
    {
        if (any(params < 0))
            return(-Inf)
        prep_list <- stan_funs$prep_multi_cens_log_lik(yo,
                                                       to, tc,
                                                       to_is, tc_is,
                                                       Jo, Jc,
                                                       params[[1]], params[[2]],
                                                       params[[3]], params[[4]],
                                                       params[[5]])
        prep_list[[1]] + log(mvtnorm::pmvnorm(lower = yc,
                                              mean = prep_list[[2]],
                                              sigma = prep_list[[3]]))
    }

    marg_log_lik(c(1, 1, 1, 1, 0.1))

    optim_res <- optim(c(magnitude_mu, length_scale_mu,
                         magnitude_eta, length_scale_eta,
                         sigma), marg_log_lik, control = list(fnscale = -1))

    param_max <- optim_res$par
    print(param_max)

    ## Jpred <- 200
    ## tpred <- seq(-10, 10, length.out = Jpred)

    post_list <-
        stan_funs$prep_multi_cond_post(
                      yo, yc,
                      to, tc,
                      tpred,
                      to_is, tc_is,
                      Jo, Jc,
                      Jpred,
                      param_max[1], param_max[2],
                      param_max[3], param_max[4],
                      param_max[5]
                  )

    ## Note: this must sometimes be run twice for some weird reason, to avoid a
    ## getting the wrong stuff back. It seems to be some bug with the R
    ## interface to tuples in Stan. Sometimes a segfault arises, which is why I
    ## save all intermediate results.
    post_list <-
        stan_funs$prep_multi_cond_post(
                      yo, yc,
                      to, tc,
                      tpred,
                      to_is, tc_is,
                      Jo, Jc,
                      Jpred,
                      param_max[1], param_max[2],
                      param_max[3], param_max[4],
                      param_max[5]
                  )

    names(post_list) <- c(
        "mean_c_cond_o",
        "cov_c_cond_o",
        "mean_mu_cond_o",
        "p_factor_mu",
        "cov_q_mu",
        "mean_eta_cond_o",
        "p_factor_eta",
        "cov_q_eta"
    )

    q_draws_mu <- t(mvtnorm::rmvnorm(n_draw, mean = rep(0, Jpred),
                                     sigma = post_list$cov_q_mu))
    p_draws_mu <- t(TruncatedNormal::rtmvnorm(n = n_draw,
                                              mu = rep(0, nrow(post_list$cov_c_cond_o)),
                                              sigma = post_list$cov_c_cond_o,
                                              lb = yc - post_list$mean_c_cond_o))
    mu_draws <- post_list$mean_mu_cond_o +
        post_list$p_factor_mu %*% p_draws_mu +
        q_draws_mu
    post_mean_mu <- rowMeans(mu_draws)
    post_95_mu <- apply(mu_draws, 1, \(x) quantile(x, c(0.05, 0.95)))

    Npred_except_1_group <- (n-1)*Jpred

    q_draws_eta <- t(mvtnorm::rmvnorm(n_draw, mean = rep(0, Npred_except_1_group),
                                      sigma = post_list$cov_q_eta))

    p_draws_eta <- t(TruncatedNormal::rtmvnorm(n = n_draw,
                                               mu = rep(0, nrow(post_list$cov_c_cond_o)),
                                               sigma = post_list$cov_c_cond_o,
                                               lb = yc - post_list$mean_c_cond_o))

    eta_draws <- post_list$mean_eta_cond_o +
        post_list$p_factor_eta %*% p_draws_eta +
        q_draws_eta

    eta_draws <- array(eta_draws, dim = c(Jpred, n-1, n_draw))
    eta_draws_last <- -apply(eta_draws, c(1, 3), sum)
    eta_draws_full <- abind::abind(eta_draws, eta_draws_last, along = 2)
    post_means_eta <- apply(eta_draws_full, c(1,2), mean)
    post_95_eta <- apply(eta_draws_full, c(1,2), \(x) quantile(x, c(0.05, 0.95)))

    f_draws <- array(NA, dim(eta_draws_full))
    for (i in 1:n) {
        f_draws[, i, ] <- mu_draws
    }
    f_draws <- f_draws + eta_draws_full
    post_means_f <- apply(f_draws, c(1,2), mean)
    post_95_f <- apply(f_draws, c(1,2), \(x) quantile(x, c(0.05, 0.95)))

    str(mu_true)
    str(mu_draws[, 1])



    t_obs <- filter(sim_tib, group == 1)$t
    J_obs <- length(t_obs)

    oracle_list <- list(
        n_obs = J_obs,
        n_group = n,
        t_obs = t_obs,
        y_obs = matrix(sim_tib$y, ncol = n),
        n_pred = Jpred,
        t_pred = tpred
    )

    log_lik_oracle <- function(params)
    {
        if (any(params < 0))
            return(-Inf)
        stan_funs$log_lik_reg_smart(
                      as.vector(oracle_list$y_obs),
                      oracle_list$y_obs,
                      oracle_list$t_obs,
                      oracle_list$n_obs,
                      oracle_list$n_group,
                      oracle_list$n_obs * oracle_list$n_group,
                      params[1], params[2],
                      params[3], params[4],
                      params[5],
                      matrix(1, nrow = oracle_list$n_group, ncol = oracle_list$n_group)
                  )
    }

    log_lik_oracle_alt <- function(params)
    {
        if (any(params < 0))
            return(-Inf)
        stan_funs$log_lik_completely_irreg(
                      as.vector(oracle_list$y_obs),
                      rep(oracle_list$t_obs, oracle_list$n_group),
                      oracle_list$n_group,
                      rep(oracle_list$n_obs, oracle_list$n_group),
                      c(0, cumsum(rep(oracle_list$n_obs, oracle_list$n_group))),
                      params[1], params[2],
                      params[3], params[4],
                      params[5]
                  )
    }

    log_lik_oracle(c(1, 2, 3, 4, 5))
    log_lik_oracle_alt(c(1, 2, 3, 4, 5))

    ## The completely irregular implementation gives the same, so it clearly works!
    ## We will use it for the other settings below.

    optim_res_oracle <- optim(c(magnitude_mu, length_scale_mu,
                                magnitude_eta, length_scale_eta,
                                sigma), log_lik_oracle, control = list(fnscale = -1))

    param_max_oracle <- optim_res_oracle$par

    oracle_list$magnitude_mu <- param_max_oracle[1]
    oracle_list$length_scale_mu <- param_max_oracle[2]
    oracle_list$magnitude_eta <- param_max_oracle[3]
    oracle_list$length_scale_eta <- param_max_oracle[4]
    oracle_list$sigma <- param_max_oracle[5]

    oracle_fit <- fit_mod$sample(oracle_list, chains = 1, iter_sampling = n_draw,
                                 fixed_param = TRUE, refresh = 25)

    mu_mat_oracle <- t(matrix(oracle_fit$draws("mu_pred"), nrow = n_draw))

    f_array_oracle <- array(NA, c(Jpred, n, n_draw))
    f_draws_oracle <- oracle_fit$draws("f_pred")
    for (i in 1:n) {
        f_array_oracle[, i, ] <- t(matrix(f_draws_oracle[, 1, ((i-1)*Jpred + 1):(i*Jpred)], nrow = n_draw))
    }



    exclude_list <- list(
        y_vec = yo,
        t = to,
        t_pred = tpred,
        n = n,
        J_pred = Jpred,
        J = Jo,
        t_is = to_is
    )

    log_lik_exclude <- function(params)
    {
        if (any(params < 0))
            return(-Inf)
        stan_funs$log_lik_completely_irreg(
                      exclude_list$y_vec,
                      exclude_list$t,
                      exclude_list$n,
                      exclude_list$J,
                      exclude_list$t_is,
                      params[1], params[2],
                      params[3], params[4],
                      params[5]
                  )
    }

    log_lik_exclude(c(1, 2, 3, 4, 5))

    ## The completely irregular implementation gives the same, so it clearly works!
    ## We will use it for the other settings below.

    optim_res_exclude <- optim(c(magnitude_mu, length_scale_mu,
                                 magnitude_eta, length_scale_eta,
                                 sigma), log_lik_exclude, control = list(fnscale = -1))

    param_max_exclude <- optim_res_exclude$par

    exclude_list$magnitude_mu <- param_max_exclude[1]
    exclude_list$length_scale_mu <- param_max_exclude[2]
    exclude_list$magnitude_eta <- param_max_exclude[3]
    exclude_list$length_scale_eta <- param_max_exclude[4]
    exclude_list$sigma <- param_max_exclude[5]

    n_iter_exclude <- n_draw
    exclude_fit <- fit_mod_irreg$sample(exclude_list, chains = 1,
                                        iter_sampling = n_iter_exclude,
                                        fixed_param = TRUE, refresh = 25)

    mu_mat_exclude <- t(matrix(exclude_fit$draws("mu_pred"), nrow = n_draw))

    f_array_exclude <- array(NA, c(Jpred, n, n_draw))
    f_draws_exclude <- exclude_fit$draws("f_pred")
    for (i in 1:n) {
        f_array_exclude[, i, ] <- t(matrix(f_draws_exclude[, 1, ((i-1)*Jpred + 1):(i*Jpred)], nrow = n_draw))
    }

    include_list <- list(
        y_vec = sim_tib$y_obs,
        t = sim_tib$t,
        t_pred = tpred,
        n = n,
        J_pred = Jpred,
        J = rep(J_obs, n),
        t_is = c(0, cumsum(rep(J_obs, n))),
        magnitude_mu = magnitude_mu,
        length_scale_mu = length_scale_mu,
        magnitude_eta = magnitude_eta,
        length_scale_eta = length_scale_eta,
        sigma = sigma
    )

    log_lik_include <- function(params)
    {
        if (any(params < 0))
            return(-Inf)
        stan_funs$log_lik_completely_irreg(
                      include_list$y_vec,
                      include_list$t,
                      include_list$n,
                      include_list$J,
                      include_list$t_is,
                      params[1], params[2],
                      params[3], params[4],
                      params[5]
                  )
    }

    log_lik_include(c(1, 2, 3, 4, 5))

    ## The completely irregular implementation gives the same, so it clearly works!
    ## We will use it for the other settings below.

    optim_res_include <- optim(c(magnitude_mu, length_scale_mu,
                                 magnitude_eta, length_scale_eta,
                                 sigma), log_lik_include, control = list(fnscale = -1))

    param_max_include <- optim_res_include$par

    include_list$magnitude_mu <- param_max_include[1]
    include_list$length_scale_mu <- param_max_include[2]
    include_list$magnitude_eta <- param_max_include[3]
    include_list$length_scale_eta <- param_max_include[4]
    include_list$sigma <- param_max_include[5]

    n_iter_include <- n_draw
    include_fit <- fit_mod_irreg$sample(include_list, chains = 1,
                                        iter_sampling = n_iter_include,
                                        fixed_param = TRUE, refresh = 25)

    mu_mat_include <- t(matrix(include_fit$draws("mu_pred"), nrow = n_draw))

    f_array_include <- array(NA, c(Jpred, n, n_draw))
    f_draws_include <- include_fit$draws("f_pred")
    for (i in 1:n) {
        f_array_include[, i, ] <- t(matrix(f_draws_include[, 1, ((i-1)*Jpred + 1):(i*Jpred)], nrow = n_draw))
    }

    mu_ise <- c(mu_ise, get_ise_mu(mu_draws, mu_true))
    mu_ise_oracle <- c(mu_ise_oracle, get_ise_mu(mu_mat_oracle, mu_true))
    mu_ise_exclude <- c(mu_ise_exclude, get_ise_mu(mu_mat_exclude, mu_true))
    mu_ise_include <- c(mu_ise_include, get_ise_mu(mu_mat_include, mu_true))
    f_ise <- c(f_ise, get_ise_f(f_draws, f_true))
    f_ise_oracle <- c(f_ise_oracle, get_ise_f(f_array_oracle, f_true))
    f_ise_exclude <- c(f_ise_exclude, get_ise_f(f_array_exclude, f_true))
    f_ise_include <- c(f_ise_include, get_ise_f(f_array_include, f_true))
    mu_iae <- c(mu_iae, get_iae_mu(mu_draws, mu_true))
    mu_iae_oracle <- c(mu_iae_oracle, get_iae_mu(mu_mat_oracle, mu_true))
    mu_iae_exclude <- c(mu_iae_exclude, get_iae_mu(mu_mat_exclude, mu_true))
    mu_iae_include <- c(mu_iae_include, get_iae_mu(mu_mat_include, mu_true))
    f_iae <- c(f_iae, get_iae_f(f_draws, f_true))
    f_iae_oracle <- c(f_iae_oracle, get_iae_f(f_array_oracle, f_true))
    f_iae_exclude <- c(f_iae_exclude, get_iae_f(f_array_exclude, f_true))
    f_iae_include <- c(f_iae_include, get_iae_f(f_array_include, f_true))

    write.csv(mu_ise, file = "code/R/IAE_ISE_results/mu_ise.csv", row.names = FALSE)
    write.csv(mu_iae, file = "code/R/IAE_ISE_results/mu_iae.csv", row.names = FALSE)
    write.csv(f_ise, file = "code/R/IAE_ISE_results/f_ise.csv", row.names = FALSE)
    write.csv(f_iae, file = "code/R/IAE_ISE_results/f_iae.csv", row.names = FALSE)
    write.csv(mu_ise_oracle, file = "code/R/IAE_ISE_results/mu_ise_oracle.csv", row.names = FALSE)
    write.csv(mu_iae_oracle, file = "code/R/IAE_ISE_results/mu_iae_oracle.csv", row.names = FALSE)
    write.csv(f_ise_oracle, file = "code/R/IAE_ISE_results/f_ise_oracle.csv", row.names = FALSE)
    write.csv(f_iae_oracle, file = "code/R/IAE_ISE_results/f_iae_oracle.csv", row.names = FALSE)
    write.csv(mu_ise_exclude, file = "code/R/IAE_ISE_results/mu_ise_exclude.csv", row.names = FALSE)
    write.csv(mu_iae_exclude, file = "code/R/IAE_ISE_results/mu_iae_exclude.csv", row.names = FALSE)
    write.csv(f_ise_exclude, file = "code/R/IAE_ISE_results/f_ise_exclude.csv", row.names = FALSE)
    write.csv(f_iae_exclude, file = "code/R/IAE_ISE_results/f_iae_exclude.csv", row.names = FALSE)
    write.csv(mu_ise_include, file = "code/R/IAE_ISE_results/mu_ise_include.csv", row.names = FALSE)
    write.csv(mu_iae_include, file = "code/R/IAE_ISE_results/mu_iae_include.csv", row.names = FALSE)
    write.csv(f_ise_include, file = "code/R/IAE_ISE_results/f_ise_include.csv", row.names = FALSE)
    write.csv(f_iae_include, file = "code/R/IAE_ISE_results/f_iae_include.csv", row.names = FALSE)
} 

