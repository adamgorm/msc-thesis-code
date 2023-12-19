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
sim_res <- sim_irregular(n_a, J_a, t_a,
                         n_b, J_b, t_b,
                         magnitude_mu, length_scale_mu,
                         magnitude_eta, length_scale_eta,
                         sigma, sim_mod, stanseed = 1533)
sim_tib_full <- sim_list_to_tib(sim_res[[2]])
sim_tib_full <- filter(sim_tib_full, type == "a") |>
    select(-type, -group_id, -observed)

sim_res[[2]]$f[1, , 1:n] |>
    apply(1, mean) ->
    mu_true

f1_true <- sim_res[[2]]$f[1, , 1]

t_obs <- t_a[(0:J) * scale_factor_f + 1]

sim_tib <- filter(sim_tib_full, t %in% t_obs)

censoring_val <- 0.75
sim_tib |>
    mutate(censored = y > censoring_val,
           y_obs = ifelse(censored, censoring_val, y),
           group = as.factor(group)) ->
    sim_tib

ggplot(sim_tib) +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.2) +
    geom_line(aes(x = t, y = y_obs, col = group)) +
    geom_point(aes(x = t, y = y, col = group, shape = censored)) +
    geom_line(aes(x = t_a, y = mu_true, col = "mu true"), data = tibble(), linewidth = 2) +
    geom_line(aes(x = t_a, y = f1_true, col = "f1 true"), data = tibble(), linewidth = 1)
    
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

Jpred <- 200
tpred <- seq(-10, 10, length.out = Jpred)

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

n_draw <- 1000

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
post_95_mu <- apply(mu_draws, 1, \(x) quantile(x, c(0.025, 0.975)))

ggplot() +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.2) +
    geom_line(aes(x = t, y = y_obs, col = group),
              sim_tib) +
    geom_point(aes(x = t, y = y, col = group, shape = censored),
               sim_tib) +
    geom_line(aes(x = tpred, y = post_mean_mu, col = "mu post mean"), linewidth = 2) +
    geom_ribbon(aes(x = tpred, ymin = post_95_mu[1,], ymax = post_95_mu[2,]),
                alpha = 0.2) +
    geom_line(aes(x = tpred, y = mu_draws[, sample(1:n_draw, 1)], col = "mu single draw")) + # A single mu draw
    geom_line(aes(x = t_a, y = mu_true, col = "mu true"), linewidth = 2)

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
post_95_eta <- apply(eta_draws_full, c(1,2), \(x) quantile(x, c(0.025, 0.975)))

f_draws <- array(NA, dim(eta_draws_full))
for (i in 1:n) {
    f_draws[, i, ] <- mu_draws
}
f_draws <- f_draws + eta_draws_full
post_means_f <- apply(f_draws, c(1,2), mean)
post_95_f <- apply(f_draws, c(1,2), \(x) quantile(x, c(0.025, 0.975)))

geom_random_fs <- function(num_f_to_plot)
{
    shuffle_f <- sample(1:n)
    lapply(1:num_f_to_plot,
           \(i)
           list(
               geom_line(aes(x = tpred, y = post_means_f[, shuffle_f[i]],
                             col = paste0("f", shuffle_f[i]))),
               geom_ribbon(aes(x = tpred,
                               ymin = post_95_f[1, , shuffle_f[i]],
                               ymax = post_95_f[2, , shuffle_f[i]],
                               fill = paste0("f", shuffle_f[i])),
                           alpha = 0.1),
               geom_line(aes(x = t_a, y = sim_res[[2]]$f[1, , shuffle_f[i]],
                             col = paste0("f", shuffle_f[i])),
                         lty = 2)
           )
           ) |>
        unlist()
}

ggplot() +
    geom_random_fs(2) +
    bayesplot::hline_at(censoring_val, lty = 3)

ggplot() +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.5) +
    geom_point(aes(x = t, y = y, shape = censored),
               data = filter(sim_tib, group != 1)) +
    geom_point(aes(x = t, y = y, shape = censored, col = "f1"),
               data = filter(sim_tib, group == 1),
               size = 2) +    
    geom_line(aes(x = tpred, y = post_means_f[, 1], col = "f1", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = post_95_f[1, , 1], ymax = post_95_f[2, , 1],
                    fill = "f1"), alpha = 0.2) +    
    geom_line(aes(x = t_a, y = f1_true, col = "f1", lty = "true")) +
    geom_line(aes(x = tpred, y = post_mean_mu, col = "mu", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = post_95_mu[1, ], ymax = post_95_mu[2, ],
                    fill = "mu"), alpha = 0.2) +
    geom_line(aes(x = t_a, y = mu_true, col = "mu", lty = "true")) +
    xlab("t") +
    ylab("")

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

set.seed(1202)
optim_res_oracle <- optim(c(magnitude_mu, length_scale_mu,
                            magnitude_eta, length_scale_eta,
                            sigma), log_lik_oracle, control = list(fnscale = -1))

param_max_oracle <- optim_res_oracle$par

oracle_list$magnitude_mu <- param_max_oracle[1]
oracle_list$length_scale_mu <- param_max_oracle[2]
oracle_list$magnitude_eta <- param_max_oracle[3]
oracle_list$length_scale_eta <- param_max_oracle[4]
oracle_list$sigma <- param_max_oracle[5]

oracle_fit <- fit_mod$sample(oracle_list, chains = 1, iter_sampling = 1000,
                             fixed_param = TRUE, seed = 1813)

mu_mat_oracle <- matrix(oracle_fit$draws("mu_pred"), nrow = 1000)
mu_post_mean_oracle <- colMeans(mu_mat_oracle)
mu_post_95_oracle <- apply(mu_mat_oracle, 2, \(x) quantile(x, c(0.025, 0.975)))

f1_mat_oracle <- matrix(oracle_fit$draws("f_pred")[, 1, 1:Jpred], nrow = 1000)
f1_post_mean_oracle <- colMeans(f1_mat_oracle)
f1_post_95_oracle <- apply(f1_mat_oracle, 2, \(x) quantile(x, c(0.025, 0.975)))

ggplot() +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.5) +
    geom_point(aes(x = t, y = y, shape = censored),
               data = filter(sim_tib, group != 1)) +
    geom_point(aes(x = t, y = y, shape = censored, col = "f1"),
               data = filter(sim_tib, group == 1),
               size = 2) +    
    geom_line(aes(x = tpred, y = f1_post_mean_oracle, col = "f1", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = f1_post_95_oracle[1, ],
                    ymax = f1_post_95_oracle[2, ],
                    fill = "f1"), alpha = 0.2) +    
    geom_line(aes(x = t_a, y = f1_true, col = "f1", lty = "true")) +
    geom_line(aes(x = tpred, y = mu_post_mean_oracle, col = "mu", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = mu_post_95_oracle[1, ],
                    ymax = mu_post_95_oracle[2, ],
                    fill = "mu"), alpha = 0.2) +
    geom_line(aes(x = t_a, y = mu_true, col = "mu", lty = "true")) +
    xlab("t") +
    ylab("")

## TODO: Modify for completely irregular version

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

set.seed(1431)
optim_res_exclude <- optim(c(magnitude_mu, length_scale_mu,
                            magnitude_eta, length_scale_eta,
                            sigma), log_lik_exclude, control = list(fnscale = -1))

param_max_exclude <- optim_res_exclude$par

exclude_list$magnitude_mu <- param_max_exclude[1]
exclude_list$length_scale_mu <- param_max_exclude[2]
exclude_list$magnitude_eta <- param_max_exclude[3]
exclude_list$length_scale_eta <- param_max_exclude[4]
exclude_list$sigma <- param_max_exclude[5]

n_iter_exclude <- 1000
exclude_fit <- fit_mod_irreg$sample(exclude_list, chains = 1,
                                    iter_sampling = n_iter_exclude,
                                    fixed_param = TRUE, seed = 949, refresh = 5)

mu_mat_exclude <- matrix(exclude_fit$draws("mu_pred"),
                         nrow = n_iter_exclude)
mu_post_mean_exclude <- colMeans(mu_mat_exclude)
mu_post_95_exclude <- apply(mu_mat_exclude, 2, \(x) quantile(x, c(0.025, 0.975)))

f1_mat_exclude <- matrix(exclude_fit$draws("f_pred")[, 1, 1:Jpred],
                         nrow = n_iter_exclude)
f1_post_mean_exclude <- colMeans(f1_mat_exclude)
f1_post_95_exclude <- apply(f1_mat_exclude, 2, \(x) quantile(x, c(0.025, 0.975)))

ggplot() +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.5) +
    geom_point(aes(x = t, y = y, shape = censored),
               data = filter(sim_tib, group != 1)) +
    geom_point(aes(x = t, y = y, shape = censored, col = "f1"),
               data = filter(sim_tib, group == 1),
               size = 2) +    
    geom_line(aes(x = tpred, y = f1_post_mean_exclude, col = "f1", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = f1_post_95_exclude[1, ],
                    ymax = f1_post_95_exclude[2, ],
                    fill = "f1"), alpha = 0.2) +    
    geom_line(aes(x = t_a, y = f1_true, col = "f1", lty = "true")) +
    geom_line(aes(x = tpred, y = mu_post_mean_exclude, col = "mu", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = mu_post_95_exclude[1, ],
                    ymax = mu_post_95_exclude[2, ],
                    fill = "mu"), alpha = 0.2) +
    geom_line(aes(x = t_a, y = mu_true, col = "mu", lty = "true")) +
    xlab("t") +
    ylab("")

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

set.seed(1433)
optim_res_include <- optim(c(magnitude_mu, length_scale_mu,
                            magnitude_eta, length_scale_eta,
                            sigma), log_lik_include, control = list(fnscale = -1))

param_max_include <- optim_res_include$par

include_list$magnitude_mu <- param_max_include[1]
include_list$length_scale_mu <- param_max_include[2]
include_list$magnitude_eta <- param_max_include[3]
include_list$length_scale_eta <- param_max_include[4]
include_list$sigma <- param_max_include[5]

n_iter_include <- 1000
include_fit <- fit_mod_irreg$sample(include_list, chains = 1,
                                    iter_sampling = n_iter_include,
                                    fixed_param = TRUE, seed = 1009, refresh = 5)

mu_mat_include <- matrix(include_fit$draws("mu_pred"), nrow = n_iter_include)
mu_post_mean_include <- colMeans(mu_mat_include)
mu_post_95_include <- apply(mu_mat_include, 2, \(x) quantile(x, c(0.025, 0.975)))

f1_mat_include <- matrix(include_fit$draws("f_pred")[, 1, 1:Jpred], nrow = n_iter_include)
f1_post_mean_include <- colMeans(f1_mat_include)
f1_post_95_include <- apply(f1_mat_include, 2, \(x) quantile(x, c(0.025, 0.975)))

ggplot() +
    bayesplot::hline_at(censoring_val, lty = 2, alpha = 0.5) +
    geom_point(aes(x = t, y = y, shape = censored),
               data = filter(sim_tib, group != 1)) +
    geom_point(aes(x = t, y = y, shape = censored, col = "f1"),
               data = filter(sim_tib, group == 1),
               size = 2) +    
    geom_line(aes(x = tpred, y = f1_post_mean_include, col = "f1", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = f1_post_95_include[1, ],
                    ymax = f1_post_95_include[2, ],
                    fill = "f1"), alpha = 0.2) +    
    geom_line(aes(x = t_a, y = f1_true, col = "f1", lty = "true")) +
    geom_line(aes(x = tpred, y = mu_post_mean_include, col = "mu", lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = mu_post_95_include[1, ],
                    ymax = mu_post_95_include[2, ],
                    fill = "mu"), alpha = 0.2) +
    geom_line(aes(x = t_a, y = mu_true, col = "mu", lty = "true")) +
    xlab("t") +
    ylab("")

get_post_tib <- function(mu_post_mean, mu_post_95,
                         f1_post_mean, f1_post_95,
                         model)
{
    bind_rows(tibble(tpred,
                     mean = mu_post_mean,
                     lwr = mu_post_95[1,],
                     upr = mu_post_95[2,],
                     parameter = "mu",
                     model = model),
              tibble(tpred,
                     mean = f1_post_mean,
                     lwr = f1_post_95[1,],
                     upr = f1_post_95[2,],
                     parameter = "f1",
                     model = model))
}

combined_tib <- bind_rows(
    get_post_tib(post_mean_mu,
                 post_95_mu,
                 post_means_f[, 1],
                 post_95_f[, , 1],
                 "censoring model"),
    get_post_tib(mu_post_mean_oracle,
                 mu_post_95_oracle,
                 f1_post_mean_oracle,
                 f1_post_95_oracle,
                 "oracle"),
    get_post_tib(mu_post_mean_exclude,
                 mu_post_95_exclude,
                 f1_post_mean_exclude,
                 f1_post_95_exclude,
                 "exclude"),
    get_post_tib(mu_post_mean_include,
                 mu_post_95_include,
                 f1_post_mean_include,
                 f1_post_95_include,
                 "include")
)

true_tib <- bind_rows(tibble(t = t_a, y = f1_true, parameter = "f1"),
                      tibble(t = t_a, y = mu_true, parameter = "mu"))

ggplot(combined_tib) +
    bayesplot::hline_at(censoring_val, lty = 3, alpha = 0.5) +
    ## geom_point(aes(x = t, y = y, shape = censored),
    ##            data = filter(sim_tib, group != 1),
    ##            alpha = 0.1) +
    geom_point(aes(x = t, y = y, shape = censored, col = "f1"),
               data = filter(sim_tib, group == 1),
               size = 2, alpha = 0.5) +
    geom_line(aes(x = tpred, y = mean, col = parameter, lty = "post mean")) +
    geom_ribbon(aes(x = tpred, ymin = lwr, ymax = upr, fill = parameter), alpha = 0.2) +
    geom_line(aes(x = t, y = y, col = parameter, lty = "true"), data = true_tib) +
    facet_wrap(~ model, ncol = 1) +
    ylab("") +
    xlab("t") +
    guides(fill = "none") +
    labs(color = "",
         linetype = "") +
    scale_shape_manual(values = c(20, 4),
                       labels = c("observed y1", "censored y1"),
                       name = "") +
    theme(legend.position = "bottom")

fig_dir <- "thesis/tex/figures/"
ggsave(paste0(fig_dir, "censoring_multi.pdf"),
       width = 17,
       height = 18,
       units = "cm")


