library(tidyverse)
options(mc.cores = parallel::detectCores())
source("code/rao/code/irregular/sim_irregular.R")
source("code//R/sim_list_to_tib.R")
sim_mod <- cmdstanr::cmdstan_model("code//rao/code//sim.stan")
stan_funs <- cmdstanr::cmdstan_model("code/rao//code/functions.stan",
                                     force_recompile = TRUE,
                                     compile_standalone = TRUE)$functions

## Settings for whether to include slow tests and plots
slow_tests <- FALSE
do_plots <- FALSE

if (slow_tests)
    log_lik_base_no_chol <- cmdstanr::cmdstan_model("code/stan//other_log_lik_implementations/regular//log_lik_base_no_chol.stan",
                                                    force_recompile = TRUE,
                                                    compile_standalone = TRUE)$functions

### Test all log-likelihoods in regular setup

cat("\n\n>>> Regular setting\n\n")

## Simulate and plot

n_a <- 5
n_b <- 7
J_a <- 13
J_b <- rep(J_a, n_b)
t_a <- seq(-5, 5, length.out = J_a)
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
                         sigma, sim_mod)
sim_tib <- sim_list_to_tib(sim_res[[2]])

if (do_plots)
    ggplot(sim_tib) +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed))

## Calculate log_lik using all implementations

test_dat_irreg <- sim_res[[1]]
n_group <- n_a + n_b
test_dat_reg <- with(test_dat_irreg, {
    y_obs_vec <- c(y_a_vec, y_b_vec)
    list(
        y_obs_vec = y_obs_vec,
        t_obs = t_a,
        n_obs = J_a,
        n_group = n_group,
        n_obs_total = J_a * n_group,
        magnitude_mu = magnitude_mu,
        length_scale_mu = length_scale_mu,
        magnitude_eta = magnitude_eta,
        length_scale_eta = length_scale_eta,
        sigma = sigma,
        zero_n_obs_total = rep(0, J_a * n_group),
        one_mat_n_group = matrix(1, ncol = n_group, nrow = n_group),
        y_obs = matrix(y_obs_vec, ncol = n_group)
    )})

cat("\nLog-likelihoods (should be equal)\n")

c(
    do.call(stan_funs$log_lik_reg_base,
            test_dat_reg[! names(test_dat_reg) %in% c("one_mat_n_group", "y_obs")]),
    do.call(stan_funs$log_lik_reg_smart,
            test_dat_reg[! names(test_dat_reg) == "zero_n_obs_total"]),
    do.call(stan_funs$log_lik_irreg_base,
            test_dat_irreg[! names(test_dat_irreg) %in% c("one_mat_n_a")]),
    do.call(stan_funs$log_lik_irreg_smart, test_dat_irreg)
)

if (slow_tests)
    do.call(log_lik_base_no_chol$log_lik_base_no_chol,
            test_dat_reg[! names(test_dat_reg) %in% c("one_mat_n_group", "y_obs")])

## Setup data for testing posterior draws

test_dat_reg_rng <- c(
    test_dat_reg[c("y_obs_vec", "y_obs", "t_obs", "n_obs",
                   "n_group", "n_obs_total",
                   "magnitude_mu", "length_scale_mu",
                   "magnitude_eta", "length_scale_eta",
                   "sigma", "one_mat_n_group")],
    list(t_pred = t_a,
         n_pred = J_a)
)

test_dat_irreg_rng <-
    c(
        test_dat_irreg,
        list(t_pred = test_dat_irreg$t_a,
             J_pred = J_a)
    )

if (do_plots) {

    ## Test base posterior draws

    mu_eta_pred_base <- do.call(stan_funs$mu_eta_reg_base_rng,
                                test_dat_reg_rng[! names(test_dat_reg_rng) %in% c("y_obs", "one_mat_n_group")])
    f_pred_base <- stan_funs$f_draws(mu_eta_pred_base)

    sim_tib |>
        mutate(f_pred = as.vector(f_pred_base)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)

    ## Test smart posterior

    mu_eta_pred_smart <- do.call(stan_funs$mu_eta_reg_smart_rng, test_dat_reg_rng)
    f_pred_smart <- stan_funs$f_draws(mu_eta_pred_smart)
    sim_tib |>
        mutate(f_pred = as.vector(f_pred_smart)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)

    ## Test iter posterior

    mu_eta_pred_iter <- do.call(stan_funs$mu_eta_reg_iter_rng, test_dat_reg_rng)
    f_pred_iter <- stan_funs$f_draws(mu_eta_pred_iter)
    sim_tib |>
        mutate(f_pred = as.vector(f_pred_iter)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)

    ## Test smart iter posterior

    mu_eta_pred_smart_iter <- do.call(stan_funs$mu_eta_reg_smart_iter_rng, test_dat_reg_rng)
    f_pred_smart_iter <- stan_funs$f_draws(mu_eta_pred_smart_iter)
    sim_tib |>
        mutate(f_pred = as.vector(f_pred_smart_iter)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)

    ## Test base irregular posterior

    mu_eta_pred_irreg_base <- do.call(stan_funs$mu_eta_irreg_base_rng,
                                      test_dat_irreg_rng[! names(test_dat_irreg_rng) %in% c("y_a", "one_mat_n_a")])
    f_pred_irreg_base <- stan_funs$f_draws(mu_eta_pred_irreg_base)
    sim_tib |>
        mutate(f_pred = as.vector(f_pred_irreg_base)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)

    ## Test smart irregular posterior

    mu_eta_pred_irreg_smart <- do.call(stan_funs$mu_eta_irreg_smart_rng,
                                       test_dat_irreg_rng)
    f_pred_irreg_smart <- stan_funs$f_draws(mu_eta_pred_irreg_smart)
    sim_tib |>
        mutate(f_pred = as.vector(f_pred_irreg_smart)) |>
        ggplot() +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = f_pred, col = group_id), lty = 2)
}

## Compare posteriors

cat("\nPosterior differences (should be zero).\n")

base_post <- do.call(stan_funs$mu_eta_reg_base_post,
                     test_dat_reg_rng[! names(test_dat_reg_rng) %in% c("y_obs", "one_mat_n_group")])
smart_post <- do.call(stan_funs$mu_eta_reg_smart_post, test_dat_reg_rng)
iter_post_chol <- do.call(stan_funs$mu_eta_reg_iter_post_chol, test_dat_reg_rng)
smart_iter_post_chol <- do.call(stan_funs$mu_eta_reg_smart_iter_post_chol, test_dat_reg_rng)
irreg_base_post <- do.call(stan_funs$mu_eta_irreg_base_post,
                           test_dat_irreg_rng[! names(test_dat_irreg_rng) %in% c("y_a", "one_mat_n_a")])
irreg_smart_post <- do.call(stan_funs$mu_eta_irreg_smart_post, test_dat_irreg_rng)

## Comparison of mean(mu), cov(mu) and mean(eta)
description <- c("mu_mean", "mu_cov", "eta_mean", "eta_cov")
for (i in 1:3) {
    print(description[i])
    print(range(base_post[[i]] - smart_post[[i]]))
    print(range(base_post[[i]] - iter_post_chol[[i]]))
    print(range(base_post[[i]] - smart_iter_post_chol[[i]]))
    print(range(base_post[[i]] - irreg_base_post[[i]]))
    print(range(base_post[[i]] - irreg_smart_post[[i]]))
}

## Comparison of cov(eta)
print("eta_cov")
print(range(base_post[[4]] - smart_post[[4]]))
print(range(base_post[[4]] - tcrossprod(iter_post_chol[[4]])))
print(range(base_post[[4]] - tcrossprod(smart_iter_post_chol[[4]])))
print(range(base_post[[4]] - irreg_base_post[[4]]))
print(range(base_post[[4]] - tcrossprod(irreg_smart_post[[4]])))

if (slow_tests && do_plots) {
    ## Sanity check that the actual simulations also seem to be from same
    ## distribution. See that they have approximately same mean and variance in each
    ## coordinate
    rep_smart <- replicate(100, do.call(stan_funs$mu_eta_reg_smart_rng, test_dat_reg_rng_with_mat))
    rep_iter <- replicate(100, do.call(stan_funs$mu_eta_reg_iter_rng, test_dat_reg_rng))
    rep_smart_iter <- replicate(100, do.call(stan_funs$mu_eta_reg_smart_iter_rng, test_dat_reg_rng))
    rep_base <- replicate(100, do.call(stan_funs$mu_eta_reg_base_rng,
                                       test_dat_reg_rng[! names(test_dat_reg_rng) == "y_obs"]))
    rep_irreg_base <- replicate(100,
                                do.call(stan_funs$mu_eta_irreg_base_rng,
                                        test_dat_irreg_rng[! names(test_dat_irreg_rng) %in% c("y_a", "one_mat_n_a")]))

    means_smart <- apply(rep_smart, c(1, 2), mean)
    means_iter <- apply(rep_iter, c(1, 2), mean)
    means_smart_iter <- apply(rep_smart_iter, c(1, 2), mean)
    means_base <- apply(rep_base, c(1, 2), mean)
    means_irreg_base <- apply(rep_irreg_base, c(1, 2), mean)
    vars_smart <- apply(rep_smart, c(1, 2), var)
    vars_iter <- apply(rep_iter, c(1, 2), var)
    vars_smart_iter <- apply(rep_smart_iter, c(1, 2), var)
    vars_base <- apply(rep_base, c(1, 2), var)
    vars_irreg_base <- apply(rep_irreg_base, c(1, 2), var)

    hist(means_base - means_smart)
    hist(means_base - means_iter)
    hist(means_base - means_smart_iter)
    hist(means_base - means_smart_iter)
    hist(vars_base - vars_smart)
    hist(vars_base - vars_iter)
    hist(vars_base - vars_smart_iter)
    hist(vars_base - vars_smart_iter)
}

### Example of fit, to show that we can get correct parameter values

if (slow_tests) {

    mod <- cmdstanr::cmdstan_model("code/rao/code/regular/smart_no_gq.stan")
    fit <- mod$sample(data = c(test_dat_reg_rng))
    fit

    fit$cmdstan_summary()

}

## We do between 3 and 15 leapfrog steps per draw (median 7).

### Check fit with rng generated quantities included

if (slow_tests && do_plots) {

    ## smart generated quantities
    mod <- cmdstanr::cmdstan_model("code/rao/code/regular/smart_smart_gq.stan")
    fit <- mod$sample(data = c(test_dat_reg_rng),
                      chains = 1,
                      iter_warmup = 1000,
                      iter_sampling = 1000)
    draws <- posterior::as_draws_rvars(fit$draws())
    mu <- posterior::draws_of(draws$mu)
    ggplot(sim_tib) +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = mu), tibble(t = t_a, mu = colMeans(mu)),
                  size = 2)

}

if (slow_tests) {

    ## base generated quantities
    mod <- cmdstanr::cmdstan_model("code/rao/code/regular/smart_base_gq.stan")
    fit <- mod$sample(data = c(test_dat_reg_rng))
    draws <- posterior::as_draws_rvars(fit$draws())
    mu <- posterior::draws_of(draws$mu)
    ggplot(sim_tib) +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed)) +
        geom_line(aes(x = t, y = mu), tibble(t = t_a, mu = colMeans(mu)),
                  size = 2)

}

### Test in irregular setup

cat("\n\n>>> Irregular setting\n\n")

## Simulate and plot to test simulation function
n_a <- 5
n_b <- 7
n <- n_a + n_b
J_a <- 13
J_b <- sample(3:17, n_b, replace = TRUE)

t_a <- runif(J_a, -5, 5)
t_b <- list()
for (i in 1:n_b) {
    t_b <- c(t_b, list(runif(J_b[i], -5, 5)))
}

## t_a <- seq(-5, 5, length.out = J_a)
## t_b <- list()
## for (i in 1:n_b) {
##     t_b <- c(t_b, list(seq(-5, 5, length.out = J_b[i])))
## }

magnitude_mu <- 1
length_scale_mu <- 1
magnitude_eta <- 0.5
length_scale_eta <- 0.5
sigma <- 0.2
sim_res <- sim_irregular(n_a, J_a, t_a,
                         n_b, J_b, t_b,
                         magnitude_mu, length_scale_mu,
                         magnitude_eta, length_scale_eta,
                         sigma, sim_mod)
sim_tib <- sim_list_to_tib(sim_res[[2]])

if (do_plots)
    ggplot(sim_tib) +
        geom_line(aes(x = t, y = f, col = group_id)) +
        geom_point(aes(x = t, y = f, col = group_id),
                   filter(sim_tib, observed))

## Calculate log_lik to test all log_lik implementations
cat("\nLog-likelihoods (should be equal)\n")
test_dat <- sim_res[[1]]
c(do.call(stan_funs$log_lik_irreg_smart, test_dat),
  do.call(stan_funs$log_lik_irreg_base,
          test_dat[! names(test_dat) %in% c("one_mat_n_a")])
)

## rng test
J_pred <- 11

test_dat_pred <-
    c(
        test_dat,
        list(t_pred = runif(J_pred, -5, 5),
             J_pred = J_pred)
    )

cat("\nPosterior differences (should be zero)\n")
post_smart <- do.call(stan_funs$mu_eta_irreg_smart_post, test_dat_pred)
post_base <- do.call(stan_funs$mu_eta_irreg_base_post,
                     test_dat_pred[! names(test_dat_pred) %in% c("y_a", "one_mat_n_a")])
post_completely_irreg <-
    stan_funs$mu_eta_completely_irreg_post(
                  test_dat_pred$y_ab_vec,
                  c(rep(test_dat_pred$t_a, test_dat_pred$n_a),
                    test_dat_pred$t_b),
                  test_dat_pred$t_pred,
                  test_dat_pred$n_a + test_dat_pred$n_b,
                  test_dat_pred$J_pred,
                  c(rep(test_dat_pred$J_a, test_dat_pred$n_a),
                    test_dat_pred$J_b),
                  c(0, cumsum(c(rep(test_dat_pred$J_a, test_dat_pred$n_a),
                                test_dat_pred$J_b))),
                  test_dat_pred$magnitude_mu,
                  test_dat_pred$length_scale_mu,
                  test_dat_pred$magnitude_eta,
                  test_dat_pred$length_scale_eta,
                  test_dat_pred$sigma
              )

for (i in 1:3) {
    print(range(post_base[[i]] - post_smart[[i]]))
    print(range(post_base[[i]] - post_completely_irreg[[i]]))
}

print(range(post_base[[4]] - tcrossprod(post_smart[[4]])))
print(range(post_base[[4]] - post_completely_irreg[[4]]))

### Example of fit, to show that we can get correct parameter values

if (slow_tests) {

    mod <- cmdstanr::cmdstan_model("code/rao/code//irregular/smart_no_gq.stan")
    fit <- mod$sample(data = test_dat)
    fit

}

## The estimates fit nicely with the true values.
