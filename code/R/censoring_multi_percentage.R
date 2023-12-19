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

cens_vec <- replicate(1000, {
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

    f1_true <- sim_res[[2]]$f[1, , 1]

    t_obs <- t_a[(0:J) * scale_factor_f + 1]

    sim_tib <- filter(sim_tib_full, t %in% t_obs)

    censoring_val <- 0.75
    sim_tib |>
        mutate(censored = y > censoring_val,
               y_obs = ifelse(censored, censoring_val, y),
               group = as.factor(group)) ->
        sim_tib

    mean(sim_tib$censored)
})

mean(cens_vec)
quantile(cens_vec, c(0.1, 0.9))
hist(cens_vec)
