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

n <- 3
n_a <- n-1
n_b <- 1
J <- 200
J_a <- J
J_b <- rep(J, n_b)
t <- t_a <- seq(-10, 10, length.out = J_a)
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

f_true <- sim_res[[2]]$f[1, , 1:n]

f_tib <- tibble()
for (i in 1:n) {
    f_tib <- bind_rows(f_tib,
                       tibble(f = f_true[, i], i, t))
}
f_tib <- mutate(f_tib, i = as.factor(i))

ggplot() +
    geom_line(aes(x = t, y = f, color = i), f_tib,
              alpha = 0.7) +
    geom_line(aes(x = t, y = mu_true), linewidth = 1) +
    theme(legend.position = "none") +
    ylab("")

test_dat_irreg <- sim_res[[1]]
n_group <- n_a + n_b

i_obs <- seq(1, length(t), by = 4)
t_obs <- t[i_obs]
J_obs <- length(t_obs)
y_obs_vec <- c(as.vector(test_dat_irreg$y_a[i_obs, ]), test_dat_irreg$y_b_vec[i_obs])

test_dat_reg <- with(test_dat_irreg, {
    list(
        y_obs_vec = y_obs_vec,
        t_obs = t_obs,
        n_obs = J_obs,
        n_group = n_group,
        n_obs_total = J_obs * n_group,
        magnitude_mu = magnitude_mu,
        length_scale_mu = length_scale_mu,
        magnitude_eta = magnitude_eta,
        length_scale_eta = length_scale_eta,
        sigma = sigma,
        zero_n_obs_total = rep(0, J_obs * n_group),
        one_mat_n_group = matrix(1, ncol = n_group, nrow = n_group),
        y_obs = matrix(y_obs_vec, ncol = n_group)
    )})
test_dat_reg_rng <- c(
    test_dat_reg[c("y_obs_vec", "y_obs", "t_obs", "n_obs",
                   "n_group", "n_obs_total",
                   "magnitude_mu", "length_scale_mu",
                   "magnitude_eta", "length_scale_eta",
                   "sigma", "one_mat_n_group")],
    list(t_pred = t_a,
         n_pred = J_a)
)

sim_tib <- sim_list_to_tib(sim_res[[2]])
mu_eta_pred_smart_iter <- do.call(stan_funs$mu_eta_reg_smart_iter_rng, test_dat_reg_rng)
f_pred_smart_iter <- stan_funs$f_draws(mu_eta_pred_smart_iter)

## smart generated quantities
mod <- cmdstanr::cmdstan_model("code/rao/code/regular/smart_smart_gq.stan")
fit <- mod$sample(data = c(test_dat_reg_rng),
                  chains = 1,
                  iter_warmup = 1000,
                  iter_sampling = 1000,
                  seed = 1723)
draws <- posterior::as_draws_rvars(fit$draws())
f_pred <- posterior::draws_of(draws$f_pred)
mu <- posterior::draws_of(draws$mu)

mu_ci <- apply(mu, 2, \(x) quantile(x, c(0.05, 0.95)))
mu_mean <- apply(mu, 2, mean)

f_mean <- apply(f_pred, c(2, 3), mean)
f_ci <- apply(f_pred, c(2, 3), \(x) quantile(x, c(0.05, 0.95)))

f_pred_tib <- tibble()
for (i in 1:n) {
    f_pred_tib <- bind_rows(f_pred_tib,
                            tibble(t, group = i,
                                   mean = f_mean[, i],
                                   lwr = f_ci[1, , i],
                                   upr = f_ci[2, , i]))
}
f_pred_tib <- mutate(f_pred_tib, group = as.factor(group))

ggplot() +
    geom_line(aes(x = t, y = mean, col = paste0("f", group), lty = "post mean"), f_pred_tib) +
    geom_ribbon(aes(x = t, ymin = lwr, ymax = upr, fill = paste0("f", group)), f_pred_tib,
                alpha = 0.1) +
    geom_line(aes(x = t, y = f, col = paste0("f", i), lty = "truth"), f_tib) +
    geom_line(aes(x = t, y = mu_mean, col = "mu", lty = "post mean")) +
    geom_ribbon(aes(x = t, ymin = mu_ci[1, ], ymax = mu_ci[2, ], fill = "mu"),
                alpha = 0.2) +
    geom_line(aes(x = t, mu_true, col = "mu", lty = "truth")) +
    theme(axis.title.y = element_blank()) +
    scale_color_discrete(name = "") +
    scale_fill_discrete(name = "") +
    scale_linetype_discrete(name = "")
    
fig_dir <- "thesis/tex/figures/"
ggsave(paste0(fig_dir, "model_fit.pdf"),
       width = 18,
       height = 9,
       units = "cm")

