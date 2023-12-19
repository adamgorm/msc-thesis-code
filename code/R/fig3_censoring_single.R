library(tidyverse)
options(mc.cores = parallel::detectCores())

theme_set(theme_bw())

stan_funs <- cmdstanr::cmdstan_model("code/rao//code/functions.stan",
                                     force_recompile = TRUE,
                                     compile_standalone = TRUE)$functions

log_lik_gp <- cmdstanr::cmdstan_model("code/stan//log_lik_single.stan",
                                      force_recompile = TRUE,
                                      compile_standalone = TRUE)$functions

sim_gp <- cmdstanr::cmdstan_model("code/stan/sim_gp_single.stan")
fit_gp <- cmdstanr::cmdstan_model("code/stan/fit_gp_single.stan")

N <- 50
t <- seq(-10, 10, length.out = N)
Nf <- 200
tf <- seq(-10, 10, length.out = Nf)
magnitude <- 1
length_scale <- 1
sigma <- 0.2
sim_dat_list <- list(
    Ny = N,
    ty = t,
    Nf = Nf,
    tf = tf,
    magnitude = magnitude,
    length_scale = length_scale,
    sigma = sigma
)

gp_samples <- sim_gp$sample(sim_dat_list,
                            fixed_param = TRUE,
                            seed = 824)
gp_draw <- as.vector(gp_samples$draws("y", format = "matrix")[1, ])
gp_draw_f <- as.vector(gp_samples$draws("f", format = "matrix")[1, ])
censor_limit <- rep_along(gp_draw, 0.75)
censored <- gp_draw > censor_limit
to <- t[!censored]
yo <- gp_draw[!censored]
tc <- t[censored]
yc <- censor_limit[censored]
ggplot() +
    geom_point(aes(x = t, y = gp_draw, col = censored)) +
    bayesplot::hline_at(censor_limit, lty = 2) +
    geom_line(aes(x = tf, y = gp_draw_f))

marg_log_lik <- function(params)
{
    prep_list <- stan_funs$prep_cens_log_lik(yo, to, tc,
                                             params[[1]], params[[2]], params[[3]])
    prep_list[[1]] + log(mvtnorm::pmvnorm(lower = yc,
                                          mean = prep_list[[2]],
                                          sigma = prep_list[[3]]))
}

set.seed(1311)
marg_log_lik(c(magnitude, length_scale, sigma))
optim_res <- optim(c(magnitude, length_scale, sigma),
                   marg_log_lik,
                   control = list(fnscale = -1))
param_max <- optim_res$par

## print(param_max)

n_pred <- 200
t_pred <- seq(-10, 10, length.out = n_pred)

post_list <- stan_funs$prep_cond_post(yo, yc,
                                      to, tc, t_pred,
                                      param_max[1], param_max[2], param_max[3])

n_draw <- 1000
set.seed(903)
q_draws <- t(mvtnorm::rmvnorm(n_draw, mean = rep(0, n_pred), sigma = post_list[[5]]))
p_draws <- t(TruncatedNormal::rtmvnorm(n_draw, rep(0, nrow(post_list[[4]])),
                                       post_list[[4]], yc - post_list[[3]]))
f_draws <- post_list[[1]] +
    post_list[[2]] %*% p_draws +
    as.vector(q_draws)
f_draws <- t(f_draws)
post_mean <- colMeans(f_draws)
post_95 <- apply(f_draws, 2, \(x) quantile(x, c(0.025, 0.975)))

log_lik_oracle <- function(params)
    log_lik_gp$log_lik_single(gp_draw, N, t, params[1], params[2], params[3])

log_lik_oracle(c(1, 1, 0.2))

set.seed(1418)
optim_res_oracle <- optim(c(magnitude, length_scale, sigma),
                          log_lik_oracle,
                          control = list(fnscale = -1))
param_max_oracle <- optim_res_oracle$par

fit_dat_oracle <- list(
    n_obs = N,
    n_pred = n_pred,
    t_obs = t,
    y_obs = gp_draw,
    t_pred = t_pred,
    magnitude = param_max_oracle[1],
    length_scale = param_max_oracle[2],
    sigma = param_max_oracle[3]
)

post_oracle <- fit_gp$sample(fit_dat_oracle, seed = 904, fixed_param = TRUE,
                             chains = 1, iter_sampling = n_draw)
draws_oracle <- matrix(post_oracle$draws(), nrow = n_draw, ncol = n_pred)
post_mean_oracle <- colMeans(draws_oracle)
post_95_oracle <- apply(draws_oracle, 2, \(x) quantile(x, c(0.025, 0.975)))

fit_dat_include <- list(
    n_obs = N,
    n_pred = n_pred,
    t_obs = c(to, tc),
    y_obs = c(yo, yc),
    t_pred = t_pred
)

log_lik_include <- function(params)
    log_lik_gp$log_lik_single(fit_dat_include$y_obs,
                              fit_dat_include$n_obs,
                              fit_dat_include$t_obs,
                              params[1], params[2], params[3])

log_lik_include(c(1, 1, 0.2))

set.seed(1422)
optim_res_include <- optim(c(magnitude, length_scale, sigma),
                          log_lik_include,
                          control = list(fnscale = -1))
param_max_include <- optim_res_include$par

fit_dat_include$magnitude <- param_max_include[1]
fit_dat_include$length_scale <- param_max_include[2]
fit_dat_include$sigma <- param_max_include[3]

post_include <- fit_gp$sample(fit_dat_include, seed = 918, fixed_param = TRUE,
                             chains = 1, iter_sampling = n_draw)
draws_include <- matrix(post_include$draws(), nrow = n_draw, ncol = n_pred)
post_mean_include <- colMeans(draws_include)
post_95_include <- apply(draws_include, 2, \(x) quantile(x, c(0.025, 0.975)))

fit_dat_exclude <- list(
    n_obs = length(to),
    n_pred = n_pred,
    t_obs = to,
    y_obs = yo,
    t_pred = t_pred
)

log_lik_exclude <- function(params)
    log_lik_gp$log_lik_single(fit_dat_exclude$y_obs,
                              fit_dat_exclude$n_obs,
                              fit_dat_exclude$t_obs,
                              params[1], params[2], params[3])

log_lik_exclude(c(1, 1, 0.2))

set.seed(1423)
optim_res_exclude <- optim(c(magnitude, length_scale, sigma),
                          log_lik_exclude,
                          control = list(fnscale = -1))
param_max_exclude <- optim_res_exclude$par

fit_dat_exclude$magnitude <- param_max_exclude[1]
fit_dat_exclude$length_scale <- param_max_exclude[2]
fit_dat_exclude$sigma <- param_max_exclude[3]

post_exclude <- fit_gp$sample(fit_dat_exclude, seed = 938, fixed_param = TRUE,
                             chains = 1, iter_sampling = n_draw)
draws_exclude <- matrix(post_exclude$draws(), nrow = n_draw, ncol = n_pred)
post_mean_exclude <- colMeans(draws_exclude)
post_95_exclude <- apply(draws_exclude, 2, \(x) quantile(x, c(0.025, 0.975)))

get_post_tib <- function(post_mean, post_95, model)
{
    tibble(t_pred,
           mean = post_mean,
           lwr = post_95[1,],
           upr = post_95[2,],
           model = model)
}

points_tib <- tibble(t = t,
                     y = gp_draw,
                     censored = censored)

f_tib <- tibble(t = tf,
                f = gp_draw_f)

bind_rows(get_post_tib(post_mean, post_95, "censoring model"),
          get_post_tib(post_mean_oracle, post_95_oracle, "oracle"),
          get_post_tib(post_mean_include, post_95_include, "include"),
          get_post_tib(post_mean_exclude, post_95_exclude, "exclude")) |>
    mutate(model = factor(model, levels = c("oracle",
                                            "censoring model",
                                            "include",
                                            "exclude"))) |>
    ggplot() +
    geom_hline(aes(yintercept = 0.75), lty = 2, alpha = 0.5) +
    geom_line(aes(x = t, y = f, col = "truth"), data = f_tib) +
    geom_point(aes(x = t, y = y, shape = censored), data = points_tib) +
    geom_line(aes(x = t_pred, y = mean, col = "post mean")) +
    geom_ribbon(aes(x = t_pred, ymin = lwr, ymax = upr, fill = "post mean"),
                alpha = 0.2) +
    facet_wrap(~ model, ncol = 1) +
    ylab("") +
    xlab("t") +
    guides(fill = "none") +
    labs(linetype = "") +
    scale_color_manual(values = c(2, 1),
                          labels = c("post mean", "truth"),
                          name = "") +
    scale_shape_manual(values = c(20, 4),
                       labels = c("observed", "censored"),
                       name = "") +
    theme(legend.position = "bottom")

fig_dir <- "thesis/tex/figures/"
ggsave(paste0(fig_dir, "censoring_single.pdf"),
       width = 17,
       height = 18,
       units = "cm")

