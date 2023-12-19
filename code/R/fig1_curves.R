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

n <- 50
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
    theme(legend.position = "none",
          axis.title.y = element_blank())

fig_dir <- "thesis/tex/figures/"
ggsave(paste0(fig_dir, "model_example.pdf"),
       width = 18,
       height = 9,
       units = "cm")

ggplot() +
    geom_line(aes(x = t, y = f, color = i), f_tib,
              alpha = 0.7) +
#    geom_line(aes(x = t, y = mu_true), linewidth = 0.2) +
    theme(legend.position = "none") +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )

ggsave(paste0(fig_dir, "front_page_curves.pdf"),
       width = 36,
       height = 20,
       units = "cm")

