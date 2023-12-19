library(tidyverse)
library(patchwork)
theme_set(theme_bw())
fig_width <- 25
fig_height <- 10
fig_units <- "cm"
fig_dir <- "thesis/tex/benchmarks/"
fig_format <- "pdf"
save_fig <- function(fig_name) {
    ggsave(paste0(fig_dir, fig_name, ".", fig_format),
           width = fig_width,
           height = fig_height,
           units = fig_units)
}
fig_width_full_page <- 17
fig_height_full_page <- 18
save_fig_full_page <- function(fig_name) {
    ggsave(paste0(fig_dir, fig_name, ".", fig_format),
           width = fig_width_full_page,
           height = fig_height_full_page,
           units = fig_units)
}
time_lab <- ylab("time (minutes)")

##### Regular

time_dir <- "code/rao/results/regular/time_mins/"
time_file_path <- function(n_obs, n_group, implementation)
    sprintf("%s/%s/%d_%d.txt", time_dir, implementation, n_obs, n_group)

res_tib_loglik <- dplyr::tibble()
for (imp in c("smart", "base")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     imp)))
        res_tib_loglik <- dplyr::bind_rows(
                              res_tib_loglik,
                              dplyr::tibble(n_obs = sim_setting[1],
                                            n_group = sim_setting[2],
                                            loglik_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_loglik |>
    dplyr::group_by(n_obs, n_group, loglik_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "loglik", gq_implementation = "none") ->
    res_tib_loglik

time_dir <- "code/rao/results/regular/time_mins_rng/"
time_file_path <- function(n_obs, n_group, n_pred, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, n_obs, n_group, n_pred)

res_tib_gq <- dplyr::tibble()
for (imp in c("smart", "base", "smart_iter", "iter")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     imp)))
        res_tib_gq <- dplyr::bind_rows(
                              res_tib_gq,
                              dplyr::tibble(n_obs = sim_setting[1],
                                            n_group = sim_setting[2],
                                            n_pred = sim_setting[3],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_gq |>
    dplyr::group_by(n_obs, n_group, n_pred, gq_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "generated quantities", loglik_implementation = "none") ->
    res_tib_gq

time_dir <- "code/rao/results/regular/time_mins_full_fit/"
time_file_path <- function(n_obs, n_group, n_pred, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, n_obs, n_group, n_pred)

res_tib_full <- dplyr::tibble()
for (imp in c("no_gq", "smart", "base", "smart_iter", "iter")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     imp)))
        res_tib_full <- dplyr::bind_rows(
                              res_tib_full,
                              dplyr::tibble(n_obs = sim_setting[1],
                                            n_group = sim_setting[2],
                                            n_pred = sim_setting[3],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_full |>
    dplyr::group_by(n_obs, n_group, n_pred, gq_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "full fit", loglik_implementation = "smart") ->
    res_tib_full

time_dir <- "code/rao/pre_tuple_results/regular/time_mins_full_fit_base/"
time_file_path <- function(n_obs, n_group, n_pred, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, n_obs, n_group, n_pred)

res_tib_fullbase <- dplyr::tibble()
## Note: We actually only use no_gq in this setting but it doesn't matter here
for (imp in c("no_gq", "smart", "base")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     imp)))
        res_tib_fullbase <- dplyr::bind_rows(
                              res_tib_fullbase,
                              dplyr::tibble(n_obs = sim_setting[1],
                                            n_group = sim_setting[2],
                                            n_pred = sim_setting[3],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}

time_dir <- "code/rao/results/regular/time_mins_full_fit_base/"
time_file_path <- function(n_obs, n_group, n_pred, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, n_obs, n_group, n_pred)

for (imp in c("no_gq", "smart", "base", "smart_iter", "iter")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     imp)))
        res_tib_fullbase <- dplyr::bind_rows(
                              res_tib_fullbase,
                              dplyr::tibble(n_obs = sim_setting[1],
                                            n_group = sim_setting[2],
                                            n_pred = sim_setting[3],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}


res_tib_fullbase |>
    dplyr::group_by(n_obs, n_group, n_pred, gq_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "full fit", loglik_implementation = "base") ->
    res_tib_fullbase

res_tib <- bind_rows(res_tib_loglik,
                     res_tib_gq,
                     res_tib_full,
                     res_tib_fullbase) |>
    mutate(gq_implementation = ifelse(gq_implementation == "no_gq",
                                      "none",
                                      gq_implementation))

res_tib |>
    mutate(type = factor(type, levels = c("loglik",
                                          "generated quantities",
                                          "full fit"))) |>
    mutate(type = fct_recode(type,
                             "log-likelihood" = "loglik",
                             "full fit with HMC" = "full fit")) ->
    res_tib

lty_manual <- 
    scale_linetype_manual(values = c("base" = "dashed",
                                     "smart" = "solid",
                                     "none" = "dotdash"))
col_manual_vec <- c("base" = "black",
                    "none" = "chartreuse3",
                    "smart_iter" = "coral2",
                    "smart" = "blue")

col_manual <- scale_color_manual(values = col_manual_vec)
reg_linewidth <- 1

reg_n <- filter(res_tib, n_obs == 100, is.na(n_pred) | n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_group, y = time_mins,
                  lty = loglik_implementation, col = gq_implementation),
              linewidth = reg_linewidth) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    xlab("n") +
    ggtitle(expression(paste("J = 100, ", J[p], " = 100"))) +
    theme(legend.position = "none") +
    lty_manual +
    col_manual +
    theme(plot.title = element_text(size=11))

## save_fig("reg_n.pdf")

reg_J <- filter(res_tib, n_group == 100, is.na(n_pred) | n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_obs, y = time_mins,
                  lty = loglik_implementation, col = gq_implementation),
              linewidth = reg_linewidth) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    xlab("J") +
    ggtitle(expression(paste("n = 100, ", J[p], " = 100"))) +
    theme(legend.position = "none") +
    lty_manual +
    col_manual +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(plot.title = element_text(size=11))

## save_fig("reg_J.pdf")

tib_flat_loglik <- filter(res_tib, n_obs == 100, n_group == 100, type == "log-likelihood")
tib_flat_loglik <- bind_rows(tib_flat_loglik, tib_flat_loglik)
tib_flat_loglik$n_pred <- rep(c(25, 300), each = 2)

reg_Jp <- filter(res_tib, n_group == 100, n_obs == 100, !is.na(n_pred)) |>
    bind_rows(tib_flat_loglik) |>
    ggplot() +
    geom_line(aes(x = n_pred, y = time_mins,
                  lty = loglik_implementation, col = gq_implementation),
              linewidth = reg_linewidth) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    xlab(expression(J[p])) +
    ggtitle("n = 100, J = 100") +
    theme(legend.position = "none") +
    lty_manual +
    col_manual +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(plot.title = element_text(size=11))

## save_fig("reg_Jp.pdf")

(reg_n +
 guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
 labs(col = "gq implementation",
      lty = "loglik implementation")) /
 (reg_J + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  labs(col = "gq implementation",
       lty = "loglik implementation")) /
 (reg_Jp + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  labs(col = "gq implementation",
       lty = "loglik implementation")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "vertical")

save_fig_full_page("regular")

##### Irregular

time_dir <- "code/rao/results/irregular/time_mins"
time_file_path <- function(J, n_a, n_b, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, J, n_a, n_b)

res_tib_loglik <- dplyr::tibble()
for (imp in c("smart", "base")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     imp)))
        res_tib_loglik <- dplyr::bind_rows(
                              res_tib_loglik,
                              dplyr::tibble(J = sim_setting[1],
                                            n_a = sim_setting[2],
                                            n_b = sim_setting[3],
                                            loglik_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_loglik |>
    dplyr::group_by(J, n_a, n_b, loglik_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "loglik", gq_implementation = "none") ->
    res_tib_loglik

time_dir <- "code/rao/results/irregular/time_mins_rng"
time_file_path <- function(J, n_a, n_b, J_pred, implementation)
    sprintf("%s/%s/%d_%d_%d_%d.txt", time_dir, implementation, J, n_a, n_b, J_pred)

res_tib_gq <- dplyr::tibble()
for (imp in c("no_gq", "smart", "base")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     sim_setting[4],
                                                     imp)))
        res_tib_gq <- dplyr::bind_rows(
                              res_tib_gq,
                              dplyr::tibble(J = sim_setting[1],
                                            n_a = sim_setting[2],
                                            n_b = sim_setting[3],
                                            J_pred = sim_setting[4],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_gq |>
    dplyr::group_by(J, n_a, n_b, J_pred, gq_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "generated quantities", loglik_implementation = "none") ->
    res_tib_gq

time_dir <- "code/rao/results/irregular/time_mins_full_fit"
time_file_path <- function(J, n_a, n_b, J_pred, implementation)
    sprintf("%s/%s/%d_%d_%d_%d.txt", time_dir, implementation, J, n_a, n_b, J_pred)

res_tib_full <- dplyr::tibble()
for (imp in c("no_gq", "smart", "base")) {
    for (filename in list.files(paste0(time_dir, "/", imp))) {
        sim_setting <- as.numeric(stringr::str_match(filename, "^([0-9]*)_([0-9]*)_([0-9]*)_([0-9]*).txt")[-1])
        times <- as.numeric(readLines(time_file_path(sim_setting[1],
                                                     sim_setting[2],
                                                     sim_setting[3],
                                                     sim_setting[4],
                                                     imp)))
        res_tib_full <- dplyr::bind_rows(
                              res_tib_full,
                              dplyr::tibble(J = sim_setting[1],
                                            n_a = sim_setting[2],
                                            n_b = sim_setting[3],
                                            J_pred = sim_setting[4],
                                            gq_implementation = imp,
                                            time_mins = times)
        )
    }
}

res_tib_full |>
    dplyr::group_by(J, n_a, n_b, J_pred, gq_implementation) |>
    dplyr::summarize(time_mins = mean(time_mins)) |>
    mutate(type = "full fit", loglik_implementation = "smart") ->
    res_tib_full

res_tib <- bind_rows(res_tib_loglik,
                     res_tib_gq,
                     res_tib_full) |>
    mutate(gq_implementation = ifelse(gq_implementation == "no_gq",
                                      "none",
                                      gq_implementation))

res_tib <- mutate(res_tib,
                  gq_implementation = ifelse(
                      gq_implementation == "smart",
                      "smart_iter",
                      gq_implementation
                  )
                  )

res_tib |>
    mutate(type = factor(type, levels = c("loglik",
                                          "generated quantities",
                                          "full fit"))) |>
    mutate(type = fct_recode(type,
                             "log-likelihood" = "loglik",
                             "full fit with HMC" = "full fit")) ->
    res_tib

lty_manual_irreg <- lty_manual

col_manual_irreg <- scale_color_manual(values = c("base" = col_manual_vec[["base"]],
                                                  "none" = col_manual_vec[["none"]],
                                                  "smart_iter" = col_manual_vec[["smart_iter"]]))

irreg_linewidth <- 1.2
irreg_alpha <- 0.8

### Smaller scale for full fit

irreg_nanb <- filter(res_tib, J == 100, n_a + n_b == 50, is.na(J_pred) | J_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_b, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation),
              linewidth = irreg_linewidth,
              alpha = irreg_alpha) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste(n[a], " + ", n[b], " = 50, J = 100, ", J[p], " = 100"))) +
    lty_manual_irreg +
    col_manual_irreg +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size=11)) +
    xlab(expression(n[b]))

irreg_nb <- filter(res_tib, J == 100, n_a == 50, is.na(J_pred) | J_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_b, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation),
              linewidth = irreg_linewidth,
              alpha = irreg_alpha) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste(n[a], " = 50, J = 100, ", J[p], " = 100"))) +
    lty_manual_irreg +
    col_manual_irreg +
    theme(legend.position = "none") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(plot.title = element_text(size=11)) +
    xlab(expression(n[b]))

irreg_J <- filter(res_tib, n_b == 10, n_a == 90, is.na(J_pred) | J_pred == 50) |>
    ggplot() +
    geom_line(aes(x = J, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation),
              linewidth = irreg_linewidth,
              alpha = irreg_alpha) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste(n[a], " = 90, ", n[b], " = 10, ", J[p], " = 50"))) +
    lty_manual_irreg +
    col_manual_irreg +
    theme(legend.position = "none") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(plot.title = element_text(size=11))

tib_flat_loglik_irreg <- filter(res_tib, n_b == 10, n_a == 90, J == 50, type == "log-likelihood")
tib_flat_loglik_irreg <- bind_rows(tib_flat_loglik_irreg, tib_flat_loglik_irreg)
tib_flat_loglik_irreg$J_pred <- rep(c(10, 150), each = 2)

irreg_Jp <- filter(res_tib, n_b == 10, n_a == 90, J == 50, !is.na(J_pred)) |>
    bind_rows(tib_flat_loglik_irreg) |>
    ggplot() +
    geom_line(aes(x = J_pred, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation),
              linewidth = irreg_linewidth,
              alpha = irreg_alpha) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste(n[a], " = 90, ", n[b], " = 10, J = 50"))) +
    lty_manual_irreg +
    col_manual_irreg +
    theme(legend.position = "none") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(plot.title = element_text(size=11)) +
    xlab(expression(J[p]))

(irreg_nanb + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
 labs(col = "gq implementation",
      lty = "loglik implementation")) /
    (irreg_nb + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
     labs(col = "gq implementation",
      lty = "loglik implementation")) /
    (irreg_J + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
     labs(col = "gq implementation",
      lty = "loglik implementation")) /
    (irreg_Jp + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
     labs(col = "gq implementation",
      lty = "loglik implementation")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box="vertical")

save_fig_full_page("irregular")

### Larger scale, not involving full fit

filter(res_tib, J == 100, n_a + n_b == 100, is.na(J_pred) | J_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_b, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation)) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste("n_a + n_b = 100, J = 100, ", J[p], " = 100"))) +
    xlab(expression(n[b]))

filter(res_tib, J == 100, n_a == 90, is.na(J_pred) | J_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_b, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation)) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste("n_a = 90, J = 100, ", J[p], " = 100"))) +
    xlab(expression(n[b]))

filter(res_tib, n_b == 10, n_a == 90, is.na(J_pred) | J_pred == 100) |>
    ggplot() +
    geom_line(aes(x = J, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation)) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle(expression(paste("n_a = 90, n_b = 10, ", J[p], " = 100")))

filter(res_tib, n_b == 10, n_a == 90, J == 100, !is.na(J_pred)) |>
    ggplot() +
    geom_line(aes(x = J_pred, y = time_mins,
                  col = gq_implementation, lty = loglik_implementation)) +
    facet_wrap(~ type) +
    scale_y_log10() +
    time_lab +
    ggtitle("n_a = 90, n_b = 10, J = 100")
