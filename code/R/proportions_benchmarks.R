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

res_tib_loglik |>
    pivot_wider(names_from = loglik_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    ggplot() +
    geom_histogram(aes(x = ord_mag))

res_tib_loglik |>
    pivot_wider(names_from = loglik_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_obs == 100) |>
    ggplot() +
    geom_line(aes(x = n_group, y = ord_mag))

res_tib_loglik |>
    pivot_wider(names_from = loglik_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_group == 100) |>
    ggplot() +
    geom_line(aes(x = n_obs, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / base)) |>
    ggplot() +
    geom_histogram(aes(x = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / base)) |>
    filter(n_group == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_obs, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / base)) |>
    filter(n_obs == 100, n_group == 100) |>
    ggplot() +
    geom_line(aes(x = n_pred, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / base)) |>
    filter(n_obs == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_group, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / smart)) |>
    filter(n_obs == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_group, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / smart)) |>
    filter(n_group == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_obs, y = ord_mag))

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart_iter / smart)) |>
    filter(n_group == 100, n_obs == 100) |>
    ggplot() +
    geom_line(aes(x = n_pred, y = ord_mag))

res_tib |>
    filter(type == "full fit") |>
    pivot_wider(names_from = c(loglik_implementation,
                               gq_implementation),
                names_sep = "__",
                values_from = time_mins) |>
    mutate(ord_mag_loglik = log10(smart__none / base__none)) |>
    filter(!is.na(base__none)) |>
    select(base__none, smart__none, ord_mag_loglik)

res_tib |>
    filter(type == "full fit") |>
    pivot_wider(names_from = c(loglik_implementation,
                               gq_implementation),
                names_sep = "__",
                values_from = time_mins) |>
    mutate(ord_mag_gq = log10(smart__smart_iter / smart__base)) |>
    filter(n_obs == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_group, y = ord_mag_gq))

res_tib |>
    filter(type == "full fit") |>
    pivot_wider(names_from = c(loglik_implementation,
                               gq_implementation),
                names_sep = "__",
                values_from = time_mins) |>
    mutate(ord_mag_gq = log10(smart__smart_iter / smart__base)) |>
    filter(n_group == 100, n_pred == 100) |>
    ggplot() +
    geom_line(aes(x = n_obs, y = ord_mag_gq))

res_tib |>
    filter(type == "full fit") |>
    pivot_wider(names_from = c(loglik_implementation,
                               gq_implementation),
                names_sep = "__",
                values_from = time_mins) |>
    mutate(ord_mag_gq = log10(smart__smart_iter / smart__base)) |>
    filter(n_group == 100, n_obs == 100) |>
    ggplot() +
    geom_line(aes(x = n_pred, y = ord_mag_gq))

### Irregular

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

res_tib_loglik |>
    pivot_wider(names_from = loglik_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_a == 90, n_b == 10)

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_a == 90, n_b == 10, J_pred == 50)

res_tib_full |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_a == 90, n_b == 10, J_pred == 50)

res_tib_gq |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_a == 90, n_b == 10, J == 50)

res_tib_full |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / base)) |>
    filter(n_a == 90, n_b == 10, J_pred == 50)

res_tib_full |>
    pivot_wider(names_from = gq_implementation,
                values_from = time_mins) |>
    mutate(ord_mag = log10(smart / no_gq)) |>
    filter(n_a == 90, n_b == 10, J_pred == 50)

