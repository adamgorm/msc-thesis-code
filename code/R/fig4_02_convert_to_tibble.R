library(tidyverse)

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

res_tib <- bind_rows(
    do.call(rbind,
            lapply(list(mu_ise, mu_ise_oracle, mu_ise_exclude, mu_ise_include,
                        f_ise, f_ise_oracle, f_ise_exclude, f_ise_include),
                   \(v) c(mean(v), quantile(v, c(0.05, 0.95))))) |>
    as_tibble() |>
    set_names("mean", "q5", "q95") |>
    mutate(method = rep(c("censoring", "oracle", "exclude", "include"), 2),
           func = rep(c("mu", "f"), each = 4)) |>
    mutate(type = "ISE"),
    do.call(rbind,
            lapply(list(mu_iae, mu_iae_oracle, mu_iae_exclude, mu_iae_include,
                        f_iae, f_iae_oracle, f_iae_exclude, f_iae_include),
                   \(v) c(mean(v), quantile(v, c(0.05, 0.95))))) |>
    as_tibble() |>
    set_names("mean", "q5", "q95") |>
    mutate(method = rep(c("censoring", "oracle", "exclude", "include"), 2),
           func = rep(c("mu", "f"), each = 4),
           type = "IAE"),
    do.call(rbind,
            lapply(list(-mu_ise_oracle + mu_ise,
                        -mu_ise_oracle + mu_ise_exclude,
                        -mu_ise_oracle + mu_ise_include,
                        -f_ise_oracle + f_ise,
                        -f_ise_oracle + f_ise_exclude,
                        -f_ise_oracle + f_ise_include),
                   \(v) c(mean(v), quantile(v, c(0.05, 0.95))))) |>
    as_tibble() |>
    set_names("mean", "q5", "q95") |>
    mutate(method = rep(c("censoring", "exclude", "include"), 2),
           func = rep(c("mu", "f"), each = 3),
           type = "ISE_contrast"),
    do.call(rbind,
            lapply(list(-mu_iae_oracle + mu_iae,
                        -mu_iae_oracle + mu_iae_exclude,
                        -mu_iae_oracle + mu_iae_include,
                        -f_iae_oracle + f_iae,
                        -f_iae_oracle + f_iae_exclude,
                        -f_iae_oracle + f_iae_include),
                   \(v) c(mean(v), quantile(v, c(0.05, 0.95))))) |>
    as_tibble() |>
    set_names("mean", "q5", "q95") |>
    mutate(method = rep(c("censoring", "exclude", "include"), 2),
           func = rep(c("mu", "f"), each = 3),
           type = "IAE_contrast")
)

write_csv(res_tib, "code/R/fig4_02_tibble.csv")
