library(tidyverse)
theme_set(theme_bw())

res_tib <- read_csv("code/R/fig4_02_tibble.csv")

res_tib |>
    mutate(method = factor(method, levels = c("oracle",
                                              "censoring",
                                              "include",
                                              "exclude") |> rev())) |>
    mutate(type = ifelse(type == "IAE_contrast", "IAE contrast", type)) |>
    mutate(type = ifelse(type == "ISE_contrast", "ISE contrast", type)) |>
    ggplot() +
    geom_pointrange(aes(x = mean, xmin = q5, xmax = q95, y = method, color = method)) +
    facet_grid(func ~ type) +
    bayesplot::vline_0(alpha = 0.2) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "bottom")

fig_dir <- "thesis/tex/figures/"
ggsave(paste0(fig_dir, "censoring_multi_forest.pdf"),
       width = 17,
       height = 10,
       units = "cm")
