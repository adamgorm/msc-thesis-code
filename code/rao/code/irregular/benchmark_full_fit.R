options(mc.cores = 4)
source("sim_irregular.R")
sim_mod <- cmdstanr::cmdstan_model("../sim.stan")
base_dir <- "../../results/irregular/"
sim_dir <- paste0(base_dir, "/sim_data/")
time_dir <- paste0(base_dir, "/time_mins_full_fit/")
time_smart_dir <- paste0(time_dir, "/smart/")
time_base_dir <- paste0(time_dir, "/base/")
time_no_gq_dir <- paste0(time_dir, "/no_gq/")
lapply(c(sim_dir, time_dir, time_smart_dir, time_base_dir, time_no_gq_dir),
       \(x) dir.create(x, showWarnings = FALSE))

fit_mods <- list(
    smart = cmdstanr::cmdstan_model("smart_smart_gq.stan"),
    base = cmdstanr::cmdstan_model("smart_base_gq.stan"),
    no_gq = cmdstanr::cmdstan_model("smart_no_gq.stan")
)

n_tib <- read.table(paste0(base_dir, "sim_args_full_fit.txt"), header = TRUE)

time_file_path <- function(J, n_a, n_b, J_pred, implementation)
    sprintf("%s/%s/%d_%d_%d_%d.txt", time_dir, implementation, J, n_a, n_b, J_pred)

sim_data_file_path <- function(J, n_a, n_b)
    sprintf("%s/%d_%d_%d.rds", sim_dir, J, n_a, n_b)

sim_data_set <- function(n_tib_row)
{
    J <- n_tib_row$J
    n_a <- n_tib_row$n_a
    n_b <- n_tib_row$n_b

    save_location <- sim_data_file_path(J, n_a, n_b)

    if (file.exists(save_location))
        return()
    
    J_b <- rep(J, n_b)
    t_a <- seq(-J/5, J/5, length.out = J)
    t_b <- list()
    for (i in 1:n_b)
        t_b <- c(t_b, list(t_a))
    sim_dat <- list(
        n_a = n_a,
        n_b = n_b,
        J_a = J,
        J_b = J_b,
        t_b_is = c(0, cumsum(J_b)),
        t_a = t_a,
        t_b = t_b,
        magnitude_mu = 1,
        length_scale_mu = 1,
        magnitude_eta = 0.5,
        length_scale_eta = 0.5,
        sigma = 0.2
    )

    sim <- sim_irregular(sim_dat$n_a,
                         sim_dat$J_a,
                         sim_dat$t_a,
                         sim_dat$n_b,
                         sim_dat$J_b,
                         sim_dat$t_b,
                         sim_dat$magnitude_mu,
                         sim_dat$length_scale_mu,
                         sim_dat$magnitude_eta,
                         sim_dat$length_scale_eta,
                         sim_dat$sigma,
                         sim_mod)
    
    saveRDS(sim[[1]], save_location)
}

split_n_tib <- split(n_tib, seq(nrow(n_tib)))

## Simulate data set for each n_obs and n_group
parallel::mclapply(split_n_tib, sim_data_set)

complete_tib <- do.call(rbind,
                        lapply(c("no_gq", "smart", "base"),
                               \(x) transform(n_tib, implementation = x)))

single_benchmark <- function(arg_tib)
{
    J_pred <- arg_tib$J_pred
    t_pred <- seq(-J_pred/5, J_pred/5, length.out = J_pred)

    ## Load simulated data set
    sim_dat <- readRDS(sim_data_file_path(arg_tib$J, arg_tib$n_a, arg_tib$n_b))

    stan_dat <- my_list <- list(
                    y_a = sim_dat$y_a,
                    y_a_vec = sim_dat$y_a_vec,
                    y_b_vec = sim_dat$y_b_vec,
                    y_ab_vec = sim_dat$y_ab_vec,
                    t_a = sim_dat$t_a,
                    t_b = sim_dat$t_b,
                    t_pred = t_pred,
                    n_a = sim_dat$n_a,
                    n_b = sim_dat$n_b,
                    J_pred = J_pred,
                    J_a = sim_dat$J_a,
                    J_b = sim_dat$J_b,
                    t_b_is = sim_dat$t_b_is,
                    magnitude_mu = sim_dat$magnitude_mu,
                    length_scale_mu = sim_dat$length_scale_mu,
                    magnitude_eta = sim_dat$magnitude_eta,
                    length_scale_eta = sim_dat$length_scale_eta,
                    sigma = sim_dat$sigma,
                    one_mat_n_a = sim_dat$one_mat_n_a
                )
    
    t1 <- Sys.time()
    fit_mods[[arg_tib$implementation]]$sample(data = stan_dat,
                                              refresh = 250,
                                              chains = 1,
                                              iter_warmup = 1000,
                                              iter_sampling = 1000,
                                              show_exceptions = FALSE)
    t2 <- Sys.time()
    
    time_mins <- as.numeric(difftime(t2, t1, units = "mins"))
    write(as.character(time_mins),
          time_file_path(arg_tib$J, arg_tib$n_a, arg_tib$n_b, J_pred, arg_tib$implementation),
          append = TRUE)
    
    return(time_mins)
}

## Converts tib to a list of tibbles. The i'the element of the list is the
## i'th row of tib.
split_complete_tib <- split(complete_tib, seq(nrow(complete_tib)))

for (i in 1:length(split_complete_tib)) {
    cat(sprintf("\n%d / %d : ", i, length(split_complete_tib)))
    cat(sprintf("J=%d, n_a=%d, n_b=%d, J_pred=%d, imp=%s\n",
                split_complete_tib[[i]]$J,
                split_complete_tib[[i]]$n_a,
                split_complete_tib[[i]]$n_b,
                split_complete_tib[[i]]$J_pred,                
                split_complete_tib[[i]]$implementation))
    cat(sprintf("\n===== %f min\n", single_benchmark(split_complete_tib[[i]])))
}

## parallel::mclapply(split_complete_tib, single_benchmark)
