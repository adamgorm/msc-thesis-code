options(mc.cores = 4)
base_dir <- "../../results/regular/"
n_tib <- read.table(paste0(base_dir, "sim_args_full_fit.txt"), header = TRUE)
sim_dir <- paste0(base_dir, "/sim_data/")
time_dir <- paste0(base_dir, "/time_mins_full_fit/")
time_no_gq_dir <- paste0(time_dir, "/no_gq/")
time_smart_dir <- paste0(time_dir, "/smart/")
time_iter_dir <- paste0(time_dir, "/iter/")
time_smart_iter_dir <- paste0(time_dir, "/smart_iter/")
time_base_dir <- paste0(time_dir, "/base/")
lapply(c(sim_dir, time_dir, time_smart_dir, time_base_dir, time_no_gq_dir,
         time_iter_dir, time_smart_iter_dir),
       \(x) dir.create(x, showWarnings = FALSE))

sim_mod <- cmdstanr::cmdstan_model("../sim.stan")
fit_mods <- list(
    smart = cmdstanr::cmdstan_model("smart_smart_gq.stan"),
    base = cmdstanr::cmdstan_model("smart_base_gq.stan"),
    no_gq = cmdstanr::cmdstan_model("smart_no_gq.stan"),
    iter = cmdstanr::cmdstan_model("smart_iter_gq.stan"),
    smart_iter = cmdstanr::cmdstan_model("smart_smart_iter_gq.stan")
)

time_file_path <- function(n_obs, n_group, n_pred, implementation)
    sprintf("%s/%s/%d_%d_%d.txt", time_dir, implementation, n_obs, n_group, n_pred)

sim_data_file_path <- function(n_obs, n_group)
    sprintf("%s/%d_%d.rds", sim_dir, n_obs, n_group)

sim_data_set <- function(n_tib_row)
{
    n_obs <- n_tib_row$n_obs
    n_group <- n_tib_row$n_group

    save_location <- sim_data_file_path(n_obs, n_group)

    if (file.exists(save_location))
        return()
        
    sim_dat <- list(
        n_obs = n_obs,
        n_group = n_group,
        t = seq(-n_obs/5, n_obs/5, length.out = n_obs),
        magnitude_mu = 1,
        length_scale_mu = 1,
        magnitude_eta = 0.5,
        length_scale_eta = 0.5,
        sigma = 0.2
    )

    sim <- sim_mod$sample(data = sim_dat,
                          fixed_param = TRUE,
                          chains = 1,
                          iter_sampling = 1)
    
    draws <- posterior::as_draws_rvars(sim$draws())
    y <- posterior::draws_of(draws$y)

    saveRDS(list(t = sim_dat$t,
                 n_obs = sim_dat$n_obs,
                 n_group = sim_dat$n_group,
                 n_obs_total = sim_dat$n_obs * sim_dat$n_group,
                 y_obs = y[1, , ],
                 y_obs_vec = as.vector(y[1, , ]),
                 magnitude_mu = sim_dat$magnitude_mu,
                 length_scale_mu = sim_dat$length_scale_mu,
                 magnitude_eta = sim_dat$magnitude_eta,
                 length_scale_eta = sim_dat$length_scale_eta,
                 zero_n_obs_total = rep(0, sim_dat$n_obs * sim_dat$n_group),
                 one_mat_n_group = matrix(1, nrow = sim_dat$n_group, ncol = sim_dat$n_group),
                 sigma = sim_dat$sigma),
            save_location)
}

single_benchmark <- function(arg_tib)
{
    ## Load simulated data set
    sim_dat <- readRDS(sim_data_file_path(arg_tib$n_obs, arg_tib$n_group))
    stan_dat <- list(
        n_obs = sim_dat$n_obs,
        n_group = sim_dat$n_group,
        t_obs = sim_dat$t,
        y_obs = sim_dat$y_obs,
        n_pred = arg_tib$n_pred,
        t_pred = seq(-arg_tib$n_obs/5, arg_tib$n_obs/5, length.out = arg_tib$n_pred)
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
          time_file_path(arg_tib$n_obs, arg_tib$n_group, arg_tib$n_pred, arg_tib$implementation),
          append = TRUE)

    return(time_mins)
}

split_n_tib <- split(n_tib, seq(nrow(n_tib)))

## Simulate data set for each n_obs and n_group
parallel::mclapply(split_n_tib, sim_data_set)

complete_tib <- do.call(rbind,
                        lapply(c("no_gq", "smart_iter", "smart", "base"),
                               \(x) transform(n_tib, implementation = x)))

## Converts tib to a list of tibbles. The i'the element of the list is the
## i'th row of tib.
split_complete_tib <- split(complete_tib, seq(nrow(complete_tib)))

for (i in 1:length(split_complete_tib)) {
    cat(sprintf("\n\n>>> %d / %d : ", i, length(split_complete_tib)))
    cat(sprintf("n_obs=%d, n_group=%d, n_pred=%d, imp=%s",
                split_complete_tib[[i]]$n_obs,
                split_complete_tib[[i]]$n_group,
                split_complete_tib[[i]]$n_pred,
                split_complete_tib[[i]]$implementation))
    cat(sprintf("\n====== %f min\n", single_benchmark(split_complete_tib[[i]])))
}

## parallel::mclapply(split_complete_tib, single_benchmark)
