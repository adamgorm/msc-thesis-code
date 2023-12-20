# Gaussian process regression for longitudinal data

This repository contains code for my MSc thesis in statistics with the title
*Gaussian process regression for longitudinal data: computational improvements
and value-censored data*.

The main code file with all the important implementations is
`/code/rao/code/functions.stan`. The rest of the files are mainly calling
functions from this file in order to benchmark, plot and test.

## Overview of directories

- `/code/rao/code`: Contains code that was kept in sync with code run on the
  `rao` server.
  - `/code/rao/code/functions.stan`: Main code file with all important functions
    implemented for the thesis.
  - `/code/rao/code/regular`:
    - R scripts for benchmarking the different implementations. These must be
      run from the `/code/rao/code/regular` directory.
    - Stan scripts for fitting the model in the regular case with different
      combinations of log-likelihood implementation (first part of filename) and
      generated quantities implementation (second part of filename).
  - `/code/rao/code/irregular`: Same as `/code/rao/code/regular` but in the
    irregular case.
- `/code/rao/results`: Results from running the benchmark scripts are saved
  here. The `sim_args` text files are read by the benchmark scripts to get the
  settings to simulate at.
- `/code/R`: R scripts, for instance for creating figures. These must be run
  from the base directory (i.e., outside of the `code` directory).
  - `/code/R/tests.R`: Tests that run all implementations on the same simulated
    data, to make sure that they give the same.
- `/code/stan`: Stan scripts for fitting and simulating from GPs. Called by
  scripts in `/code/R`.
