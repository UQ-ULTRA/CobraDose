#' Write scripts for marginal prior tuning
#'
#' This function writes two R scripts for marginal prior tuning:
#' `marg_fun.R`, which sets up and compiles the NIMBLE model, and
#' `sim_marg.R`, which runs the simulation study over marginal distribution
#' scenarios.
#'
#' @param output_dir Directory where the scripts should be written.
#' @param n Total sample size.
#' @param group_sizes Numeric vector giving the sample size in each dose group.
#'   Must have length >= 2 (one entry per dose group). Each element should
#'   typically be >= 3 to allow stable estimation of group-specific parameters.
#'   The default `c(6, 4, 4, 4)` corresponds to the four-dose example from
#'   COBRA-DOSE paper.
#' @param data_gen_hyper_mu Hyper-mean parameters used to generate simulated data.
#' @param data_gen_hyper_prec Hyper-precision parameters used to generate simulated data.
#' @param prior_alpha_marg Dirichlet prior parameters for the two marginal model probabilities.
#' @param prior_hyper_mu Hyper-mean parameters used as priors in the fitted model.
#' @param prior_hyper_prec Hyper-precision parameters used as priors in the fitted model.
#' @param n_sim Number of simulation repetitions per true marginal scenario.
#' @param n_cores Number of parallel cores to use in `sim_marg.R`.
#' @param mcmc_niter Number of MCMC iterations for each simulated dataset.
#' @param mcmc_nburnin Number of burn-in iterations for each simulated dataset.
#' @param mcmc_nchains Number of MCMC chains.
#' @param output_prefix Prefix used for output CSV files.
#' @param seed Random seed base used for data generation.
#' @param overwrite Logical; if `FALSE`, existing files will not be overwritten.
#'
#' @details
#' The generated `marg_fun.R` script constructs a list of initial values for
#' MCMC, with one initial value list per chain. Initial values are designed
#' so that each of the four (`z1`, `z2`) starting combinations is assigned an
#' equal number of chains; this requires `mcmc_nchains` to be a positive
#' multiple of 4. The default `mcmc_nchains = 16` assigns 4 chains per
#' (`z1`, `z2`) combination, matching the configuration used in COBRA-DOSE paper.
#' Other valid values are 4, 8, 12, 20, etc. Smaller values run faster but
#' may give noisier estimates of the posterior model selection probabilities.
#'
#' @return Invisibly returns a list with the paths to the generated scripts.
#'
#' @examples
#' \dontrun{
#' # Generate the two scripts:
#' writeMarginalTuningScripts(
#'   output_dir = tempdir(),
#'   n = 18,
#'   group_sizes = c(6, 4, 4, 4),
#'   overwrite = TRUE
#' )
#'
#' # Then source the simulation script:
#' source(out$sim_marg, chdir = TRUE)
#'
#' # Inspect a result file:
#' csv <- read.csv(file.path(tempdir(), "model_select_marg1_1_marg2_1.csv"))
#' summarizeModelProbs(csv[, 2:5])
#' }
#'
#'@section When to re-run prior tuning:
#' The default `prior_alpha_marg = c(1.625, 1)` was tuned in COBRA-DOSE paper
#' for a specific setting: total sample size 18, four dose groups of sizes
#' c(6, 4, 4, 4), and the default `prior_hyper_mu` / `prior_hyper_prec`.
#' Under that setting, `c(1.625, 1)` is the smallest Dirichlet parameter
#' producing well-calibrated marginal model selection (see Table 2 of COBRA-DOSE paper).
#'
#' Users applying CobraDose to a substantially different setting --- different
#' total sample size, different group structure, or different prior
#' hyperparameters --- should re-run `writeMarginalTuningScripts()` with
#' their own setting and verify that the diagonal of the resulting model
#' probability matrix is approximately equal across rows and clearly higher
#' than the off-diagonal entries. If there is bias toward the lognormal distribution,
#' increase the first Dirichlet parameter until this holds. If there is bias toward the
#' gamma distribution, increase the second Dirichlet parameter until this holds.
#'
#' @export
writeMarginalTuningScripts <- function(
    output_dir = ".",
    n = 18,
    group_sizes = c(6, 4, 4, 4),

    data_gen_hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
    data_gen_hyper_prec = c(400, 400, 20, 800, 800, 20),

    prior_alpha_marg = c(1.625, 1),
    prior_hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
    prior_hyper_prec = c(5, 5, 2, 5, 5, 2),

    n_sim = 500,
    n_cores = max(1, parallel::detectCores() - 1),
    mcmc_niter = 3500,
    mcmc_nburnin = 1000,
    mcmc_nchains = 16,

    output_prefix = "model_select",
    seed = 1,
    overwrite = FALSE
) {

  # ---- input checks ----

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("`output_dir` must be a single character string.")
  }

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a single positive number.")
  }

  if (!is.numeric(group_sizes) || any(group_sizes <= 0)) {
    stop("`group_sizes` must be a numeric vector of positive values.")
  }

  if (length(group_sizes) < 2) {
    stop("`group_sizes` must have length >= 2 (a hierarchical model requires at least two dose groups).")
  }

  if (sum(group_sizes) != n) {
    stop("`n` must equal the sum of `group_sizes`.")
  }

  if (length(data_gen_hyper_mu) != 6 || length(data_gen_hyper_prec) != 6) {
    stop("`data_gen_hyper_mu` and `data_gen_hyper_prec` must both have length 6.")
  }

  if (length(prior_hyper_mu) != 6 || length(prior_hyper_prec) != 6) {
    stop("`prior_hyper_mu` and `prior_hyper_prec` must both have length 6.")
  }

  if (length(prior_alpha_marg) != 2 || any(prior_alpha_marg <= 0)) {
    stop("`prior_alpha_marg` must be a positive numeric vector of length 2.")
  }

  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim <= 0) {
    stop("`n_sim` must be a single positive number.")
  }

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    stop("`n_cores` must be a single positive number.")
  }

  if (!is.numeric(mcmc_niter) || length(mcmc_niter) != 1 || mcmc_niter <= 0) {
    stop("`mcmc_niter` must be a single positive number.")
  }

  if (!is.numeric(mcmc_nburnin) || length(mcmc_nburnin) != 1 || mcmc_nburnin < 0) {
    stop("`mcmc_nburnin` must be a single non-negative number.")
  }

  if (mcmc_nburnin >= mcmc_niter) {
    stop("`mcmc_nburnin` must be smaller than `mcmc_niter`.")
  }

  if (!is.numeric(mcmc_nchains) || length(mcmc_nchains) != 1 || mcmc_nchains <= 0) {
    stop("`mcmc_nchains` must be a single positive number.")
  }

  if (mcmc_nchains %% 4 != 0) {
    stop("`mcmc_nchains` must be a positive multiple of 4 (one or more chains per (z1, z2) initial value combination).")
  }

  if (!is.character(output_prefix) || length(output_prefix) != 1) {
    stop("`output_prefix` must be a single character string.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  marg_fun_file <- file.path(output_dir, "marg_fun.R")
  sim_marg_file <- file.path(output_dir, "sim_marg.R")

  if (!overwrite && (file.exists(marg_fun_file) || file.exists(sim_marg_file))) {
    stop("Output files already exist. Use `overwrite = TRUE` to replace them.")
  }

  # ---- helper for writing R vectors into generated scripts ----

  vec_to_r <- function(x) {
    paste0("c(", paste(x, collapse = ", "), ")")
  }

  G <- length(group_sizes)

  # ---- marg_fun.R ----

  marg_fun_code <- c(
    "library(nimble)",
    "library(CobraDose)",
    "",
    "# Auto-generated by CobraDose::writeMarginalTuningScripts()",
    "# This script sets up the NIMBLE model for marginal prior tuning.",
    "",
    paste0("n <- ", as.integer(n)),
    paste0("group_sizes <- ", vec_to_r(group_sizes)),
    paste0("mcmc_nchains <- ", as.integer(mcmc_nchains)),
    paste0("G <- ", as.integer(G)),
    paste0("prior_alpha_marg <- ", vec_to_r(prior_alpha_marg)),
    paste0("prior_hyper_mu <- ", vec_to_r(prior_hyper_mu)),
    paste0("prior_hyper_prec <- ", vec_to_r(prior_hyper_prec)),
    "",
    "# Register the custom bivariate copula distribution.",
    "CobraDose::registerBivarCopulaDistribution()",
    "",
    "# NIMBLE model for marginal prior tuning.",
    "# The copula is fixed to independence because this script focuses on marginal priors.",
    "modelCode <- nimbleCode({",
    "",
    "  z1 ~ dcat(pi_marg1[1:2])",
    "  z2 ~ dcat(pi_marg2[1:2])",
    "",
    "  pi_marg1[1:2] ~ ddirch(alpha_marg[1:2])",
    "  pi_marg2[1:2] ~ ddirch(alpha_marg[1:2])",
    "",
    "  mu_m1 ~ dnorm(hyper_mu[1], hyper_prec[1])",
    "  mu_m2 ~ dnorm(hyper_mu[1], hyper_prec[1])",
    "  mu_c1 ~ dnorm(hyper_mu[2], hyper_prec[2])",
    "  mu_c2 ~ dnorm(hyper_mu[2], hyper_prec[2])",
    "",
    "  sigma_m1 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)",
    "  sigma_m2 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)",
    "  sigma_c1 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)",
    "  sigma_c2 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)",
    "",
    "  for(g in 1:G){",
    "",
    "    lower_m1[g] <- -mu_m1 / sigma_m1",
    "    lower_m2[g] <- -mu_m2 / sigma_m2",
    "    lower_c1[g] <- -mu_c1 / sigma_c1",
    "    lower_c2[g] <- -mu_c2 / sigma_c2",
    "",
    "    z_m1[g] ~ T(dnorm(0, 1), lower_m1[g], Inf)",
    "    z_m2[g] ~ T(dnorm(0, 1), lower_m2[g], Inf)",
    "    z_c1[g] ~ T(dnorm(0, 1), lower_c1[g], Inf)",
    "    z_c2[g] ~ T(dnorm(0, 1), lower_c2[g], Inf)",
    "",
    "    m1[g] <- mu_m1 + sigma_m1 * z_m1[g]",
    "    m2[g] <- mu_m2 + sigma_m2 * z_m2[g]",
    "    c1[g] <- mu_c1 + sigma_c1 * z_c1[g]",
    "    c2[g] <- mu_c2 + sigma_c2 * z_c2[g]",
    "  }",
    "",
    "  for(i in 1:N){",
    "    y[i, 1:2] ~ dBivarCopula(",
    "      marg1 = z1,",
    "      m1 = m1[group[i]], c1 = c1[group[i]],",
    "      marg2 = z2,",
    "      m2 = m2[group[i]], c2 = c2[group[i]],",
    "      cop = 1,",
    "      tau = 0.5",
    "    )",
    "  }",
    "})",
    "",
    "# Initial values.",
    "# G is derived from group_sizes, so this works for any number of groups.",
    "inits1 <- list(",
    "  mu_m1 = 1, mu_c1 = 1, mu_m2 = 1, mu_c2 = 1,",
    "  sigma_m1 = 1, sigma_c1 = 1, sigma_m2 = 1, sigma_c2 = 1,",
    "  z_m1 = rep(0, G),",
    "  z_m2 = rep(0, G),",
    "  z_c1 = rep(0, G),",
    "  z_c2 = rep(0, G),",
    "  pi_marg1 = rep(0.5, 2),",
    "  pi_marg2 = rep(0.5, 2),",
    "  z1 = 1,",
    "  z2 = 1",
    ")",
    "",
    "# Initial values for all marginal-distribution combinations.",
    "# Each of the four (z1, z2) combinations is repeated 4 times,",
    "# producing 16 inits to match the default mcmc_nchains = 16.",
    "zz_grid <- expand.grid(z1 = 1:2, z2 = 1:2)",
    "n_chains_per_combo <- mcmc_nchains %/% nrow(zz_grid)",
    "inits <- vector(\"list\", mcmc_nchains)",
    "for (j in seq_len(nrow(zz_grid))) {",
    "  for (m in seq_len(n_chains_per_combo)) {",
    "    idx <- (j - 1) * n_chains_per_combo + m",
    "    inits[[idx]] <- inits1",
    "    inits[[idx]]$z1 <- zz_grid$z1[j]",
    "    inits[[idx]]$z2 <- zz_grid$z2[j]",
    "  }",
    "}",
    "",
    "# Build a dummy dataset to initialise and compile the model.",
    "# The simulation script will replace y for each simulation repetition.",
    "data1 <- CobraDose::dataGen(",
    "  n = n,",
    "  gs = group_sizes,",
    "  marg1 = 1,",
    "  marg2 = 1,",
    "  cop = 1,",
    "  hyper_mu = prior_hyper_mu,",
    "  hyper_prec = prior_hyper_prec,",
    "  seed = 1",
    ")",
    "",
    "model <- nimbleModel(",
    "  modelCode,",
    "  data = data1[[2]],",
    "  inits = inits[[1]],",
    "  constants = list(",
    "    N = data1[[1]]$N,",
    "    G = data1[[1]]$G,",
    "    group = data1[[1]]$group,",
    "    alpha_marg = prior_alpha_marg,",
    "    hyper_mu = prior_hyper_mu,",
    "    hyper_prec = prior_hyper_prec",
    "  ),",
    "  check = FALSE",
    ")",
    "",
    "conf <- configureMCMC(model)",
    "",
    "conf$addMonitors(c(\"z1\", \"z2\", \"m1\", \"c1\", \"m2\", \"c2\"))",
    "",
    "conf$removeSamplers(c(",
    "  \"z_m1\", \"z_c1\", \"z_m2\", \"z_c2\",",
    "  \"sigma_m1\", \"sigma_c1\", \"sigma_m2\", \"sigma_c2\"",
    "))",
    "",
    "for(i in 1:data1[[1]]$G) {",
    "  conf$addSampler(target = paste0(\"z_m1[\", i, \"]\"), type = \"slice\")",
    "  conf$addSampler(target = paste0(\"z_c1[\", i, \"]\"), type = \"slice\")",
    "  conf$addSampler(target = paste0(\"z_m2[\", i, \"]\"), type = \"slice\")",
    "  conf$addSampler(target = paste0(\"z_c2[\", i, \"]\"), type = \"slice\")",
    "}",
    "",
    "conf$addSampler(target = \"sigma_m1\", type = \"slice\")",
    "conf$addSampler(target = \"sigma_c1\", type = \"slice\")",
    "conf$addSampler(target = \"sigma_m2\", type = \"slice\")",
    "conf$addSampler(target = \"sigma_c2\", type = \"slice\")",
    "",
    "Rmcmc <- buildMCMC(conf)",
    "Cmodel <- compileNimble(model)",
    "Cmcmc <- compileNimble(Rmcmc, project = Cmodel)"
  )

  # ---- sim_marg.R ----

  sim_marg_code <- c(
    "library(nimble)",
    "library(foreach)",
    "library(doSNOW)",
    #"library(parallel)",
    "library(CobraDose)",
    "",
    "# Auto-generated by CobraDose::writeMarginalTuningScripts()",
    "# This script runs marginal prior tuning simulations.",
    "#",
    "# Note on runtime:",
    "# Each parallel worker compiles its own NIMBLE model on startup",
    "# (via source(marg_fun.R)). This means the initial setup time",
    "# scales with n_cores rather than n_sim, so you will see the",
    "# progress bar sit at 0% for some time before simulations begin.",
    "# Larger n_sim amortises this fixed cost; larger n_cores does not.",
    "",
    paste0("output_dir <- ", deparse(normalizePath(output_dir, winslash = "/", mustWork = TRUE))),
    paste0("n <- ", as.integer(n)),
    paste0("group_sizes <- ", vec_to_r(group_sizes)),
    paste0("data_gen_hyper_mu <- ", vec_to_r(data_gen_hyper_mu)),
    paste0("data_gen_hyper_prec <- ", vec_to_r(data_gen_hyper_prec)),
    paste0("n_sim <- ", as.integer(n_sim)),
    paste0("n_cores <- ", as.integer(n_cores)),
    paste0("mcmc_niter <- ", as.integer(mcmc_niter)),
    paste0("mcmc_nburnin <- ", as.integer(mcmc_nburnin)),
    paste0("mcmc_nchains <- ", as.integer(mcmc_nchains)),
    paste0("output_prefix <- \"", output_prefix, "\""),
    paste0("seed_base <- ", as.integer(seed)),
    "",
    "true_marginal_grid <- expand.grid(marg1 = 1:2, marg2 = 1:2)",
    "",
    "# Set up parallel backend.",
    "cl <- parallel::makeCluster(n_cores, type = \"SOCK\")",
    "doSNOW::registerDoSNOW(cl)",
    "",
    "pb <- utils::txtProgressBar(max = n_sim * nrow(true_marginal_grid), style = 3)",
    "progress <- function(n) utils::setTxtProgressBar(pb, n)",
    "opts <- list(progress = progress)",
    "",
    "# Compile NIMBLE model once on each worker.",
    "# Use an absolute path to reliably find marg_fun.R.",
    "marg_fun_path <- file.path(output_dir, \"marg_fun.R\")",
    "",
    "parallel::clusterExport(",
    "  cl,",
    "  varlist = \"marg_fun_path\",",
    "  envir = environment()",
    ")",
    "",
    "# Compile NIMBLE on workers serially to avoid Cygwin/RTools fork conflicts",
    "for (worker_idx in seq_along(cl)) {",
    "  print(paste0('Initializing worker ', worker_idx, ' of ', max(seq_along(cl))))",
    "  parallel::clusterCall(cl[worker_idx], function(path) {",
    "    source(path)",
    "  }, marg_fun_path)",
    "}",
    "",
    "scenario_counter <- 0",
    "",
    "for (scenario_index in seq_len(nrow(true_marginal_grid))) {",
    "",
    "  ii <- true_marginal_grid$marg1[scenario_index]",
    "  jj <- true_marginal_grid$marg2[scenario_index]",
    "",
    "  sim_res <- foreach::foreach(",
    "    k = seq_len(n_sim),",
    "    .packages = c(\"nimble\", \"CobraDose\"),",
    "    .combine = rbind,",
    "    .errorhandling = \"remove\",",
    "    .options.snow = opts,",
    "    .noexport = c(",
    "      \"model\", \"Cmodel\", \"Cmcmc\", \"inits\",",
    "      \"modelCode\", \"conf\", \"Rmcmc\", \"data1\"",
    "    )",
    "  ) %dopar% {",
    "",
    "    sim_seed <- seed_base + k + (jj - 1) * n_sim + 2 * n_sim * (ii - 1)",
    "",
    "    data_i <- CobraDose::dataGen(",
    "      n = n,",
    "      gs = group_sizes,",
    "      marg1 = ii,",
    "      marg2 = jj,",
    "      cop = 1,",
    "      hyper_mu = data_gen_hyper_mu,",
    "      hyper_prec = data_gen_hyper_prec,",
    "      seed = sim_seed",
    "    )",
    "",
    "    Cmodel$setData(list(y = data_i[[2]]$y))",
    "",
    "    Cmodel$initializeInfo()",
    "    model$calculate()",
    "    Cmcmc$run(1, reset = TRUE)",
    "",
    "    samples <- nimble::runMCMC(",
    "      Cmcmc,",
    "      niter = mcmc_niter,",
    "      nburnin = mcmc_nburnin,",
    "      nchains = mcmc_nchains,",
    "      inits = inits",
    "    )",
    "",
    "    samples_all <- do.call(rbind, samples)",
    "",
    "    c(",
    "      simulation = k,",
    "      gamma_gamma = mean(samples_all[, \"z1\"] == 1 & samples_all[, \"z2\"] == 1),",
    "      gamma_lognormal = mean(samples_all[, \"z1\"] == 1 & samples_all[, \"z2\"] == 2),",
    "      lognormal_gamma = mean(samples_all[, \"z1\"] == 2 & samples_all[, \"z2\"] == 1),",
    "      lognormal_lognormal = mean(samples_all[, \"z1\"] == 2 & samples_all[, \"z2\"] == 2)",
    "    )",
    "  }",
    "",
    "  out_file <- file.path(output_dir, paste0(output_prefix, \"_marg1_\", ii, \"_marg2_\", jj, \".csv\"))",
    "  utils::write.csv(sim_res, out_file, row.names = FALSE)",
    "",
    "  print(paste(\"Saved\", out_file))",
    "  print(CobraDose::summarizeModelProbs(sim_res[, 2:5]))",
    "}",
    "",
    "close(pb)",
    "parallel::stopCluster(cl)"
  )

  writeLines(marg_fun_code, marg_fun_file)
  writeLines(sim_marg_code, sim_marg_file)

  invisible(list(
    marg_fun = marg_fun_file,
    sim_marg = sim_marg_file
  ))
}


