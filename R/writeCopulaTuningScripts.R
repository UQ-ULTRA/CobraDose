#' Write scripts for copula prior tuning
#'
#' This function writes two R scripts for copula prior tuning:
#' `cop_fun.R`, which sets up and compiles the NIMBLE model, and
#' `sim_cop.R`, which runs the simulation study over copula scenarios.
#'
#' @param output_dir Directory where the scripts should be written.
#' @param n Total sample size.
#' @param group_sizes Numeric vector giving the sample size in each dose group.
#'   Must have length >= 2 (one entry per dose group). Each element should
#'   typically be >= 3 to allow stable estimation of group-specific parameters.
#'   The default `c(6, 4, 4, 4)` corresponds to the four-dose example from
#'   Cobra-Dose paper.
#' @param data_gen_hyper_mu Hyper-mean parameters used to generate simulated data.
#' @param data_gen_hyper_prec Hyper-precision parameters used to generate simulated data.
#'   These are typically more concentrated (higher precision) than the analysis
#'   priors, ensuring the simulated data has appreciable dependence to allow
#'   meaningful discrimination among copula families.
#' @param prior_alpha_cop Dirichlet prior parameters for the four copula model
#'   probabilities (independence, Clayton, Gaussian, Gumbel). Must be a positive
#'   numeric vector of length 4.
#' @param prior_hyper_mu Hyper-mean parameters used as priors in the fitted model.
#' @param prior_hyper_prec Hyper-precision parameters used as priors in the fitted model.
#' @param n_sim Number of simulation repetitions per true copula scenario.
#' @param n_cores Number of parallel cores to use in `sim_cop.R`.
#' @param mcmc_niter Number of MCMC iterations for each simulated dataset.
#' @param mcmc_nburnin Number of burn-in iterations for each simulated dataset.
#' @param mcmc_nchains Number of MCMC chains. Must be a positive multiple of 4.
#' @param output_prefix Prefix used for output CSV files.
#' @param seed Random seed base used for data generation. The default
#'   `seed = 8001` is offset from the marginal tuning default (`seed = 1`)
#'   so that simulated datasets in marginal and copula tuning workflows
#'   are disjoint.
#' @param overwrite Logical; if `FALSE`, existing files will not be overwritten.
#'
#' @details
#' The generated `cop_fun.R` script fits a Bayesian model in which the
#' marginal distributions are fixed to gamma for both biomarkers, while
#' the copula family is treated as a discrete unknown selected via the
#' latent variable `zc` (1 = independence, 2 = Clayton, 3 = Gaussian,
#' 4 = Gumbel). This follows the second simulation study in Cobra-Dose
#' paper, in which copula priors are tuned under the assumption of
#' gamma marginals.
#'
#' Initial values cover the four `zc` starting values (one per candidate
#' copula); this requires `mcmc_nchains` to be a positive multiple of 4.
#' The default `mcmc_nchains = 16` assigns 4 chains per copula starting
#' value, matching the configuration used in Cobra-Dose paper. Other
#' valid values are 4, 8, 12, 20, etc. Smaller values run faster but may
#' give noisier estimates of the posterior model selection probabilities.
#'
#' @return Invisibly returns a list with the paths to the generated scripts.
#'
#' @examples
#' \dontrun{
#' # Generate the two scripts:
#' out <- writeCopulaTuningScripts(
#'   output_dir = tempdir(),
#'   n = 18,
#'   group_sizes = c(6, 4, 4, 4),
#'   overwrite = TRUE
#' )
#'
#' # Then source the simulation script:
#' source(out$sim_cop, chdir = TRUE)
#'
#' # Inspect a result file:
#' csv <- read.csv(file.path(tempdir(), "model_select_cop_1.csv"))
#' summarizeModelProbs(csv[, 2:5])
#' }
#'
#' @section When to re-run prior tuning:
#' The default `prior_alpha_cop = c(1, 1, 1, 1)` was used in Cobra-Dose
#' paper for a specific setting: total sample size 18, four dose groups
#' of sizes c(6, 4, 4, 4), gamma marginals, and the default
#' `prior_hyper_mu` / `prior_hyper_prec`. Under that setting, a uniform
#' Dirichlet prior produced well-calibrated copula model selection (see
#' Table 3 of Cobra-Dose paper).
#'
#' Users applying CobraDose to a substantially different setting ---
#' different total sample size, different group structure, or different
#' prior hyperparameters --- should re-run `writeCopulaTuningScripts()`
#' with their own setting and verify that the diagonal of the resulting
#' model probability matrix is approximately equal across rows and
#' clearly higher than the off-diagonal entries. If a particular copula
#' family is systematically under-selected, increase the corresponding
#' entry of `alpha_cop`.
#'
#' If gamma is not an appropriate marginal family for the user's data,
#' the calibrated `alpha_cop` may not transfer directly; run
#' `writeMarginalTuningScripts()` first to confirm.
#'
#' @export
writeCopulaTuningScripts <- function(
    output_dir = ".",
    n = 18,
    group_sizes = c(6, 4, 4, 4),

    data_gen_hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
    data_gen_hyper_prec = c(400, 400, 20, 800, 800, 20),

    prior_alpha_cop = c(1, 1, 1, 1),
    prior_hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
    prior_hyper_prec = c(5, 5, 2, 5, 5, 2),

    n_sim = 500,
    n_cores = max(1, parallel::detectCores() - 1),
    mcmc_niter = 3500,
    mcmc_nburnin = 1000,
    mcmc_nchains = 16,

    output_prefix = "model_select",
    seed = 8001,
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

  if (length(prior_alpha_cop) != 4 || any(prior_alpha_cop <= 0)) {
    stop("`prior_alpha_cop` must be a positive numeric vector of length 4 (one entry per candidate copula family).")
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
    stop("`mcmc_nchains` must be a positive multiple of 4 (one or more chains per copula starting value).")
  }

  if (!is.character(output_prefix) || length(output_prefix) != 1) {
    stop("`output_prefix` must be a single character string.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cop_fun_file <- file.path(output_dir, "cop_fun.R")
  sim_cop_file <- file.path(output_dir, "sim_cop.R")

  if (!overwrite && (file.exists(cop_fun_file) || file.exists(sim_cop_file))) {
    stop("Output files already exist. Use `overwrite = TRUE` to replace them.")
  }

  # ---- helper for writing R vectors into generated scripts ----

  vec_to_r <- function(x) {
    paste0("c(", paste(x, collapse = ", "), ")")
  }

  G <- length(group_sizes)

  # ---- cop_fun.R ----

  cop_fun_code <- c(
    "library(nimble)",
    "library(CobraDose)",
    "",
    "# Auto-generated by CobraDose::writeCopulaTuningScripts()",
    "# This script sets up the NIMBLE model for copula prior tuning.",
    "# Marginal distributions are fixed to gamma for both biomarkers; the",
    "# copula family (zc in {1, 2, 3, 4}) is the unknown to be selected.",
    "",
    paste0("n <- ", as.integer(n)),
    paste0("group_sizes <- ", vec_to_r(group_sizes)),
    paste0("mcmc_nchains <- ", as.integer(mcmc_nchains)),
    paste0("G <- ", as.integer(G)),
    paste0("prior_alpha_cop <- ", vec_to_r(prior_alpha_cop)),
    paste0("prior_hyper_mu <- ", vec_to_r(prior_hyper_mu)),
    paste0("prior_hyper_prec <- ", vec_to_r(prior_hyper_prec)),
    "",
    "# Register the custom bivariate copula distribution.",
    "CobraDose::registerBivarCopulaDistribution()",
    "",
    "# NIMBLE model for copula prior tuning.",
    "# Marginals are fixed to gamma (marg1 = marg2 = 1); only the copula",
    "# family zc is treated as a discrete unknown.",
    "modelCode <- nimbleCode({",
    "",
    "  zc ~ dcat(pi_cop[1:4])",
    "",
    "  pi_cop[1:4] ~ ddirch(alpha_cop[1:4])",
    "",
    "  mu_m1 ~ dnorm(hyper_mu[1], hyper_prec[1])",
    "  mu_m2 ~ dnorm(hyper_mu[1], hyper_prec[1])",
    "  mu_c1 ~ dnorm(hyper_mu[2], hyper_prec[2])",
    "  mu_c2 ~ dnorm(hyper_mu[2], hyper_prec[2])",
    "  mu_logittau ~ dnorm(hyper_mu[3], hyper_prec[3])",
    "",
    "  sigma_m1 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)",
    "  sigma_m2 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)",
    "  sigma_c1 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)",
    "  sigma_c2 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)",
    "  sigma_logittau ~ T(dnorm(hyper_mu[6], hyper_prec[6]), 1e-6, Inf)",
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
    "    z_tau[g] ~ dnorm(0, 1)",
    "",
    "    m1[g] <- mu_m1 + sigma_m1 * z_m1[g]",
    "    m2[g] <- mu_m2 + sigma_m2 * z_m2[g]",
    "    c1[g] <- mu_c1 + sigma_c1 * z_c1[g]",
    "    c2[g] <- mu_c2 + sigma_c2 * z_c2[g]",
    "",
    "    logit_tau[g] <- mu_logittau + sigma_logittau * z_tau[g]",
    "    tau[g] <- ilogit(logit_tau[g])",
    "  }",
    "",
    "  for(i in 1:N){",
    "    y[i, 1:2] ~ dBivarCopula(",
    "      marg1 = 1,",
    "      m1 = m1[group[i]], c1 = c1[group[i]],",
    "      marg2 = 1,",
    "      m2 = m2[group[i]], c2 = c2[group[i]],",
    "      cop = zc,",
    "      tau = tau[group[i]]",
    "    )",
    "  }",
    "})",
    "",
    "# Initial values.",
    "# G is derived from group_sizes, so this works for any number of groups.",
    "inits1 <- list(",
    "  mu_m1 = 1, mu_c1 = 1, mu_m2 = 1, mu_c2 = 1,",
    "  sigma_m1 = 1, sigma_c1 = 1, sigma_m2 = 1, sigma_c2 = 1,",
    "  mu_logittau = 0, sigma_logittau = 1,",
    "  z_m1 = rep(0, G),",
    "  z_m2 = rep(0, G),",
    "  z_c1 = rep(0, G),",
    "  z_c2 = rep(0, G),",
    "  z_tau = rep(0, G),",
    "  pi_cop = rep(0.25, 4),",
    "  zc = 1",
    ")",
    "",
    "# Initial values for all copula starting values.",
    "# Each of the four zc starting values is repeated multiple times,",
    "# producing mcmc_nchains inits in total.",
    "zc_grid <- data.frame(zc = 1:4)",
    "n_chains_per_combo <- mcmc_nchains %/% nrow(zc_grid)",
    "inits <- vector(\"list\", mcmc_nchains)",
    "for (j in seq_len(nrow(zc_grid))) {",
    "  for (m in seq_len(n_chains_per_combo)) {",
    "    idx <- (j - 1) * n_chains_per_combo + m",
    "    inits[[idx]] <- inits1",
    "    inits[[idx]]$zc <- zc_grid$zc[j]",
    "  }",
    "}",
    "",
    "# Build a dummy dataset to initialise and compile the model.",
    "# Marginals are fixed to gamma (marg1 = marg2 = 1); copula is independence",
    "# (cop = 1) for the dummy data. The simulation script will replace y for",
    "# each simulation repetition with data generated under different copulas.",
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
    "    alpha_cop = prior_alpha_cop,",
    "    hyper_mu = prior_hyper_mu,",
    "    hyper_prec = prior_hyper_prec",
    "  ),",
    "  check = FALSE",
    ")",
    "",
    "conf <- configureMCMC(model)",
    "",
    "conf$addMonitors(c(\"zc\", \"m1\", \"c1\", \"m2\", \"c2\"))",
    "",
    "conf$removeSamplers(c(",
    "  \"z_m1\", \"z_c1\", \"z_m2\", \"z_c2\",",
    "  \"sigma_m1\", \"sigma_c1\", \"sigma_m2\", \"sigma_c2\",",
    "  \"sigma_logittau\"",
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
    "conf$addSampler(target = \"sigma_logittau\", type = \"slice\")",
    "",
    "Rmcmc <- buildMCMC(conf)",
    "Cmodel <- compileNimble(model)",
    "Cmcmc <- compileNimble(Rmcmc, project = Cmodel)"
  )

  # ---- sim_cop.R ----

  sim_cop_code <- c(
    "library(nimble)",
    "library(foreach)",
    "library(doSNOW)",
    "library(CobraDose)",
    "",
    "# Auto-generated by CobraDose::writeCopulaTuningScripts()",
    "# This script runs copula prior tuning simulations.",
    "#",
    "# Note on runtime:",
    "# Each parallel worker compiles its own NIMBLE model on startup",
    "# (via source(cop_fun.R)). This means the initial setup time",
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
    "# Four copula scenarios: 1 = independence, 2 = Clayton, 3 = Gaussian, 4 = Gumbel.",
    "true_cop_grid <- data.frame(cop = 1:4)",
    "",
    "# Set up parallel backend.",
    "cl <- parallel::makeCluster(n_cores, type = \"SOCK\")",
    "doSNOW::registerDoSNOW(cl)",
    "",
    "pb <- utils::txtProgressBar(max = n_sim * nrow(true_cop_grid), style = 3)",
    "progress <- function(n) utils::setTxtProgressBar(pb, n)",
    "opts <- list(progress = progress)",
    "",
    "# Compile NIMBLE model once on each worker.",
    "# Use an absolute path to reliably find cop_fun.R.",
    "cop_fun_path <- file.path(output_dir, \"cop_fun.R\")",
    "",
    "parallel::clusterExport(",
    "  cl,",
    "  varlist = \"cop_fun_path\",",
    "  envir = environment()",
    ")",
    "",
    "# Compile NIMBLE on workers serially to avoid Cygwin/RTools fork conflicts.",
    "for (worker_idx in seq_along(cl)) {",
    "  parallel::clusterCall(cl[worker_idx], function(path) {",
    "    source(path)",
    "  }, cop_fun_path)",
    "}",
    "",
    "for (scenario_index in seq_len(nrow(true_cop_grid))) {",
    "",
    "  ii <- true_cop_grid$cop[scenario_index]",
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
    "    sim_seed <- seed_base + k + (ii - 1) * n_sim",
    "",
    "    data_i <- CobraDose::dataGen(",
    "      n = n,",
    "      gs = group_sizes,",
    "      marg1 = 1,",
    "      marg2 = 1,",
    "      cop = ii,",
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
    "      independence = mean(samples_all[, \"zc\"] == 1),",
    "      clayton      = mean(samples_all[, \"zc\"] == 2),",
    "      gaussian     = mean(samples_all[, \"zc\"] == 3),",
    "      gumbel       = mean(samples_all[, \"zc\"] == 4)",
    "    )",
    "  }",
    "",
    "  out_file <- file.path(output_dir, paste0(output_prefix, \"_cop_\", ii, \".csv\"))",
    "  utils::write.csv(sim_res, out_file, row.names = FALSE)",
    "",
    "  print(paste(\"Saved\", out_file))",
    "  print(CobraDose::summarizeModelProbs(sim_res[, 2:5]))",
    "}",
    "",
    "close(pb)",
    "parallel::stopCluster(cl)"
  )

  writeLines(cop_fun_code, cop_fun_file)
  writeLines(sim_cop_code, sim_cop_file)

  invisible(list(
    cop_fun = cop_fun_file,
    sim_cop = sim_cop_file
  ))
}
