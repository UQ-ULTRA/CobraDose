#' Approximate posterior using Bayesian model averaging
#'
#' This function takes in a bivariate sample, prior distributions, and
#' parameters for various MCMC settings. The function outputs a data frame
#' with a sample from the posterior distribution, R hat statistics to
#' assess Markov chain convergence, and summaries of posterior probabilities
#' for candidate model selection
#'
#' @param output_dir Directory where the scripts should be written.
#' @param data Matrix with two columns (one column per biomarker).
#' @param group a vector of group assignments (length must be the column length of `data`).
#' @param hyper_mu Hyper-mean parameters for data analysis.
#' @param hyper_prec Hyper-precision parameters for data analysis.
#' @param alpha_marg Dirichlet prior parameters for the two marginal model probabilities.
#' @param alpha_cop Dirichlet prior parameters for the four copula models.
#' @param mcmc_niter Number of MCMC iterations per chain (16 chains used).
#' @param mcmc_nburnin Number of burn-in iterations.
#' @param mcmc_nthin Level of thinning for each MCMC chain.
#' @param output_prefix Prefix used for output CSV files.
#' @param overwrite Logical; if `FALSE`, existing files will not be overwritten.
#'
#' @details
#' It may take a few minutes for the MCMC model to compile using NIMBLE. Once compiled,
#' a progress bar will populate in the console to display progress for each MCMC chain.
#'
#' @importFrom nimble nimbleModel nimbleCode configureMCMC buildMCMC
#' @importFrom nimble compileNimble runMCMC registerDistributions
#' @importFrom nimble getNimbleOption
#'
#' @return Returns a list with the posterior sample, R hat statistics, and posterior probabilities for
#' candidate model selection; the posterior sample may be output to a .csv file
#'
#' @examples
#' \dontrun{
#' ## get example data set
#' data.temp <- cbind(rgamma(18, 2, 1), rgamma(18, 2, 1))
#'
#' ## save posterior approximation results
#' approximatePosterior(data = data.temp,
#'        group = rep(1:4, c(6, 4, 4, 4)),
#'        hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
#         hyper_prec = c(400, 400, 20, 800, 800, 20))
#' }
#'
#' @export

approximatePosterior <- function(
    output_dir = ".",
    data = NULL,
    group = NULL,

    hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
    hyper_prec = c(400, 400, 20, 800, 800, 20),

    alpha_marg = c(1.625, 1),
    alpha_cop = c(1, 1, 1, 1),

    mcmc_niter = 3500,
    mcmc_nburnin = 1000,
    mcmc_nthin = 1,

    output_prefix = "bma_post",
    overwrite = FALSE
) {

  # ---- input checks ----

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("`output_dir` must be a single character string.")
  }

  if (!is.matrix(data) || ncol(data) != 2 || !is.numeric(data[,1]) ||
      !is.numeric(data[,2]) || nrow(data) < 2) {
    stop("`data` must be a two-column matrix.")
  }

  if (any(data[,1] <= 0) || any(data[,2] <= 0)) {
    stop("All entries in `data` must be positive.")
  }

  if (sum(1:max(group) %in% group) != max(group)) {
    stop("the `group` numbers must start at 1 and all numbers between 1 and `max(group)` should be used.")
  }

  if (length(hyper_mu) != 6 || length(hyper_prec) != 6) {
    stop("`hyper_mu` and `hyper_prec` must both have length 6.")
  }

  if (length(alpha_marg) != 2 || any(alpha_marg <= 0)) {
    stop("`alpha_marg` must be a positive numeric vector of length 2.")
  }

  if (length(alpha_cop) != 4 || any(alpha_cop <= 0)) {
    stop("`alpha_cop` must be a positive numeric vector of length 4.")
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

  if (!is.numeric(mcmc_nthin) || length(mcmc_nthin) != 1 || mcmc_nthin <= 0) {
    stop("`mcmc_nthin` must be a single positive number.")
  }

  if (!is.character(output_prefix) || length(output_prefix) != 1) {
    stop("`output_prefix` must be a single character string.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  out_file <- paste0(file.path(output_dir, output_prefix), "_dose1.csv")

  if (!overwrite && (file.exists(out_file))) {
    stop("Output files already exist. Use `overwrite = TRUE` to replace them.")
  }

  ## set one chain per candidate model
  mcmc_nchains <- 16

  ## register bivariate copula distribution
  nimble::registerDistributions(
    list(
      dBivarCopula = list(
        BUGSdist = "dBivarCopula(marg1, m1, c1, marg2, m2, c2, cop, tau)",
        types = c("value = double(1)")
      )
    )
  )

  ## code for nimble model
  modelCode <- nimble::nimbleCode({

    z1 ~ dcat(pi_marg1[1:2]) ## categorical variable for marg1
    z2 ~ dcat(pi_marg2[1:2]) ## categorical variable for marg2
    zc ~ dcat(pi_cop[1:4]) ## categorical value for cop

    pi_marg1[1:2] ~ ddirch(alpha_marg[1:2]) ## prior on marginal distributions
    pi_marg2[1:2] ~ ddirch(alpha_marg[1:2]) ## prior on marginal distributions
    pi_cop[1:4]  ~ ddirch(alpha_cop[1:4]) ## prior on copula models

    mu_m1 ~ dnorm(hyper_mu[1], hyper_prec[1])
    mu_m2 ~ dnorm(hyper_mu[1], hyper_prec[1])
    mu_c1 ~ dnorm(hyper_mu[2], hyper_prec[2])
    mu_c2 ~ dnorm(hyper_mu[2], hyper_prec[2])
    mu_logittau ~ dnorm(hyper_mu[3], hyper_prec[3])

    sigma_m1 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)
    sigma_m2 ~ T(dnorm(hyper_mu[4], hyper_prec[4]), 1e-6, Inf)
    sigma_c1 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)
    sigma_c2 ~ T(dnorm(hyper_mu[5], hyper_prec[5]), 1e-6, Inf)
    sigma_logittau ~ T(dnorm(hyper_mu[6], hyper_prec[6]), 1e-6, Inf)

    for(g in 1:G){

      lower_m1[g] <- -mu_m1/sigma_m1
      lower_m2[g] <- -mu_m2/sigma_m2
      lower_c1[g] <- -mu_c1/sigma_c1
      lower_c2[g] <- -mu_c2/sigma_c2

      z_m1[g] ~ T(dnorm(0,1), lower_m1[g], Inf)
      z_m2[g] ~ T(dnorm(0,1), lower_m2[g], Inf)
      z_c1[g] ~ T(dnorm(0,1), lower_c1[g], Inf)
      z_c2[g] ~ T(dnorm(0,1), lower_c2[g], Inf)
      z_tau[g] ~ dnorm(0,1)

      m1[g] <- mu_m1 + sigma_m1*z_m1[g]
      m2[g] <- mu_m2 + sigma_m2*z_m2[g]
      c1[g] <- mu_c1 + sigma_c1*z_c1[g]
      c2[g] <- mu_c2 + sigma_c2*z_c2[g]

      logit_tau[g] <- mu_logittau + sigma_logittau*z_tau[g]
      tau[g] <- ilogit(logit_tau[g])

    }

    for(i in 1:N){

      y[i,1:2] ~ dBivarCopula(
        marg1 = z1,
        m1 = m1[group[i]], c1 = c1[group[i]],
        marg2 = z2,
        m2 = m2[group[i]], c2 = c2[group[i]],
        cop = zc,
        tau = tau[group[i]]  ## interpreted inside
      )
    }
  })

  ## check whether data are negatively dependent
  cor_data <- stats::cor(data[,1], data[,2])
  invert <- FALSE
  if (cor_data < 0){
    data[,1] <- 1/data[,1]
    invert <- TRUE
  }

  ## get initial values for each group
  inits.end1 <- initFromData(data[,1], group = group)
  inits.end2 <- initFromData(data[,2], group = group)

  ## get sample size and # of group
  N <- length(group)
  G <- max(group)

  ## define initial values for random variables used in MCMC (for NIMBLE)
  inits1 <- list(mu_m1 = inits.end1$mu_m, mu_c1 = inits.end1$mu_c,
                 mu_m2 = inits.end2$mu_m, mu_c2 = inits.end2$mu_c,
                 sigma_m1 = inits.end1$sigma_m, sigma_c1 = inits.end1$sigma_c,
                 sigma_m2 = inits.end2$sigma_m, sigma_c2 = inits.end2$sigma_c,
                 mu_logittau = 0, sigma_logittau = 1,
                 z_m1 = inits.end1$z_m, z_m2 = inits.end2$z_m,
                 z_c1 = inits.end1$z_c, z_c2 = inits.end2$z_c,
                 z_tau = rep(0,G), pi_marg1 = c(0.5, 0.5), pi_marg2 = c(0.5, 0.5),
                 pi_cop = rep(0.25, 4), z1 = 1, z2 = 1, zc = 1)

  ## create list of initial values (one for each chain)
  zz1 <- rep(1:2, each = 8)
  zz2 <- rep(c(1,1, 1, 1, 2, 2, 2, 2), 2)
  zzc <- rep(1:4, 4)
  for (j in 2:16){
    assign(paste0("inits", j), inits1)
    tmp <- get(paste0("inits", j))
    tmp$z1 <- zz1[j]
    tmp$z2 <- zz2[j]
    tmp$zc <- zzc[j]
    assign(paste0("inits", j), tmp)
  }

  inits <- list(rep(NULL, 16))
  for (j in 1:16){
    inits[[j]] <- get(paste0("inits", j))
  }

  model <- nimble::nimbleModel(modelCode, data = list(y = data), inits = inits[[1]],
                               constants = list(N = N, G = G,
                                        group = group, alpha_marg = alpha_marg,
                                        alpha_cop = alpha_cop,
                                        hyper_mu = hyper_mu,
                                        hyper_prec = hyper_prec), check = FALSE)

  conf <- nimble::configureMCMC(model)

  conf$addMonitors(c("z1", "z2", "zc", "m1","c1","m2", "c2", "tau"))

 conf$removeSamplers(c(
    "z_m1", "z_c1", "z_m2", "z_c2",
    "sigma_m1", "sigma_c1", "sigma_m2",
    "sigma_c2", "sigma_logittau"
  ))

  for(i in 1:G) {
    conf$addSampler(target = paste0("z_m1[", i, "]"), type = "slice")
    conf$addSampler(target = paste0("z_c1[", i, "]"), type = "slice")
    conf$addSampler(target = paste0("z_m2[", i, "]"), type = "slice")
    conf$addSampler(target = paste0("z_c2[", i, "]"), type = "slice")
  }

  conf$addSampler(target = "sigma_m1", type = "slice")
  conf$addSampler(target = "sigma_c1", type = "slice")
  conf$addSampler(target = "sigma_m2", type = "slice")
  conf$addSampler(target = "sigma_c2", type = "slice")
  conf$addSampler(target = "sigma_logittau", type = "slice")

  Rmcmc <- nimble::buildMCMC(conf)

  Cmodel <- nimble::compileNimble(model)
  Cmcmc  <- nimble::compileNimble(Rmcmc, project = Cmodel)

  samples <- nimble::runMCMC(Cmcmc, niter = mcmc_niter, nburnin = mcmc_nburnin,
                             nchains = mcmc_nchains, thin = mcmc_nthin,
                             inits = inits)

  ## combine the MCMC draws from all chains
  samplesAll <- do.call(rbind, samples)

  ## save group-specific parameters in smaller files
  dose_list <- NULL
  for (j in 1:G){
    dose_list[[j]] <- data.frame(
      z1 = samplesAll[,"z1"],
      z2 = samplesAll[,"z2"],
      zc = samplesAll[, "zc"],
      m1 = samplesAll[, paste0("m1[", j, "]")],
      c1 = samplesAll[, paste0("c1[", j, "]")],
      m2 = samplesAll[, paste0("m2[", j, "]")],
      c2 = samplesAll[, paste0("c2[", j, "]")],
      tau = samplesAll[, paste0("tau[", j, "]")]
    )
    utils::write.csv(dose_list[[j]], paste0(file.path(output_dir, output_prefix), "_dose", j, ".csv"),
              row.names = FALSE)
  }

  ## compute R hat statistics
  samplesArray <- array(NA, dim = c(nrow(samples[[1]]), length(samples),
                                    ncol(samples[[1]])))
  for (i in 1:length(samples)){
    samplesArray[,i,] <- samples[[i]]
  }

  rHat <- NULL
  for (i in 1:ncol(samples[[1]])){
    rHat[i] <- posterior::rhat(samplesArray[, , i])
  }
  names(rHat) <- colnames(samples[[1]])
  rHat <- rHat[!grepl("^(mu_|sigma_)", names(rHat))]

  ## compute the posterior probabilities for each model

  ## marginal distributions
  marg_summary <- c(round(mean(samplesAll[, "z1"] == 1 & samplesAll[, "z2"] == 1), 3),
                    round(mean(samplesAll[, "z1"] == 1 & samplesAll[, "z2"] == 2), 3),
                    round(mean(samplesAll[, "z1"] == 2 & samplesAll[, "z2"] == 1), 3),
                    round(mean(samplesAll[, "z1"] == 2 & samplesAll[, "z2"] == 2), 3))

  ## copula models
  cop_summary <- c(round(mean(samplesAll[, "zc"] == 1), 3),
                   round(mean(samplesAll[, "zc"] == 2), 3),
                   round(mean(samplesAll[, "zc"] == 3), 3),
                   round(mean(samplesAll[, "zc"] == 4), 3))

  ## both marginals and copula
  for (j in 1:4){
    assign(paste0("both", j),
           c(round(mean(samplesAll[, "z1"] == 1 & samplesAll[, "z2"] == 1 & samplesAll[, "zc"] == j),3),
             round(mean(samplesAll[, "z1"] == 1 & samplesAll[, "z2"] == 2 & samplesAll[, "zc"] == j),3),
             round(mean(samplesAll[, "z1"] == 2 & samplesAll[, "z2"] == 1 & samplesAll[, "zc"] == j),3),
             round(mean(samplesAll[, "z1"] == 2 & samplesAll[, "z2"] == 2 & samplesAll[, "zc"] == j),3)))
  }

  both_summary <- rbind(both1, both2, both3, both4)

  names(marg_summary) <- c("gamma_gamma", "gamma_lognormal", "lognormal_gamma", "lognormal_lognormal")
  names(cop_summary) <- c("independence", "clayton", "gaussian", "gumbel")

  colnames(both_summary) <- c("gamma_gamma", "gamma_lognormal", "lognormal_gamma", "lognormal_lognormal")
  rownames(both_summary) <- c("independence", "clayton", "gaussian", "gumbel")

  return(structure(list(samples_dose = dose_list, rHat = rHat,
                        post_probs = list(marg_summary, cop_summary, both_summary), invert = invert),
                   class = "bma.sample"))

}

#' @title Printing summaries of posterior approximations
#' @name bma.sample-methods
#'
#' @description Helper function to print bma.sample objects
#'
#' @param x an bma.sample object
#'
#' @rdname bma.sample-methods
#' @method print bma.sample
#' @keywords internal
#' @export

print.bma.sample <- function(x, ...) {

  ## extract rHat
  rhat <- x$rHat

  ## build R-hat table
  tab <- data.frame(Variable = names(rhat),
    R.hat = as.numeric(rhat),row.names = NULL)

  ## Sort by R-hat descending
  tab <- tab[order(tab$R.hat, decreasing = TRUE), , drop = FALSE]
  rownames(tab) <- NULL

  ## print header for R-hat statistics
  cat("The R-hat statistics for the posterior samples are as follows:\n\n")

  ## print R-hat table
  print(tab)

  ## Return invisibly
  invisible(x)

  ## extract posterior probability summaries
  both_summary <- x$post_probs[[3]]

  cat("\nThe posterior probabilities for the candidate models are as follows:\n\n")
  tab_probs <- cbind(rbind(both_summary, colSums(both_summary)), c(rowSums(both_summary), 1))
  rownames(tab_probs) <- c("Independence", "Clayton", "Gaussian", "Gumbel", "All")
  colnames(tab_probs) <- c("Gamma-Gamma", "Gamma-Lognormal", "Lognormal-Gamma", "Lognormal-Lognormal", "All")

  print(tab_probs)

  invert <- x[[4]]
  if (invert == TRUE){
    cat("\nNote: the ratios for the first endpoint were inverted to account for negative dependence between the endpoints.")
  }
}
