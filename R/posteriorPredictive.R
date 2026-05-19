#' Posterior predictive inference using Bayesian model averaging
#'
#' This function takes in the output from the function `approximatePosterior()`.
#' The function outputs a table that summarizes the posterior predictive
#' probabilities for the various events input into the function with corresponding
#' credible intervals.
#'
#' @param output_dir Directory where the scripts should be written.
#' @param bma_post Output saved from a previous `approximatePosterior()` function call.
#' @param inner_rep Number of samples per posterior draw.
#' @param bounds1 A two-column matrix of lower and upper bounds for the first
#' endpoint that define joint events. The first column is for lower bounds, and the
#' second one is for upper bounds. Use `0`/`Inf` if there is no lower or upper bound.
#' These lower and upper bounds should not account for inversion of the first endpoint.
#' @param bounds2 A two-column matrix of lower and upper bounds for the second
#' endpoint that define joint events. The rows must correspond to the rows from `bounds1`. The
#' joint event in row \eqn{k} of the matrices is such that \eqn{\tilde{Y}_1 \in} `c(bounds1[k, 1], bounds1[k, 2])`
#' and \eqn{\tilde{Y}_2 \in} `c(bounds2[k, 1], bounds2[k, 2])`.
#' @param n_cores Number of parallel cores to use with parallelization.
#' @param credible credibility level for the credible intervals (0.95 by default).
#' @param output_prefix Prefix used for output CSV files.
#' @param overwrite Logical; if `FALSE`, existing files will not be overwritten.
#'
#' @importFrom foreach %dopar%
#'
#' @return Returns a table with posterior predictive probabilities and credible intervals
#' for each dose; this table will be output to a .csv file
#'
#' @examples
#' \dontrun{
#' ## save posterior approximation results
#' bma.post.temp <- approximatePosterior(...)
#'
#' ## obtain posterior predictive probabilities
#' posteriorPredictive(bma_post = bma.post.temp,
#'                     inner_rep = 500,
#'                     bounds1 = cbind(c(0, 0, 1, 1), c(1, 1, Inf, Inf)),
#'                     bounds2 = cbind(c(0, 1, 0, 1), c(1, Inf, 1, Inf)),
#'                     n_cores = max(1, parallel::detectCores() - 1),
#'                     credible = 0.95,
#'                     output_prefix = "post_pred",
#'                     overwrite = TRUE)
#' }
#'
#' @export

posteriorPredictive <- function(
    output_dir = ".",
    bma_post = NULL,
    inner_rep = 500,

    bounds1 = cbind(c(0, 0, 1, 1), c(1, 1, Inf, Inf)),
    bounds2 = cbind(c(0, 1, 0, 1), c(1, Inf, 1, Inf)),

    n_cores = max(1, parallel::detectCores() - 1),
    credible = 0.95,

    output_prefix = "post_pred",
    overwrite = TRUE
) {

  # ---- input checks ----

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("`output_dir` must be a single character string.")
  }

  if(!(methods::is(bma_post, "bma.sample"))){
    stop("Please input valid bma.sample object.")
  }

  if (!is.matrix(bounds1) || ncol(bounds1) != 2 || !is.numeric(bounds1[,1]) ||
      !is.numeric(bounds1[,2]) || nrow(bounds1) < 2) {
    stop("`bounds1` must be a two-column matrix with at least two rows.")
  }

  if (any(bounds1[,1] < 0) || any(bounds1[,2] < 0)) {
    stop("All entries in `bounds1` must be non-negative.")
  }

  if (!is.matrix(bounds2) || ncol(bounds2) != 2 || !is.numeric(bounds2[,1]) ||
      !is.numeric(bounds2[,2]) || nrow(bounds2) < 2) {
    stop("`bounds2` must be a two-column matrix with at least two rows.")
  }

  if (any(bounds2[,1] < 0) || any(bounds2[,2] < 0)) {
    stop("All entries in `bounds2` must be non-negative.")
  }

  if (nrow(bounds1) != nrow(bounds2)) {
    stop("`bounds1` and `bounds2` must have the same number of rows.")
  }

  if (!is.numeric(credible) || length(credible) != 1 || credible < 0.5 || credible >= 1) {
    stop("`credible` must be a single number between 0.5 and 1.")
  }

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    stop("`n_cores` must be a single positive number.")
  }

  if (!is.character(output_prefix) || length(output_prefix) != 1) {
    stop("`output_prefix` must be a single character string.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  out_file <- paste0(file.path(output_dir, output_prefix), "_probs.csv")

  if (!overwrite && (file.exists(out_file))) {
    stop("Output file already exists. Use `overwrite = TRUE` to replace it.")
  }

  ## load internal function for simulating from bivariate copula distribution
  ## n sample size (can be greater than n here)
  ## marg1 The marginal distribution for endpoint 1 (1 = gamma, 2 = lognormal)
  ## m1 The mean of the marginal distribution for endpoint 1
  ## c1 The coefficient of variation of the marginal distribution for endpoint 1
  ## marg2 The marginal distribution for endpoint 2 (1 = gamma, 2 = lognormal)
  ## m2 The mean of the marginal distribution for endpoint 2
  ## c2 The coefficient of variation of the marginal distribution for endpoint 2
  ## cop The copula model (1 = independence, 2 = Clayton, 3 = Gaussian, 4 = Gumbel)
  ## tau Kendall's tau (must be between 0 and 1 for this example)
  rBivarCop <- function(n, marg1, m1, c1, marg2, m2, c2, cop, tau) {

    ## convert from mean and variance to model parameters
    if (marg1 == 1L){ ## gamma shape and rate
      a1 <- 1/c1^2; b1 <- 1/(m1*c1^2)
    } else if (marg1 == 2L){ ## lognormal mean and sd
      s1 <- sqrt(log(1 + c1^2)); mu1  <- log(m1) - 0.5*s1^2
    }

    if (marg2 == 1L){
      a2 <- 1/c2^2; b2 <- 1/(m2*c2^2)
    } else if (marg2 == 2L){
      s2 <- sqrt(log(1 + c2^2)); mu2  <- log(m2) - 0.5*s2^2
    }

    ## sample from a copula
    if(cop == 1L) { ## independence
      u1 <- runif(n); u2 <- runif(n)
    } else if(cop == 2L) { ## Clayton
      theta <- 2*tau/(1-tau) ## convert from tau to theta
      u.temp <- copula::rCopula(n, copula::claytonCopula(theta)) ## use copula package
      u1 <- as.numeric(u.temp[,1]); u2 <- as.numeric(u.temp[,2])
    } else if (cop == 3L){  ## Gaussian
      theta <- sin(pi*tau/2)
      z1 <- stats::rnorm(n) ## convert Gaussian to unit scale
      z2 <- stats::rnorm(n, theta*z1, sqrt(1-theta^2))
      u1 <- stats::pnorm(z1); u2 <- stats::pnorm(z2)
    } else if (cop == 4L){  ## Gumbel
      theta <- 1/(1-tau) ## use copula package to generate
      u.temp <- copula::rCopula(n, copula::gumbelCopula(theta))
      u1 <- as.numeric(u.temp[,1]); u2 <- as.numeric(u.temp[,2])
    }

    x <- matrix(0, nrow = n, ncol = 2)
    ## copula observation to marignals
    if(marg1 == 1L) {
      x[,1] <- stats::qgamma(u1, shape=a1, rate=b1)
    } else {
      x[,1] <- stats::qlnorm(u1, meanlog=mu1, sdlog=s1)
    }

    if(marg2 == 1L) {
      x[,2] <- stats::qgamma(u2, shape=a2, rate=b2)
    } else {
      x[,2] <- stats::qlnorm(u2, meanlog=mu2, sdlog=s2)
    }

    return(x)
  }

  ## function to obtain posterior predictive probability and credible intervals
  ## mat a matrix of posterior predictive probability realizations (one per summary/event)
  ## cred the credibility level (taken from main function input)
  get_post_pred <- function(mat, cred){
    ppp <- colMeans(mat) ## posterior predictive probability is the mean
    lb <- apply(mat, 2, stats::quantile, probs = (1-cred)/2, names = FALSE) ## use quantile to get credible interval
    ub <- apply(mat, 2, stats::quantile, probs = 1 - (1-cred)/2, names = FALSE)

    res_list <- list(round(ppp, 3), round(lb, 3), round(ub, 3)) ## one probability and CI per summary
    return(res_list)
  }

  cores <- parallel::detectCores()
  if (n_cores > cores[1]-1){n_cores <- cores[1] - 1}
  cl <- snow::makeSOCKcluster(n_cores)

  samples_dose <- bma_post[[1]]
  invert <- bma_post[[4]]
  m <- nrow(samples_dose[[1]])
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(max = m, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  tab_probs <- NULL
  tab_ci_low <- NULL
  tab_ci_up <- NULL
  for (kk in 1:length(samples_dose)){
    cat(paste0("\nPosterior predictive sampling for dose ", kk))
    params.kk <- samples_dose[[kk]]
    z1 <- params.kk$z1; z2 <- params.kk$z2; zc <- params.kk$zc
    m1 <- params.kk$m1; m2 <- params.kk$m2
    c1 <- params.kk$c1; c2 <- params.kk$c2
    tau <- params.kk$tau

    ## explore posterior probability of model selection under different marginals
    sim_res <- foreach::foreach(k=1:m, .packages=c('copula'), .combine=rbind,
                                .errorhandling = "remove", .options.snow=opts) %dopar% {

                                set.seed(k + m*(kk-1))
                                d.temp <- rBivarCop(inner_rep, marg1 = z1[k], m1 = m1[k], c1 = c1[k],
                                                    marg2 = z2[k], m2 = m2[k], c2 = c2[k],
                                                    cop = zc[k], tau = tau[k])

                                if (invert == TRUE){
                                  d.temp[,1] <- 1/d.temp[,1]
                                }

                                ## return realizations of posterior predictive probability
                                res.temp <- NULL
                                for (ii in 1:nrow(bounds1)){
                                  res.temp <- c(res.temp, mean(bounds1[ii, 1] < d.temp[,1] & d.temp[,1] < bounds1[ii, 2] &
                                                                 bounds2[ii, 1] < d.temp[,2] & d.temp[,2] < bounds2[ii, 2]))
                                }

                                res.temp

                                }

    list.temp <- get_post_pred(sim_res, cred = credible)

    tab_probs <- rbind(tab_probs, list.temp[[1]])
    tab_ci_low <- rbind(tab_ci_low, list.temp[[2]])
    tab_ci_up <- rbind(tab_ci_up, list.temp[[3]])
  }

  rownames(tab_probs) <- sapply(1:length(samples_dose), function(x){paste0("dose_", x)})
  colnames(tab_probs) <- sapply(1:nrow(bounds1), function(x){paste0("event_", x)})

  rownames(tab_ci_low) <- sapply(1:length(samples_dose), function(x){paste0("dose_", x)})
  colnames(tab_ci_low) <- sapply(1:nrow(bounds1), function(x){paste0("event_", x)})

  rownames(tab_ci_up) <- sapply(1:length(samples_dose), function(x){paste0("dose_", x)})
  colnames(tab_ci_up) <- sapply(1:nrow(bounds1), function(x){paste0("event_", x)})

  utils::write.csv(tab_probs,
                   paste0(file.path(output_dir, output_prefix), "_probs.csv"),
                   row.names = TRUE)

  utils::write.csv(tab_ci_low,
                   paste0(file.path(output_dir, output_prefix), "_lbs.csv"),
                   row.names = TRUE)

  utils::write.csv(tab_ci_up,
                   paste0(file.path(output_dir, output_prefix), "_ubs.csv"),
                   row.names = TRUE)

  parallel::stopCluster(cl)
  close(pb)

  return(structure(list(post_pred_probs = tab_probs, ci_lb = tab_ci_low, ci_ub = tab_ci_up),
                   class = "post.pred"))

}

#' @title Printing summaries of posterior approximations
#' @name post.pred-methods
#'
#' @description Helper function to print post.pred objects
#'
#' @param x a post.pred object
#'
#' @rdname post.pred-methods
#' @method print post.pred
#' @keywords internal
#' @export

print.post.pred <- function(x, ...) {

  ## extract all three tables
  post_pred_probs <- x$post_pred_probs
  ci_lb <- x$ci_lb
  ci_ub <- x$ci_ub

  output_mat <- matrix("", nrow = nrow(post_pred_probs),
                       ncol = ncol(post_pred_probs))
  for (ii in 1:nrow(post_pred_probs)){
    for (jj in 1:ncol(post_pred_probs)){
      output_mat[ii,jj] <- paste0(post_pred_probs[ii,jj], " (",
                                  ci_lb[ii,jj], ", ", ci_ub[ii,jj], ")")
    }
  }

  ## print header for posterior predictive summaries
  cat("The posterior predictive probabilities and their credible intervals (in parentheses) are\n\n")

  rownames(output_mat) <- sapply(1:nrow(output_mat), function(x){paste0("Dose ", x)})
  colnames(output_mat) <- sapply(1:ncol(output_mat), function(x){paste0("Event ", x)})

  ## print summary table
  print(output_mat)

  ## Return invisibly
  invisible(x)
}
