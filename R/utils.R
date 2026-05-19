#' Register the custom bivariate copula distribution for NIMBLE
#'
#' This function registers the custom `dBivarCopula` distribution so that it
#' can be used inside NIMBLE model code.
#'
#' @return Invisibly returns `NULL`.
#'
#' @examples
#' registerBivarCopulaDistribution()
#'
#' @export

registerBivarCopulaDistribution <- function() {
  nimble::registerDistributions(
    list(
      dBivarCopula = list(
        BUGSdist = "dBivarCopula(marg1, m1, c1, marg2, m2, c2, cop, tau)",
        types = c("value = double(1)")
      )
    )
  )

  invisible(NULL)
}

#' Helper function to get initial values for the means and CVs for each dose
#'
#' This function takes in the data corresponding to one endpoint and
#' returns a list of hyperparameters for MCMC initialization
#'
#' @param y a vector of observations (one endpoint)
#' @param group a vector denoting group membership
#'
#' @return A list of hyperparameters
#' @keywords internal

initFromData <- function(y, group) {

  G <- max(group) ## extract the number of group

  m <- numeric(G) ## initialize group-specific vectors for mean
  v <- numeric(G) ## variance
  c <- numeric(G) ## and cv

  z_m <- numeric(G) ## initialize variables on the standard normal scale
  z_c <- numeric(G)

  for(g in 1:G) {
    idx <- which(group == g) ## get observations in group
    y_g <- y[idx]

    m[g] <- mean(y_g) ## compute group-specific summaries
    v[g]  <- max(var(y_g), 1e-3)

    c[g] <- sqrt(v[g])/m[g]
  }

  mu_m <- mean(m) ## get initial hyperparameters for truncated normal on group means
  sigma_m <- max(sd(m), 1e-2)

  mu_c <- mean(c) ## initial hyperparameters for truncated normal on group CVs
  sigma_c <- max(sd(c), 1e-2)

  lower_m <- -mu_m/sigma_m
  lower_c <- -mu_c/sigma_c

  for (g in 1:G){  ## get initial normal variable values (truncate a value corresponding to 0)
    z_m[g] <- max((m[g] - mu_m)/sigma_m, lower_m + 1e-4)
    z_c[g] <- max((c[g]   - mu_c)/sigma_c, lower_c + 1e-4)
  }

  ## return list of hyperparameters
  list(z_m = z_m, z_c = z_c, mu_m = mu_m, mu_c = mu_c,
       sigma_m = sigma_m, sigma_c = sigma_c)
}

#' Determine the probability of each candidate model having largest posterior probability
#'
#' This function takes in a matrix of posterior probabilities, with one column
#' per candidate model and returns the probability of each candidate model having the largest
#' posterior probability
#'
#' @param mat a matrix of posterior probabilities
#'
#' @return A vector of probabilities
#' @export

summarizeModelProbs <- function(mat){

  # Coerce to matrix early so we have a consistent ncol() to refer to
  if (is.null(mat)) {
    warning("summarizeModelProbs received NULL input. Likely all simulations failed; check the foreach .errorhandling setting and worker errors.")
    return(NA_real_)   # NA (number of columns is unknown)
  }
  if (is.null(dim(mat))) {
    mat <- matrix(mat, nrow = 1)
  }
  if (nrow(mat) == 0) {
    warning("summarizeModelProbs received a matrix with zero rows. Likely all simulations failed; check the foreach .errorhandling setting and worker errors.")
    return(rep(NA_real_, ncol(mat)))   # NA based on number of columns
  }

  probs_max <- apply(mat, 1, max)

  probs <- numeric(ncol(mat))
  for (i in seq_len(ncol(mat))) { ## what percentage of time does model i have largest post prob
    probs[i] <- mean(mat[, i] >= probs_max)
  }
  return(probs)
}


#' Helper function to get initial values for the means and CVs for each dose
#'
#' This function takes in the data corresponding to one endpoint and
#' returns a list of hyperparameters for MCMC initialization. This function
#' is only used internally.
#'
#' @param y a vector of observations (one endpoint)
#' @param group a vector denoting group/dose membership
#'
#' @return A list of hyperparameters
#' @export

initFromData <- function(y, group) {

  G <- max(group) ## extract the number of group

  m <- numeric(G) ## initialize group-specific vectors for mean
  v <- numeric(G) ## variance
  c <- numeric(G) ## and cv

  z_m <- numeric(G) ## initialize variables on the standard normal scale
  z_c <- numeric(G)

  for(g in 1:G) {
    idx <- which(group == g) ## get observations in group
    y_g <- y[idx]

    m[g] <- mean(y_g) ## compute group-specific summaries
    v[g]  <- max(stats::var(y_g), 1e-3)

    c[g] <- sqrt(v[g])/m[g]
  }

  mu_m <- mean(m) ## get initial hyperparameters for truncated normal on group means
  sigma_m <- max(stats::sd(m), 1e-2)

  mu_c <- mean(c) ## initial hyperparameters for truncated normal on group CVs
  sigma_c <- max(stats::sd(c), 1e-2)

  lower_m <- -mu_m/sigma_m
  lower_c <- -mu_c/sigma_c

  for (g in 1:G){  ## get initial normal variable values (truncate a value corresponding to 0)
    z_m[g] <- max((m[g] - mu_m)/sigma_m, lower_m + 1e-4)
    z_c[g] <- max((c[g]   - mu_c)/sigma_c, lower_c + 1e-4)
  }

  ## return list of hyperparameters
  list(z_m = z_m, z_c = z_c, mu_m = mu_m, mu_c = mu_c,
       sigma_m = sigma_m, sigma_c = sigma_c)
}
