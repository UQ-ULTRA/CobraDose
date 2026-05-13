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



#' Determine probability of each submodel having largest posterior probability
#'
#' This function takes in a matrix of posterior probabilities, with one column
#' per submodel and returns the probability of each submodel having the largest
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



