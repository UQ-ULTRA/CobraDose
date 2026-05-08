#' Compute (log) density of the bivariate distribution
#'
#' This function takes in a bivariate x as well as inputs for the
#' marginal distributions, means, coefficients of variation, copula,
#' and Kendall's tau and returns the density function evaluated at x
#'
#' @param x Bivarate observation
#' @param marg1 The marginal distribution for endpoint 1 (1 = gamma, 2 = lognormal)
#' @param m1 The mean of the marginal distribution for endpoint 1
#' @param c1 The coefficient of variation of the marginal distribution for endpoint 1
#' @param marg2 The marginal distribution for endpoint 2 (1 = gamma, 2 = lognormal)
#' @param m2 The mean of the marginal distribution for endpoint 2
#' @param c2 The coefficient of variation of the marginal distribution for endpoint 2
#' @param cop The copula model (1 = independence, 2 = Clayton, 3 = Gaussian, 4 = Gumbel)
#' @param tau Kendall's tau (must be between 0 and 1 for this example)
#' @param log Boolean that determines whether (1) or not (0) to return density on log scale
#' @importFrom stats dgamma dlnorm pgamma plnorm pnorm qgamma qlnorm qnorm rexp rgamma rnorm runif
#'
#' @return A scalar density or log-density value
#' @export

dBivarCopula <- nimble::nimbleFunction(
  run = function(x = double(1),
                 marg1 = integer(0), m1 = double(0), c1 = double(0),
                 marg2 = integer(0), m2 = double(0), c2 = double(0),
                 cop = integer(0), tau = double(0),
                 log = integer(0, default = 0)) {

    returnType(double(0))
    y1 <- x[1]; y2 <- x[2] ## extract data

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

    ## compute marginal log-densities
    if(marg1 == 1L) logf1 <- dgamma(y1, shape=a1, rate=b1, log=TRUE)
    if(marg1 == 2L) logf1 <- dlnorm(y1, meanlog=mu1, sdlog=s1, log=TRUE)

    if(marg2 == 1L) logf2 <- dgamma(y2, shape=a2, rate=b2, log=TRUE)
    if(marg2 == 2L) logf2 <- dlnorm(y2, meanlog=mu2, sdlog=s2, log=TRUE)

    ## map y1 and y2 to [0,1] via the marginal CDF
    if(marg1 == 1L) u1 <- pgamma(y1, shape=a1, rate=b1)
    if(marg1 == 2L) u1 <- plnorm(y1, meanlog=mu1, sdlog=s1)

    if(marg2 == 1L) u2 <- pgamma(y2, shape=a2, rate=b2)
    if(marg2 == 2L) u2 <- plnorm(y2, meanlog=mu2, sdlog=s2)

    ## guard against numerical boundary issues
    eps <- 1e-12
    u1 <- max(eps, min(1 - eps, u1))
    u2 <- max(eps, min(1 - eps, u2))

    ## compute the copula log-density
    if(cop == 1L) { ## independence
      logc <- 0.0
    } else if(cop == 2L) { ## clayton
      theta <- 2*tau/(1-tau)
      if(theta < 1e-12) {
        logc <- 0} else {
          logc <- log(theta + 1) - (theta + 1) * (log(u1) + log(u2)) -
            (2 + 1/theta) * log(u1^(-theta) + u2^(-theta) - 1)}

    } else if(cop == 3L) { ## gaussian
      theta <- sin(pi*tau/2)
      a <- qnorm(u1); b <- qnorm(u2)
      logc <- -0.5*log(1 - theta^2) -
        0.5*((a^2 + b^2)*theta^2 - 2*a*b*theta)/(1 - theta^2)
    } else if(cop == 4L) { ## gumbel
      theta <- 1/(1-tau)
      if(theta < 1e-12) {
        logc <- 0} else {
          a <- -log(u1); b <- -log(u2)
          log_a <- log(a); log_b <- log(b)

          M <- max(theta*log_a, theta*log_b)
          logS <- M + log(exp(theta*log_a - M) + exp(theta*log_b - M))

          t <- exp(logS/theta)

          logc <- -t + (theta - 1)*(log_a + log_b) -
            log(u1) - log(u2) + (2/theta - 2)*logS +
            log1p((theta - 1)*exp(-logS/theta))
        }
    }

    ## get log-density by adding contributions from marginals and copula
    logdens <- logf1 + logf2 + logc
    if(log) return(logdens) else return(exp(logdens))
  }
)

#' Generate from the bivariate distribution (only used to stop error messages)
#'
#' This function takes in a sample size of n = 1 as well as inputs for the
#' marginal distributions, means, variances, copula, and Kendall's tau
#' and returns a bivariate observation x
#'
#' @param n sample size (must be 1)
#' @param marg1 The marginal distribution for endpoint 1 (1 = gamma, 2 = lognormal)
#' @param m1 The mean of the marginal distribution for endpoint 1
#' @param c1 The coefficient of variation of the marginal distribution for endpoint 1
#' @param marg2 The marginal distribution for endpoint 2 (1 = gamma, 2 = lognormal)
#' @param m2 The mean of the marginal distribution for endpoint 2
#' @param c2 The coefficient of variation of the marginal distribution for endpoint 2
#' @param cop The copula model (1 = independence, 2 = Clayton, 3 = Gaussian, 4 = Gumbel)
#' @param tau Kendall's tau (must be between 0 and 1 for this example)
#' @importFrom nimble nimStop nimNumeric
#'
#' @return A bivariate observation
#' @export

rBivarCopula <- nimble::nimbleFunction(
  run = function(n = integer(0), marg1 = integer(0), m1 = double(0), c1 = double(0),
                 marg2 = integer(0), m2 = double(0), c2 = double(0),
                 cop = integer(0), tau = double(0)) {

    returnType(double(1))  ## length-2 vector

    if(n != 1) stop("Only n=1 supported for NIMBLE r-function")

    x <- numeric(2, init = FALSE)

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

    ## sample (u1, u2) from the desired copula

    if (cop == 1L) {
      ## independence
      u1 <- runif(1); u2 <- runif(1)

    } else if (cop == 2L) {
      ## Clayton via gamma mixture
      theta <- 2 * tau / (1 - tau)

      if (theta < 1e-8) {
        u1 <- runif(1)
        u2 <- runif(1)
      } else {
        W <- rgamma(1, shape = 1/theta, rate = 1)
        E1 <- rexp(1)
        E2 <- rexp(1)
        u1 <- (1 + E1/W)^(-1/theta)
        u2 <- (1 + E2/W)^(-1/theta)
      }

    } else if (cop == 3L) {
      ## Gaussian copula
      rho <- sin(pi * tau / 2)

      z1 <- rnorm(1)
      z2 <- rnorm(1, mean = rho*z1, sd = sqrt(1 - rho^2))

      u1 <- pnorm(z1); u2 <- pnorm(z2)

    } else if (cop == 4L) {
      ## Gumbel copula via positive stable approximation
      theta <- 1 / (1 - tau)

      if (theta <= 1 + 1e-8) {
        u1 <- runif(1); u2 <- runif(1)
      } else {
        ## simulate positive stable V (alpha = 1/theta)
        alpha <- 1/theta
        U <- runif(1, 0, pi)
        E <- rexp(1)

        V <- (sin(alpha*U)/(sin(U))^(alpha))*
          ((sin((1 - alpha)*U)/E)^((1 - alpha)))

        E1 <- rexp(1)
        E2 <- rexp(1)

        u1 <- exp(-E1/V); u2 <- exp(-E2/V)
      }
    }

    ## transform copulas to marginal distribution
    if (marg1 == 1L) {
      x[1] <- qgamma(u1, shape = a1, rate = b1)
    } else {
      x[1] <- qlnorm(u1, meanlog = mu1, sdlog = s1)
    }

    if (marg2 == 1L) {
      x[2] <- qgamma(u2, shape = a2, rate = b2)
    } else {
      x[2] <- qlnorm(u2, meanlog = mu2, sdlog = s2)
    }

    return(x)
  }
)
