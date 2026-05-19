#' Contour plots of the posterior predictive distribution
#'
#' This function takes in the output from the function `approximatePosterior()`.
#' The function outputs contour plots for each dose with the data input into
#' `approximatePosterior()` superimposed.
#'
#' @param output_dir Directory where the scripts should be written.
#' @param bma_post Output saved from a previous `approximatePosterior()` function call.
#' @param n_grid Number of grid points in each dimension for the contour plot (at most 500).
#' @param n_samp Number of observations simulated from the posterior predictive distribution for plotting
#' @param n_cores Number of parallel cores to use with parallelization.
#' @param output_prefix Prefix used for output PDF and CSV files.
#' @param overwrite Logical; if `FALSE`, existing CSV files will not be overwritten.
#'
#' @importFrom foreach %dopar%
#'
#' @return Returns a contour plot for each dose as output to a PDF file; the samples from
#' the posterior predictive distribution used to make the contour plot for each dose will
#' be output to CSV files. Those samples can be used to create more tailored contour plots.
#'
#' @examples
#' \dontrun{
#' ## save posterior approximation results
#' bma.post.temp <- approximatePosterior(...)
#'
#' ## obtain contour plots
#' contourDoses(bma_post = bma.post.temp,
#'              n_grid = 300,
#'              n_samp = 25000,
#'              n_cores = max(1, parallel::detectCores() - 1),
#'              output_prefix = "contour",
#'              overwrite = TRUE)
#' }
#'
#' @export

contourDoses <- function(
    output_dir = ".",
    bma_post = NULL,
    n_grid = 300,

    n_samp = 25000,
    n_cores = max(1, parallel::detectCores() - 1),

    output_prefix = "contour",
    overwrite = TRUE
) {

  # ---- input checks ----

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("`output_dir` must be a single character string.")
  }

  if(!(methods::is(bma_post, "bma.sample"))){
    stop("Please input valid bma.sample object.")
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid < 100 || n_grid > 500) {
    stop("`n_grid` must be a single integer between 100 and 500.")
  }

  if (!is.numeric(n_samp) || length(n_samp) != 1 || n_samp < 10000 || n_samp > 100000) {
    stop("`n_samp` must be a single integer between 10000 and 100000.")
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

  out_file <- paste0(file.path(output_dir, output_prefix), "_plot_dose1.pdf")

  if (!overwrite && (file.exists(out_file))) {
    stop("Output files already exist. Use `overwrite = TRUE` to replace them.")
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

  cores <- parallel::detectCores()
  if (n_cores > cores[1]-1){n_cores <- cores[1] - 1}
  cl <- snow::makeSOCKcluster(n_cores)

  samples_dose <- bma_post[[1]]
  invert <- bma_post[[4]]
  m <- n_samp
  doSNOW::registerDoSNOW(cl)
  pb <- utils::txtProgressBar(max = m, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  ppp_samp <- vector("list", length(samples_dose))
  for (kk in 1:length(samples_dose)){
    cat(paste0("\nPosterior predictive sampling for dose ", kk))
    params.kk <- samples_dose[[kk]]
    params.kk <- params.kk[sample(1:nrow(params.kk), n_samp, replace = TRUE),]
    z1 <- params.kk$z1; z2 <- params.kk$z2; zc <- params.kk$zc
    m1 <- params.kk$m1; m2 <- params.kk$m2
    c1 <- params.kk$c1; c2 <- params.kk$c2
    tau <- params.kk$tau

    ## explore posterior probability of model selection under different marginals
    sim_res <- foreach::foreach(k=1:m, .packages=c('copula'), .combine=rbind,
                                .errorhandling = "remove", .options.snow=opts) %dopar% {

                                set.seed(k + m*(kk-1))
                                d.temp <- rBivarCop(1, marg1 = z1[k], m1 = m1[k], c1 = c1[k],
                                                    marg2 = z2[k], m2 = m2[k], c2 = c2[k],
                                                    cop = zc[k], tau = tau[k])

                                if (invert == TRUE){
                                  d.temp[,1] <- 1/d.temp[,1]
                                }

                                as.numeric(d.temp)

                                }

    ppp_samp[[kk]] <- data.frame(endpoint1 = sim_res[,1], endpoint2 = sim_res[,2])
  }

  parallel::stopCluster(cl)
  close(pb)

  names(ppp_samp) <- sapply(1:length(samples_dose), function(x){paste0("dose", x)})

  for (kk in 1:length(samples_dose)){
    utils::write.csv(ppp_samp[[kk]],
                     paste0(file.path(output_dir, output_prefix), "_sample_dose", kk, ".csv"),
                     row.names = FALSE)
  }


  ## now obtain and output the contour plots
  data <- bma_post[[5]]
  group <- bma_post[[6]]

  lim_x <- range(data[,1])
  lim_y <- range(data[,2])
  for (kk in 1:length(samples_dose)){
    samp.kk <- ppp_samp[[kk]]
    lim_x_temp <- as.numeric(stats::quantile(samp.kk[,1], probs = c(0.025, 0.975)))
    if (lim_x_temp[1] < lim_x[1]){lim_x[1] <- lim_x_temp[1]}
    if (lim_x_temp[2] > lim_x[2]){lim_x[2] <- lim_x_temp[2]}

    lim_y_temp <- as.numeric(stats::quantile(samp.kk[,2], probs = c(0.025, 0.975)))
    if (lim_y_temp[1] < lim_y[1]){lim_y[1] <- lim_y_temp[1]}
    if (lim_y_temp[2] > lim_y[2]){lim_y[2] <- lim_y_temp[2]}
  }

  for (kk in 1:length(samples_dose)){
    samp.kk <- ppp_samp[[kk]]
    kde.temp = MASS::kde2d(samp.kk[,1], samp.kk[,2], n = n_grid, c(lim_x, lim_y))

    data.kk <- data[group == kk, ]

    grDevices::pdf(file = paste0(file.path(output_dir, output_prefix), "_plot_dose", kk, ".pdf"),
                   width = 6, height = 6)

    graphics::par(mar = c(3.75, 4.5, 1.75, 0.35) + 0.1, mgp=c(2.7,1,0))

    graphics::contour(kde.temp$x, kde.temp$y, kde.temp$z, xlim = lim_x, ylim = lim_y,
                      main = paste0("Dose ", kk), cex = 1.2, cex.axis = 1.2, cex.main = 1.2, cex.lab = 1.2, labcex = 1,
                      xlab = bquote(italic(tilde(Y)[1])), ylab = bquote(italic(tilde(Y)[2])),
                      method = "flattest")
    graphics::points(data.kk[,1], data.kk[,2],
                     pch = 19, col = grDevices::adjustcolor("firebrick", 0.6), cex = 1.2)

    grDevices::dev.off()

    ## repeat plot to the console
    graphics::par(mar = c(3.75, 4.5, 1.75, 0.35) + 0.1, mgp=c(2.7,1,0))

    graphics::contour(kde.temp$x, kde.temp$y, kde.temp$z, xlim = lim_x, ylim = lim_y,
                      main = paste0("Dose ", kk), cex = 1.2, cex.axis = 1.2, cex.main = 1.2, cex.lab = 1.2, labcex = 1,
                      xlab = bquote(italic(tilde(Y)[1])), ylab = bquote(italic(tilde(Y)[2])),
                      method = "flattest")
    graphics::points(data.kk[,1], data.kk[,2],
                     pch = 19, col = grDevices::adjustcolor("firebrick", 0.6), cex = 1.2)
  }

  return(ppp_samp)

}
