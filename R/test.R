#' Generate a cubic coefficient search grid
#' @param eigenvalues Residual covariance eigenvalues, including zeros.
#' @param resolution Number of grid points per dimension (at least two).
#' @return Matrix with columns l0,l1,l2,l3 satisfying the supplementary
#'   material's monotonicity and discriminant restrictions.
#' @references Li, H., Aue, A. and Paul, D. (2020).
#'   High-dimensional general linear hypothesis tests via non-linear spectral
#'   shrinkage. Bernoulli, 26(4), 2541--2571. \doi{10.3150/19-BEJ1186}.
#' @export
cubic_grid <- function(eigenvalues, resolution = 5L) {
  if (!is.numeric(eigenvalues) || !length(eigenvalues) ||
      any(!is.finite(eigenvalues)) || any(eigenvalues < 0) || mean(eigenvalues) <= 0)
    stop("Supply nonnegative eigenvalues with positive mean.")
  if (length(resolution) != 1 || !is.finite(resolution) || resolution < 2 ||
      resolution != as.integer(resolution)) stop("resolution must be an integer >= 2.")
  mu <- mean(eigenvalues)
  upper <- ceiling(max(eigenvalues)) + .01 * mu
  g <- expand.grid(l3 = c(-1, 1) %o%
      (10^seq(-3, 2, length.out = resolution) / (upper^2 * mu)),
      inflection = seq(-.1, 1.1, length.out = resolution) * upper,
      slope = seq(0, 2, length.out = resolution))
  g$l2 <- -3 * g$l3 * g$inflection
  g$l1 <- g$slope - g$l3 * mu^2 - g$l2 * mu
  keep <- vapply(seq_len(nrow(g)), function(i) {
    a <- g[i, ]; x <- c(0, upper)
    if (a$inflection > 0 && a$inflection < upper) x <- c(x, a$inflection)
    discriminant <- 18*a$l3*a$l2*a$l1 - 4*a$l2^3 +
      a$l2^2*a$l1^2 - 4*a$l3*a$l1^3 - 27*a$l3^2
    min(a$l1 + 2*a$l2*x + 3*a$l3*x^2) >= 0 &&
      abs(discriminant) >= 1e-6 / upper^2
  }, logical(1))
  out <- cbind(l0 = 1, as.matrix(g[keep, c("l1", "l2", "l3")]))
  if (!nrow(out)) stop("No admissible cubic candidates; increase resolution.")
  out
}

#' Spectral shrinkage general linear hypothesis test
#'
#' Tests BC = 0 in Y = BX + error. Columns of Y are observations.
#' Selection maximizes the estimated local power numerator divided by the
#' square root of the estimated variance on a finite candidate grid.
#' @details Calibration is asymptotic. A small Gaussian MANOVA study found
#'   mild finite-sample LH size inflation (up to 7.6 percent at nominal
#'   5 percent). See the package README and installed VALIDATION.md.
#' @param Y Numeric p by N response matrix.
#' @param X Numeric k by N full row rank design matrix.
#' @param C Numeric k by q full column rank contrast matrix.
#' @param family One of linear (ridge), cubic (reciprocal polynomial), general.
#' @param mode Single prior test or maximum across multiple priors.
#' @param prior Nonnegative length-three weights for a single test.
#' @param priors Matrix of nonnegative prior rows for a composite test.
#' @param parameters Candidate parameters; see shrinkage_moments. Defaults to
#'   60 logarithmically spaced negative ridge parameters or cubic_grid.
#' @param regularizer Vectorized function(z, parameter) for general shrinkage.
#' @param contour Circular integration contour; see shrinkage_moments.
#' @param nquad Contour quadrature node count.
#' @param nsim Number of Gaussian simulations for composite calibration.
#' @param alpha Significance level.
#' @return A list with statistic, p.value and reject for LR, LH and BNP;
#'   component statistics, selected parameters, prior weights, moments,
#'   covariance before and after PSD projection, and residual dimensions.
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(10 * 40), 10, 40)
#' X <- rbind(1, rep(c(0, 1), each = 20))
#' C <- matrix(c(0, 1), 2, 1)
#' glh_test(Y, X, C)
#' glh_test(Y, X, C, mode = "composite", nsim = 1000)
#' @references Li, H., Aue, A. and Paul, D. (2020).
#'   High-dimensional general linear hypothesis tests via non-linear spectral
#'   shrinkage. Bernoulli, 26(4), 2541--2571. \doi{10.3150/19-BEJ1186}.
#' @export
glh_test <- function(Y, X, C, family = c("linear", "cubic", "general"),
    mode = c("single", "composite"), prior = c(1, 0, 0), priors = diag(3),
    parameters = NULL, regularizer = NULL, contour = NULL,
    nquad = 256L, nsim = 10000L, alpha = .05) {
  family <- match.arg(family); mode <- match.arg(mode)
  for (a in list(Y, X, C)) if (!is.matrix(a) || !is.numeric(a) ||
      any(!is.finite(a)) || !all(dim(a) > 0)) stop("Y, X, C must be finite numeric matrices.")
  p <- nrow(Y); N <- ncol(Y); k <- nrow(X); q <- ncol(C); n <- N-k
  if (ncol(X) != N || nrow(C) != k || n <= 0 ||
      qr(t(X))$rank != k || qr(C)$rank != q)
    stop("Incompatible dimensions or rank-deficient design/contrasts.")
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must lie strictly between zero and one.")
  if (length(nsim) != 1 || !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim))
    stop("nsim must be a positive integer.")
  prior_matrix <- if (mode == "single") matrix(prior, nrow = 1) else as.matrix(priors)
  if (!is.numeric(prior_matrix) || ncol(prior_matrix) != 3 || !nrow(prior_matrix) ||
      any(!is.finite(prior_matrix)) || any(prior_matrix < 0) ||
      any(rowSums(prior_matrix) <= 0)) stop("Each prior must have three nonnegative weights and positive sum.")
  xx <- solve(tcrossprod(X))
  cx <- eigen(crossprod(C, xx %*% C), symmetric = TRUE)
  invsqrt <- cx$vectors %*% diag(1/sqrt(cx$values), nrow = q) %*% t(cx$vectors)
  Q <- t(X) %*% xx %*% C %*% invsqrt
  residual <- t(qr.resid(qr(t(X)), t(Y)))
  # A full eigenbasis also captures the nullspace when p > n.
  ee <- eigen(tcrossprod(residual)/n, symmetric = TRUE)
  eigenvalues <- pmax(ee$values, 0)
  eigenvalues[eigenvalues < max(eigenvalues) * 1e-12] <- 0
  if (mean(eigenvalues) <= 0) stop("Residual covariance is zero.")
  projected <- crossprod(ee$vectors, Y %*% Q) / sqrt(n)
  if (is.null(parameters)) {
    parameters <- switch(family,
      linear = -10^seq(-2, 2, length.out = 60),
      cubic = cubic_grid(eigenvalues),
      general = stop("Specify a parameter grid for the general family."))
  }
  moments <- shrinkage_moments(eigenvalues, p/n, family, parameters,
                              regularizer, contour, nquad)
  variance <- diag(moments$delta)
  if (any(!is.finite(variance)) || any(variance <= 0))
    stop("Nonpositive variance estimate; change grid or refine contour.")
  objective <- sweep(moments$numerator %*% t(prior_matrix), 1, sqrt(variance), "/")
  selected <- apply(objective, 2, which.max)
  components <- t(vapply(selected, function(i) {
    a <- moments$parameters[i, ]
    f <- switch(family, linear = 1/(eigenvalues - a[1]),
      cubic = 1/(a[1] + a[2]*eigenvalues + a[3]*eigenvalues^2 + a[4]*eigenvalues^3),
      general = regularizer(eigenvalues, a))
    if (length(f) != p) stop("regularizer must return one value per eigenvalue.")
    f <- .real(f)
    if (any(f <= 0)) stop("Shrinkage must be positive on the observed spectrum.")
    roots <- eigen(crossprod(projected, f * projected), symmetric = TRUE, only.values = TRUE)$values
    o <- moments$omega[i]; d <- variance[i]
    if (o <= -1) stop("Invalid mean for LR and BNP transformations.")
    sqrt(n/q/d) * c(LR = (1+o)*(sum(log1p(roots))-q*log1p(o)),
      LH = sum(roots)-q*o,
      BNP = (1+o)^2*(sum(roots/(1+roots))-q*o/(1+o)))
  }, numeric(3)))
  colnames(components) <- c("LR", "LH", "BNP")
  covariance <- moments$delta[selected, selected, drop = FALSE] /
    sqrt(outer(variance[selected], variance[selected]))
  covariance <- (covariance + t(covariance))/2
  diag(covariance) <- 1
  ev <- eigen(covariance, symmetric = TRUE)
  factor <- sweep(ev$vectors, 2, sqrt(pmax(ev$values, 0)), "*")
  covariance_psd <- tcrossprod(factor)
  statistic <- apply(components, 2, max)
  if (mode == "single" || length(unique(selected)) == 1) {
    pvalue <- stats::pnorm(statistic, lower.tail = FALSE)
  } else {
    maxima <- apply(factor %*% matrix(stats::rnorm(length(selected)*nsim),
                                    nrow = length(selected)), 2, max)
    pvalue <- vapply(statistic, function(s) (1+sum(maxima >= s))/(nsim+1), numeric(1))
  }
  list(statistic = statistic, p.value = pvalue, reject = pvalue < alpha,
    family = family, mode = mode, alpha = alpha,
    selected = selected, parameters = moments$parameters[selected, , drop = FALSE],
    priors = prior_matrix, components = components, objective = objective,
    moments = moments, covariance = covariance, covariance_psd = covariance_psd,
    n = n, p = p, q = q, nsim = if (mode == "composite") nsim else 0L)
}
