# Empirical resolvent kernels; sums include the zero eigenvalues.
.kernels <- function(eigenvalues, gamma) {
  at <- function(z) {
    m <- vapply(z, function(a) mean(1 / (eigenvalues - a)), complex(1))
    mp <- vapply(z, function(a) mean(1 / (eigenvalues - a)^2), complex(1))
    th <- 1 / (1 - gamma - gamma * z * m)
    rho1 <- th * (1 + z * m)
    list(theta = th, rho = cbind(m, rho1,
         th * (mean(eigenvalues) + z * rho1)),
         diagonal = gamma * (1 + z * m) * th^3 +
           gamma * z * (m + z * mp) * th^4)
  }
  delta <- function(a, b) {
    aa <- at(a); bb <- at(b)
    out <- outer(a * aa$theta, b * bb$theta, "-")
    den <- outer(a, b, "-")
    same <- abs(den) < 1e-10 * pmax(1, abs(outer(a, b, "+")))
    den[same] <- 1
    out <- outer(aa$theta, bb$theta) * (out / den - 1)
    ij <- which(same, arr.ind = TRUE)
    if (nrow(ij)) out[ij] <- aa$diagonal[ij[, 1]]
    out
  }
  list(at = at, delta = delta)
}

.real <- function(x) {
  if (any(!is.finite(x))) stop("Non-finite numerical integral or kernel.")
  if (max(abs(Im(x))) > 1e-6 * max(1, max(abs(Re(x)))))
    stop("Unexpected imaginary remainder; refine or change the contour.")
  Re(x)
}

#' Estimate null moments and prior power numerators
#' @param eigenvalues Nonnegative eigenvalues of the residual covariance,
#'   including its zero eigenvalues.
#' @param gamma Dimension divided by residual degrees of freedom.
#' @param family One of linear, cubic, general.
#' @param parameters Negative ridge parameters; cubic coefficient matrix with
#'   columns l0,l1,l2,l3; or matrix of user-defined parameter rows.
#' @param regularizer Vectorized function(z, parameter) supporting complex z.
#' @param contour List with center and radius for a counterclockwise circle.
#'   Must enclose the empirical spectrum and kernel singularities, but no
#'   singularities of the regularizers. Required for the general family.
#' @param nquad Number of equally spaced midpoint contour nodes.
#' @return List containing omega, delta, numerator (three prior basis columns),
#'   parameters, and complex numerical estimates for the general family.
#' @references Li, H., Aue, A. and Paul, D. (2020).
#'   High-dimensional general linear hypothesis tests via non-linear spectral
#'   shrinkage. Bernoulli, 26(4), 2541--2571. \doi{10.3150/19-BEJ1186}.
#' @export
shrinkage_moments <- function(eigenvalues, gamma,
    family = c("linear", "cubic", "general"), parameters,
    regularizer = NULL, contour = NULL, nquad = 256L) {
  family <- match.arg(family)
  if (!is.numeric(eigenvalues) || !length(eigenvalues) ||
      any(!is.finite(eigenvalues)) || any(eigenvalues < 0) ||
      length(gamma) != 1L || !is.finite(gamma) || gamma <= 0)
    stop("Supply nonnegative finite eigenvalues and positive gamma.")
  k <- .kernels(eigenvalues, gamma)
  if (family == "linear") {
    parameters <- matrix(as.numeric(parameters), ncol = 1)
    if (!nrow(parameters) || any(!is.finite(parameters)) || any(parameters >= 0))
      stop("Linear parameters must be finite and negative.")
    poles <- lapply(parameters[, 1], identity)
    weights <- lapply(poles, function(z) 1)
  } else if (family == "cubic") {
    parameters <- as.matrix(parameters)
    if (!is.numeric(parameters) || ncol(parameters) != 4L ||
        !nrow(parameters) || any(!is.finite(parameters)))
      stop("Cubic parameters need four finite columns: l0,l1,l2,l3.")
    poles <- weights <- vector("list", nrow(parameters))
    for (i in seq_len(nrow(parameters))) {
      a <- parameters[i, ]
      if (a[1] <= 0 || a[4] == 0) stop("Require l0 > 0 and l3 != 0.")
      upper <- ceiling(max(eigenvalues)) + 0.01 * mean(eigenvalues)
      xx <- c(0, upper)
      vertex <- -a[3] / (3 * a[4])
      if (vertex > 0 && vertex < upper) xx <- c(xx, vertex)
      if (min(a[2] + 2*a[3]*xx + 3*a[4]*xx^2) < -1e-12)
        stop("Cubic denominator must be nondecreasing on the spectral interval.")
      r <- polyroot(a)
      if (min(abs(r[1] - r[2]), abs(r[1] - r[3]), abs(r[2] - r[3])) <
          1e-6 * max(1, abs(r))) stop("Cubic roots must be distinct and separated.")
      poles[[i]] <- r
      weights[[i]] <- 1 / (a[2] + 2*a[3]*r + 3*a[4]*r^2)
    }
  } else {
    parameters <- as.matrix(parameters)
    if (!nrow(parameters) || any(!is.finite(parameters)) || !is.function(regularizer))
      stop("General family requires finite parameter rows and a regularizer.")
    if (is.null(contour) || length(contour$center) != 1L ||
        length(contour$radius) != 1L || !is.finite(contour$center) ||
        !is.finite(contour$radius) || contour$radius <= 0 ||
        any(abs(eigenvalues - contour$center) >= contour$radius))
      stop("Supply a circular contour enclosing the spectrum.")
    if (length(nquad) != 1L || !is.finite(nquad) || nquad < 32 || nquad != as.integer(nquad))
      stop("nquad must be an integer >= 32.")
    angle <- 2*pi*((seq_len(nquad)-0.5)/nquad)
    z <- contour$center + contour$radius * exp(1i*angle)
    # dz/(2*pi*i), including the quadrature spacing.
    dz <- contour$radius * exp(1i*angle) / nquad
    f <- vapply(seq_len(nrow(parameters)), function(i) {
      value <- regularizer(z, parameters[i, ])
      if (length(value) != length(z)) stop("regularizer must be vectorized.")
      as.complex(value)
    }, complex(nquad))
    w <- f * dz
    at <- k$at(z)
    omega <- -colSums(w * (at$theta - 1))
    numerator <- -t(w) %*% at$rho
    delta <- 2 * t(w) %*% k$delta(z, z) %*% w
    raw <- list(omega = omega, delta = delta, numerator = numerator)
    return(list(omega = .real(omega), delta = .real(delta),
      numerator = .real(numerator), parameters = parameters, complex = raw))
  }
  omega <- vapply(seq_along(poles), function(i)
    .real(sum(weights[[i]] * (k$at(poles[[i]])$theta - 1))), numeric(1))
  numerator <- t(vapply(seq_along(poles), function(i)
    .real(colSums(weights[[i]] * k$at(poles[[i]])$rho)), numeric(3)))
  delta <- matrix(0, length(poles), length(poles))
  for (i in seq_along(poles)) for (j in seq_len(i)) {
    delta[i,j] <- delta[j,i] <- .real(2 * sum(
      outer(weights[[i]], weights[[j]]) * k$delta(poles[[i]], poles[[j]])))
  }
  list(omega = omega, delta = delta, numerator = numerator, parameters = parameters)
}
