\name{glh_test}
\alias{glh_test}
\title{Spectral shrinkage general linear hypothesis test}
\description{Tests BC = 0 in Y = BX + error. Columns of Y are observations. Selection maximizes the estimated local power numerator divided by the square root of the estimated variance on a finite candidate grid.}
\usage{
glh_test(Y, X, C, family = c("linear", "cubic", "general"),
  mode = c("single", "composite"), prior = c(1, 0, 0),
  priors = diag(3), parameters = NULL, regularizer = NULL,
  contour = NULL, nquad = 256L, nsim = 10000L, alpha = 0.05)
}
\arguments{
\item{Y}{Numeric p by N response matrix.}
\item{X}{Numeric k by N full row rank design matrix.}
\item{C}{Numeric k by q full column rank contrast matrix.}
\item{family}{One of linear (ridge), cubic (reciprocal polynomial), general.}
\item{mode}{Single prior test or maximum across multiple priors.}
\item{prior}{Nonnegative length-three weights for a single test.}
\item{priors}{Matrix of nonnegative prior rows for a composite test.}
\item{parameters}{Candidate parameters; see shrinkage_moments. Defaults to 60 logarithmically spaced negative ridge parameters or cubic_grid.}
\item{regularizer}{Vectorized function(z, parameter) for general shrinkage.}
\item{contour}{Circular integration contour; see shrinkage_moments.}
\item{nquad}{Contour quadrature node count.}
\item{nsim}{Number of Gaussian simulations for composite calibration.}
\item{alpha}{Significance level.}
}
\value{A list with statistic, p.value and reject for LR, LH and BNP; component statistics, selected parameters, prior weights, moments, covariance before and after PSD projection, and residual dimensions.}
\examples{
set.seed(1)
Y <- matrix(rnorm(10 * 40), 10, 40)
X <- rbind(1, rep(c(0, 1), each = 20))
C <- matrix(c(0, 1), 2, 1)
glh_test(Y, X, C)
glh_test(Y, X, C, mode = "composite", nsim = 1000)
}
\references{
Li, H., Aue, A. and Paul, D. (2020).
High-dimensional general linear hypothesis tests via non-linear spectral shrinkage.
\emph{Bernoulli}, 26(4), 2541--2571.
\doi{10.3150/19-BEJ1186}.
}
\details{Calibration is asymptotic. A small Gaussian MANOVA study found
mild finite-sample LH size inflation (up to 7.6 percent at nominal 5 percent).
See the package README and installed VALIDATION.md.}
