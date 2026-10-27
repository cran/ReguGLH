\name{ReguGLH-package}
\alias{ReguGLH}
\alias{ReguGLH-package}
\docType{package}
\title{General Linear Hypothesis Tests via Spectral Shrinkage}
\description{Implements the methods of Li, Aue and Paul (2020), including ridge,
reciprocal cubic polynomial, and general analytic regularization.}
\details{P-values use asymptotic approximations. A small Gaussian MANOVA study
found mild finite-sample LH size inflation; see the package README and installed
VALIDATION.md for the design and results. General analytic functions require an
admissible contour and numerical convergence checks.}
\references{Li, H., Aue, A. and Paul, D. (2020).
High-dimensional general linear hypothesis tests via non-linear spectral shrinkage.
\emph{Bernoulli}, 26(4), 2541--2571. \doi{10.3150/19-BEJ1186}.}
\seealso{\code{\link{glh_test}}, \code{\link{shrinkage_moments}}, \code{\link{cubic_grid}}}
