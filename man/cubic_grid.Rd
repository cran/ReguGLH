\name{cubic_grid}
\alias{cubic_grid}
\title{Generate a cubic coefficient search grid}
\description{Generate a cubic coefficient search grid}
\usage{
cubic_grid(eigenvalues, resolution = 5L)
}
\arguments{
\item{eigenvalues}{Residual covariance eigenvalues, including zeros.}
\item{resolution}{Number of grid points per dimension (at least two).}
}
\value{Matrix with columns l0,l1,l2,l3 satisfying the supplementary material's monotonicity and discriminant restrictions.}
\references{
Li, H., Aue, A. and Paul, D. (2020).
High-dimensional general linear hypothesis tests via non-linear spectral shrinkage.
\emph{Bernoulli}, 26(4), 2541--2571.
\doi{10.3150/19-BEJ1186}.
}
