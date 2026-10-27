\name{shrinkage_moments}
\alias{shrinkage_moments}
\title{Estimate null moments and prior power numerators}
\description{Estimate null moments and prior power numerators}
\usage{
shrinkage_moments(eigenvalues, gamma,
  family = c("linear", "cubic", "general"), parameters,
  regularizer = NULL, contour = NULL, nquad = 256L)
}
\arguments{
\item{eigenvalues}{Nonnegative eigenvalues of the residual covariance, including its zero eigenvalues.}
\item{gamma}{Dimension divided by residual degrees of freedom.}
\item{family}{One of linear, cubic, general.}
\item{parameters}{Negative ridge parameters; cubic coefficient matrix with columns l0,l1,l2,l3; or matrix of user-defined parameter rows.}
\item{regularizer}{Vectorized function(z, parameter) supporting complex z.}
\item{contour}{List with center and radius for a counterclockwise circle. Must enclose the empirical spectrum and kernel singularities, but no singularities of the regularizers. Required for the general family.}
\item{nquad}{Number of equally spaced midpoint contour nodes.}
}
\value{List containing omega, delta, numerator (three prior basis columns), parameters, and complex numerical estimates for the general family.}
\references{
Li, H., Aue, A. and Paul, D. (2020).
High-dimensional general linear hypothesis tests via non-linear spectral shrinkage.
\emph{Bernoulli}, 26(4), 2541--2571.
\doi{10.3150/19-BEJ1186}.
}
