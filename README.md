# ReguGLH

An R implementation of **High-dimensional general linear hypothesis tests via
non-linear spectral shrinkage** by Haoran Li, Alexander Aue, and Debashis Paul
(2020), published in *Bernoulli*. Includes the higher-order regularizers described
in the supplementary material.

## Publication and citation

Li, H., Aue, A., and Paul, D. (2020).
High-dimensional general linear hypothesis tests via non-linear spectral shrinkage.
*Bernoulli*, **26**(4), 2541–2571.
[doi:10.3150/19-BEJ1186](https://doi.org/10.3150/19-BEJ1186).

To retrieve the publication citation and BibTeX entry in R:

```r
citation("ReguGLH")
```

## Install and run

```r
install.packages("ReguGLH_0.1.0.tar.gz", repos = NULL, type = "source")
library(ReguGLH)
set.seed(1)
Y <- matrix(rnorm(10 * 40), 10, 40)
X <- rbind(1, rep(c(0, 1), each = 20))
C <- matrix(c(0, 1), 2, 1)
fit <- glh_test(Y, X, C, family = "linear", mode = "single",
                prior = c(1, 0, 0))
fit$p.value
fit <- glh_test(Y, X, C, family = "cubic", mode = "composite",
                priors = diag(3), nsim = 10000)
fit$parameters
fit$p.value
```

The two independent choices are `family = "linear"`, `"cubic"`, or `"general"`,
and `mode = "single"` or `"composite"`. All calls return LR, LH, and BNP results.
Y is p by N, X is k by N, and C is k by q. The hypothesis is BC = 0.
The residual covariance is YY-transpose after projection, divided by N-k.

## Families and selection

* Linear: f(x) = 1/(x - ell), ell < 0. Pass a negative vector in `parameters`.
  Default: 60 log-spaced candidates from -0.01 to -100.
* Cubic: f(x) = 1/(l0 + l1*x + l2*x^2 + l3*x^3). Each `parameters` row
  contains l0,l1,l2,l3 in increasing degree order. Closed forms use simple
  roots and partial-fraction weights. Repeated/near-repeated roots are rejected,
  following the supplement. The default `cubic_grid()` implements the suggested
  ranges with a modest resolution of five per dimension; raise `resolution`
  and supply the grid to study sensitivity. Monotonicity is checked analytically.
* General: supply a finite matrix of parameter rows and a vectorized
  `regularizer(z, parameter)` accepting complex z. Null means, covariances,
  and prior power numerators use midpoint trapezoidal contour quadrature.

A single test selects the grid candidate maximizing the estimated power
numerator divided by sqrt(Delta) under `prior`. Composite tests select once
per row of `priors`, take componentwise maxima, and calibrate against correlated
Gaussian maxima. Negative covariance eigenvalues are set to zero as in the
manuscript. Monte Carlo p-values use (exceedances+1)/(nsim+1); coincident selections
use the exact univariate Gaussian tail. Set the R random seed before calling to
reproduce composite results. Priors are nonnegative weights of I, Sigma, Sigma^2.

Selection is optimal **over the supplied finite grid**, rather than a guarantee
of a continuous global optimum. The implementation currently computes the full
candidate covariance matrix, so memory and time grow quadratically with grid size.
It also uses a full p by p eigen decomposition; large-p optimization remains future work.

## General-family example and integration controls

```r
# Reproduce two ridge candidates using the general integration path.
# Use an upper bound enclosing poles of the empirical Theta kernel as well.
N <- ncol(Y); n <- N - nrow(X)
residual <- t(qr.resid(qr(t(X)), t(Y)))
e <- pmax(eigen(tcrossprod(residual)/n, symmetric = TRUE)$values, 0)
edge <- max(e) + (nrow(Y)/n)*mean(e) + 1
contour <- list(center = edge/2, radius = edge/2 + 0.1)
general <- glh_test(Y, X, C, family = "general", mode = "composite",
  parameters = matrix(c(-1, -2), ncol = 1),
  regularizer = function(z, parameter) 1/(z - parameter[1]),
  contour = contour, nquad = 512)
```

A user-supplied circle must enclose the spectrum and empirical kernel poles,
while excluding singularities of every regularizer. Analyticity of arbitrary
user functions cannot be established automatically. Check stability after doubling
`nquad` and changing the contour within its admissible region. Returned
`moments$complex` retains the unmodified complex quadrature approximations;
real parts enter the manuscript formulas only after checking imaginary remainders.
No absolute value is substituted for signed contour estimates. Nonpositive variance
estimates fail explicitly. The implementation requires positive shrinkage at the
observed eigenvalues for LR/LH/BNP.

## Validation and provenance

`tests/numerics.R` checks contour integrals against ridge and cubic reduced
forms (including conjugate complex roots), p below and above residual sample
size, one-dimensional responses/contrasts, coincident selections, and invalid inputs.
`tests/legacy.R` compares ridge statistics and calibrated p-values against the original
RGLHComlin implementation using the same grid and residual eigendecomposition.
The theoretical basis is the main manuscript's null limit, empirical counterparts,
prior rho recursion, and composite calibration, and the supplement's cubic algorithm.
Numerical agreement is not a finite-sample validity guarantee: the manuscript's
asymptotic assumptions still apply. A broader simulation study is recommended
before drawing conclusions about settings beyond those evaluated below.

License: GPL-3. Package author and maintainer: Haoran Li.
Alexander Aue and Debashis Paul are credited as coauthors of the cited publication.
Maintainer: Haoran Li <hzl0152@auburn.edu>.
Secondary contact: <haoranli1990@gmail.com>.

## Finite-sample validation and limitations

A reproducible small simulation used balanced three-group Gaussian MANOVA,
AR(1) error covariance with correlation 0.4, N = 90, p = 45 or 135, and 500
null replicates per setting. It evaluated all three families, three canonical
single priors, and the composite test on small explicit candidate grids.

Across these settings, rejection rates at nominal 5% were 3.8--6.6% for LR,
5.4--7.6% for LH, and 2.8--6.2% for BNP. The Monte Carlo standard error near
5% is approximately one percentage point. Four of 72 pointwise exact 95%
binomial confidence intervals excluded 5%; these comparisons are correlated
and no simultaneous size-control conclusion is claimed. LH shows mild
finite-sample inflation and BNP is conservative in some settings.

Power increased at every tested signal and sample-size step. At fixed signal
norm 0.8, composite power increased from 42.0--63.0% at (N,p) = (90,45) to
97.0--99.5% at (360,180), with 200 replicates per alternative setting.
This supports the expected trend but does not prove asymptotic consistency.

No method evaluation failed. General-family integration estimates at 256
and 512 nodes differed by at most 8.22e-15 in the checked exponential-family
examples. Such agreement does not guarantee accuracy for arbitrary user functions.
The study does not cover heavy tails or all covariance and alternative structures.
See the installed VALIDATION.md for numerical formula checks and method details.
