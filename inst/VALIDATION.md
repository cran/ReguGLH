# Numerical implementation review

Reference: Li, Aue and Paul (2020), Bernoulli 26(4), 2541--2571,
doi:10.3150/19-BEJ1186, and supplementary material Section S.1.

The supplied manuscript TeX sources were reviewed against the implementation.
This is an implementation review, not an independent peer review of the theory.

## Formula mapping

* Model and hypothesis: Y is p by N, X is k by N, C is k by q; H0: BC=0.
  Residual degrees of freedom are n=N-k, not N-2k. Residual covariance is
  Y(I-X' (XX')^-1 X)Y'/n. The implementation computes residuals with QR.
* Q=X'(XX')^-1 C [C'(XX')^-1 C]^-1/2. The normalized projected responses
  yield M(f)=Q'Y'f(S)YQ/n. A complete eigenbasis retains nullspace components
  when p exceeds residual sample size. Zero eigenvalues are included in m(z).
* m(z)=mean(1/(e-z)), m'(z)=mean(1/(e-z)^2),
  Theta(z)=[1-gamma-gamma*z*m(z)]^-1. Both the diagonal derivative formula
  and off-diagonal divided-difference formula for delta match the null theorem.
* Omega=-integral f(z)(Theta(z)-1) dz/(2*pi*i),
  Delta=2*double-integral f(z1)f(z2)delta(z1,z2) dz1 dz2/(2*pi*i)^2.
  Counterclockwise midpoint quadrature uses dz/(2*pi*i)=r*exp(i*t)/nquad.
  Complex calculations retain their signs. Real estimates are returned after
  checking imaginary remainders; complex absolute values are not substituted.
* The three prior kernels are m(z), Theta(z)(1+z*m(z)), and
  Theta(z)(mean(e)+z*rho1(z)). Selection maximizes their weighted integral
  divided by sqrt(Delta), on the supplied finite grid.
* Ridge moments use Omega=Theta(ell)-1 and Delta=2*delta(ell,ell).
* Cubic moments use simple roots r_j and weights 1/P'(r_j), giving
  Omega=sum w_j(Theta(r_j)-1) and
  Delta=2*sum w_j*w_k*delta(r_j,r_k). No conjugation is needed in this
  bilinear form; the conjugated form in the supplement is equivalent because
  roots and weights for a real polynomial occur in conjugate pairs.
* LR, LH and BNP use the manuscript's respective centering and scaling,
  including sqrt(n/q/Delta), (1+Omega), and (1+Omega)^2.
* Composite calibration uses standardized cross-Delta, symmetrization, and
  clipping negative eigenvalues to zero before Gaussian maximum simulation.
  Exactly coincident selected functions use the univariate normal tail.
  Monte Carlo p-values use (exceedances+1)/(nsim+1), a documented small
  difference from the simulation code's uncorrected exceedance fraction.

## Meaningful numerical checks

The package tests compare reduced ridge and cubic moments with contour
integration, including conjugate cubic roots, and compare ridge statistics
with the supplied simulation code on a shared parameter grid.

Independent identities for f(x)=1 check Omega=gamma*mean(e),
Delta=2*gamma*(mean(e^2)-gamma*mean(e)^2), and the three power numerators
(1, mean(e), mean(e^2)-gamma*mean(e)^2). A direct response-space calculation
checks all three standardized statistics without using the implementation's
projected eigenbasis. Additional checks cover contrast-basis invariance,
exponential regularization with two contours and resolutions, p>n,
one-dimensional responses/contrasts, and coincident selections.

## Scope and finite-sample results

Optimization is over finite candidate grids; no continuous global optimum
is guaranteed. Priors are restricted to nonnegative weights of I, Sigma,
and Sigma^2. Cubic denominators must be positive and nondecreasing on the
observed spectral interval, with separated simple roots. General functions
must accept complex arguments and produce positive shrinkage on observed
eigenvalues. The user must choose a contour enclosing the empirical spectrum
and kernel poles while excluding regularizer singularities, and check
quadrature convergence. Imaginary-remainder checks alone cannot establish
quadrature accuracy or analyticity.

A small simulation used three balanced groups, Gaussian errors with AR(1)
correlation 0.4, N=90 and p=45 or 135, with 500 null replicates per setting.
It tested three families, three single priors and composite calibration.
Nominal 5% rejection rates ranged from 3.8--6.6% (LR), 5.4--7.6% (LH),
and 2.8--6.2% (BNP). Four of 72 pointwise exact 95% intervals excluded 5%.
The null Monte Carlo standard error near 5% was about one percentage point.
LH displayed mild size inflation; BNP was conservative in some settings.

With 200 replicates per alternative, power increased at every tested signal
and sample-size step. Composite power at fixed signal norm 0.8 increased
from 42.0--63.0% at (N,p)=(90,45) to 97.0--99.5% at (360,180).
No evaluations failed. Checked exponential-family moment estimates at 256
and 512 nodes differed by at most 8.22e-15. These results support numerical
behavior in the evaluated settings, not universal size control or a proof
of consistency. Heavy tails and other covariance/alternative structures
require further study.
