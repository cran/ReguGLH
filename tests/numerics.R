library(ReguGLH)
set.seed(47)
e <- c(0, 0, .3, .8, 1.2, 2)
a <- c(-1, -2)
ridge <- shrinkage_moments(e, 1.5, parameters = a)
edge <- max(e) + 1.5*mean(e) + 1
contour <- list(center = edge/2, radius = edge/2 + .2)
numerical <- shrinkage_moments(e, 1.5, "general", matrix(a),
  function(z, a) 1/(z-a[1]), contour, 512)
stopifnot(max(abs(ridge$omega-numerical$omega)) < 1e-7,
          max(abs(ridge$delta-numerical$delta)) < 1e-7,
          max(abs(ridge$numerator-numerical$numerator)) < 1e-7)
# Complex conjugate poles, with simple roots: (x+1)(x^2+x+1).
cubic <- shrinkage_moments(e, 1.5, "cubic", matrix(c(1,2,2,1), 1))
numeric_cubic <- shrinkage_moments(e, 1.5, "general", matrix(1),
  function(z, a) 1/(1+2*z+2*z^2+z^3),
  list(center = edge/2, radius = edge/2+.1), 1024)
stopifnot(max(abs(cubic$omega-numeric_cubic$omega)) < 1e-7,
          max(abs(cubic$delta-numeric_cubic$delta)) < 1e-7,
          max(abs(cubic$numerator-numeric_cubic$numerator)) < 1e-7)
# Both sides of the p/n = 1 boundary, including q = 1 and p = 1.
for (p in c(1, 10, 55)) {
  Y <- matrix(rnorm(p*40), p, 40)
  X <- rbind(1, rep(0:1, each=20)); C <- matrix(c(0,1),2)
  single <- glh_test(Y,X,C,parameters=a)
  composite <- glh_test(Y,X,C,mode="composite",priors=rbind(c(1,0,0),c(1,0,0)),parameters=a)
  stopifnot(max(abs(single$p.value-composite$p.value)) < 1e-12,
            all(is.finite(single$statistic)))
  fit <- glh_test(Y,X,C,mode="composite",parameters=a,nsim=1000)
  stopifnot(all(fit$p.value > 0 & fit$p.value <= 1),
            min(eigen(fit$covariance_psd,symmetric=TRUE)$values) > -1e-10)
  for (mode in c("single", "composite")) {
    fit <- glh_test(Y,X,C,family="cubic",mode=mode,
      parameters=rbind(c(1,2,2,1),c(1,3,3,1.1)),nsim=1000)
    stopifnot(all(is.finite(fit$statistic)))
  }
}
# Invalid inputs must fail early.
fails <- function(expr) inherits(tryCatch(force(expr), error=identity), "error")
stopifnot(fails(shrinkage_moments(e,1,parameters=0)),
 fails(shrinkage_moments(e,1,"cubic",matrix(c(1,3,3,1),1))),
 fails(glh_test(Y,X,C,prior=c(-1,0,0))),
 fails(glh_test(Y,rbind(X,X[1,]),C)))
# Exercise both test modes through numerical integration end to end.
e <- pmax(eigen(tcrossprod(t(qr.resid(qr(t(X)),t(Y))))/(ncol(Y)-nrow(X)),
                symmetric=TRUE)$values,0)
edge <- max(e)+(nrow(Y)/(ncol(Y)-nrow(X)))*mean(e)+1
for (mode in c("single","composite")) {
  exact <- glh_test(Y,X,C,mode=mode,parameters=a,nsim=500)
  integrated <- glh_test(Y,X,C,mode=mode,family="general",parameters=matrix(a),
    regularizer=function(z,a) 1/(z-a[1]),
    contour=list(center=edge/2,radius=edge/2+.5),nquad=1024,nsim=500)
  stopifnot(max(abs(exact$statistic-integrated$statistic)) < 1e-7,
            is.numeric(integrated$moments$omega),
            !is.complex(integrated$moments$delta))
}
