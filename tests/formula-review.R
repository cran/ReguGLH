library(ReguGLH)
# Independent large-z identities for f(x)=1 verify contour signs,
# the factor of two in Delta, and the rho recursion without ridge residues.
e <- c(.1,.3,.7,1.2,2,3); gamma <- .7
mu <- mean(e); pop2 <- mean(e^2)-gamma*mu^2
edge <- max(e)+gamma*mu+1
circle <- list(center=edge/2,radius=edge/2+1)
constant <- function(z,a) rep(1+0i,length(z))
m <- shrinkage_moments(e,gamma,'general',matrix(1),constant,circle,512)
stopifnot(abs(m$omega-gamma*mu)<1e-9,
          abs(m$delta-2*gamma*pop2)<1e-9,
          max(abs(m$numerator-c(1,mu,pop2)))<1e-9)
# A genuinely non-rational regularizer: contour and resolution invariance.
exponential <- function(z,a) exp(-a[1]*z)
a <- matrix(c(.1,.5),ncol=1)
u <- shrinkage_moments(e,gamma,'general',a,exponential,circle,128)
v <- shrinkage_moments(e,gamma,'general',a,exponential,
                      list(center=edge/2,radius=edge/2+1.3),512)
stopifnot(max(abs(u$omega-v$omega))<1e-9,
          max(abs(u$delta-v$delta))<1e-9,
          max(abs(u$numerator-v$numerator))<1e-9)
# Compare the complete statistic with a direct response-space construction,
# independent of the package's projected eigenbasis calculation.
set.seed(123)
Y <- matrix(rnorm(7*48),7,48)
X <- t(model.matrix(~factor(rep(1:3,each=16))-1))
C <- cbind(c(1,-1,0),c(0,1,-1)); n <- ncol(Y)-nrow(X); q <- ncol(C)
xx <- solve(tcrossprod(X)); ec <- eigen(crossprod(C,xx%*%C),symmetric=TRUE)
Q <- t(X)%*%xx%*%C%*%ec$vectors%*%diag(1/sqrt(ec$values))%*%t(ec$vectors)
residual <- Y-Y%*%t(X)%*%xx%*%X
S <- tcrossprod(residual)/n; e <- eigen(S,symmetric=TRUE)$values
edge <- max(e)+(nrow(Y)/n)*mean(e)+1
fit <- glh_test(Y,X,C,family='general',parameters=matrix(1),regularizer=constant,
               contour=list(center=edge/2,radius=edge/2+1),nquad=256)
roots <- eigen(crossprod(Y%*%Q)/n,symmetric=TRUE)$values
o <- (nrow(Y)/n)*mean(e)
d <- 2*(nrow(Y)/n)*(mean(e^2)-(nrow(Y)/n)*mean(e)^2)
expected <- sqrt(n/q/d)*c((1+o)*(sum(log1p(roots))-q*log1p(o)),
                         sum(roots)-q*o,
                         (1+o)^2*(sum(roots/(1+roots))-q*o/(1+o)))
stopifnot(max(abs(fit$statistic-expected))<1e-8)
# Reparameterizing a contrast basis leaves the tested hypothesis unchanged.
fit1 <- glh_test(Y,X,C,parameters=c(-.1,-1,-10))
fit2 <- glh_test(Y,X,C%*%matrix(c(1,2,0,1),2),parameters=c(-.1,-1,-10))
stopifnot(max(abs(fit1$statistic-fit2$statistic))<1e-8)
