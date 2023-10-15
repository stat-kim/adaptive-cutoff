set.seed(1)
library(abind)
library(mvnfast)

N <- 60 # Two dimensional unit square regular grid of size 60×60
nu <- 0.5 # Smoothness parameter: 0.5 or 1.0 

cov.matern=function(x, nu=1, beta=1, vars=1)
{
  if(nu == 0.5)
    return(vars*exp( - x / beta)) 
  ismatrix <- is.matrix(x)
  if(ismatrix){nr=nrow(x); nl=ncol(x)}
  x <- c(x / beta) 
  output <- rep(1, length(x))
  n <- sum(x > 0)
  if(n > 0) {
    x1 <- x[x > 0]
    output[x > 0] <-
      (1/((2^(nu - 1)) * gamma(nu))) * (x1^nu) * besselK(x1, nu)
  }
  if(ismatrix){
    output <- matrix(output, nr, nl)
  }
  return(vars*output)
}

signed.power <- function(x, p){
  return(abs(x)^p * sign(x))
}

gen.matern.data <- function(n.samples, N, nu, beta, exponent=1){
  locs=cbind(rep(0:(N-1), N)/(N-1), rep(0:(N-1), each=N)/(N-1))
  if(beta == 0){
    V=diag(N*N)
  } else{
    V=cov.matern(as.matrix(dist(locs)),nu=nu, beta=beta)  
  }
  X <- array(dim = c(N, N, n.samples))
  x <- rmvn(n = n.samples, mu=rep(0, N^2), sigma=V)

  if(exponent !=1){
    x <- signed.power(x, exponent)
  }
  for(i in 1:n.samples){
    z <- x[i,]
    z <- matrix(z, ncol=N)
    z <- (z -mean(z))/sd(z)
    X[,,i] <- z
  }
  return(X)
}

gen.matern.many <- function(n.samples, N, nu, betas, exponents = 1){
  X <- c()
  for(beta in betas){
    for(exponent in exponents){
      X <- abind(X, gen.matern.data(n.samples, N, nu, beta, exponent))
    }
  }
  return(X)
}

samples.per.setting <- 200 # number of sample points for each setup

### Parameters for training data
betas <- seq(0.0, 0.234, length.out=30) # when nu is equal to 0.5
# betas <- seq(0.0, 0.175, length.out=30) # when nu is equal to 1.0
powers <- c(1.2, 1.4, 1.6, 1.8) # transformation for non-normal data

### Parameters for testing data
betas.test <- seq(0.0, 0.234, length.out=50) # when nu is equal to 0.5
# betas.test <- seq(0.0, 0.175, length.out=50) # when nu is equal to 1.0
powers.test <- seq(from = 1.1, to = 2.0, by = .1) # transformation for non-normal data

train.nonnormal <- gen.matern.many(samples.per.setting, N, nu, betas, powers)
train.normal <- gen.matern.many(samples.per.setting*length(powers), N, nu, betas, 1)

test.nonnormal <- gen.matern.many(samples.per.setting, N, nu, betas.test, powers.test)
test.normal <- gen.matern.many(samples.per.setting*length(powers.test), N, nu, betas.test, 1)

save(N, nu, betas, powers,
     samples.per.setting,
     train.nonnormal, train.normal,
     file = paste("Train-data-beta", length(betas),"-p", length(powers),
                  "-nu", nu, "-gridsize", N, ".Rdata", sep=""))

save(N, nu, betas.test, powers.test,
     samples.per.setting,
     test.nonnormal, test.normal,
     file = paste("Test-data-beta", length(betas.test),"-p", length(powers.test),
                  "-nu", nu, "-gridsize", N, ".Rdata", sep=""))
