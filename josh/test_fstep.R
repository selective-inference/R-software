
#library(ggplot2)

library(MASS)
library(intervals)
source("fstep.R")

set.seed(1)

instance <- function(n, p, sparsity, snr, index) {
    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- y + x %*% beta
    }

    fit <- fstep(x, y, index, steps = 4)

    pvals <- interval.fstep(fit, x, y, index)
    return(pvals)
}

n <- 20
p <- 10
index <- sort(rep(1:(p/2),2 ))
sparsity <- 0
snr <- 1

pvals <- replicate(5, instance(n, p, sparsity, snr, index))



##  rob's script
set.seed(30)
n=20
p=10
index <- sort(rep(1:(p/2),2 ))

    x <- matrix(rnorm(n*p), nrow=n)
x=scale(X,T,T)
nsteps=4
pvals=matrix(NA,nsim,nsteps)
nsim=100
for(ii in 1:nsim){
    cat(ii)
    y=rnorm(n)
 fit <- fstep(x, y, index, steps = nsteps)
  pvals[ii,] <- interval.fstep(fit, x, y, index)
}
par(mfrow=c(2,2))

for(k in 1:4)
    plot((1:nsim)/nsim,sort(pvals[,k]))
abline(0,1)
}
