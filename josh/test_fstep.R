
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
