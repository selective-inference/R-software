
library(MASS)
library(intervals)
source("../josh/quadratic.R")
source("../josh/cv.R")
source("../josh/groupfs.R")

set.seed(1)
niters <- 100
n <- 50
p <- 100
steps <- 10
sparsity <- 5
snr <- 3

instance <- function(n, p, sparsity, snr, index, steps, nfolds) {

    x <- matrix(rnorm(n*p), nrow=n)
    y <- rnorm(n)
    
    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
      y <- y + x %*% beta
    }

    fit <- cv_fs(x, y, steps=steps, nfolds=nfolds)    
    x <- x[fit$cvperm, ]
    y <- y[fit$cvperm]
    pvals <- interval.fstep(fit, x, y, index = 1:ncol(x))

    return(list(variable = fit$variable, pvals = pvals))
}

output <- replicate(niter, instance(n, p, sparsity, snr, index, steps))

pvals <- do.call(c, output[2,])
vars <- do.call(c, output[1,])

save(pvals, vars, file = paste0(
                      "results_cv_full_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      ".RData"))
