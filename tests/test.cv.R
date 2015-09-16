
library(MASS)
library(intervals)
source("../josh/quadratic.R")
source("../josh/cv.R")
source("../josh/groupfs.R")

set.seed(1)
niters <- 500
n <- 40
p <- 20
steps <- 10
sparsity <- 5
snr <- 2
nfolds <- 5

instance <- function(n, p, sparsity, snr, steps, nfolds) {

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
    
    pvals_naive <- interval.groupfs(fit, x, y, index = 1:ncol(x))
    
    fit$projections <- c(fit$projections, fit$rssprojections)
    pvals_reduced <- interval.groupfs(fit, x, y, index = 1:ncol(x))
    
    fit$projections <- c(fit$projections, fit$foldprojections)
    pvals <- interval.groupfs(fit, x, y, index = 1:ncol(x))
    
    return(list(variable = fit$variable, pvals = pvals,
                pvals_naive = pvals_naive, pvals_reduced = pvals_reduced))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, steps, nfolds))
})

pvals_reduced <- do.call(c, list(output[4,]))
pvals_naive <- do.call(c, list(output[3,]))
pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, pvals_reduced, pvals_naive, vars, file = paste0(
                      "results_cv_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      ".RData"))

print(time)
