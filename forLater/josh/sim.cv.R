library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 500
n <- 50
p <- 100
maxsteps <- 10
sparsity <- 5
snr <- 1
nfolds <- 5

instance <- function(n, p, sparsity, snr, maxsteps, nfolds) {

    x <- matrix(rnorm(n*p), nrow=n)
    y <- rnorm(n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
      y <- y + x %*% beta
    }

    fit <- cvfs(x, y, maxsteps=maxsteps, nfolds=nfolds)
    pvals <- groupfsInf(fit, sigma = 1, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps, nfolds))
})

pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, vars, file = paste0(
                      "results_cv_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      ".RData"))

print(time)
