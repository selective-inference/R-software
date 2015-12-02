library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 4
n <- 50
p <- 100
maxsteps <- 10
sparsity <- 5
snr <- 2
nfolds <- 5

instance <- function(n, p, sparsity, snr, maxsteps, nfolds) {

    x <- matrix(rnorm(n*p), nrow=n)
    y <- rnorm(n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
      y <- y + x %*% beta
    }

    fit <- cvfs(x, y, maxsteps=maxsteps, sigma=1, nfolds=nfolds)
    fit2 <- groupfs(x, y, index = 1:p, maxsteps = attr(fit, "maxsteps"), sigma = 1)
    pvals <- groupfsInf(fit, verbose=T)
    pv2 <- groupfsInf(fit2, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv, var2 = fit2$action, pvals2 = pv2$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps, nfolds))
})

pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))
vars2 <- do.call(c, list(output[3,]))
pvals2 <- do.call(c, list(output[4,]))

save(pvals, vars, vars2, pvals2, file = paste0(
                      "results_cv_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      "_comparison.RData"))

print(time)
