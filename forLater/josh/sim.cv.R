library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 50
n <- 100
p <- 200
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

    fit <- cvfs(x, y, maxsteps=maxsteps, sigma=1, nfolds=nfolds)
    vars <- fit$action
    #fit2 <- groupfs(x, y, index = 1:p, maxsteps = attr(fit, "maxsteps"), sigma = 1)
    pvals <- groupfsInf(fit, verbose=T)
    fit$cvobj <- NULL
    nocvpv <- groupfsInf(fit, verbose=T)
    Y <- y - mean(y)
    cols <- which(1:p %in% vars)
    noselpv <- summary(lm(Y~x[, cols]-1))$coefficients[,4]
    names(noselpv) <- as.character(sort(vars))
    return(list(vars = vars, pvals = pvals$pv,
                nocvvars = vars, nocvpvals = nocvpv$pv,
                noselvars = sort(vars), noselpvals = noselpv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps, nfolds))
})

vars <- do.call(c, list(output[1,]))
pvals <- do.call(c, list(output[2,]))
nocvvars <- do.call(c, list(output[3,]))
nocvpvals <- do.call(c, list(output[4,]))
noselvars <- do.call(c, list(output[5,]))
noselpvals <- do.call(c, list(output[6,]))

save(vars, pvals, nocvvars, nocvpvals, noselvars, noselpvals, file = paste0(
                      "results_cv_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      "_comparison.RData"))

print(time)
