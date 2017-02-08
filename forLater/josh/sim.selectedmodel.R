library(selectiveInference)
library(intervals)
setwd("~/Dropbox/work/R-software/forLater/josh")
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")
source("../../selectiveInference/R/funs.fs.R")
source("../../selectiveInference/R/funs.lar.R")
source("../../selectiveInference/R/funs.inf.R")
library(MASS)
pinv = ginv

set.seed(19)
niters <- 500
known <- TRUE
n <- 50
p <- 100
maxsteps <- 8
sparsity <- 5
snr <- 2
index <- 1:p

x <- matrix(rnorm(n*p), nrow=n)

instance <- function(n, p, sparsity, snr, maxsteps) {
    y <- rnorm(n)
    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
      y <- y + x %*% beta
    }
    y <- y - mean(y)
    fit <- groupfs(x, y, index, maxsteps=maxsteps, sigma=1, intercept=F, center=F, normalize=F)
    fitfs <- fs(x, y, maxsteps=maxsteps, intercept=F, normalize=F)
    if (any(fit$action != fitfs$action)) stop("Model paths did not agree")
    pvfs <- fsInf(fitfs, sigma=1, k = maxsteps, type = "all")
    pv <- groupfsInf(fit)
    return(list(vars = fit$action, pvals = pv$pv, selpvals = pvfs$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps))
})

vars <- do.call(c, list(output[1,]))
pvals <- do.call(c, list(output[2,]))
selpvals <- do.call(c, list(output[3,]))

save(vars, pvals, selpvals,
     file = paste0("results/selected",
         "_", ifelse(known, "TC", "TF"),
         "_n", n,
         "_p", p,
         "_sparsity", sparsity,
         "_snr", as.character(snr),
         ".RData"))

print(time)

