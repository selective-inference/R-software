library(intervals)
source("../selectiveInference/R/cv.R")
source("../../../selectiveInference/R/funs.groupfs.R")
source("../../../selectiveInference/R/funs.quadratic.R")
source("../../../selectiveInference/R/funs.common.R")

set.seed(1)
n <- 50
p <- 100
maxsteps <- 10
sparsity <- 5
snr <- 1
nfolds <- 5
x <- matrix(rnorm(n*p), nrow=n)
y <- rnorm(n)
beta <- rep(0, p)
beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
y <- y + x %*% beta
fit <- cvfs(x, y, maxsteps=maxsteps, nfolds=nfolds)
pvals <- groupfsInf(fit, verbose=T)

