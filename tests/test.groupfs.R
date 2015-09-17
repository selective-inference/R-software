
library(intervals)
source("../selectiveInference/R/funs.common.R")
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
source("../selectiveInference/R/funs.fs.R")
source("../selectiveInference/R/funs.lar.R")
#library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

set.seed(1)
n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
maxsteps <- 10
sparsity <- 5
snr <- 3

y <- rnorm(n)
x <- matrix(rnorm(n*p), nrow=n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta

system.time({
    fit <- groupfs(x, y, index, maxsteps = maxsteps)
    pvals <- groupfsInf(fit)
})

# Compare to step function in R
index <- 1:ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta
df <- data.frame(y = y, x = x)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = 10)
fit <- groupfs(x, y, index, maxsteps)

names(fsfit$coefficients)[-1]
paste0("x.", fit$action)
# They all match
