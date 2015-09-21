
library(intervals)
source("../selectiveInference/R/funs.common.R")
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
source("../selectiveInference/R/funs.fs.R")
source("../selectiveInference/R/funs.lar.R")
library(selectiveInference)
library(lars)

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


# Compare to step function in R
index <- 1:ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta
df <- data.frame(y = y, x = x)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = 20)
fit <- groupfs(x, y, index, maxsteps)

names(fsfit$coefficients)[-1]
paste0("x.", fit$action)

junk=fs(x,y)
# They all match
junk2=lars(x,y,type="step")
cbind(names(fsfit$coefficients)[-1],unlist(junk2$act)[1:20])


