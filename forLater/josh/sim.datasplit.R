library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 1000
n <- 100
p <- 50
maxsteps <- 20
sparsity <- 5
snr <- 1
train <- 1:(7*n/10)
test <- setdiff(1:n, train)
index <- 1:p

instance <- function(n, p, sparsity, snr, maxsteps, nfolds) {

    x <- matrix(rnorm(n*p), nrow=n)
    y <- rnorm(n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[1:sparsity] <- snr * sample(c(-1,1), sparsity, replace=T)
      y <- y + x %*% beta
    }

    ytr <- y[train]
    xtr <- x[train, ]
    yte <- y[test]
    xte <- x[test, ]

    #log(length(train))
    trfit <- groupfs(xtr, ytr, index, maxsteps=maxsteps, aicstop=1, k = log(length(train)))
    fit <- groupfs(x, y, index, maxsteps=maxsteps, aicstop=1, k = log(n))

    trcols <- which(1:p %in% trfit$action)
    tepv <- summary(lm(yte~xte[, trcols]-1))$coefficients[,4]
    names(tepv) <- as.character(sort(trfit$action))
    pv <- groupfsInf(fit)
    return(list(vars = fit$action, pvals = pv$pv,
                splitvars = sort(trfit$action), splitpvals = tepv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps, nfolds))
})

vars <- do.call(c, list(output[1,]))
pvals <- do.call(c, list(output[2,]))
splitvars <- do.call(c, list(output[3,]))
splitpvals <- do.call(c, list(output[4,]))


save(vars, pvals, splitvars, splitpvals, file = paste0(
                      "results_datasplit_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", as.character(snr),
                      "_bic.RData"))

print(time)

