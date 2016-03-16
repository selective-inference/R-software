library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 400
known <- FALSE
n <- 100
p <- 50
maxsteps <- 20
sparsity <- 10
snr <- 1
rho <- 0.1
ratio <- 0.75
train <- 1:(ratio*n)
test <- setdiff(1:n, train)
index <- 1:p
nfolds <- 5

instance <- function(n, p, sparsity, snr, maxsteps, rho) {

    x <- matrix(rnorm(n*p), nrow=n)
    if (rho != 0) {
        z <- matrix(rep(t(rnorm(n)), p), nrow = n)
        x <- sqrt(1-rho)*x + sqrt(rho)*z
    }
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

    if (known) {
        trfit <- cvfs(xtr, ytr, maxsteps=maxsteps, sigma = 1, nfolds=nfolds)
        fit <- cvfs(x, y, maxsteps=maxsteps, sigma = 1, nfolds=nfolds)
    } else {
        trfit <- cvfs(xtr, ytr, maxsteps=maxsteps, nfolds=nfolds)
        fit <- cvfs(x, y, maxsteps=maxsteps, nfolds=nfolds)
    }

    trcols <- which(1:p %in% trfit$action)
    tepv <- summary(lm(yte~xte[, trcols]-1))$coefficients[,4]
    names(tepv) <- as.character(sort(trfit$action))
    pv <- groupfsInf(fit)
    trpv <- groupfsInf(trfit)
    return(list(vars = fit$action, pvals = pv$pv,
                splitvars = sort(trfit$action), splitpvals = tepv,
                trpvals = trpv$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, sparsity, snr, maxsteps, rho))
})

vars <- do.call(c, list(output[1,]))
pvals <- do.call(c, list(output[2,]))
splitvars <- do.call(c, list(output[3,]))
splitpvals <- do.call(c, list(output[4,]))
trpvals <- do.call(c, list(output[5,]))

save(vars, pvals, splitvars, splitpvals, trpvals,
     file = paste0("results/carvecv",
         "_", ifelse(known, "TC", "TF"),
         "_n", n,
         "_p", p,
         "_rho", gsub(".", "pt", rho, fixed=T),
         "_sparsity", sparsity,
         "_ratio", gsub(".", "pt", round(ratio, 2), fixed=T),
         "_snr", as.character(snr),
         ".RData"))

print(time)

