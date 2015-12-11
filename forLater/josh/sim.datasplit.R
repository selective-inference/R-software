library(intervals)
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 400
known <- FALSE
n <- 100
p <- 100
maxsteps <- 10
sparsity <- 5
snr <- 1
rho <- 0.1
ratio <- 0.7
ratio2 <- 0.85
train <- 1:(ratio*n)
test <- setdiff(1:n, train)
train2 <- 1:(ratio2*n)
test <- setdiff(1:n, train2)
index <- 1:p

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
    
    ytr2 <- y[train2]
    xtr2 <- x[train2, ]
    yte2 <- y[test2]
    xte2 <- x[test2, ]

    if (known) {
        trfit <- groupfs(xtr, ytr, index, maxsteps=maxsteps, sigma=1, aicstop=1, k = 2*log(p))
        fit <- groupfs(xtr2, ytr2, index, maxsteps=maxsteps, sigma=1, aicstop=1, k = 2*log(p))
    } else {
        trfit <- groupfs(xtr, ytr, index, maxsteps=maxsteps, aicstop=1, k = log(length(train)))
        fit <- groupfs(xtr2, ytr2, index, maxsteps=maxsteps, aicstop=1, k = log(length(train2)))
    }

    trcols <- which(1:p %in% trfit$action)
    tr2cols <- which(1:p %in% fit$action)
    tepv <- summary(lm(yte~xte[, trcols]-1))$coefficients[,4]
    tepv2 <- summary(lm(yte2~xte2[, tr2cols]-1))$coefficients[,4]
    names(tepv) <- as.character(sort(trfit$action))
    names(tepv2) <- as.character(sort(trfit$action))
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
     file = paste0("results/datasplit",
         "_", ifelse(known, "TC", "TF"),
         "_n", n,
         "_p", p,
         "_rho", gsub(".", "pt", rho, fixed=T),
         "_sparsity", sparsity,
         "_ratio", gsub(".", "pt", round(ratio, 2), fixed=T),
         "_snr", as.character(snr),
         "_bic.RData"))

print(time)

