#library(selectiveInference)
#library(lars)
library(intervals)
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
source("../selectiveInference/R/funs.common.R")

set.seed(1)
n <- 100
G <- 10
maxsteps <- 10
snr <- 1
niter <- 100

print("Comparing step with groupfs on random categorical designs")
aicdiffs <- numeric(niter)
mismatchcount <- 0
for (iter in 1:niter) {
    rles <- 2 + rpois(G, 2)
    df <- data.frame(do.call(cbind, lapply(rles, function(g) {
        sample(LETTERS[1:g], n, replace = TRUE, prob = runif(g))
    })), stringsAsFactors = TRUE)
    if (any(apply(df, 2, function(col) length(unique(col))) == 1)) next
    fd <- factorDesign(df)
    if (any(duplicated(fd$x, MARGIN = 2))) next
    y <- rnorm(n)
    x1 <- fd$x[, fd$index == 1]
    y <- y + x1 %*% c(snr, rep(0, ncol(x1) - 2), -snr)
    y <- y - mean(y)
    df$y <- y
    capture.output(fsfit <- step(lm(y ~ 0, df), direction="forward", scope = formula(lm(y~.-1, df)), steps = maxsteps, trace = 1000), file = "/dev/null")
    fit <- groupfs(fd$x, df$y, fd$index, maxsteps = 10, intercept = F, center = F, normalize = F, aicstop = 1)
    fsnames <- names(fsfit$coefficients)
    if (length(fsnames) > 0) {
        fsnames <- unique(substr(fsnames, 1, nchar(fsnames) - 1))
        k <- length(fsnames)
        fitnames <- attr(fit, "varnames")[fit$action][1:(length(fit$action)-attr(fit, "aicstop"))]
        aicdiffs[iter] <- AIC(fsfit) - fit$log$AIC[k]
        if (length(fitnames) !=k || any(fsnames != fitnames)) {
            print(paste("Mismatch at iteration", iter))
            print(fsnames)
            print(fitnames)
            mismatchcount <- mismatchcount + 1
        }
    }
}
print(paste("Total mismatches:", mismatchcount, "out of", niter))
print("Summary of differences in AIC")
summary(aicdiffs)
