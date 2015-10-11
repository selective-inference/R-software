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

check.mismatch <- function(fsfit, fit) {
    fsnames <- names(fsfit$coefficients)
    if (length(fsnames) > 0) {
        fsnames <- unique(substr(fsnames, 1, nchar(fsnames) - 1))
        k <- length(fsnames)
        fitnames <- attr(fit, "varnames")[fit$action][1:(length(fit$action)-attr(fit, "aicstop"))]
        if (is.null(fit$sigma)) {
            aicdiff <- AIC(fsfit) - fit$log$AIC[k]
        } else {
            aicdiff <- extractAIC(fsfit, scale = fit$sigma)[2] - fit$log$AIC[k]
        }
        if (length(fitnames) !=k || any(fsnames != fitnames)) {
            print(paste("Mismatch at iteration", iter, ifelse(is.null(fit$sigma), "unknown", "known")))
            print(fsnames)
            print(fitnames)
            return(list(count = 1, aicdiff = aicdiff))
        }
        return(list(count = 0, aicdiff = aicdiff))
    }
    return(list(count = 0, aicdiff = 0))
}

print("Comparing step with groupfs on random categorical designs")
umismatchcount <- kmismatchcount <- 0
uaicdiffs <- kaicdiffs <- numeric(niter)
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
    capture.output(fsfit <- step(lm(y ~ 0, df), direction="forward", scope = formula(lm(y~.-1, df)), steps = maxsteps), file = "/dev/null")
    fit <- groupfs(fd$x, df$y, fd$index, maxsteps = 10, intercept = F, center = F, normalize = F, aicstop = 1)
    mismatches <- check.mismatch(fsfit, fit)
    umismatchcount <- umismatchcount + mismatches$count
    uaicdiffs[iter] <- mismatches$aicdiff
    capture.output(fsfit <- step(lm(y ~ 0, df), scale = 1, direction="forward", scope = formula(lm(y~.-1, df)), steps = maxsteps), file = "/dev/null")
    fit <- groupfs(fd$x, df$y, fd$index, maxsteps = 10, sigma = 1, intercept = F, center = F, normalize = F, aicstop = 1)
    mismatches <- check.mismatch(fsfit, fit)
    kmismatchcount <- kmismatchcount + mismatches$count
    kaicdiffs[iter] <- mismatches$aicdiff
}
print(paste("Mismatches:", umismatchcount, "for unknown sigma and", kmismatchcount, "for known"))
summary(uaicdiffs)
summary(kaicdiffs)
