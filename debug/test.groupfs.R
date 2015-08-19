
library(MASS)
library(intervals)
source("../josh/quadratic.R")
source("../josh/groupfs.R")

set.seed(1)

instance <- function(n, p, sparsity, snr, index, steps) {
    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- y + x %*% beta
    }

    fit <- fstep(x, y, index, steps = steps)

    pvals <- interval.fstep(fit, x, y, index)
    return(list(variable = fit$variable, pvals = pvals))
}

n <- 15
p <- 80
index <- sort(rep(1:(p/2), 2))
steps <- 4
sparsity <- 2
snr <- 20

output <- replicate(50, instance(n, p, sparsity, snr, index, steps))

pvals <- do.call(rbind, output[2,])
vars <- do.call(rbind, output[1,])
colnames(pvals) <- paste0("step", 1:steps)
colnames(vars) <- paste0("step", 1:steps)

save(pvals, vars, file = paste0(
                      "results_size4_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      ".RData"))

