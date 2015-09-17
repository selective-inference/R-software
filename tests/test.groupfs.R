
library(intervals)
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")


set.seed(1)

instance <- function(n, p, sparsity, snr, index, steps) {
    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)

    if (sparsity > 0) {
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- y + x %*% beta
    }

    fit <- groupfs(x, y, index, maxsteps = steps)
    pvals <- groupfsInf(fit, sigma = 1)
    
    return(list(variable = fit$variable, pvals = pvals))
}

n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
steps <- 10
sparsity <- 5
snr <- 3

system.time({
output <- replicate(500, instance(n, p, sparsity, snr, index, steps))
})

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

