library(intervals)
source("funs.sims.R")
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 200
n <- 100
p <- 400
G <- 200
maxsteps <- 5
sparsity <- 3
snr <- 2
rho <- .1

instance <- function(n, p, G, sparsity, snr, rho, maxsteps) {
    simd <- randomGaussianFixedP(n, p, G, sparsity, snr, sigma = 1, rho)
    x <- simd$x
    y <- simd$y
    index <- simd$index
    fit <- groupfs(x, y, index, maxsteps, k = log(n))
    pvals <- groupfsInf(fit, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, G, sparsity, snr, rho, maxsteps))
})

pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, vars, file = paste0(
                      "results_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      "_F.RData"))

print(time)
