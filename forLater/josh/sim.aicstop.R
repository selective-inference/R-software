library(intervals)
source("funs.sims.R")
source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
niters <- 400
n <- 100
p <- 100
G <- 50
maxsteps <- 15
sparsity <- 5
snr <- 1
sigma <- 1
rho <- 0
aicstop <- 1

instance <- function(n, p, G, sparsity, snr, rho, maxsteps, aicstop) {

    simd <- randomGaussianFixedP(n, p, G, sparsity, snr, sigma, rho)
    x <- simd$x
    y <- simd$y
    index <- simd$index

    fit <- groupfs(x, y, index, maxsteps, k = log(n), aicstop = aicstop)
    pvals <- groupfsInf(fit, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv, stopped = attr(fit, "stopped")))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, G, sparsity, snr, rho, maxsteps, aicstop))
})

stopped <- do.call(c, list(output[3,]))
pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, vars, stopped, file = paste0(
                      "results_aic", aicstop, "_n", n,
                      "_p", p,
                      "_sparsity", sparsity,
                      "_snr", snr,
                      "_F_rho0.RData"))

print(time)
