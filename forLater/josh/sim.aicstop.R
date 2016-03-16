library(intervals)
source("funs.sims.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
known <- FALSE
niters <- 500
n <- 50
p <- 150
G <- 75
maxsteps <- 10
sparsity <- 4
snr <- 3
rho <- 0
aicstop <- 1

instance <- function(n, p, G, sparsity, snr, rho, maxsteps, aicstop) {
    simd <- randomGaussianFixedP(n, p, G, sparsity, snr, sigma = 1, rho)
    x <- simd$x
    y <- simd$y
    index <- simd$index
    if (known) {
        fit <- groupfs(x, y, index, maxsteps, sigma = 1, k = 2*log(G), aicstop = aicstop, verbose = T)
    } else {
        fit <- groupfs(x, y, index, maxsteps, k = 2*log(G), aicstop = aicstop, verbose = T)
    }
    pvals <- groupfsInf(fit, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv, stopped = attr(fit, "stopped")))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, G, sparsity, snr, rho, maxsteps, aicstop))
})

stopped <- do.call(c, list(output[3,]))
pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, vars, stopped,
     file = paste0(
         "results/aic", 
         "_", ifelse(known, "TC", "TF"),
         "_n", n,
         "_p", p,
         "_g", G,
         "_rho", gsub(".", "pt", rho, fixed=T),
         "_maxsteps", maxsteps,
         "_sparsity", sparsity,
         "_snr", round(snr),
         ".RData"))

print(time)
