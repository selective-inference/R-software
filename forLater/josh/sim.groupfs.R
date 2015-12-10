library(intervals)
source("funs.sims.R")
#source("selectiveInference/R/cv.R")
source("../../selectiveInference/R/funs.groupfs.R")
source("../../selectiveInference/R/funs.quadratic.R")
source("../../selectiveInference/R/funs.common.R")

set.seed(1)
known <- TRUE
niters <- 300
n <- 50
p <- 150
G <- 75
maxsteps <- 8
sparsity <- 4
snr <- 2
rho <- 0

instance <- function(n, p, G, sparsity, snr, rho, maxsteps) {
    simd <- randomGaussianFixedP(n, p, G, sparsity, snr, sigma = 1, rho)
    x <- simd$x
    y <- simd$y
    index <- simd$index
    if (known) {
        fit <- groupfs(x, y, index, maxsteps, sigma = 1, k = log(n))
    } else {
        fit <- groupfs(x, y, index, maxsteps, k = log(n))
    }
    pvals <- groupfsInf(fit, verbose=T)
    return(list(variable = fit$action, pvals = pvals$pv))
}

time <- system.time({
          output <- replicate(niters, instance(n, p, G, sparsity, snr, rho, maxsteps))
})

pvals <- do.call(c, list(output[2,]))
vars <- do.call(c, list(output[1,]))

save(pvals, vars,
     file = paste0("results/",
         ifelse(known, "TC", "TF"),
         "_n", n,
         "_p", p,
         "_g", G,
         "_rho", gsub(".", "pt", rho, fixed=T),
         "_maxsteps", maxsteps,
         "_sparsity", sparsity,
         "_snr", round(snr),
         ".RData"))

print(time)
