
library(intervals)
source("../josh/quadratic.R")
source("../josh/grouptests.R")
source("../josh/groupfs.R")
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")
set.seed(1)
n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
steps <- 10
sparsity <- 5
snr <- 3


    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)

  
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- y + x %*% beta
  

    fit <- groupfs(x, y, index, maxsteps = steps)
    pvals <- groupfsInf(fit)
   

