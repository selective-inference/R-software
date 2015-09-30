library(selectiveInference)

source("../selectiveInference/R/cv.R")

n = 100
p = 20
x <- matrix(rnorm(n*p), nrow=n)
y <- rnorm(n)
fit <- cvfs(x, y, maxsteps = 5)
