# Choices

# RSS: least-squares or penalized beta?
# depends on final model. Go with least-squares for now

# fixed vs lar? (lar, apparently)
# fixed probably slower, but advantage of same lambda grid?
# is same lambda grid necessary? -- doesn't lar algorithm give all possible models anyway?
# i.e. for non-knot lambda just find where it is in lar path, take corresponding model

# groups? later

# TODO

# copy larInf or groupfsInf?
# larInf: add CV quadratic constraints* & break/fix p-value computation
# -------- *but can we even use the ydecomp we use for quadratic?
# groupfsInf: some ugly rewriting, no cumprojs etc, but straightforward
# -------- downside: need to implement larInf basically

# larInf
# [ ] is.null(sigma) don't estimate it

# plan:
# expand Gamma for [-fold] indices?
# stack all the Gammas? or iterate through them?
# work backward from poly.pval <- larInf


# big picture / long term
# what OOP kind of design would lend itself to easily implementing more cv things?

# Gamma: something x n
# Gamma %*% y >= 0

# pass 0-padded x[-fold] and y[-fold] to lar?

library(selectiveInference)
setwd("/Users/joftius/Dropbox/work/R-software/forLater/josh")
source("selectiveInference/R/cv.R")

set.seed(1)
n <- 100
p <- 50
maxsteps <- 10
sparsity <- 3
snr <- 2
rho <- 0.1
nfolds <- 5

x <- matrix(rnorm(n*p), nrow=n)
y <- rnorm(n)
beta <- rep(0, p)
beta[1:sparsity] <- 2* sqrt(2*log(p)/n) * sample(c(-1,1), sparsity, replace=T)
y <- y + x %*% beta
my <- mean(y)
y <- y - my

