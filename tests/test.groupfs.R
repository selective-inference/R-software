#library(selectiveInference)
#library(lars)
library(intervals)
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
source("../selectiveInference/R/funs.common.R")

set.seed(1)
n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
maxsteps <- 8
sparsity <- 4
snr <- 2

system.time({
for (iter in 1:100) {
    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)
    beta <- rep(0, p)
    beta[which(index %in% 1:sparsity)] <- snr
    y <- y + x %*% beta
    fit <- groupfs(x, y, index, maxsteps = maxsteps)
    pvals <- groupfsInf(fit, verbose = T)
}
})

# Compare to step function in R
index <- 1:ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta
df <- data.frame(y = y, x = x)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = maxsteps)
fit <- groupfs(x, y, index, maxsteps)

names(fsfit$coefficients)[-1]
paste0("x.", fit$action)

# They all match

n <- nrow(state.x77)
ndiv <- length(levels(state.division))
gsizes <- c(6,8,3,4,2,5,3,4, ndiv)
cnames <- c(colnames(state.x77), "state.division")
cnames <- gsub(" ", ".", cnames)
index <- rep(1:(ncol(state.x77)+1), gsizes)
labels <- unique(index)
maxsteps <- max(labels)-1
sparsity <- 3
snr <- 5
states <- data.frame(matrix(NA, nrow=n, ncol=ncol(state.x77)))
colnames(states) <- colnames(state.x77)
for (j in 1:ncol(state.x77)) {
    var <- state.x77[,j]
    qs <- quantile(var, probs = seq(0, 1, length.out = gsizes[j]+1))
    qs <- qs[-length(qs)]
    for (i in 1:n) {
        var[i] <- sum(var[i] >= qs)
    }
    var <- as.factor(var)
    states[,j] <- var
}
states <- cbind(states, state.division)
x <- factorDesign(states)$x
X <- scaleGroups(x, index)$x

p <- ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
nz <- sample(labels,sparsity)
nzinds <- which(index %in% nz)
beta[nzinds] <- snr
y <- y + x %*% beta
y <- y-mean(y)
df <- data.frame(y = y, states)
fsfit <- step(lm(y ~ 0, df), direction="forward", scope = formula(lm(y~., df)), steps = maxsteps, k = 2)
fit <- groupfs(x, y, index, maxsteps, k = 2, intercept = F, center = F, normalize = T)
# names(fsfit$coefficients)[-1]
if (length(fsfit$coefficients) > 0) {
    fsnames <- cnames[which(!is.na(charmatch(cnames,names(fsfit$coefficients)[-1])))][order(unlist(lapply(cnames, function(cn) {
    matches = grep(cn, names(fsfit$coefficients)[-1])
    if (length(matches) > 0) min(matches)
    else NULL
})))]
fsnames
cnames[fit$action]#[1:length(fsnames)]
} else {
    print("empty")
}

n = 100
p = 120
maxsteps = 9
niter = 500
# 10 groups of size 10, 10 groups of size 2
index = sort(c(c(1, 1), rep(2:11, 10), rep(12:20, 2)))
pvalm = pvalmk = matrix(NA, nrow=niter, ncol=maxsteps)

for (iter in 1:niter) {
    x = matrix(rnorm(n*p), nrow=n)
    y = rnorm(n)
    fit = groupfs(x, y, index, maxsteps)
    pvals = groupfsInf(fit)
    pvalm[iter, ] = pvals$pv
    fitk = groupfs(x, y, index, maxsteps, sigma = 1)
    pvalsk = groupfsInf(fitk)
    pvalmk[iter, ] = pvalsk$pv
    cat(paste("Iteration", iter, "\n"))
}

library(ggplot2)
library(reshape2)
df <- melt(pvalm)[,2:3]
colnames(df) <- c("step", "Pvalue")
df$step <- as.factor(df$step)
p <- ggplot(df, aes(Pvalue, colour = step, linetype = step)) + stat_ecdf()
ggsave(filename = "test_unknown.pdf", plot = p)
df <- melt(pvalmk)[,2:3]
colnames(df) <- c("step", "Pvalue")
df$step <- as.factor(df$step)
p <- ggplot(df, aes(Pvalue, colour = step, linetype = step)) + stat_ecdf()
ggsave(filename = "test_known.pdf", plot = p)

print(colMeans(pvalm))
print(colMeans(pvalmk))

print(mean(pvalm))
print(mean(pvalmk))

