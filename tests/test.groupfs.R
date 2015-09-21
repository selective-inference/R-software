
library(intervals)
source("../selectiveInference/R/funs.common.R")
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
source("../selectiveInference/R/funs.fs.R")
source("../selectiveInference/R/funs.lar.R")
#library(selectiveInference)
#library(lars)

set.seed(1)
n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
maxsteps <- 10
sparsity <- 5
snr <- 3

y <- rnorm(n)
x <- matrix(rnorm(n*p), nrow=n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta

system.time({
    fit <- groupfs(x, y, index, maxsteps = maxsteps, sigma = 1)
    pvals <- groupfsInf(fit, sigma = 1)
})

# Compare to step function in R
index <- 1:ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta
df <- data.frame(y = y, x = x)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = 20)
fit <- groupfs(x, y, index, maxsteps)

names(fsfit$coefficients)[-1]
paste0("x.", fit$action)

# They all match


n <- nrow(state.x77)
index <- rep(1:(ncol(state.x77)+1), gsizes)
labels <- unique(index)
maxsteps <- max(labels)-1
sparsity <- 3
snr <- 5
ndiv <- length(levels(state.division))
states <- data.frame(matrix(NA, nrow=n, ncol=ncol(state.x77)))
colnames(states) <- colnames(state.x77)
gsizes <- c(6,8,3,4,2,5,3,4)
x <- matrix(NA, nrow = n, ncol = sum(gsizes) + ndiv)
colnames(x) <- 1:ncol(x)
ind <- 1
for (j in 1:ncol(state.x77)) {
    var <- state.x77[,j]
    qs <- quantile(var, probs = seq(0, 1, length.out = gsizes[j]+1))
    qs <- qs[-length(qs)]
    for (i in 1:n) {
        var[i] <- sum(var[i] >= qs)
    }
    var <- as.factor(var)
    submat <- model.matrix(~ var - 1)
    inds <- ind:(ind+gsizes[j]-1)
    x[, inds] <- submat
    colnames(x)[inds] <- paste(colnames(state.x77)[j], 1:ncol(submat), sep=":")        
    ind <- ind+gsizes[j]
    states[,j] <- var
}
states <- cbind(states, state.division)
submat <- model.matrix(~ state.division - 1)
inds <- ind:(ind+ndiv-1)
x[, inds] <- submat
colnames(x)[inds] <- paste("state.division", 1:ncol(submat), sep=":")
X <- scale_groups(x, index)$x
gsizes <- c(gsizes, ndiv)



p <- ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
nzinds <- which(index %in% sample(labels,sparsity))
beta[nzinds] <- snr 
y <- y + x %*% beta
#y <- y-mean(y)
df <- data.frame(y = y, states)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = 10, k = 2)
fit <- groupfs(x, y, index, maxsteps, k = 2)

names(fsfit$coefficients)
c(colnames(state.x77), "state.division")[fit$action]

