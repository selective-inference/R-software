#library(selectiveInference)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

set.seed(12345)

n = 100 # sample size
signal = 3 # signal size
mu = c(rep(signal, floor (n/5)), rep(0, n-floor(n/5))) # 20% of elements get the signal; rest 0
mu=sample(mu)
y = mu + rnorm (n, 0, 1)

mmObj = manyMeans(y, k=10, verbose=T)

