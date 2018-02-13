library(selectiveInference)
library(glmnet)
set.seed(32)

n=6000; p=50000

coverage = c()
for (i in 1:100) {
print("creating X")
print(system.time(x <- matrix(rnorm(n*p),n,p)))
print(system.time(x <- scale(x,T,T)/sqrt(n-1)))
print('running glmnet')
y=rnorm(n)
print(system.time(G <- glmnet(x,y,intercept=F,standardize=F)))
lam=G$lam[5] #for p< 1000
lam=G$lam[10] #for p=1000
lam=G$lam[2]  #for p=50,000
print("solving randomized LASSO")
print(system.time(rand_lasso <- randomizedLasso(x, y, n*lam)))
print(rand_lasso$active_set)
print("inference for randomizedLasso")
tim=system.time(rand_lasso_inf<-randomizedLassoInf(rand_lasso, sampler='norejection', nsample=20000, burnin=10000))[1]
print(rand_lasso_inf$pvalues)
print(names(rand_lasso_inf$pvalues))
print(rand_lasso_inf$ci)
coverage = c(coverage, (rand_lasso_inf$ci[,1] <= 0) * (rand_lasso_inf$ci[,2] >= 0))
print(mean(coverage))
}