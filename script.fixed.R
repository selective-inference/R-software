library(glmnet)

 library(selectiveInference,lib.loc="mylib")
library(truncnorm)
library(MASS)


set.seed(33)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(3,3,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=glmnet(x,y,standardize=F)
lam = .1
bhat = coef(a, s=lam/n)[-1]


# compute fixed lambda p-values

a4=fixedLassoInf(x,y,bhat,lam,sigma,compute.ci=T)

