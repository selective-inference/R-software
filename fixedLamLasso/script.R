library(glmnet)


source("newfuns.R")
source("ryanjon.R")

set.seed(33)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)

#generate y
beta=c(3,3,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)



a=glmnet(x,y,standardize=F)
lam = .5
bhat = coef(a, s=lam)[-1]

# compute fixed lambda p-values

a4=mylasso.pv(x,y,bhat,lam,sigma)

