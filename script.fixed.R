library(glmnet)

 library(selectiveInference,lib.loc="mylib")
library(truncnorm)
library(MASS)


set.seed(133)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(-20,3,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=glmnet(x,y,standardize=F)
lambda = .1
bhat = coef(a, s=lambda/n)[-1]


# compute fixed lambda p-values

a4=fixedLassoInf(x,y,bhat,lambda,sigma,compute.si=T)


## check of numerics for one case
 fun = function(x,etay,vm,vp,sigma.eta) return(1-ptruncnorm(etay,vm,vp,x,sigma.eta))
xx=seq(-62,63,length=50000)
vm=18.98
etay=19.81
vp=23.37
sigma.eta=1.4
u=0

 1-ptruncnorm(etay, vm, vp, u, sigma.eta)

w=fun(xx,etay,vm,vp,sigma.eta) 


plot(xx,w)
