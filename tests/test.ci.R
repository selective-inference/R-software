library(selectiveInference)
#library(selectiveInference,lib.loc="mylib")

#
##check coverage
set.seed(3)

n=50
p=10
sigma=2

x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)    #try with and without standardization

beta=c(5,4,3,2,1,rep(0,p-5))

nsim=100
seeds=sample(1:9999,size=nsim)
pv=rep(NA,nsim)
ci=matrix(NA,nsim,2)
btrue=rep(NA,nsim)
lambda=20

for(ii in 1:nsim){
    cat(ii)
    set.seed(seeds[ii])
    mu=x%*%beta
   y=mu+sigma*rnorm(n)
    y=y-mean(y)  
   gfit=glmnet(x,y,standardize=F,lambda.min.ratio=1e-9)
#   ilam=trunc(length(gfit$lam)/4)
#    lambda=gfit$lam[ilam]*n
     bhat = predict(gfit, s=lambda/n,type="coef",exact=TRUE)[-1]

    junk= fixedLassoInf(x,y,bhat,lambda,sigma=sigma)
    pv[ii]=junk$pv[1]
   ## oo=junk$pred # for old package
    oo=junk$var   # for new package
     btrue[ii]=lsfit(x[,oo],mu,intercept=F)$coef[1]
     ci[ii,]=junk$ci[1,]
}

mean(ci[,1]> btrue)
mean(ci[,2]< btrue)
