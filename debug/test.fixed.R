#library(selectiveInference)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

#options(error=dump.frames)
attach("/Users/tibs/dropbox/PAPERS/lasso/lasso3/.RData")

set.seed(133)
n=100
p=200
sigma=1

x = matrix(rnorm(n*p),n,p)
x=scale(x,T,F)/sqrt(n-1)
beta=c(-6,3,2,-1,rep(0,p-4))
y=x%*%beta+sigma*rnorm(n)

### RJT COMMENTS: this seems to work with every combination of
### intercept and standardize, EXCEPT intercept=F and standardize=T.
### In this case, glmnet simply refuses to fit very many lambdas
### along the path, and so lambda=10 is way to small for it
a = glmnet(x,y,intercept=F,standardize=F,lambda.min.ratio=1e-6,thresh=1e-10)
lambda = 1
bhat = (coef(a, s=lambda/n, exact=TRUE))[-1]
out= fixedLassoInf(x,y,bhat,lambda,sigma=sigma)
out


critf=function(b,lam,x,y){
     yhat=x%*%b
     .5*sum( (y-yhat)^2) + lam*sum(abs(b))
   }

##check coverage
set.seed(3)

n=15
p=10
sigma=2.3

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

beta=c(rep(2,5),rep(0,p-5))
     lambda = 2
nsim=100
seeds=sample(1:9999,size=nsim)

ci=matrix(NA,nsim,2)
btrue=rep(NA,nsim)
for(ii in 1:nsim){
    cat(ii)
    set.seed(seeds[ii])
    mu=x%*%beta
y=mu+sigma*rnorm(n)

y=y-mean(y)
    

    
gfit=glmnet(x,y,standardize=F,lambda.min.ratio=1e-9)

     bhat = predict(gfit, s=lambda/n,type="coef",exact=T)[-1]
bhat2=lasso2lam(x,y,lambda,int=F,stand=F)
 junk= lassoInf(x,y,bhat,lambda,sigma=sigma)
ci[ii,]=junk$ci[3,]
    xx=x[,junk$pred]
    btrue[ii]=(solve(t(xx)%*%xx)%*%t(xx)%*%mu)[3]
}

sum(ci[,1]< btrue & ci[,]> btrue)/nsim


## BIG example
 library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

options(error=dump.frames)
attach("/Users/tibs/dropbox/PAPERS/lasso/lasso3/.RData")
critf=function(b,lam,x,y){
     yhat=x%*%b
     .5*sum( (y-yhat)^2) + lam*sum(abs(b))
 }
set.seed(4)
n=100
p=500
sigma=1
x=matrix(rnorm(n*p),ncol=p)
x=scale(x,T,F)
beta=c(rep(2.5,10),rep(0,p-10))
y=x%*%beta+sigma*rnorm(n)
y=y-mean(y)
    

    
gfit=glmnet(x,y,standardize=F)
cvf=cv.glmnet(x,y)

lambda=n*cvf$lambda.min
#lambda=10
     bhat = as.numeric(predict(gfit, s=lambda/n,type="coef",exact=T))[-1]

bhat2=lasso2lam(x,y,lambda,int=F,stand=F)$coef
plot(bhat,bhat2)

critf(bhat,lambda,x,y)
critf(bhat2,lambda,x,y)
 junk= fixedLassoInf(x,y,bhat,lambda,sigma=sigma)


 # check of KKT
ch=function(bhat,tol.beta=1e-5,tol.kkt=0.1){
    xx=cbind(1,x)
    bhatt=c(0,bhat)
   g0=t(xx)%*%(y-xx%*%bhatt)
   g=g0-lambda*sign(bhatt)
    gg=g0/lambda
    oo=abs(bhatt)>tol.beta
    cat(c(max(abs(g[oo]))>tol.kkt,min(gg[!oo])< -1-tol.kkt,max(gg[!oo])>1 +tol.kkt),fill=T)
}


#
x=read.table("/Users/tibs/dropbox/PAPERS/FourOfUs/data64.txt")
x=as.matrix(x)
x=scale(x,T,T)
n=length(y)
nams=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/data64.names",what="")
y=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/diab.y")
y=y-mean(y)

cvf=cv.glmnet(x,y)
sigmahat=estimateSigma(x,y,stand=F)$si

lambda=n*cvf$lambda.min

gfit=glmnet(x,y,standardize=F)
 
bhat=coef(gfit, s=lambda/n, exact=TRUE)[-1]
bhat2=lasso2lam(x,y,lambda,int=F,stand=F)$coef

critf(bhat,lambda,x,y)
critf(bhat2,lambda,x,y)


fixedLassoInf(x,y,bhat,lambda,sigma=sigmahat)
##
set.seed(44)
n=50
p=10
sigma=.7
x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)
beta=c(3,2,0,0,rep(0,p-4))
y=x%*%beta+sigma*rnorm(n)
y=y-mean(y)
# first run  glmnet
gfit=glmnet(x,y,standardize=F)
lambda = .1
#extract coef for a given lambda; Note the 1/n factor!
beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
fixedLassoInf(x,y,beta,lambda,sigma=sigma)

# as above, but use lar function to get initial lasso fit
 fit=lar(x,y,normalize=F)
beta=coef(fit,s=lambda,mode="lambda")
fixedLassoInf(x,y,beta,lambda,sigma=sigma)

###
x=state.x77[,-4]
y=state.x77[,4]
x=scale(x,T,T)
n=nrow(x)

cvf=cv.glmnet(x,y)
sigmahat=estimateSigma(x,y,stand=F)$si

lambda=n*cvf$lambda.min

gfit=glmnet(x,y,standardize=F)
 
bhat=coef(gfit, s=lambda/n, exact=TRUE)[-1]

fixedLassoInf(x,y,bhat,lambda,sigma=sigmahat)
