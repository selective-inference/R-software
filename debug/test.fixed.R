#library(selectiveInference)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

#options(error=dump.frames)
#attach("/Users/tibs/dropbox/PAPERS/lasso/lasso3/.RData")

set.seed(133)
n=100
p=10
sigma=3

x = matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
beta=c(-6,3,2,-1,rep(0,p-4))
y=x%*%beta+sigma*rnorm(n)

### RJT COMMENTS: this seems to work with every combination of
### intercept and standardize, EXCEPT intercept=F and standardize=T.
### In this case, glmnet simply refuses to fit very many lambdas
### along the path, and so lambda=10 is way to small for it
a = glmnet(x,y,intercept=F,standardize=F,lambda.min.ratio=1e-3,thresh=1e-10)
nlam2=trunc(length(a$lam)/2)
lambda = n*(a$lam[nlam2])
bhat = (coef(a, s=lambda/n, exact=TRUE))[-1]
out= fixedLassoInf(x,y,bhat,lambda,sigma=sigma,intercept=F)
out

##
a=fs(x,y)
aa=fsInf(a)

a=lar(x,y)
aa=larInf(a)

critf=function(b,lam,x,y){
     yhat=x%*%b
     .5*sum( (y-yhat)^2) + lam*sum(abs(b))
   }

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
for(ii in 1:nsim){
    cat(ii)
    set.seed(seeds[ii])
    mu=x%*%beta
   y=mu+sigma*rnorm(n)
    y=y-mean(y)  
   gfit=glmnet(x,y,standardize=F,lambda.min.ratio=1e-9)
   ilam=trunc(length(gfit$lam)/4)
    lambda=gfit$lam[ilam]*n
     bhat = predict(gfit, s=lambda/n,type="coef",exact=F)[-1]

     junk= fixedLassoInf(x,y,bhat,lambda,sigma=sigma)
    pv[ii]=junk$pv[1]
   # oo=junk$pred # for old package
    oo=junk$var   # for new package
     btrue[ii]=lsfit(x[,oo],mu)$coef[2]
     ci[ii,]=junk$ci[1,]
}

sum(ci[,1]> btrue)
sum(ci[,2]< btrue)




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


# lucas example
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

set.seed(44)

p <- 300
n <- 200
s0 <- 2
b <- 1
b0 <- 0
sigma <- .5
alpha = 0.05

X = matrix(rnorm(n*p),n,p)
X = scale(X,center=T,scale=T)

m = 1000
eps = matrix(rnorm(m*n),n,m)
lam = 2*mean(apply(t(X)%*%eps,2,max))


theta0 <- c(rep(b,s0),rep(0,p-s0));
w <- sigma*rnorm(n);
y <- (b0+X%*%theta0+w);



tic = proc.time()
gfit = glmnet(X,y,standardize=F)
coef = coef(gfit, s=lam/n, exact=T)[-1]
sint = fixedLassoInf(X,y,coef,lam,sigma=sigma,alpha=alpha)

### lucas example with sims

set.seed(44)

p <- 300
n <- 200
s0 <- 2
b <- 1
b0 <- 0
sigma <- .5
alpha = 0.05
#set.seed('1')

#X <- rbinom(p*n,1,prob=0.15);
#dim(X) <- c(n,p);
#X <- X %*% diag(1+9*runif(p))
X = matrix(rnorm(n*p),n,p)
X = scale(X,center=T,scale=T)  # original
#X = scale(X,center=T,scale=T)/sqrt(n-1)   #CHANGED

m = 1000
eps = matrix(rnorm(m*n),n,m)

lam = 2*mean(apply(t(X)%*%eps,2,max))   #original
theta0 <- c(rep(b,s0),rep(0,p-s0))  #original
#theta0 <- c(rep(b,s0),rep(0,p-s0))*sqrt(n-1)  #CHANGED
  mu=b0+X%*%theta0
nsim=100
int=matrix(NA,nsim,2)
btrue=rep(NA,nsim)
for(ii in 1:nsim){
    cat(ii)
w <- sigma*rnorm(n);
       
y <- (mu+w);

tic = proc.time()
gfit = glmnet(X,y,standardize=F)
    nz=colSums(gfit$beta!=0)
 #   lam=gfit$lambda[nz>=2]/1.1 # CHANGED
 #   lam=lam[1]*n               #CHANGED
coef = coef(gfit, s=lam/n, exact=T)[-1]
oo=which(coef!=0)
btrue[ii]=lsfit(X[,oo],mu)$coef[2]
sint = fixedLassoInf(X,y,coef,lam,sigma=sigma,alpha=alpha)
int[ii,]=sint$ci[1,]
}


areaf=function(tt,mean,sigma.eta,a,b,nd=20){
 #compute Prob_mean (W<tt| W in [a,b]
  val=NA
  if(tt>=a & tt <=b){
#      val0=pnorm(mpfr((tt-mean)/sigma.eta,nd),log.p=TRUE)
 #     val1=pnorm(mpfr((a-mean)/sigma.eta,nd),log.p=TRUE)
 #     val2=pnorm(mpfr((b-mean)/sigma.eta,nd),log.p=TRUE)

#  val=(1-mpfr(exp(val2-val0),nd))/(mpfr(exp(val1-val0),nd)-mpfr(exp(val2-val0),nd))
      val0=pnorm(mpfr((tt-mean)/sigma.eta,nd),log.p=F)
      val1=pnorm(mpfr((a-mean)/sigma.eta,nd),log.p=F)
      val2=pnorm(mpfr((b-mean)/sigma.eta,nd),log.p=F)
 if(val1!=val2)  val=(val0-val2)/(val1-val2)
}
return(val)
}

tnorm.surv <- function(z, mean, sd, a, b) {
  z = max(min(z,b),a)

  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1

  o = is.finite(mean)
  p[o] = bryc.tnorm.surv(z,mean[o],sd,a,b)
  #p[o] = gsell.tnorm.surv(z,mean[o],sd,a,b)
  return(p)
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)

  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)

  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}
ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
         (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

areaf2=function(x,a,b,mean=0,sigma=1,nd=500){
   x = max(x, a)
    x = min(x, b)
    if (a > 0 & b > 0){
        Fx= pnorm(mpfr(pnorm(-x,mean=mean,sd=sigma),nd)); Fa=pnorm(mpfr(pnorm(-a,mean=mean,sd=sigma),nd)); Fb= pnorm(mpfr(pnorm(-b,mean=mean,sd=sigma),nd))
        return (1- ( Fa - Fx ) / ( Fa - Fb ) )
    }
    else{
        Fx= pnorm(mpfr(pnorm(x,mean=mean,sd=sigma),nd));Fa=pnorm(mpfr(pnorm(a,mean=mean,sd=sigma),nd));Fb= pnorm(mpfr(pnorm(b,mean=mean,sd=sigma),nd))
        return ( 1-( Fx - Fa ) / ( Fb - Fa ) )
      }
    }

library(Rmpfr)
a=10
b=14
x=12
mean=0
sigma=1.4
tnorm.surv(x,mean,sigma,a,b)
areaf(x,mean,sigma,a,b)
areaf2(x,a,b,mean=mean,sigma=sigma)
 
