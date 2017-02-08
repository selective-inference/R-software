library(selectiveInference)
#library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

library(glmnet)
library(MASS)
library(scalreg)
#options(error=dump.frames)
#attach("/Users/tibs/dropbox/PAPERS/lasso/lasso3/.RData")

#####
#gaussian
n=50
p=10
sigma=.7
beta=c(3,2,0,0,rep(0,p-4))
set.seed(43)
nsim = 200
pvals <- matrix(NA, nrow=nsim, ncol=p)
x = matrix(rnorm(n*p),n,p)
x = scale(x,T,T)/sqrt(n-1)
mu = x%*%beta
for (i in 1:nsim) {
    cat(i)
y=mu+sigma*rnorm(n)
#y=y-mean(y)
# first run  glmnet
gfit=glmnet(x,y,intercept=F,standardize=F,thresh=1e-8)
lambda = 1
#extract coef for a given lambda; Note the 1/n factor!
beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]
# compute fixed lambda p-values and selection intervals
aa = fixedLassoInf(x,y,beta,lambda,intercept=F,sigma=sigma)
pvals[i, which(beta != 0)] <- aa$pv
}
nulls = which(!is.na(pvals[,1]) & !is.na(pvals[,2]))
np = pvals[nulls,-(1:2)]
mean(np[!is.na(np)] < 0.1)
o=!is.na(np)
plot((1:sum(o))/sum(o),sort(np))
abline(0,1)
#####


S <- diag(10)
n <- 100
p <- 10
pval <- matrix(1, nrow = 100, ncol = p)
for(i in 1:100){
    cat(i)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = S)
  Y <- X[, 1] + X[, 2] + rnorm(n)
  sig.L <- scalreg(X, Y)$hsigma

  lam <- cv.glmnet(X, Y, standardize = FALSE, intercept = FALSE)$lambda.min
  bl <- glmnet(X, Y, lambda = lam, standardize = FALSE, intercept = FALSE)$beta[, 1]
  idx <- which(bl != 0)
  pval[i, idx] <- fixedLassoInf(X, Y, beta = bl, lambda = lam * n, intercept = FALSE, sigma = sig.L, alpha = 0.05)$pv
}

p <- pval[, -(1:2)]
mean(p[p < 1] < 0.05)

##logistic

n=50
p=10
beta=c(3,2,0,0,rep(0,p-4))
beta=rep(0,p)
set.seed(3)
nsim = 200
pvals=matrix(NA, nrow=nsim, ncol=p)
ci=array(NA,c(nsim,p,2))
x = matrix(rnorm(n*p),n,p)
x = scale(x,T,T)/sqrt(n-1)
mu = x%*%beta
for (ii in 1:nsim) {
    cat(ii)
y=mu+rnorm(n)
y=1*(y>mean(y))
# first run  glmnet
gfit=glmnet(x,y,standardize=F,thresh=1e-8,family="binomial")
lambda = .25
#extract coef for a given lambda; Note the 1/n factor!
beta = as.numeric(coef(gfit, s=lambda/n, exact=TRUE))
# compute fixed lambda p-values and selection intervals
    aa = fixedLassoInf(x,y,beta,lambda,family="binomial")
    pvals[ii, which(beta[-1] != 0)] <- aa$pv
   ci[ii,which(beta[-1] != 0),]=aa$ci
}

o=!is.na(pvals)
plot((1:sum(o))/sum(o),sort(pvals))
abline(0,1)
o=ci[,1,1]>0 | ci[,1,2]<0
mean(o,na.rm=T)


## cox

n=50
p=10
#beta=c(6,6,0,0,rep(0,p-4))
beta=rep(0,p)
set.seed(3)
nsim = 200
pvals=matrix(NA, nrow=nsim, ncol=p)
ci=array(NA,c(nsim,p,2))
x = matrix(rnorm(n*p),n,p)
x = scale(x,T,T)/sqrt(n-1)
mu = x%*%beta
for (ii in 1:nsim) {
    cat(ii)
tim=as.vector(mu+rnorm(n))+10
status=sample(c(0,1),size=n,replace=T)
    lambda=0.2
    y=cbind(time=tim,status=status)
    gfit=glmnet(x,y,family="cox",standardize=FALSE)
     b=as.numeric(coef(gfit,s=lambda/n,exact=TRUE))
 
   aa= fixedLassoInf(x,tim,b,lambda,status=status,family="cox")
      
pvals[ii, which(b != 0)] <- aa$pv[1:sum(!is.na(aa$pv))]
  ci[ii,which(b != 0),]=aa$ci
}

o=!is.na(pvals)
plot((1:sum(o))/sum(o),sort(pvals))
abline(0,1)


#####more Gaussian

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
nsim=100

x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)    #try with and without standardization

beta=c(5,4,3,2,1,rep(0,p-5))


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
    pvals[ii, which(bhat != 0)] <- aa$pv[1:sum(!is.na(aa$pv))]
  ci[ii,which(bhat != 0),]=aa$ci
   
}
o=!is.na(pvals)
plot((1:sum(o))/sum(o),sort(pvals))

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
junk= fixedLassoInf(x,y,bhat,lambda,sigma=sigma,bits=200)

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
x=scale(x,T,F)
#x=scale(x,T,T)
n=length(y)
nams=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/data64.names",what="")
y=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/diab.y")
y=y-mean(y)

cvf=cv.glmnet(x,y)
sigmahat=estimateSigma(x,y,stand=F)$si

lambda=n*cvf$lambda.min

lambda=estimateLambda(x,sigma=sigmahat)/2

gfit=glmnet(x,y,standardize=F)
 
bhat=coef(gfit, s=lambda/n, exact=TRUE)[-1]
bhat2=lasso2lam(x,y,lambda,int=F,stand=F)$coef

plot(bhat,bhat2)

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

# now try penalty factors
  set.seed(43)
     n = 50
     p = 10
     sigma = 1
     
     x = matrix(rnorm(n*p),n,p)
     x=scale(x,T,T)
     
     beta = c(3,2,rep(0,p-2))
     y = x%*%beta + sigma*rnorm(n)
    
     pf=c(rep(1,7),rep(.1,p-7))
     pf=p*pf/sum(pf)   # penalty factors should be rescaled so they sum to p
     xs=scale(x,F,pf) #scale cols of x
     # first run glmnet
     gfit = glmnet(xs,y,standardize=F)
     
     # extract coef for a given lambda; note the 1/n factor!
     # (and we don't save the intercept term)
     lambda = .8
     beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]
     
     # compute fixed lambda p-values and selection intervals
     out = fixedLassoInf(xs,y,beta,lambda,sigma=sigma)
     #rescale conf points to undo the penalty factor
     out$ci=t(scale(t(out$ci),F,pf[out$vars]))


###
x=state.x77[,-4]
y=state.x77[,4]
x=scale(x,T,T)
n=nrow(x)


sigmahat=estimateSigma(x,y,stand=F)$si

lambda=65

gfit=glmnet(x,y,standardize=F, thresh=1e-9)
 
bhat=coef(gfit, s=lambda/n, exact=TRUE)[-1]
#bhat2=lasso2lam(x,y,lambda,int=F,stand=F)$coef

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
#X = scale(X,center=T,scale=T)  # original
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
## new bugs from lucas
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

set.seed(1)
p <- 500
n <- 400
s0 <- 10
b <- 1.2
b0 <- 0
sigma <- 5
alpha = 0.05

X = matrix(rnorm(n*p),n,p)
X = scale(X,center=T,scale=T)
theta0 <- c(rep(b,s0),rep(0,p-s0))
w <- sigma*rnorm(n)
y <- b0 + X%*%theta0 + w

# Pick lambda as in Negahban et al. (2012), as done in Lee et al. (2015)
m = 1000
eps = sigma*matrix(rnorm(m*n),n,m)
lam = 2*mean(apply(abs(t(X)%*%eps),2,max))

gfit = glmnet(X,y,standardize=F)
coef = coef(gfit, s=lam/n, exact=T)[-1]
sint = fixedLassoInf(X,y,coef,lam,sigma=sigma,alpha=alpha,type="partial")
# Error in v %*% diag(d) : non-conformable arguments

## lucas again

load("params_for_Rob.rdata") #variables: X, y, coef, lam, sigma, alpha

sint = fixedLassoInf(X,y,coef,lam,sigma=sigma,alpha=alpha,type="partial")


#### bug from Sen at UW

library(MASS)
library(scalreg)

S <- diag(10)
n <- 100
p <- 10
pval <- matrix(1, nrow = 100, ncol = p)
for(i in 1:100){
    cat(i)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = S)
  Y <- X[, 1] + X[, 2] + rnorm(n)
  sig.L <- scalreg(X, Y)$hsigma

  lam <- cv.glmnet(X, Y, standardize = FALSE, intercept = FALSE)$lambda.min
  bl <- glmnet(X, Y, lambda = lam, standardize = FALSE, intercept = FALSE)$beta[, 1]
  idx <- which(bl != 0)
  pval[i, idx] <- fixedLassoInf(X, Y, beta = bl, lambda = lam * n, intercept = FALSE, sigma = sig.L, alpha = 0.05)$pv
}

p <- pval[, -(1:2)]
mean(p[p < 1] < 0.05)


#test from Chong


library(selectiveInference)

library(glmnet);library(MASS);#library(grplasso);library(gvlma);library(grpreg)
library(penalized)
load("fooXY.RData")

#d=read.csv("DesignMatrixX_and_y.csv");dim(d); head(d)

#source("temp.R")
n=length(Y)
p=ncol(X)
#X=scale(X,T,F)
X=X+.01*matrix(rnorm(n*p),n,p) # I added noise to avoid collinearity
#X=scale(X,T,T)/sqrt(n-1)

Y=Y-mean(Y)

#X=as.matrix(d[,1:192]); ### design matrix, no intercept
#Y=d$y; ### Response variable Y.
fit = glmnet(x=X, y=Y, family="gaussian",alpha = 1, thresh = 1e-9, standardize=F)
set.seed(39)
lam= fit$lambda[30];
#lam= fit$lambda[15]## Try getting coefficient at this lambda
beta = coef(fit, s=lam, exact=TRUE)[-1];length(beta);table(beta!=0)

 aa=penalized(Y~X,lambda1=lam*n,model="linear",standardize=F)
  b=coef(aa,which="all")[-1]

lam2=n*lam

g=t(X)%*%(Y-X%*%beta)/lam2

g[beta!=0]

g=t(X)%*%(Y-X%*%b)/lam2
out = fixedLassoInf(X,Y,beta,lam*n)





#

#gaussian
n=50
p=10
sigma=.7
beta=c(0,0,0,0,rep(0,p-4))
set.seed(43)
nsim = 1000
pvals <- matrix(NA, nrow=nsim, ncol=p)
x = matrix(rnorm(n*p),n,p)
x = scale(x,T,T)/sqrt(n-1)
mu = x%*%beta
for (i in 1:nsim) {
    cat(i)
y=mu+sigma*rnorm(n)
#y=y-mean(y)
# first run  glmnet
    pf=c(rep(.001,4),rep(1,p-4))
     xs=scale(x,FALSE,pf) #scale cols of x by penalty factors
     # first run glmnet
     gfit = glmnet(xs,y,standardize=FALSE)
     
    
     lambda = .8
     beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]
     
     # compute fixed lambda p-values and selection intervals
     aa = fixedLassoInf(xs,y,beta,lambda,sigma=sigma)
     
pvals[i, which(beta != 0)] <- aa$pv
}
nulls = 1:nsim
np = pvals[nulls,-(1:4)]
mean(np[!is.na(np)] < 0.1)
o=!is.na(np)
plot((1:sum(o))/sum(o),sort(np))
abline(0,1)
#####


