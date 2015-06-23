

 library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

options(error=dump.frames)
#attach("/Users/tibs/dropbox/PAPERS/lasso/lasso3/.RData")

set.seed(133)
n=45
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(-6,3,2,-1,rep(0,p-4))
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=glmnet(x,y,standardize=F)
lambda = 1
bhat = coef(a, s=lambda/n)[-1]


# compute fixed lambda p-values

a4=fixedLassoInf(x,y,bhat,lambda,compute.si=T)

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
x=scale(x,T,F)

beta=c(rep(2,5),rep(0,p-5))
     lambda = 2
nsim=500
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
#bhat2=lasso2lam(x,y,lambda,int=F,stand=F)
 junk= fixedLassoInf(x,y,bhat,lambda,sigma=sigma)
ci[ii,]=junk$ci[1,]
    xx=x[,junk$pred]
    btrue[ii]=(solve(t(xx)%*%xx)%*%t(xx)%*%mu)[1]
}

sum(ci[,1]< btrue & ci[,]> btrue)/nsim

# another example
set.seed(13)

n=15
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)

#generate y

beta=c(5,1,-.5,-4,2,rep(0,p-5))
   nsim=10
seeds=sample(1:9999,size=nsim)
for(ii in 1:nsim){
    set.seed(seeds[ii])
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)
    

    
gfit=glmnet(x,y,standardize=F,lambda.min.ratio=1e-9)
     lambda = 9
     bhat = predict(gfit, s=lambda/n,type="coef",exact=T)[-1]
 


  #  a=lasso2lam(x,y,lambda,int=F,stand=F)
  #  critf(bhat,lambda,x,y)
# critf(a$coef,lambda,x,y)

   print(fixedLassoInf(x,y,bhat,lambda,sigma))
}

junk=tf.jonab(y,x,bhat,lambda)
junk$A%*%y-junk$b
    o=bhat!=0
    b=lsfit(x[,o],y)$coef[-1]
    c(aa$vm[1],b[1],aa$vp[1])
    eta=(solve(t(x[,o])%*%x[,o])%*%t(x[,o]))[1,]
vs=list(vm=aa$vm[1],vp=aa$vp[1])
alpha=.1

##
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
