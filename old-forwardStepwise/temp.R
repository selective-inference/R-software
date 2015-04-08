source("funs.R")

source("../fixedLamLasso/funs.R")  # to get pv+CI funcs
set.seed(33)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(4,2,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)
y=y-mean(y)

a=forwardStep(x,y)
aic=rss=val=rep(NA,p)

rss[1]=sum( (y-mean(y))^2)
aic[1]=rss[1]+2*sigma^2
for(j in 2:p){
      xx=x[,a$pred[1:j],drop=F]
  rss[j]=sum(lsfit(xx,y)$res^2)
  aic[j]=rss[j]+2*j*sigma^2
  b=lsfit(xx,y)
  b=b$coef[length(b$coef)]
  v=solve(t(xx)%*%xx)
  v=v[ncol(xx),ncol(xx)]
 val[j]=b/(sqrt(v)*sigma)
}
d=diff(aic)
o=which(d<0)
o=o[length(o)]

aa=forwardStepInf(a,x,y,sigma,compute.ci=F)


n=100
p=10
nsim=500
sigma=1
tt=pv=matrix(NA,nsim,p)
for(ii in 1:nsim){
    cat(ii)
   x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)
    beta=c(4,4,4,rep(0,p-3))
y=as.numeric(x%*%beta)+rnorm(n)
        y=y-mean(y)
a=forwardStep(x,y)
aa=forwardStepInf(a,x,y,sigma,compute.ci=F)
    
    pv[ii,]=aa$pv
}
        par(mfrow=c(2,2))
        for(i in 1:4){
            plot((1:nsim)/nsim,sort(pv[,i]))
             abline(0,1)
        }
