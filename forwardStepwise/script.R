source("morefuns.R")

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


a=myfs(x,y)

aa=myfs.pval(a,x,y,sigma,compute.ci=F)


n=100
p=10
nsim=200
sigma=1
tt=pv=matrix(NA,nsim,p)
for(ii in 1:nsim){
    cat(ii)
   x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)
    beta=c(4,4,4,rep(0,p-3))
y=as.numeric(x%*%beta)+rnorm(n)
        y=y-mean(y)
a=myfs(x,y)
aa=myfs.pval(a,x,y,sigma,compute.ci=T)
    
    pv[ii,]=aa$pv
}
        par(mfrow=c(2,2))
        for(i in 1:4){
            plot((1:nsim)/nsim,sort(pv[,i]))
             abline(0,1)
        }
