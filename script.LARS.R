





library(selectiveInference,lib.loc="mylib")

library(polypath,lib.loc="/Users/tibs/dropbox/gamma/mylib")
library(truncnorm)

options(error=dump.frames)


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

larfit=lar(x,y,verbose=TRUE)
                                      

aa2=larInf(x,y,larfit,sigma)

plot(larfit)




#check covtest or lar
n=100
p=10
nsim=200
sigma=1
tt=pv2=matrix(NA,nsim,p)
for(ii in 1:nsim){
    cat(ii)
   x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)
    beta=c(3,3,0,rep(0,p-3))
y=as.numeric(x%*%beta)+rnorm(n)
        y=y-mean(y)
larfit=lar(y,x)

#tt[ii,]=covtest(larfit,x,y,sigma=sigma)
    pv2[ii,]=larInf(x,y,larfit,sigma,compute.ci=F)$pv
}
pv=1-pexp(tt,1)
        par(mfrow=c(2,2))
        for(i in 1:4){
            plot((1:nsim)/nsim,sort(pv2[,i]))
             abline(0,1)
        }
