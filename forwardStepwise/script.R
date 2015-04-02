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

aa=forwardStepInf(a,x,y,sigma,compute.ci=T,nsteps=2)

aa2=forwardStepInf(a,x,y,sigma,compute.ci=T,fixed.step=4)

aa3=myfs.pval(a,x,y,sigma,compute.ci=F)

n=100
p=10
nsim=500
sigma=1
tt=pv=pred=matrix(NA,nsim,p)
pred=matrix(NA,nsim,4)
for(ii in 1:nsim){
    cat(ii)
   x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)
    beta=c(4,4,4,rep(0,p-3))
y=as.numeric(x%*%beta)+rnorm(n)
        y=y-mean(y)
a=forwardStep(x,y)
#aa=forwardStepInf(a,x,y,sigma,compute.ci=F)
 aa=forwardStepInf(a,x,y,sigma,compute.ci=F,fixed.step=4)   
    pv[ii,]=aa$pv
    pred[ii,]=a$pred[1:4]
}
        par(mfrow=c(2,2))
        for(i in 1:4){
            if(i<4) o=pred[,i]<4
            if(i==4) o=pred[,i]>3
            plot((1:sum(o))/sum(o),sort(pv[o,i]))
             abline(0,1)
        }
