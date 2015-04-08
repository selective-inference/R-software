source("funs.R")
options(error=dump.frames)
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

#aa2=forwardStepInf(a,x,y,sigma,compute.ci=T,fixed.step=4)
aa3=forwardStepInf(a,x,y,sigma,compute.ci=T, aic.stop=T)




###########
# test aic case
set.seed(33)
n=20
p=5
nsim=2000


x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(0,0,rep(0,p-2))
sigma=1
pv=matrix(NA,nsim,p)
ran=matrix(NA,nsim,2)
aichat=rep(NA,nsim)
seeds=sample(1:99999,size=nsim)
for(ii in 1:nsim){
    set.seed(seeds[ii])
    cat(ii)
y=x%*%beta+sigma*rnorm(n)
y=y-mean(y)

fsfit=forwardStep(x,y)

#aa=forwardStepInf(fsfit,x,y,sigma,compute.ci=F,nsteps=2)

aa2=forwardStepInf(fsfit,x,y,sigma,compute.ci=F,fixed.step=2)
#aa3=forwardStepInf(fsfit,x,y,sigma,compute.ci=F, aic.stop=T)
pv[ii,]=aa2$pv
aichat[ii]=fsfit$aichat
    aichat[ii]=2
 #   aichat[ii]=4
    if(!is.na(aa3$A)) ran[ii,]=range(aa3$A%*%y-aa3$b)
}

pvv=NULL
for(ii in 1:nsim){
     pvv=c(pvv,pv[ii,1:aichat[ii]])
 }
 #pvv=pv[,1]
o=!is.na(pvv)
 plot((1:sum(o))/sum(o),sort(pvv[o]))
 abline(0,1)

######

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
