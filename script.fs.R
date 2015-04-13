 library(selectiveInference,lib.loc="mylib")
library(truncnorm)

options(error=dump.frames)
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

z=rnorm(n)

a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,sigma,compute.ci=T,nsteps=2)

aa2=forwardStepInf(a,x,y,sigma,compute.ci=T,fixed.step=4)

aa3=forwardStepInf(a,x,y,sigma,compute.ci=T, aic.stop=T)

a=forwardStep(x,y,z=z)
aa4=forwardStepInf(a,x,y,sigma,z=z,compute.ci=T,nsteps=2)

###########
# test aic case
library(selectiveInference,lib.loc="mylib")
library(truncnorm)
options(error=dump.frames)
set.seed(333)
n=20
p=5
nsim=500


x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(0,0,rep(0,p-2))
sigma=1
pv=matrix(NA,nsim,p)
ran=matrix(NA,nsim,2)
aichat=rep(NA,nsim)
pvz=rep(NA,nsim)
seeds=sample(1:99999,size=nsim)
for(ii in 1:nsim){
    set.seed(seeds[ii])
    cat(ii)
y=x%*%beta+sigma*rnorm(n)
y=y-mean(y)
z=rnorm(n)
    z=z-mean(z)
fsfit=forwardStep(x,y,z)
#aa=forwardStepInf(fsfit,x,y,sigma,compute.ci=F,nsteps=2)

#aa2=forwardStepInf(fsfit,x,y,sigma,compute.ci=F,fixed.step=2)
#aa3=forwardStepInf(fsfit,x,y,sigma,compute.ci=F, aic.stop=T)
aa4=forwardStepInf(fsfit,x,y,sigma,z=z,compute.ci=F,nsteps=5)
#    aa4=forwardStepInf(fsfit,x,y,sigma,z=z,compute.ci=F,fixed.step=4)
pv[ii,]=aa4$pv
  pvz[ii]=aa4$pvz
aichat[ii]=fsfit$aichat
    aichat[ii]=2

}

pvv=NULL
#for(ii in 1:nsim){
#     pvv=c(pvv,pv[ii,1:aichat[ii]])
# }
 pvv=pv[,1]
pvv=pvz
par(mfrow=c(2,3))
plot((1:nsim)/nsim,sort(pvz))
     abline(0,1)
for(k in 1:1){
 plot((1:nsim)/nsim,sort(pv[,k]))
 abline(0,1)
}

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
