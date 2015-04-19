 library(selectiveInference,lib.loc="mylib")
library(truncnorm)

options(error=dump.frames)
set.seed(3)
n=40
p=16

n=100
p=100
sigma=.7

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(-5,0,0,0,rep(0,p-4))
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,compute.si=T,alpha=.05,trace=T)

aa2=forwardStepInf(a,x,y,sigma=sigma,compute.si=T,fixed.step=4)

aa3=forwardStepInf(a,x,y,sigma=sigma,compute.si=T, aic.stop=T)

fsfit=a
sigma=a$sigma
aic.stop=F
trace=F
alpha=.1
fixed.step=NULL
nsteps=10
compute.ci=T
one.sided=T
gridfac=50


###########
# tests
library(selectiveInference,lib.loc="mylib")
library(truncnorm)
options(error=dump.frames)

setting=1  # sequential steps
#setting=2  #fixed stop
#setting=3  #AIC stop

set.seed(333)
n=20
p=10
nsim=500


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

fsfit=forwardStep(x,y,sigma=sigma)
if(setting==1) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.ci=F,nsteps=2)
if(setting==2) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.ci=F,fixed.step=2)
if(setting==3) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.ci=F, aic.stop=T)
pv[ii,]=aa$pv

aichat[ii]=fsfit$aichat

}


if(setting<3){
par(mfrow=c(2,3))
for(k in 1:4){
    o=!is.na(pv[,k])
 plot((1:sum(o))/sum(o),sort(pv[o,k]))
 abline(0,1)
}
}

if(setting==3){
pvall=NULL
for(ii in 1:nsim){
    o=1:aichat[ii]
    pvall=c(pvall,pv[ii,o])
    }

 o=!is.na(pvall)
 plot((1:sum(o))/sum(o),sort(pvall))
 abline(0,1)
 }


#
