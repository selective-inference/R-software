library(selectiveInference)#,lib.loc="/Users/tibs/dropbox/git/R/mylib")

options(error=dump.frames)
set.seed(333)
n=40
p=16

n=200
p=20
sigma=.7

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)

#generate y
#beta=c(3,-2,0,0,rep(0,p-4))
#beta=c(rep(2,10),rep(0,p-10))
beta=c(rep(100,10),rep(0,p-10))
theta = x%*%beta
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,compute.si=T,trace=T)

aa2=forwardStepInf(a,x,y,sigma=sigma,compute.si=T,fixed.step=12)

aa3=forwardStepInf(a,x,y,sigma=sigma,compute.si=T, aic.stop=T)

x=state.x77[,-4]
y=state.x77[,4]
a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,compute.si=T,trace=T)


sigmahat=estimateSigma(x,y)
aa=forwardStepInf(a,x,y,sigma=sigmahat,compute.si=T,trace=T)


n=20
p=100
sigma=.7

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,F)

#generate y
#beta=c(3,-2,0,0,rep(0,p-4))
#beta=c(rep(2,10),rep(0,p-10))
beta=c(rep(10,10),rep(0,p-10))
y=x%*%beta+sigma*rnorm(n)

y=y-mean(y)

a=forwardStep(x,y)

a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,compute.si=T,trace=T)

a=forwardStep(x,y,sigma=sigmahat)
sigmahat=estimateSigma(x,y)
aa=forwardStepInf(a,x,y,sigma=sigmahat,compute.si=T,trace=T)


fsfit=a
sigma=a$sigma
aic.stop=F
trace=F
alpha=.1
fixed.step=NULL
nsteps=NULL
compute.Si=T
one.sided=T


###########
# tests
library(selectiveInference,lib.loc="mylib")
library(truncnorm)
options(error=dump.frames)

setting=1  # sequential steps
#setting=2  #fixed stop
setting=3  #AIC stop

set.seed(333)
n=100
p=50
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
if(setting==1) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.si=F,nsteps=2)
if(setting==2) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.si=F,fixed.step=2)
if(setting==3) aa=forwardStepInf(fsfit,x,y,sigma=sigma,compute.si=F, aic.stop=T)
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


A=matrix(rnorm(100*20),ncol=20)
y=rnorm(20)
b=rep(0,100)
eta=rnorm(20)
pp=100
sigma=1

#####
set.seed(40)
 n=20
 p=8
sigma=1
 nsteps=1
 nsim=500
 
 x=matrix(rnorm(n*p),ncol=p)
x=scale(x,T,T)/sqrt(n-1)
beta=c(1,rep(0,p-1))
y=x%*%beta+sigma*rnorm(n)
x=scale(x,T,F)/sqrt(n-1)
y=y-mean(y)
a=forwardStep(x,y,sigma=sigma)

aa=forwardStepInf(a,x,y,compute.si=T,trace=T)
