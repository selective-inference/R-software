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

a=forwardStep(x,y)

aa=forwardStepInf(a,x,y,compute.ci=T,nsteps=2)

aa2=forwardStepInf(a,x,y,compute.ci=T,fixed.step=4)

aa3=forwardStepInf(a,x,y,compute.ci=T, aic.stop=T)



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
beta=c(4,0,rep(0,p-2))
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
#aa=forwardStepInf(fsfit,x,y,compute.ci=F,nsteps=2)

#aa=forwardStepInf(fsfit,x,y,compute.ci=F,fixed.step=2)
aa=forwardStepInf(fsfit,x,y,compute.ci=F, aic.stop=T)
pv[ii,]=aa$pv

aichat[ii]=fsfit$aichat


}

par(mfrow=c(2,3))
for(k in 3:4){
    o=!is.na(pv[,k])
 plot((1:sum(o))/sum(o),sort(pv[o,k]))
 abline(0,1)
}

#
