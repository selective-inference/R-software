

library(glmnet)
#library(R.utils)
#library(truncnorm)
library(selectiveInference)

#source("/Users/tibs/Dropbox/git/R-software/selectiveInference-currentCRAN/R/funs.fixed.R")
source("funs.fixed.R")

#source("/Users/tibs/Dropbox/git/R-software/selectiveInference-currentCRAN/R/funs.inf.R")
source("funs.inf.R")


source("myfuns.R")

pinv=solve


options(error=dump.frames)
set.seed(32103)
nsim=100




dorand=T

#type="full"
type="partial"

n=100
p=12

verbose=F
sign=rep(NA,p)
ci=ci0=ci2=ci4=ci5=ci6=vlims0=array(NA,c(p,2,nsim))

betaall=matrix(NA,nsim,p)
btruepart=matrix(NA,p,nsim)
betastarall=sestarall=array(NA,c(nsim,p,nboot))
aall=acterr=rep(NA,nsim)
mc2=matrix(NA,nsim,p)
x=matrix(rnorm(n*p),n,p)
#x=matrix(rnorm(n*p),n,p)+scale(matrix(rnorm(n),n,p),F,sign(rnorm(p)))
x=scale(x,T,T)/sqrt(n-1)
sigma=1
alpha=.1
zalpha=-qnorm(alpha/2)



lam=.03 / sqrt(5)  #n=20 p=12

btrue=c(1,4,-4,rep(0,p-3))  #n=20 p=12

 
seeds=sample(1:99999,size=nsim)

for(ii in 1:nsim){
    set.seed(seeds[ii])
    cat(ii)
    mutrue=x%*%btrue
     y=mutrue+sigma*rnorm(n)
    y=y-mean(y)
    if(p>1){   
     a=glmnet(x,y,standardize=F)
     beta=as.numeric(coef(a,s=lam,exact=T,x=x,y=y))[-1]
    }
    if(p==1){
    coef=lsfit(x,y)$coef[2]
    beta=sign(coef)*(abs(coef)-n*lam)*(abs(coef)>n*lam)
    }
    if(sum(beta!=0)>0){
       
        
  betaall[ii,]=beta
  a=lsfit(x[,beta!=0],y)
    aa=ls.diag(a)
   bhat=a$coef[-1]
   bhat0=a$coef[1]
  act=which(beta!=0)
  se=aa$std.err[-1]
   btruepart[,ii]=0
   btruepart[act,ii]=lsfit(x[,act,drop=F],mutrue)$coef[-1]
#naive intervals
  ci0[beta!=0,1,ii]=bhat-zalpha*se
  ci0[beta!=0,2,ii]=bhat+zalpha*se

   #bonf-adj naive
   alpha4=alpha/p
   zalpha4=-qnorm(alpha4/2)
   ci4[beta!=0,1,ii]=bhat-zalpha4*se
  ci4[beta!=0,2,ii]=bhat+zalpha4*se

#lee et al intervals
lee=foo(x,y,beta,lam*n,alpha,type=type)
        ci[,,ii]=lee$ci

 #randomized
        if(dorand){
            rand_lasso_soln = randomizedLasso(x, y, n*lam)
        junk=randomizedLassoInf(rand_lasso_soln, nsample=5000, burnin=1000)
    ci5[rand_lasso_soln$act,1,ii]=junk$ci[,1]
    ci5[rand_lasso_soln$act,2,ii]=junk$ci[,2]
}

}}
   



if(type=="partial") btrue=btruepart

mc0=mean(ci0[,1,]>btrue | ci0[,2,]<btrue,na.rm=T)
len0=mean(ci0[,2,]-ci0[,1,],na.rm=T)
len0m=median(ci0[,2,]-ci0[,1,],na.rm=T)

ninf=mean(abs(ci)==Inf,na.rm=T)
ci[abs(ci)==Inf]=NA
mc=mean(ci[,1,]>btrue | ci[,2,]<btrue,na.rm=T)
len=mean(ci[,2,]-ci[,1,],na.rm=T)
lenm=median(ci[,2,]-ci[,1,],na.rm=T)


mc4=mean(ci4[,1,]>btrue | ci4[,2,]<btrue,na.rm=T)
len4=mean(ci4[,2,]-ci4[,1,],na.rm=T)
len4m=median(ci4[,2,]-ci4[,1,],na.rm=T)

mc5=mean(ci5[,1,]>btrue | ci5[,2,]<btrue,na.rm=T)
len5=mean(ci5[,2,]-ci5[,1,],na.rm=T)
len5m=median(ci5[,2,]-ci5[,1,],na.rm=T)



cat(c("ave# nz",mean(apply(betaall!=0,1,sum,na.rm=T))),fill=T)

miscov=c(mc0,mc,mc4,mc5)

res=
matrix(
    c(len0,len0m,0,
      len,lenm,ninf,
      len4,len4m,0,
        len5,len5m,0
      ),4,3,byrow=T)
res=cbind(res,miscov)

dimnames(res)=list(c("naive","Lee","bonf-naive","rand"),c("mean","median","propInf","miscov"))
print(res)


#out=cbind(as.vector(ci0[,2,]-ci0[,1,]),
#             as.vector(ci[,2,]-ci[,1,]),
#             as.vector(ci2[,2,]-ci2[,1,]),
##             as.vector(ci3[,2,]-ci3[,1,]),
#          as.vector(ci4[,2,]-ci4[,1,]),
#                    as.vector(ci5[,2,]-ci5[,1,] ))
#print(out)
