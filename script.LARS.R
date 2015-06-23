


library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

library(lars)

options(error=dump.frames)


set.seed(33)
n=200
p=20

sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(3,-2,rep(0,p-2))
beta=c(rep(3,trunc(p/2)),rep(0,p-trunc(p/2)))
y=x%*%beta+sigma*rnorm(n)
larfit=lar(x,y,verbose=TRUE)
fit=predict.lar(larfit,x,s=2,type="fit")
                                      

aa2=larInf(x,y,larfit,nsteps=6)

plot(larfit)




#check qqplots for covtest and lar
set.seed(3)
n=100
p=10
nsim=300
sigma=2
tt=tt2=pv0=pv=pv2=sel=matrix(NA,nsim,p)
   x=matrix(rnorm(n*p),n,p)

#   x=scale(x,T,T)/sqrt(n-1)

 

beta=c(3,-4,3,rep(0,p-3))
for(ii in 1:nsim){
    cat(ii)

y=as.numeric(x%*%beta)+sigma*rnorm(n)
larfit=lar(x,y)
    pv0[ii,]=larInf(x,y,larfit,sigma,compute.si=F)$pv
sel[ii,]=larfit$act
tt[ii,]=covtest(larfit,x,y,sigma=sigma)
    b=lars(x,y)
    tt2[ii,]=covTest(b,x,y,sigma.est=sigma)$res[,2]
 
   
}

pv=1-pexp(tt,1)
        par(mfrow=c(2,2))
        for(i in 1:4){
            o=rep(T,nsim)
            if(i>3) o=rowSums(sel[,1:3])==6
            plot((1:sum(o))/sum(o),sort(pv0[o,i]))
             points((1:sum(o))/sum(o),sort(pv[o,i]),col="green")
             abline(0,1)
        }


#extra check vs covTest lib
library(covTest)
set.seed(4)
sigma=4.3
n=20
p=6

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(3,-2,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)
larfit=lar(x,y,verbose=TRUE)
                                      
#
aa2=larInf(x,y,larfit,sigma=sigma,nsteps=6)

b=lars(x,y)
bb=covTest(b,x,y,sigma.est=sigma)
