                                      
library(selectiveInference)

library(glmnet)

set.seed(424)

#n=100
#p=30

n=20
p=40
sigma=.4
beta=c(3,2,-1,4,-2,2,rep(0,p-6))
#beta=rep(0,p)

tr=beta!=0

#type="full"
type="part"

nsim = 1000
lambda=.3
nzb=0
pvals <- matrix(NA, nrow=nsim, ncol=p)
x = matrix(rnorm(n*p),n,p)
x = scale(x,T,T)/sqrt(n-1)
mu = x%*%beta

for (i in 1:nsim) {
    cat(i)
y=mu+sigma*rnorm(n)
y=y-mean(y) 
# first run  glmnet
gfit=glmnet(x,y,intercept=F,standardize=F,thresh=1e-8)
   
#extract coef for a given lambda; Note the 1/n factor!
    bhat = coef(gfit, s=lambda/n, exact=TRUE,x=x,y=y)[-1]
    nzb=nzb+sum(bhat!=0)
# compute fixed lambda p-values and selection intervals
aa = fixedLassoInf(x,y,bhat,lambda,intercept=F,sigma=sigma,type=type)
pvals[i, aa$vars] <- aa$pv
}

# summarize results

if(type=="partial"){
nulls=rowSums(is.na(pvals[,tr]))==0  # for type=partial, nonnull setting
np = pvals[nulls,-(1:sum(beta!=0))]
}

if(type=="full"){
nulls=1:nrow(pvals)   # for type=full  non null setting
np = pvals[nulls,-(1:sum(beta!=0))]
}



#np=pvals         #for null setting

o=!is.na(np)

#check uniformity 

plot((1:sum(o))/sum(o),sort(np[o]),xlab="Expected pvalue",ylab="Observed pvalue")
abline(0,1)


 # estimate and plot FDR

pvadj=pvadj.by=pvadj.holm=matrix(NA,nsim,p)
for(ii in 1:nsim){
    o=!is.na(pvals[ii,])
    pvadj[ii,o]=p.adjust(pvals[ii,o],method="BH")
    pvadj.by[ii,o]=p.adjust(pvals[ii,o],method="BY")
      pvadj.holm[ii,o]=p.adjust(pvals[ii,o],method="holm")
    }
qqlist=fdr=se=fdr.by=se.by=fdr.holm=se.holm=c(.05, .1,.15,.2,.25,.3)
jj=0
for(qq in qqlist){
    jj=jj+1

r=v=r.by=v.by=r.holm=v.holm=rep(NA,nsim)
for(ii in 1:nsim){
    v[ii]=sum( (pvadj[ii,]<qq & !tr), na.rm=T)
    r[ii]=sum( (pvadj[ii,]<qq), na.rm=T)
    v.by[ii]=sum( (pvadj.by[ii,]<qq & !tr), na.rm=T)
    r.by[ii]=sum( (pvadj.by[ii,]<qq), na.rm=T)
      v.holm[ii]=sum( (pvadj.holm[ii,]<qq & !tr), na.rm=T)
    r.holm[ii]=sum( (pvadj.holm[ii,]<qq), na.rm=T)
    
}
oo=r!=0
    fdr[jj]=mean((v/r)[oo])
    se[jj]=sqrt(var((v/r)[oo])/sum(oo))
    oo=r.by!=0
     fdr.by[jj]=mean((v.by/r.by)[oo])
    se.by[jj]=sqrt(var((v.by/r.by)[oo])/sum(oo))
     oo=r.by!=0
     fdr.holm[jj]=mean((v.holm/r.holm)[oo])
    se.holm[jj]=sqrt(var((v.holm/r.holm)[oo])/sum(oo))
}


plot(qqlist,fdr,type="b",xlab="target FDR",ylab="observed FDR",ylim=c(0,.6),xlim=c(0,.6))
lines(qqlist,fdr.by,type="b",col=3)
lines(qqlist,fdr.holm,type="b",col=4)
abline(0,1,lty=2)
title(paste("n=",as.character(n)," p=",as.character(p),"  ",as.character(type)))
legend("bottomright",c("BH","BY","Holm"),col=c(1,3,4),lty=1)
