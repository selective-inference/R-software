
myfs=function(x,y,nsteps=min(nrow(x),ncol(x))){
p=ncol(x)
# fs by minimizing scaled ip
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
pred=s=scor=bhat=rep(NA,nsteps)
   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
for(j in 2:nsteps){
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r)$res
  xr= lsfit(x[,mod],x)$res
  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
 ip[mod]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
 bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
}
return(list(pred=pred,s=s,scor=scor,bhat=bhat))
}



myfs.pval=
function(fsfit,x,y,sigma,nsteps=NULL,signed.eta=FALSE,alpha=.10,spacing.paper=TRUE,which.pred=1:nsteps){
# pvalues for forwrd stepwise
# does ryan's two-sided version
#  returns pval and interval for predictor just entered (default with which.pred=1:nsteps) 
# otherwise which.pred is a vector of length nteps, and it returns 
#  info for which.pred[k]th pred at step k
n=nrow(x)
p=ncol(x)
if(is.null(nsteps)){ nsteps=length(a$pred)}
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
SMALL=1e-7
pred=fsfit$pred
s=fsfit$s
a=NULL
proj=0

stepind=NULL
for(j in 1:nsteps){
  notin=rep(T,p)
  if(j>1) {
  mod=pred[1:(j-1)]
  notin[mod]=F
  x[,-mod]=lsfit(x[,mod],x[,-mod])$res
  proj=x[,mod]%*%solve(t(x[,mod])%*%x[,mod])%*%t(x[,mod])
  }
  jj=pred[j]
  rest=notin;rest[jj]=F
 w=sqrt(sum(x[,jj]^2))
  for(k in which(rest)){
 w2=sqrt(sum(x[,k]^2))
a=rbind(a, -(s[j]*x[,jj]/w-x[,k]/w2)%*%(diag(n)-proj))
a=rbind(a, -(s[j]*x[,jj]/w+x[,k]/w2)%*%(diag(n)-proj))
stepind=c(stepind,j,j)
}
}
# last step- constraint on sign only
mod=pred[1:nsteps]
fitter=solve(t(x[,mod])%*%x[,mod])%*%t(x[,mod])
aa=-1*fitter[nsteps,]*s[nsteps]
a=rbind(a,aa)
stepind=c(stepind,nsteps)

# pv tests
vmall=vpall=vector("list",nsteps)
pv=rep(NA,p)
ci=cov=matrix(NA,nrow=p,ncol=2)

for(k in 1:nsteps){
pp=sum(stepind<=k)
b=rep(0,pp)
mod=pred[1:k]
 temp=(solve(t(x[,mod,drop=F])%*%x[,mod,drop=F])%*%t(x[,mod,drop=F]))
# compute pvalue only for predictor which.pred[k]
for(jj in which.pred[k]){
  eta=as.vector(temp[jj,])
  if(signed.eta) eta=eta*s[which.pred[k]]
  alp=as.vector(a%*%eta/sum(eta^2))
  alp[abs(alp)<SMALL]=0
  vp=rep(Inf,pp)
  vm=rep(-Inf,pp)
  for(j in 1:pp){
   if(alp[j]<0) vm[j]=(b[j]-(a%*%y)[j]+alp[j]*sum(eta*y))/alp[j]
   if(alp[j]>0) vp[j]=(b[j]-(a%*%y)[j]+alp[j]*sum(eta*y))/alp[j]
   }
   vmm=max(vm,na.rm=T)
   vpp=min(vp,na.rm=T)
   vmall[[k]][jj]=vmm
   vpall[[k]][jj]=vpp
   tt=sum(eta*y)
   sigma.eta=sigma*sqrt(sum(eta^2))
   u=0  #null
   pv[k]=1-(pnorm((tt-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))/(pnorm((vpp-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))
    pv[k]=2*min(pv[k],1-pv[k])
junk=selint2(vmm,vpp,tt,alpha,sigma.eta=sigma.eta,spacing.paper=spacing.paper)

# Check this is correct!!

#ci[k,]=junk$int*s[k]
ci[k,]=junk$int
cov[k,]=junk$cov
   }
 }
return(list(pv=pv,vm=vmall,vp=vpall,ci=ci,cov=cov))
}


selint2=function(vmm,vpp,tt,alpha,sigma.eta=1,maxbound=20,maxdel=50,np=1000,eps=1e-6,spacing.paper=T){
# construct interval
# spacing.paper flag I added aft spacing paper was done (it used T)
cov=rep(NA,np)
low=tt-maxdel*sigma.eta
up=tt+maxdel*sigma.eta
# compute tail areas for a range of bvals
bval=sort(seq(low,up,length=np))
for(i in 2:(np-1)){
cov[i]=selint3(vmm,vpp,bval[i],tt,alpha,sigma.eta,maxbound=maxbound)
}
# find points that achieve tail areas of alpha/2
# find points that achieve tail areas of alpha/2
temp=1-cov
if(spacing.paper) temp[temp<(alpha/2)]=NA
temp[bval<tt]=NA
cup=Inf;covright=0
if(sum(!is.na(temp))>0){
o1=which.min(abs(temp-(alpha/2)))
cup=bval[o1]
covright=1-cov[o1]
}
temp2=cov
if(spacing.paper){temp2[temp2<(alpha/2)]=NA}
temp2[bval>tt]=NA
clow=-Inf;covleft=0
if(sum(!is.na(temp2))>0){
o2=which.min(abs(temp2-(alpha/2)))
clow=bval[o2]
covleft=cov[o2]
}
return(list(interval=c(clow,cup),cov=c(covleft,covright)))
}

selint3=function(vmm,vpp,betaval,tt,alpha, sigma.eta, maxbound,eps=1e-4,np=100){
# find tail area to the left of tt, assuming truth is betaval
low=max(vpp,-maxbound)
up=min(vmm,maxbound)
cov=NA
#coverage in left tail
 cov=areaf(tt,betaval=betaval,vmm=vmm,vpp=vpp,sigma.eta=sigma.eta)
return(cov)
}
areaf=function(tt,betaval,vmm,vpp,sigma.eta){
 #compute Prob_betaval (W<tt| W in [vmm,vpp]
  val=NA
  if(tt>=vmm & tt <=vpp){
      val0=pnorm((tt-betaval)/sigma.eta,log.p=TRUE)
      val1=pnorm((vmm-betaval)/sigma.eta,log.p=TRUE)
      val2=pnorm((vpp-betaval)/sigma.eta,log.p=TRUE)

  val=(1-exp(val2-val0))/(exp(val1-val0)-exp(val2-val0))
#   val0=pnorm((tt-betaval)/sigma.eta,log.p=F)
#      val1=pnorm((vmm-betaval)/sigma.eta,log.p=F)
#      val2=pnorm((vpp-betaval)/sigma.eta,log.p=F)
# if(val1!=val2)  val=(val0-val2)/(val1-val2)
}
return(val)
}


genx.pw=
function(n = 20, rr = 0, p = 8)
{
#    generate x's multivariate normal with pairwise corr rr
        inds <- 1:p
        Sigma <- matrix(rr,nrow=p,ncol=p)
Sigma[row(Sigma)==col(Sigma)]=1
        bb <- svd(Sigma)
#        hh <- bb$u %*% (sqrt(diag(bb$d))) %*% t(bb$v)
hh=scale(bb$u,center=F,scale=1/sqrt(bb$d))%*%t(bb$v)

        x <- matrix(rnorm(p * n), n, p) %*% hh
        return(list(x=x, Sigma=Sigma))
}

