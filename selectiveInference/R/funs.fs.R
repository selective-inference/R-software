
forwardStep=function(x,y,sigma=NULL,nsteps=min(nrow(x)-1,ncol(x))){
    this.call=match.call()
    checkargs(x=x,y=y,nsteps=nsteps,sigma=sigma)
    BIG=10e9
    n=nrow(x)
p=ncol(x)
    
# fs by minimizing scaled ip
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
    if(is.null(sigma) & p>=n){cat("Warning: p ge n; sigma=1 used for AIC",fill=T); sigma=1}
  if(is.null(sigma) & n>p){
      sigma=sqrt(sum(lsfit(x,y)$res^2)/(n-p))
  }
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
#compute AIC
    aic=rss=val=rep(NA,nsteps+1)
rss[1]=sum( (y-mean(y))^2)
    aic[1]=rss[1]
for(j in 1:nsteps){
      xx=x[,pred[1:j],drop=F]
  rss[j+1]=sum(lsfit(xx,y)$res^2)
  aic[j+1]=rss[j+1]+2*j*sigma^2
 }
 d=diff(aic)
o=which(d>0)
    if(length(o)>0) aichat=o[1]-1
    if(length(o)==0) aichat=length(aic)-1
    out=list(pred=pred,scor=scor,bhat=bhat,rss=rss,sigma=sigma,aic=aic,aichat=aichat,call=this.call)
    class(out)="forwardStep"
return(out)
}




print.forwardStep=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat("",fill=T)
tab=cbind(1:length(x$pred),x$pred,x$bhat,x$scor)
      dimnames(tab)=list(NULL,c("step","predictor","coef","scoreStat"))
      print(tab)
      cat("",fill=T)
            cat("",fill=T)
      cat(c("Value used for error standard deviation (sigma)=",round(x$sigma,6)),fill=T)
       cat("",fill=T)
      cat(c("AIC optimal model size=",x$aichat),fill=T)
  }

forwardStepInf=
function(fsfit,x,y,sigma=NULL,nsteps=NULL,alpha=.10,fixed.step=NULL,aic.stop=FALSE,trace=F,compute.si=TRUE,one.sided=TRUE){
# pvalues for forward stepwise
#  returns pval and interval for predictor just entered (default with which.pred=1:nsteps) 
# otherwise which.pred is a vector of length nsteps, and it returns 
#  info for which.pred[k]th pred at step k
    this.call=match.call()
    
      checkargs(x=x,y=y,nsteps=nsteps,sigma=sigma,alpha=alpha, fsfit=fsfit)
    
n=nrow(x)
p=ncol(x)
    
  if(is.null(sigma)) sigma=fsfit$sigma
    
    if(aic.stop){fixed.step=fsfit$aichat}
if(is.null(nsteps)){ nsteps=min(length(fsfit$pred),20)}
    if(!is.null(fixed.step))  nsteps=fixed.step
# first center x and y

x=scale(x,TRUE,FALSE)
y=y-mean(y)

xx=x
 
SMALL=1e-7
pred=fsfit$pred
s=sign(fsfit$scor)
a=NULL
proj=0
which.steps=1:nsteps
if(!is.null(fixed.step)){ which.steps=rep(fixed.step,fixed.step)}
if(aic.stop & fsfit$aichat==0){
    out=list(pv=rep(NA,p),ci=NA,A=NA,b=NA)
    cat("Null model picked by AIC; no  results returned",fill=T)
    return(out)
}

# construct matrices A and b
    vv=1/sum(xx[,pred[1]]^2)
stepind=NULL

 A=b=NULL


if(trace) cat("Constructing constraint matrix",fill=T)
for(j in 1:nsteps){
if(trace) cat(c("step=",j),fill=T)
  notin=rep(T,p)
  if(j>1) {
  mod=pred[1:(j-1)]
  notin[mod]=F
  xx[,-mod]=lsfit(xx[,mod],xx[,-mod])$res
  vv=solve(t(xx[,mod,drop=F])%*%xx[,mod,drop=F])
  proj=xx[,mod]%*%vv%*%t(xx[,mod])
  }
  jj=pred[j]
  rest=notin;rest[jj]=F
 w=sqrt(sum(xx[,jj]^2))
  for(k in which(rest)){
 w2=sqrt(sum(xx[,k]^2))
A=rbind(A, -(s[j]*xx[,jj]/w-xx[,k]/w2)%*%(diag(n)-proj))
A=rbind(A, -(s[j]*xx[,jj]/w+xx[,k]/w2)%*%(diag(n)-proj))
 b=c(b,0,0)
stepind=c(stepind,j,j)
}
 if(aic.stop){
     #constrain that change in AIC > 2sigma
      mod2=pred[1:j]
       vv=solve(t(xx[,mod2,drop=F])%*%xx[,mod2,drop=F])
   temp=vv%*%t(xx[,mod2,drop=F])
      temp2=-s[j]*temp[j,,drop=F]/(sigma*sqrt(vv[j,j]))
  A=rbind(A, temp2)
  b=c(b,-sqrt(2))
}
}

    
# last step- constraint on sign only (for fixed.step=NULL only)
if(is.null(fixed.step)){
mod=pred[1:nsteps]
fitter=solve(t(xx[,mod])%*%xx[,mod])%*%t(xx[,mod])
aa=-1*fitter[nsteps,]*s[nsteps]
A=rbind(A,aa)
stepind=c(stepind,nsteps)
b=c(b,0)
}
#for aic stop, add constraint that diff in AIC is < 2sigma, so that it stops
 # (equiv to bhat/se < sqrt(2))
if(aic.stop){
      mod2=pred[1:(nsteps+1)]
 vv=solve(t(xx[,mod2,drop=F])%*%xx[,mod2,drop=F])
     pp=length(mod2)
      temp=vv%*%t(xx[,mod2,drop=F])
     snew=sign(sum(temp[pp,]*y))
     temp2=snew*temp[pp,,drop=F]/(sigma*sqrt(vv[pp,pp]))
     A=rbind(A, temp2)
  b=c(b,sqrt(2))
 }
    
    
# pv tests

vmall=vpall=vector("list",nsteps)
pv=rep(NA,p)
ci=miscov=matrix(NA,nrow=p,ncol=2)


if(trace) cat("Computing p-values",fill=T)
for(kk in 1:length(which.steps)){
if(trace) cat(c("step=",kk),fill=T)
if(is.null(fixed.step)) pp=sum(stepind<=which.steps[kk])
if(!is.null(fixed.step)) pp=nrow(A)
mod=pred[1:which.steps[kk]]

 temp=(solve(t(x[,mod,drop=F])%*%x[,mod,drop=F])%*%t(x[,mod,drop=F]))

# compute pvalues and CIs  
kkk=nrow(temp)
if(!is.null(fixed.step)) kkk=kk
  eta=as.vector(temp[kkk,])

    bhat=sum(eta*y)
    if(one.sided) eta=eta*sign(bhat)
flip=(one.sided & sign(bhat)==-1)
junk=compute.vmvp(eta,A,b,pp,sigma)
vmm=junk$vm;vpp=junk$vp
   vmall[[kk]][kk]=vmm
   vpall[[kk]][kk]=vpp
   tt=sum(eta*y)
   sigma.eta=sigma*sqrt(sum(eta^2))
   u=0  #null
  # pv[kk]=1-(pnorm((tt-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))/(pnorm((vpp-u)/sigma.eta)-pnorm((vmm-u)/sigma.
  pv[kk]= 1-rob.ptruncnorm(tt, vmm, vpp, u,sigma.eta)
if(!one.sided)  pv[kk]=2*min(pv[kk],1-pv[kk])

  
  if(compute.si)
      {
          if(trace) cat("Computing selection interval",fill=T)
           vs=list(vm=vmm,vp=vpp)
          junk=selection.int(y,eta,sigma,vs,alpha,flip=flip)
#     cat(c(vs$vm,sum(eta*y),vs$vp,sigma,sigma.eta,alpha),fill=T)
          ci[kk,]=junk$ci;miscov[kk,]=junk$miscov
      }


}
    
forwardStopHat=NULL
if(is.null(fixed.step)) forwardStopHat=forwardStop(pv,alpha)
        
        
out=list(pv=pv,vm=vmall,vp=vpall,ci=ci,tailarea=miscov,pred=pred,which.steps=which.steps,stepind=stepind,forwardStopHat=forwardStopHat,alpha=alpha,sigma=sigma,one.sided=one.sided,A=A,b=b,call=this.call)

    class(out)="forwardStepInf"
return(out)
}


print.forwardStepInf=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat("",fill=T)
      cat(c("alpha=",x$alpha),fill=T)
      if(x$one.sided) cat("p-values below are one-sided",fill=T)
       if(!x$one.sided) cat("p-values below are two-sided",fill=T)
cat("",fill=T)
cat("",fill=T)
      nn=length(x$which.steps)
tab=cbind(x$which.steps,x$pred[1:nn],round(x$pv[1:nn],6),round(x$ci[1:nn,],6),round(x$tailarea[1:nn,],6))
      
      
        dimnames(tab)=list(NULL,c("step","predictor","p-value","lowerConfPt","upperConfPt","lowerArea","upperArea"))
      print(tab)

if(!is.null(x$forwardStopHat)){
cat("",fill=T)
cat(c("Estimated stopping point from forwardStop rule=", x$forwardStopHat),fill=T)
         cat("",fill=T)
}
      cat(c("Value used for error standard deviation (sigma)=",round(x$sigma,6)),fill=T)
  }
  
  




compute.vmvp=function(eta,A,b,pp,sigma,SMALL=1e-7){
 # we use first pp constraints
    A=A[1:pp,,drop=F]
    b=b[1:pp]
  alp=as.vector(A%*%eta/sum(eta^2))
  alp[abs(alp)<SMALL]=0
  vp=rep(Inf,pp)
  vm=rep(-Inf,pp)
  for(j in 1:pp){
   if(alp[j]<0) vm[j]=(b[j]-(A%*%y)[j]+alp[j]*sum(eta*y))/alp[j]
   if(alp[j]>0) vp[j]=(b[j]-(A%*%y)[j]+alp[j]*sum(eta*y))/alp[j]
   }
   vmm=max(vm,na.rm=T)
   vpp=min(vp,na.rm=T)
  return(list(vm=vmm,vp=vpp))
}

forwardStop=function(pv,alpha=.10){
 val=-(1/(1:length(pv)))*cumsum(log(1-pv))
 oo=which(val <= alpha)
 if(length(oo)==0) out=0
 if(length(oo)>0) out=oo[length(oo)]
return(out)
}
