
forwardStep=function(x,y,sigma=NULL,intercept=TRUE,nsteps=min(nrow(x)-1,ncol(x)),trace=FALSE){
    this.call=match.call()
    checkargs(x=x,y=y,nsteps=nsteps,sigma=sigma)
    BIG=10e9
    n=nrow(x)
p=ncol(x)
    
# fs by minimizing scaled ip
# first center x and y if specified
    if(intercept){
      x=scale(x,T,F)
      y=y-mean(y)
     }
    if(is.null(sigma) & p>=n){cat("Warning: p ge n; the value sigma=1 used for AIC; you may want to estimate sigma using the estimateSigma function",fill=T); sigma=1}
  if(is.null(sigma) & n>p){
      sigma=sqrt(sum(lsfit(x,y,intercept=FALSE)$res^2)/(n-p-1*intercept))
      cat("Standard deviation of noise estimated from mean squared residual",fill=T)
  }
pred=s=scor=bhat=rep(NA,nsteps)
   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
    if(trace) cat("step=1",fill=T)
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y,intercept=FALSE)$res
    
 if(nsteps>1){
for(j in 2:nsteps){
     if(trace) cat(c("step=",j),fill=T)
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r,intercept=FALSE)$res
  xr= lsfit(x[,mod],x,intercept=FALSE)$res
  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  ip[mod]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
 bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
}}
    
#compute AIC
    aic=rss=val=rep(NA,nsteps+1)
rss[1]=sum( (y-mean(y))^2)
    aic[1]=rss[1]
for(j in 1:nsteps){
      xx=x[,pred[1:j],drop=F]
  rss[j+1]=sum(lsfit(xx,y,intercept=FALSE)$res^2)
  aic[j+1]=rss[j+1]+2*j*sigma^2
 }
 d=diff(aic)
o=which(d>0)
    if(length(o)>0) aichat=o[1]-1
    if(length(o)==0) aichat=length(aic)-1
    out=list(pred=pred,scor=scor,bhat=bhat,rss=rss,sigma=sigma,intercept=intercept,aic=aic,aichat=aichat,call=this.call)
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
#  returns pval and interval for predictor just entered (default with fixed.step=NULL) 
# otherwise if fixed.step=k, returns results for all predictors at step k
#  
    this.call=match.call()
    
      checkargs(x=x,y=y,nsteps=nsteps,sigma=sigma,alpha=alpha, fsfit=fsfit)
    
n=nrow(x)
p=ncol(x)
    intercept=fsfit$intercept
    
  if(is.null(sigma)) sigma=fsfit$sigma
    if(aic.stop & !is.null(fixed.step)) stop("Only one of fixed.step and aic.stop can be specified")
    if(aic.stop){fixed.step=fsfit$aichat}
if(is.null(nsteps)){ nsteps=min(length(fsfit$pred),20)}
    if(!is.null(fixed.step))  nsteps=fixed.step
    
# first center x and y if specified in original call

    if(intercept){
    x=scale(x,TRUE,FALSE)
    y=y-mean(y)
}
xx=x
 
SMALL=1e-7
pred=fsfit$pred
s=sign(fsfit$scor)
a=NULL
proj=0
which.steps=1:nsteps
    if(!is.null(fixed.step)){
     if(!is.wholenumber(fixed.step)) stop("fixed.step must be an integer")
     if(fixed.step<1 | fixed.step>ncol(x)) stop("fixed.step must be between 1 and number of predictors")
 }
if(!is.null(fixed.step)){ which.steps=rep(fixed.step,fixed.step)}

if(aic.stop & fsfit$aichat==0){
    out=list(pv=rep(NA,p),ci=NA,A=NA,b=NA)
    cat("Null model picked by AIC; no  results returned",fill=T)
    return(out)
}

# construct matrices A and b
    vv=1/sum(xx[,pred[1]]^2)


# A=b=stepind=NULL
    nr=2*(nsteps*p-sum(1:nsteps))
    if(is.null(fixed.step)) nr=nr+1
    if(aic.stop & fsfit$aichat<ncol(x)) nr=nr+nsteps+1   
A=matrix(NA,nrow=nr,ncol=n)
    b=stepind=rep(NA,nr)

if(trace) cat("Constructing constraint matrix",fill=T)
    ii=0
for(j in 1:nsteps){
if(trace) cat(c("step=",j),fill=T)
  notin=rep(T,p)
  if(j>1) {
  mod=pred[1:(j-1)]
  notin[mod]=F
  xx[,-mod]=lsfit(xx[,mod],xx[,-mod],intercept=FALSE)$res
 # vv=solve(t(xx[,mod,drop=F])%*%xx[,mod,drop=F])
 # proj=xx[,mod]%*%vv%*%t(xx[,mod])
  }
  jj=pred[j]
  rest=notin;rest[jj]=F
 w=sqrt(sum(xx[,jj]^2))
  for(k in which(rest)){
 w2=sqrt(sum(xx[,k]^2))
  ii=ii+1
 #A[ii,]=-(s[j]*xx[,jj]/w-xx[,k]/w2)%*%(diag(n)-proj)
 A[ii,]=-(s[j]*xx[,jj]/w-xx[,k]/w2)
 if(j>1) A[ii,]=lsfit(xx[,mod],A[ii,],intercept=FALSE)$res
 b[ii]=0;  stepind[ii]=j
 ii=ii+1
 #A[ii,]=-(s[j]*xx[,jj]/w+xx[,k]/w2)%*%(diag(n)-proj)
  A[ii,]=-(s[j]*xx[,jj]/w+xx[,k]/w2)
 if(j>1) A[ii,]=lsfit(xx[,mod],A[ii,],intercept=FALSE)$res
 b[ii]=0;  stepind[ii]=j
}
 if(aic.stop){
     #constraint so that change in AIC > 2sigma
      mod2=pred[1:j]
       vv=solve(t(xx[,mod2,drop=F])%*%xx[,mod2,drop=F])
   temp=vv%*%t(xx[,mod2,drop=F])
      temp2=-s[j]*temp[j,,drop=F]/(sigma*sqrt(vv[j,j]))
      ii=ii+1
      A[ii,]=temp2
      b[ii]=-sqrt(2)
      stepind[ii]=j
}
}

    
# last step- constraint on sign only (for fixed.step=NULL only)
if(is.null(fixed.step)){
mod=pred[1:nsteps]
fitter=solve(t(xx[,mod])%*%xx[,mod])%*%t(xx[,mod])
aa=-1*fitter[nsteps,]*s[nsteps]
ii=ii+1
A[ii,]=aa
b[ii]=0
stepind[ii]=nsteps
}
    
#for aic stop, add constraint that diff in AIC is < 2sigma, so that it stops
 # (equiv to bhat/se < sqrt(2))
if(aic.stop & fsfit$aichat<ncol(x)){
      mod2=pred[1:(nsteps+1)]
 vv=solve(t(xx[,mod2,drop=F])%*%xx[,mod2,drop=F])
     pp=length(mod2)
      temp=vv%*%t(xx[,mod2,drop=F])
     snew=sign(sum(temp[pp,]*y))
     temp2=snew*temp[pp,,drop=F]/(sigma*sqrt(vv[pp,pp]))
   #  A=rbind(A, temp2)
#  b=c(b,sqrt(2))
      #stepind=nsteps
      ii=ii+1
      A[ii,]=temp2
b[ii]=sqrt(2)
stepind[ii]=nsteps
 }
    
    
# pv tests

vall = matrix(NA,n,nsteps)
vmall=vpall=rep(NA,nsteps)
pv=rep(NA,nsteps)
ci=miscov=matrix(NA,nrow=p,ncol=2)

 
if(trace & !compute.si) cat("Computing p-values",fill=T)
    if(trace & compute.si) cat("Computing p-values and selection intervals",fill=T)
for(kk in 1:length(which.steps)){
if(trace) cat(c("step=",kk),fill=T)
if(is.null(fixed.step)) pp=sum(stepind<=which.steps[kk])
if(!is.null(fixed.step)) pp=nrow(A)
  mod=pred[1:which.steps[kk]]

 contr=(solve(t(x[,mod,drop=F])%*%x[,mod,drop=F])%*%t(x[,mod,drop=F]))

# compute pvalues and CIs  
kkk=nrow(contr)
A2=A;b2=b
if(!is.null(fixed.step)) kkk=kk
  eta=as.vector(contr[kkk,])

    bhat=sum(eta*y)


    if(one.sided) eta=eta*sign(bhat)
vall[,kk] = eta
flip=(one.sided & sign(bhat)==-1)

if(!is.null(fixed.step)){  #add sign constraint for predictor being tested
  A2=rbind(A2,-1*contr[kkk,]*sign(bhat))
  b2=c(b2,0)
  }
    
junk=compute.vmvp(y,eta,A2,b2,pp)
vmm=junk$vm;vpp=junk$vp

   vmall[kk]=vmm
   vpall[kk]=vpp
   tt=sum(eta*y)
   sigma.eta=sigma*sqrt(sum(eta^2))
   u=0  #null
  # pv[kk]=1-(pnorm((tt-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))/(pnorm((vpp-u)/sigma.eta)-pnorm((vmm-u)/sigma.
  pv[kk]= 1-myptruncnorm(tt, vmm, vpp, u,sigma.eta)
if(!one.sided)  pv[kk]=2*min(pv[kk],1-pv[kk])

  
  if(compute.si)
      {
           vs=list(vm=vmm,vp=vpp)
          junk=selection.int(y,eta,sigma,vs,alpha,flip=flip)
#     cat(c(vs$vm,sum(eta*y),vs$vp,sigma,sigma.eta,alpha),fill=T)
          ci[kk,]=junk$ci;miscov[kk,]=junk$miscov
      }


}
    
forwardStopHat=NULL
#if(is.null(fixed.step)) forwardStopHat=forwardStop(pv,alpha)
        
        
out=list(pv=pv,v=vall,vm=vmall,vp=vpall,ci=ci,tailarea=miscov,pred=pred,which.steps=which.steps,stepind=stepind,forwardStopHat=forwardStopHat,alpha=alpha,sigma=sigma,one.sided=one.sided,A=A,b=b,call=this.call)

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
tab=cbind(x$which.steps,x$pred[1:nn],round(x$pv[1:nn],6),round(x$ci[1:nn,,drop=FALSE],6),round(x$tailarea[1:nn,,drop=FALSE],6))
      
      
        dimnames(tab)=list(NULL,c("step","predictor","p-value","lowerConfPt","upperConfPt","lowerArea","upperArea"))
      print(tab)

if(!is.null(x$forwardStopHat)){
cat("",fill=T)
cat(c("Estimated stopping point from forwardStop rule=", x$forwardStopHat),fill=T)
         cat("",fill=T)
}
      cat(c("Value used for error standard deviation (sigma)=",round(x$sigma,6)),fill=T)
  }
  
  





forwardStop=function(pv,alpha=.10){
    if(alpha<=0 | alpha>=1) stop("alpha must be in [0,1]")
    if(min(pv,na.rm=T)<0 | max(pv,na.rm=T)>1) stop("pvalues must be in [0,1]")
 val=-(1/(1:length(pv)))*cumsum(log(1-pv))
 oo=which(val <= alpha)
 if(length(oo)==0) out=0
 if(length(oo)>0) out=oo[length(oo)]
return(out)
}
