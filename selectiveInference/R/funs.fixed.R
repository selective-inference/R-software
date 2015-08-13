
fixedLassoInf=function(x,y,bhat,lambda, sigma=NULL,coeftype=c("partial","full"),alpha=.10,verbose=FALSE,tol.beta=1e-5,tol.kkt=0.1,one.sided=TRUE,mingap=0.05,gridfac = 25, gridpts = 1000){
    # inference for fixed lam lasso
    # careful!  lambda is for usual lasso problem; glmnet uses lambda/n
    #assumes glmnet (or lasso solver)  is run with intercept=T and standardize=F
    this.call=match.call()
    coeftype=match.arg(coeftype)
    checkargs2(x=x,y=y,bhat=bhat,lambda=lambda,sigma=sigma,alpha=alpha,tol.beta=tol.beta,tol.kkt=tol.kkt)
     if(sum(bhat!=0)==0) stop("Value of lambda too large, bhat is zero")

    n=length(y)
    p=ncol(x)
    # Center predictors and ouctome
      x=scale(x,T,F)
      y=y-mean(y)
    
    tol.poly=0.01 # tolerance for checking polyedral lemma
    
##check KKT
   g0=t(x)%*%(y-x%*%bhat)
   g=g0-lambda*sign(bhat)
    gg=g0/lambda
    oo=abs(bhat)>tol.beta
    if(max(abs(g[oo]))>tol.kkt) cat(c("Warning: Solution bhat may not satisfy the KKT conditions"),fill=T)
    if(sum(!oo)>0){
     if(min(gg[!oo])< -1-tol.kkt | max(gg[!oo])>1 +tol.kkt)
         cat(c("Warning: Solution bhat may not satisfy the KKT conditions"),fill=T)
 }
##       
junk=tf.jonab(y,x,bhat,lambda,tol=tol.beta)
a=junk$A
b=junk$b

  
      if(is.null(sigma) & p>=n){cat("Warning: p ge n; the value sigma=1 used for error standard deviation; you may want    to estimate sigma using the estimateSigma function",fill=T); sigma=1}
  if(is.null(sigma) & n>p){
      sigma=sqrt(sum(lsfit(x,y)$res^2)/(n-p))
       cat("Standard deviation of noise estimated from mean squared residual",fill=T)
  }
 

    if(max(a%*%y-b)>tol.poly*sqrt(var(y))){
     stop("Polyhedral constraints not satisfied;
      you need to recompute bhat more
     accurately. With glmnet, be sure to use exact=TRUE in coef() and also try decreasing lambda.min
       in call to glmnet. If p>N, the value of lambda specified may also be too small---
      smaller than the value yielding zero training error")}
    
e=which(bhat!=0)
xe=x[,e,drop=F]

pp=length(e)
pv=vmall=vpall=rep(NA,pp)
ci=tailarea=matrix(NA,pp,2)
etaall=matrix(NA,nrow=pp,ncol=n)
SMALL=1e-7
ttall=rep(NA,pp)
    
 
for(k in 1:pp){
   #construct contrast vectors
    if(verbose) cat(k,fill=T)
    if (coeftype == "partial") 
            eta = ginv(t(xe))[, k]
        if (coeftype == "full") 
            eta = ginv(t(x))[,e[k]]
 
    bhat0=sum(eta*y)
    if(one.sided) eta=eta*sign(bhat0)
    flip=(one.sided & sign(bhat0)==-1)
    etaall[k,]=eta
 #compute truncation limits
  #  vs=compute.vmvp(y,eta,a,b,nrow(a))    
  #  vpp=vs$vp;vmm=vs$vm
  #  vmall[k]=vmm
  #  vpall[k]=vpp
   tt=sum(eta*y)
    ttall[k]=tt
 #  sigma.eta=sigma*sqrt(sum(eta^2))
    #compute p-values
 #  u=0  #null
#  pv[k]=  1-myptruncnorm(tt, vmm, vpp, u, sigma.eta)
    junk= poly.pval(y,-a,-b,eta,sigma) 
    pv[k] = junk$pv
    vmall[k]=junk$vlo
    vpall[k]=junk$vup
    if(!one.sided)  pv[k]=2*min(pv[k],1-pv[k])
    #compute selection intervals
    #  junk2=selection.int(y,eta,sigma,vs,alpha,flip=flip,mingap=mingap)
    #   ci[k,]=junk2$ci;miscov[k,]=junk2$miscov
    eta2=eta/sqrt(sum(eta^2))
   junk2= poly.int(y,-a,-b,eta2,sigma,alpha,gridfac=gridfac,gridpts=gridpts,
        flip=flip)

       ci[k,] = junk2$int
      tailarea[k,] = junk2$tailarea
                  
                
}
    
out=list(pv=pv,ci=ci,tailarea=tailarea,eta=etaall,vlo=vmall,vup=vpall,pred=e,alpha=alpha,sigma=sigma,one.sided=one.sided,lambda=lambda,coeftype=coeftype)
class(out)="fixedLassoInf"
    out$call=this.call
return(out)
}

print.fixedLassoInf=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat(c("lambda=",x$lam,", alpha=",x$alpha),fill=T)
      cat("",fill=T)
tab=cbind(x$pred,x$pv,x$ci,x$tailarea)
      dimnames(tab)=list(NULL,c("predictor","p-value","lowerConfPt","upperConfPt","lowerArea","upperArea"))
      print(tab)
       cat("",fill=T)
    cat(c("Value used for error standard deviation (sigma)=",round(x$sigma,6)),fill=T)
  }


### functions from Ryan
#tf.jonbands = function(y,x,beta,k,lambda,alpha,sigma2,verb=FALSE) {
 # ht = tf.ht(y,x,beta,k)
 # H = ht$H
 # theta = ht$theta
 # if (verb) cat(sprintf("Approx error: %f\n",max(abs(H%*%theta-beta))))

 # H1 = H[,1:(k+1),drop=FALSE]
 # H2 = H[,(k+2):n,drop=FALSE]
  #P = diag(1,nrow(H1)) - H1 %*% solve(t(H1) %*% H1) %*% t(H1)
 # yP = P %*% y
 # HP = P %*% H2
 # thetaP = theta[(k+2):n]
  
 # ab = tf.jonab(yP,HP,thetaP,lambda)
#  A = ab$A %*% P ## Important! Multiply by P here
#  b = ab$b

 # if (verb) cat(sprintf("Poly error: %f\n",min(b - A%*%y)))

 # n = length(y)
#  bands = cover = matrix(0,n,2)
#  for (i in 1:n) {
 #   if (verb) cat(sprintf("%i...",i))
  #  eta = rep(0,n); eta[i] = 1
    #eta[pmin(pmax(i+(-2):2,1),n)] = 1/5
 #   vs = tf.jonvs(y,A,b,eta)
#    bc = tf.jonint(y,eta,sigma2,vs,alpha)
 #   bands[i,] = bc[1:2]
 #   cover[i,] = bc[3:4]
#  }

#  return(list(bands=bands,cover=cover))
#}



tf.h = function(n,k,x=1:n) {
  if (k==0) return(lower.tri(matrix(0,n,n),diag=TRUE)+0)
  H = matrix(0,n,n)
  H[,1] = rep(1,n)
  for (j in 2:n) {
    if (j<=k+1) {
      H[,j] = apply(matrix(rep(x,each=j-1),ncol=n)-
         matrix(rep(x[1:(j-1)],n),ncol=n),2,prod)
    }
    else {
      H[,j] = apply(pmax(matrix(rep(x,each=k),ncol=n)-
         matrix(rep(x[j-k-1+1:k],n),ncol=n),0),2,prod)
    }
  }
  return(H)
}

#tf.theta = function(y,x,beta,k) {
#  n = length(y)
#  D = getDtfPos(n,k,x)
#  theta2 = D %*% beta / factorial(k)

#  alpha = beta[Seq(1,k+1)]
#  theta1 = numeric(k+1)
#  theta1[1] = alpha[1]
#  for (j in Seq(2,k+1)) {
#   delta = x[Seq(j,k+1)]-x[Seq(1,k+1-j+1)]
#    alpha = diff(alpha)/delta
#   theta1[j] = alpha[1]
# }
  
#  if (FALSE) {
# alpha = beta[Seq(1,k+1)]
#  for (j in Seq(1,k)) {
#   for (l in rev(Seq(k+1,j+1))) {
#     alpha[l] = (alpha[l]-alpha[l-1])/(x[l]-x[l-j])
#    }
# }
# theta1 = alpha
#  }
  
#return(c(theta1,theta2))
#}

Seq = function(a,b,...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}

tf.jonab = function(y,X,beta,lambda,tol=1e-5) {
    # compute constraint matrices  A and b for fixed lam lasso
  E = abs(beta)>tol
  XE = X[,E,drop=F]
  XEc = X[,!E,drop=F]
 # XEi = solve(t(XE) %*% XE)
#ROB changed this
  XEi = ginv(t(XE) %*% XE)
  XEp = XEi %*% t(XE)
  zE = sign(beta[E])
  if(length(zE)>1){dz=diag(zE)}
  if(length(zE)==1){dz=matrix(zE,ncol=1,nrow=1)}

  PEo = diag(1,nrow(XE)) - XE %*% XEp
  A = rbind(
    1/lambda * t(XEc) %*% PEo,
    -1/lambda * t(XEc) %*% PEo,
    -dz %*% XEp)
  b = c(
    1 - t(XEc) %*% t(XEp) %*% zE,
    1 + t(XEc) %*% t(XEp) %*% zE,
    -lambda * dz %*% XEi %*% zE)

  return(list(A=A,b=b))
}


##NOT USED
tf.jonvs = function(y,A,b,eta) {
  g = A %*% eta/sum(eta^2)
  f = b - A%*%y + g*sum(eta*y)
  vm = suppressWarnings(max((f/g)[g<0]))
  vp = suppressWarnings(min((f/g)[g>0]))
  vz = suppressWarnings(min(f[g==0]))
  return(list(vm=vm,vp=vp,vz=vz))
}

tf.jonint = function(y,eta,sigma2,vs,alpha,del=1e-3) {
  fun = function(x) return(tnorm.surv(sum(eta*y),x,sqrt(sigma2),vs$vm,vs$vp))
  inc = sqrt(sum(eta^2)*sigma2)*del
  lo = bin.search(sum(eta*y),fun,alpha/2,inc=inc)
  hi = bin.search(sum(eta*y),fun,1-alpha/2,inc=inc)
  return(c(lo,hi))
}
######
