require(genlasso)
require(truncnorm)
require(MASS)

fixedLassoInf=function(x,y,bhat,lambda,sigma,alpha=.10,trace=F,compute.ci=F,tol=1e-5,one.sided=TRUE,gridfac=50){
    # inference for fixed lam lasso
#assumes data is centered
 # careful!  lambda is for usual lasso problem; glmnet uses lambda/n
    this.call=match.call()
junk=tf.jonab(y,x,bhat,lambda,tol=tol)
n=length(y)
a=junk$A
b=junk$b
e=which(bhat!=0)
xe=x[,e,drop=F]

pp=length(e)
pv=vmall=vpall=rep(NA,pp)
ci=miscov=matrix(NA,pp,2)
etaall=matrix(NA,nrow=pp,ncol=n)
SMALL=1e-7
for(k in 1:pp){
    if(trace) cat(k,fill=T)
eta=ginv(t(xe))[,k]
    bhat0=sum(eta*y)
    if(one.sided) eta=eta*sign(bhat0)
etaall[k,]=eta
   
vs=tf.jonvs(y,a,b,eta)
vpp=vs$vp;vmm=vs$vm
vmall[k]=vmm
vpall[k]=vpp
   tt=sum(eta*y)
   sigma.eta=sigma*sqrt(sum(eta^2))
   u=0  #null
#  pv[k]= 1-mytruncnorm(tt, vmm, vpp, sigma.eta, u)
  pv[k]=  1-ptruncnorm(tt, vmm, vpp, u, sigma.eta)
   
if(!one.sided)  pv[k]=2*min(pv[k],1-pv[k])
  if(compute.ci) { junk=selection.int(y,eta,sigma^2,vs,alpha,gridfac=gridfac)
                   ci[k,]=junk$ci;miscov[k,]=junk$miscov
               
}}
out=list(pv=pv,ci=ci,miscov=miscov,eta=etaall,vm=vmall,vp=vpall,pred=e,alpha=alpha,sigma=sigma,one.sided=one.sided,lambda=lambda)
class(out)="fixedLassoInf"
    out$call=this.call
return(out)
}



print.fixedLassoInf=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat(c("lambda=",x$lam,", alpha=",x$alpha),fill=T)
      cat("",fill=T)
tab=cbind(x$pred,x$pv,x$ci,x$miscov)
      dimnames(tab)=list(NULL,c("predictor","p-value","lowerConfPt","upperConfPt","lowerArea","upperArea"))
      print(tab)
  }


### functions from Ryan
tf.jonbands = function(y,x,beta,k,lambda,alpha,sigma2,verb=FALSE) {
  ht = tf.ht(y,x,beta,k)
  H = ht$H
  theta = ht$theta
  if (verb) cat(sprintf("Approx error: %f\n",max(abs(H%*%theta-beta))))

  H1 = H[,1:(k+1),drop=FALSE]
  H2 = H[,(k+2):n,drop=FALSE]
  P = diag(1,nrow(H1)) - H1 %*% solve(t(H1) %*% H1) %*% t(H1)
  yP = P %*% y
  HP = P %*% H2
  thetaP = theta[(k+2):n]
  
  ab = tf.jonab(yP,HP,thetaP,lambda)
  A = ab$A %*% P ## Important! Multiply by P here
  b = ab$b

  if (verb) cat(sprintf("Poly error: %f\n",min(b - A%*%y)))

  n = length(y)
  bands = cover = matrix(0,n,2)
  for (i in 1:n) {
    if (verb) cat(sprintf("%i...",i))
    eta = rep(0,n); eta[i] = 1
    #eta[pmin(pmax(i+(-2):2,1),n)] = 1/5
    vs = tf.jonvs(y,A,b,eta)
    bc = tf.jonint(y,eta,sigma2,vs,alpha)
    bands[i,] = bc[1:2]
    cover[i,] = bc[3:4]
  }

  return(list(bands=bands,cover=cover))
}



tf.ht = function(y,x,beta,k) {
  H = tf.h(length(y),k,x)
  theta = tf.theta(y,x,beta,k)
  return(list(H=H,theta=theta))
}

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

tf.theta = function(y,x,beta,k) {
  n = length(y)
  D = getDtfPos(n,k,x)
  theta2 = D %*% beta / factorial(k)

  alpha = beta[Seq(1,k+1)]
  theta1 = numeric(k+1)
  theta1[1] = alpha[1]
  for (j in Seq(2,k+1)) {
    delta = x[Seq(j,k+1)]-x[Seq(1,k+1-j+1)]
    alpha = diff(alpha)/delta
    theta1[j] = alpha[1]
  }
  
  if (FALSE) {
  alpha = beta[Seq(1,k+1)]
  for (j in Seq(1,k)) {
    for (l in rev(Seq(k+1,j+1))) {
      alpha[l] = (alpha[l]-alpha[l-1])/(x[l]-x[l-j])
    }
  }
  theta1 = alpha
  }
  
  return(c(theta1,theta2))
}

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



tnorm.cdf = function(x,mean,sd,a,b) {
  return((pnorm(x,mean,sd)-pnorm(a,mean,sd))/
         (pnorm(b,mean,sd)-pnorm(a,mean,sd)))
}

tnorm.surv = function(x,mean,sd,a,b) {
    #prob(X>x)
  return((pnorm(b,mean,sd)-pnorm(x,mean,sd))/
         (pnorm(b,mean,sd)-pnorm(a,mean,sd)))
}

bin.search = function(x, fun, val, inc=0.01, tol=1e-2) {  
  # Find left and right bounds
  xl = x
  while (fun(xl)>val) xl = xl-inc
  if (xl==x) {
    xr = x
    while (fun(xr)<val) xr = xr+inc
    xl = xr-inc
  }
  else xr = xl+inc
  
  # Grid search
  G = 100
  xg = seq(xl,xr,length=G)
  vals = numeric(G)
  for (i in 1:G) vals[i] = fun(xg[i])
  return(xg[which.min(abs(vals-val))])
 
}

grid.search=function(x, fun, val, xL,xR, inc=0.01, tol=1e-2) {
#
    #  here is used grid search instead of binary search
G = 1000
  xg = seq(xL,xR,length=G)
  vals = numeric(G)
  for (i in 1:G) vals[i] = fun(xg[i])
  return(xg[which.min(abs(vals-val))])

}



mytruncnorm = function(etay, vneg,vpos, sigma, etamu){
    # From Sam Gross- uses exp approximation in extreme tails
	# if range is too many sds away from mu, then there
	# will be numerical errors from using truncnorm
	if(max(vneg-etamu,etamu-vpos)/sigma < 7){
		     return(ptruncnorm(etay, vneg, vpos, etamu, sigma))
                 }
		   
	else{
	   
	      return(1 - pexp(vpos-etay, etamu-vpos)/ pexp(vpos-vneg, etamu-vpos))
	        
          }
    }

selection.int = function(y,eta,sigma,vs,alpha,del=1e-4,gridfac=50) {
    #Rob's version using grid search
    etay=sum(eta*y)
#    fun = function(x) return(tnorm.surv(etay,x,sigma,vs$vm,vs$vp))
    fun = function(x) return(1-ptruncnorm(etay,vs$vm,vs$vp,x,sigma))
 #   fun = function(x) return(1-mytruncnorm(etay,vs$vm,vs$vp,sigma,x))
  inc = sqrt(sum(eta^2)*sigma)*del
    sigma.eta=sqrt(sum(eta^2))*sigma
    xL=etay-gridfac*sigma.eta
    xR=etay+gridfac*sigma.eta
 #   cat(c(gridfac,etay,xL,sigma.eta),fill=T)
  lo = grid.search(etay,fun,alpha/2,xL,xR,inc=inc)
  hi = grid.search(etay,fun,1-alpha/2,xL,xR,inc=inc)
    covlo=fun(lo)
    covhi=1-fun(hi)
  return(list(ci=c(lo,hi), miscov=c(covlo,covhi)))
}

#tf.jonint.rob = function(y,eta,sigma2,vs,alpha,del=1e-3) {
#    #Rob's version using grid search
#  fun = function(x) return(tnorm.surv(sum(eta*y),x,sqrt(sigma2),vs$vm,vs$vp))
#  inc = sqrt(sum(eta^2)*sigma2)*del
#  lo = rob(sum(eta*y),fun,alpha/2,inc=inc)
#  hi = rob(sum(eta*y),fun,1-alpha/2,inc=inc)
#  return(c(lo,hi))
#}

forwardStop=function(pv,alpha=.10){
 val=-(1/(1:length(pv)))*cumsum(log(1-pv))
 oo=which(val <= alpha)
 if(length(oo)==0) out=0
 if(length(oo)>0) out=oo[length(oo)]
return(out)
}
