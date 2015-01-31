library(genlasso)

fixedLamInf=function(x,y,bhat,lam,sigma,alpha=.05){
    # inference for fixed lam lasso
#assumes data is centered
junk=tf.jonab(y,x,bhat,lam)
n=length(y)
a=junk$A
b=junk$b
e=which(bhat!=0)
xe=x[,e,drop=F]

pp=length(e)
pv=vmall=vpall=rep(NA,pp)
ci=matrix(NA,pp,2)
etaall=matrix(NA,nrow=pp,ncol=n)
SMALL=1e-7
for(k in 1:pp){
eta=ginv(t(xe))[,k]
etaall[k,]=eta
vs=tf.jonvs(y,a,b,eta)
vpp=vs$vp;vmm=vs$vm
vmall[k]=vmm
vpall[k]=vpp
   tt=sum(eta*y)
   sigma.eta=sigma*sqrt(sum(eta^2))
   u=0  #null
      val0=pnorm((tt-u)/sigma.eta,log.p=TRUE)
      val1=pnorm((vmm-u)/sigma.eta,log.p=TRUE)
      val2=pnorm((vpp-u)/sigma.eta,log.p=TRUE)
  pv[k]=(1-exp(val2-val0))/(exp(val1-val0)-exp(val2-val0))
   #pv[k]=1-(pnorm((tt-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))/(pnorm((vpp-u)/sigma.eta)-pnorm((vmm-u)/sigma.eta))
    pv[k]=2*min(pv[k],1-pv[k])
    ci[k,]=tf.jonint(y,eta,sigma^2,vs,alpha)
}
return(list(lam=lam,eta=etaall,vm=vmall,vp=vpall,pv=pv,ci=ci))
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
  #XEi = solve(t(XE) %*% XE)
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

tf.jonint.rob = function(y,eta,sigma2,vs,alpha,del=1e-3) {
    #Rob's version using grid search
  fun = function(x) return(tnorm.surv(sum(eta*y),x,sqrt(sigma2),vs$vm,vs$vp))
  inc = sqrt(sum(eta^2)*sigma2)*del
  lo = rob(sum(eta*y),fun,alpha/2,inc=inc)
  hi = rob(sum(eta*y),fun,1-alpha/2,inc=inc)
  return(c(lo,hi))
}

tnorm.cdf = function(x,mean,sd,a,b) {
  return((pnorm(x,mean,sd)-pnorm(a,mean,sd))/
         (pnorm(b,mean,sd)-pnorm(a,mean,sd)))
}

tnorm.surv = function(x,mean,sd,a,b) {
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

rob=function(x, fun, val, inc=0.01, tol=1e-2,xl=x-10,xr=x+10) {
# note - hard coded 10 above
    #  here is used grid serach instead of binary search
G = 1000
  xg = seq(xl,xr,length=G)
  vals = numeric(G)
  for (i in 1:G) vals[i] = fun(xg[i])
  return(xg[which.min(abs(vals-val))])

}


