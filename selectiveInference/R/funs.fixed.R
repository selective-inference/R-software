# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

lassoInf <- function(x, y, beta, lambda, intercept=TRUE, sigma=NULL, alpha=0.1,
                     type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                     gridrange=c(-100,100), gridpts=1000, verbose=FALSE) {
  
<<<<<<< HEAD
  this.call = match.call()
  type = match.arg(type)
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                 gridrange=gridrange,gridpts=gridpts,
                 tol.beta=tol.beta,tol.kkt=tol.kkt)
  n = nrow(x)
  p = ncol(x)
  beta = as.numeric(beta)
  if (length(beta) != p) stop("beta must have length equal to ncol(x)")
  
  # If glmnet was run with an intercept term, center x and y
  if (intercept==TRUE) {
    obj = standardize(x,y,TRUE,FALSE)
    x = obj$x
    y = obj$y
=======
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
>>>>>>> b9bfc4fb52c176d0abc6b23ebe952ccfdd2677fb
  }
  
  # Check the KKT conditions
  g = t(x)%*%(y-x%*%beta) / lambda
  if (any(abs(g) > 1+tol.kkt))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))
  vars = which(abs(beta) > tol.beta)
  if (any(sign(g[vars]) != sign(beta[vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))
  
  # Get lasso polyhedral region, of form Gy >= u
  out = lasso.poly(x,y,beta,lambda,tol.beta)
  G = out$G
  u = out$u

  # Check polyhedral region
  tol.poly = 0.01 
  if (min(G %*% y - u) < -tol.poly * sd(y))
    stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, be sure to use exact=TRUE in coef().",
               "Also, it may be the case that the value of lambda is too small",
               "(beyond the grid of values visited by glmnet)"))

  # Estimate sigma
  if (is.null(sigma)) {
    if (n >= 2*p) sigma = sqrt(sum(lsfit(x,y,intercept=F)$res^2)/(n-p))
    else {
      sigma = sd(y)
      warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                    "you may want to use the estimateSigma function"))
    }
  }
 
  k = length(vars)
  pv = vlo = vup = numeric(k) 
  vmat = matrix(0,k,n)
  ci = tailarea = matrix(0,k,2)
  sign = numeric(k)
  
  if (type=="partial") {
    xa = x[,vars,drop=F]
    M = pinv(crossprod(xa)) %*% t(xa)
  }
  else {
    M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=FALSE]
  }
  
  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
    
    vj = M[j,]
    sign[j] = sign(sum(vj*y))
    vj = vj / sqrt(sum(vj^2))
    vj = sign[j] * vj

    a = poly.pval(y,G,u,vj,sigma)
    pv[j] = a$pv
    vlo[j] = a$vlo
    vup[j] = a$vup
    vmat[j,] = vj

    a = poly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
      gridpts=gridpts,flip=(sign[j]==-1))
    ci[j,] = a$int
    tailarea[j,] = a$tailarea
  }
  
  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    call=this.call)
  class(out) = "lassoInf"
  return(out)
}

##############################

lasso.poly <- function(x, y, beta, lambda, tol.beta=1e-5) {
  a = abs(beta) > tol.beta
  xa = x[,a,drop=F]
  xac = x[,!a,drop=F]
  xai = pinv(crossprod(xa))
  xap = xai %*% t(xa)
  za = sign(beta[a])
  if (length(za)>1) dz = diag(za)
  if (length(za)==1) dz = matrix(za,1,1)

  P = diag(1,nrow(xa)) - xa %*% xap
  G = -rbind(1/lambda * t(xac) %*% P,
    -1/lambda * t(xac) %*% P,
    -dz %*% xap)
  u = -c(1 - t(xac) %*% t(xap) %*% za,
    1 + t(xac) %*% t(xap) %*% za,
    -lambda * dz %*% xai %*% za)

  return(list(G=G,u=u))
}

# Moore-Penrose pseudo inverse for symmetric matrices

pinv <- function(A, tol=.Machine$double.eps) {
  e = eigen(A)
  v = Re(e$vec)
  d = Re(e$val)
  d[d > tol] = 1/d[d > tol]
  d[d < tol] = 0
  return(v %*% diag(d) %*% t(v))
}

##############################

print.lassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$sign*x$vmat%*%x$y,3),round(x$sign*x$vmat%*%x$y/x$sigma,3),
    round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "StdzCoef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)
 
  cat(sprintf("\nNote: displayed are coefficients in the %s regression model\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}

