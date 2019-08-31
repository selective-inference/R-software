
fixedLogitLassoInf=function(x,y,beta,lambda,alpha=.1, type=c("partial,full"), tol.beta=1e-5, tol.kkt=0.1,
                            gridrange=c(-100,100), bits=NULL, verbose=FALSE,
                            linesearch.try=10, this.call=NULL){
  
  
  type = match.arg(type)
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  
  
  n=length(y)
  p=ncol(x)
  # I assume that intcpt was used
  if(length(beta)!=p+1) stop("Since family='binomial', beta must be of length ncol(x)+1, that is, it should include an intercept")
  vars = which(abs(beta[-1]) > tol.beta / sqrt(colSums(x^2)))
  nvar=length(vars)
  pv=vlo=vup=sd=rep(NA, nvar)
  ci=tailarea=matrix(NA,nvar,2)
  
  #do we need to worry about standardization?
  
  #  obj = standardize(x,y,TRUE,FALSE)
  #  x = obj$x
  #  y = obj$y
  
  bhat=c(beta[1],beta[-1][vars]) # intcpt plus active vars
  s2=sign(bhat)
  
  xm=cbind(1,x[,vars])
  xnotm=x[,-vars]
  
  etahat = xm %*% bhat
  prhat = as.vector(exp(etahat) / (1 + exp(etahat)))
  ww=prhat*(1-prhat)
  w=diag(ww)
  
  #check KKT
  z=etahat+(y-prhat)/ww
  g=scale(t(x),FALSE,1/ww)%*%(z-etahat)/lambda # negative gradient scaled by lambda
  if (any(abs(g) > 1+tol.kkt) )
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))
  
  if(length(vars)==0){
    cat("Empty model",fill=T)
    return()
  }
  if (any(sign(g[vars]) != sign(beta[-1][vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances). You might try rerunning",
                  "glmnet with a lower setting of the",
                  "'thresh' parameter, for a more accurate convergence."))
  
  #constraints for active variables             
  MM=solve(scale(t(xm),F,1/ww)%*%xm)
  gm = c(0,-g[vars]*lambda) # gradient at LASSO solution, first entry is 0 because intercept is unpenalized
  # at exact LASSO solution it should be s2[-1]
  dbeta = MM %*% gm
  
  MM = MM*n
  
  # bbar=(bhat+lam2m%*%MM%*%s2)  # JT: this is wrong, shouldn't use sign of intercept anywhere...
  bbar = (bhat - dbeta)*sqrt(n)
  
  A1= matrix(-(mydiag(s2))[-1,],nrow=length(s2)-1)
  b1= ((s2 * dbeta)[-1])*sqrt(n)
  V = (diag(length(bbar))[-1,])/sqrt(n)
  null_value = rep(0,nvar)
  
  if (type=='full') {
    
    is_wide = n < (2 * p) # somewhat arbitrary decision -- it is really for when we don't want to form with pxp matrices
    
    # Approximate inverse covariance matrix for when (n < p) from lasso_Inference.R
    if (!is_wide) {
      M = debiasingMatrix(1/n*(scale(t(x),FALSE,1/ww)%*%x), is_wide, n, vars, verbose=FALSE, max_try=linesearch.try, warn_kkt=TRUE)
    } else {
      M = debiasingMatrix(t(scale(t(x),1/sqrt(ww))), is_wide, n, vars, verbose=FALSE, max_try=linesearch.try, warn_kkt=TRUE)
    }
    
    #M <- matrix(InverseLinfty(hsigma,n,dim(xm)[2],verbose=F,max.try=10),ncol=p+1)[-1,] # remove intercept row
    I <- matrix(diag(dim(xm)[2])[-1,],nrow=dim(xm)[2]-1)
    if (is.null(dim(M))) {
      M_notE <- M[-vars]
    } else {
      M_notE <- M[,-vars]
    }
    M_notE = matrix(M_notE,nrow=nvar)
    V <- matrix(cbind(I/sqrt(n),M_notE[,-1]/n),nrow=dim(xm)[2]-1)
    
    xnotm_w = scale(t(xnotm),FALSE,1/ww)
    xnotm_w_xm = xnotm_w%*%xm
    c <- matrix(c(gm[-1],xnotm_w_xm%*%(-dbeta)),ncol=1)
    d <- -dbeta[-1]
    null_value = -(M[,-1]%*%c/n - d)
    
    A0 = matrix(0,ncol(xnotm),length(bbar))
    A0 = cbind(A0,diag(nrow(A0)))
    fill = matrix(0,nrow(A1),ncol(xnotm))
    A1 = cbind(A1,fill)
    A1 = rbind(A1,A0,-A0)
    
    b1 = matrix(c(b1,rep(lambda,2*nrow(A0))),ncol=1)
    
    # full covariance
    MMbr = (xnotm_w%*%xnotm - xnotm_w_xm%*%(MM/n)%*%t(xnotm_w_xm))*n
    MM = cbind(MM,matrix(0,nrow(MM),ncol(MMbr)))
    MMbr = cbind(matrix(0,nrow(MMbr),nrow(MM)),MMbr)
    MM = rbind(MM,MMbr)
    
    etahat_bbar = xm %*% (bbar/sqrt(n))
    gnotm = (scale(t(xnotm),FALSE,1/ww)%*%(z-etahat_bbar))*sqrt(n)
    bbar = matrix(c(bbar,gnotm),ncol=1)
  }
  
  if (is.null(dim(V))) V=matrix(V,nrow=1)
  
  tol.poly = 0.01 
  if (max((A1 %*% bbar) - b1) > tol.poly)
    stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
               "and check whether the specified value of lambda is too small",
               "(beyond the grid of values visited by glmnet).",
               "You might also try rerunning glmnet with a lower setting of the",
               "'thresh' parameter, for a more accurate convergence."))
  
  sign=numeric(nvar)
  coef0=numeric(nvar)
  
  for(j in 1:nvar){
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
    
    if (is.null(dim(V))) vj = V
    else vj = matrix(V[j,],nrow=1)
    coef0[j] = vj%*%bbar # sum(vj * bbar)
    sign[j] = sign(coef0[j])
    vj = vj * sign[j]
    
    # compute p-values
    limits.info = TG.limits(bbar, A1, b1, vj, Sigma=MM)
    # if(is.null(limits.info)) return(list(pv=NULL,MM=MM,eta=vj))
    a = TG.pvalue.base(limits.info, null_value=null_value[j], bits=bits)
    pv[j] = a$pv
    if (is.na(s2[j])) { # for variables not in the active set, report 2-sided pvalue
      pv[j] = 2 * min(pv[j], 1 - pv[j])
    }  

    vlo[j] = a$vlo # * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup # * mj # Unstandardize (mult by norm of vj)
    sd[j] = a$sd # * mj # Unstandardize (mult by norm of vj)
    if (type=='full') { # rescale because of local alternatives
      vlo[j] = vlo[j]/sqrt(n)
      vup[j] = vup[j]/sqrt(n)
      sd[j] = sd[j]/sqrt(n)
    }
    
    a = TG.interval.base(limits.info, 
                         alpha=alpha,
                         gridrange=gridrange,
                         flip=(sign[j]==-1),
                         bits=bits)
    ci[j,] = (a$int-null_value[j]) # * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }
  
  se0 = sqrt(diag(V%*%MM%*%t(V)))
  zscore0 = (coef0+null_value)/se0
  
  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
             tailarea=tailarea,vlo=vlo,vup=vup,sd=sd,
             vars=vars,alpha=alpha,coef0=coef0,zscore0=zscore0,
             call=this.call,
             info.matrix=MM) # info.matrix is output just for debugging purposes at the moment
  class(out) = "fixedLogitLassoInf"
  return(out)
  
}



print.fixedLogitLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
              round(x$coef0,3),
              round(x$zscore0,3),
              round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  
  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}
